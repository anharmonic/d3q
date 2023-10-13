!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! WIP: compute the scattering terms in two ways to enforce detailed balancee
!
MODULE variational_tk_symq
#include "mpi_thermal.h"
  USE kinds,           ONLY : DP
  USE more_constants,  ONLY : eps_vel, eps_freq
  USE q_grids,         ONLY : q_grid
  USE mpi_thermal,     ONLY : ionode
  USE posix_signal,    ONLY : check_graceful_termination
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute the diagonal matrix A_out
  SUBROUTINE compute_A_out_symq(A_out, input, basis, out_grid, in_grid, S, fc2, fc3)
    USE constants,          ONLY : RY_TO_CMM1
    USE linewidth,          ONLY : linewidth_q
    USE q_grids,            ONLY : q_basis, q_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE mpi_thermal,        ONLY : mpi_bsum
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)     :: basis
    TYPE(q_grid),INTENT(in)      :: out_grid, in_grid
    REAL(DP),INTENT(out) :: A_out(input%nconf,S%nat3, out_grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    REAL(DP) :: A_out_aux(S%nat3, input%nconf)
    REAL(DP) :: lw_casimir(S%nat3)
    !
    INTEGER  :: iq, it, ix, nu
    !
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    lw_casimir  = 0._dp
    A_out_aux = 0._dp

    QPOINT_LOOP : &
    DO iq = 1,out_grid%nq
      ! 
      CALL print_percent_wall(20._dp, 600._dp, iq, out_grid%nq, (iq==1))
      !
      timer_CALL t_lwphph%start()
      !
      ! Compute the ph-ph and the isotope scattering terms
      A_out_aux = A_out_q_symq(out_grid%xq(:,iq), input%nconf, input%T, sigma_ry, S, in_grid, fc2, fc3, &
                           input%isotopic_disorder)
      timer_CALL t_lwphph%stop()
      !
      ! Exchange the order of indexes (the new order is better for the CG algorithm)
      DO it = 1, input%nconf
      DO nu = 1, S%nat3
        A_out(it,nu,iq) = A_out_aux(nu,it)
      ENDDO
      ENDDO
      !
      !
    ENDDO QPOINT_LOOP
    !
  END SUBROUTINE compute_A_out_symq
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute 3*nat lines of the matrix A_out (those corresponding to the
  ! q-point xq0), including both ph-ph and isotopic scattering processes
  FUNCTION A_out_q_symq(xq0, nconf, T, sigma, S, grid, fc2, fc3, isotopic_disorder) &
    RESULT(A_out_q)
    USE q_grids,            ONLY : q_basis
    USE fc2_interpolate,    ONLY : forceconst2_grid, bose_phq, set_nu0, &
                                   freq_phq_safe, bose_phq
    USE fc3_interpolate,    ONLY : forceconst3, ip_cart2pat
    USE isotopes_linewidth, ONLY : sum_isotope_scattering_modes
    USE input_fc,           ONLY : ph_system_info
    USE mpi_thermal,        ONLY : mpi_bsum
    USE merge_degenerate,   ONLY : merge_degen
    USE timers
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    ! CAREFUL! grid%nq is scattered over MPI, basis%nq is the total number of q vectors
    REAL(DP),INTENT(in)  :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    LOGICAL,INTENT(in)  :: isotopic_disorder
    !
    ! FUNCTION RESULT:
    REAL(DP) :: A_out_q(S%nat3,nconf)
    !
    REAL(DP) :: A_out(S%nat3,nconf)
    REAL(DP) :: A_out_isot(S%nat3,nconf)
!     REAL(DP) :: lw(S%nat3,nconf)
    REAL(DP) :: P3_isot(S%nat3,S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:), V3Bsq(:,:,:)
    INTEGER :: iq, jq, mu, nu, it, nu0(4), ix
    !
    REAL(DP) :: freq(S%nat3,5), bose(S%nat3,5), xq(3,3)
    REAL(DP),PARAMETER :: epsq = 1.e-12_dp
    INTEGER, PARAMETER :: a=1, b=2, c=3
    !
    ALLOCATE(U(S%nat3, S%nat3,5))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(V3Bsq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    A_out = 0._dp
    A_out_isot = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
      timer_CALL t_freq%start()
    xq(:,a) = xq0
    nu0(1)  = set_nu0(xq(:,1), S%at)
    !freq(:,1) = basis%w(:,iq0)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      timer_CALL t_freq%stop()
    !
    QPOINT_INNER_LOOP : &
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
        timer_CALL t_freq%start()
        xq(:,c) = grid%xq(:,iq)
        xq(:,b) = xq(:,3)-xq(:,1)
      !
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      U(:,:,5) = CONJG(U(:,:,2))
        timer_CALL t_freq%stop()
      !
      ! Interpolate D3(q1,q2,-q1-q2)
        timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,b), xq(:,c), S%nat3, D3)
        timer_CALL t_fc3int%stop()
        timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
        timer_CALL t_fc3rot%stop()
        timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
        timer_CALL t_fc3m2%stop()
      !
      CONF_LOOP : &
      DO it = 1,nconf
          timer_CALL t_bose%start()
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
          timer_CALL t_bose%stop()
        !
          timer_CALL t_sum%start()
        A_out(:,it) =  A_out(:,it) + grid%w(iq)&
                    *sum_A_out_modes_symq( S%nat3, sigma(it), freq, bose, V3sq, V3Bsq, nu0)
          timer_CALL t_sum%stop()
        !
      ENDDO &
      CONF_LOOP
      !
    ENDDO &
    QPOINT_INNER_LOOP
    !
    A_out = A_out + A_out_isot
    !
    ! Recollect over MPI processes if necessary
    IF(grid%scattered) THEN
        timer_CALL t_mpicom%start()
      CALL mpi_bsum(S%nat3, nconf, A_out)
        timer_CALL t_mpicom%stop()
    ENDIF
    DO it = 1, nconf
      CALL merge_degen(S%nat3, A_out(:,it), freq(:,1))
    ENDDO
    !
    ! Assign to output variable (better not to pass it to mpi)
    A_out_q = A_out
    !
    CALL check_graceful_termination
    !
    DEALLOCATE(U, V3sq, V3Bsq, D3)
    !
  END FUNCTION A_out_q_symq
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the scattering probability required for the matrix A.
  ! The form is a bit more complex than in  PRB 88, 045430 (2013)
  ! because we have rearranged the indexes to only have cohalescence
  ! processes instead of mixed scattering-cohalescence. This trick ensures
  ! that everything is positive definite.
  FUNCTION sum_A_out_modes_symq(nat3, sigma, freq, bose, V3sq, V3Bsq, nu0) &
    RESULT(sum_A_out_modes)
    USE functions,        ONLY : f_gauss => f_gauss
    USE constants,        ONLY : tpi, RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(nat3,3)
    REAL(DP),INTENT(in) :: bose(nat3,3)
    REAL(DP),INTENT(in) :: V3sq(nat3,nat3,nat3)
    REAL(DP),INTENT(in) :: V3Bsq(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    REAL(DP) :: sum_A_out_modes(nat3)
    !
    REAL(DP) :: bose_a, bose_b, bose_c
    REAL(DP) :: dom_a, dom_b, dom_c
    REAL(DP) :: ctm_a, ctm_b, ctm_c
    REAL(DP) :: norm_a, norm_b, norm_bc
    REAL(DP) :: sum_a, sum_b, sum_c, sum_ac
    REAL(DP) :: freqm1(nat3,4)
    !REAL(DP),SAVE :: leftover_e = 0._dp
    !
    INTEGER :: i,j,k
    REAL(DP) :: sum_A_out(nat3)
    sum_A_out(:) = 0._dp

    freqm1 = 0._dp
    DO i = 1,nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
    !
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP          PRIVATE(i, j, k, bose_a, bose_c, dom_a, dom_c) &
!$OMP          PRIVATE(ctm_a, ctm_c, sum_a, sum_ac, norm_a, norm_c) &
!$OMP          REDUCTION(+: sum_A_out)
!$OMP DO  COLLAPSE(3)
    DO k = 1,nat3
      DO j = 1,nat3
        norm_bc = tpi*freqm1(j,2)*freqm1(k,3)
        DO i = 1,nat3
          !
          bose_b = bose(i,1) * bose(j,2) * (bose(k,3)+1)
          bose_c = (bose(i,1)+1) * bose(j,2) * bose(k,3)
          !
          dom_b =  freq(i,1) + freq(j,2) - freq(k,3)
          dom_c = -freq(i,1) + freq(j,2) + freq(k,3)
          !
          ctm_b = bose_b *  f_gauss(dom_b, sigma)
          ctm_c = bose_c *  f_gauss(dom_c, sigma)
          !
          sum_b = norm_bc*freqm1(i,1)* (ctm_b * V3sq(i,k,j) + ctm_c * V3sq(j,k,i))
          !
          sum_A_out(i) = sum_A_out(i) + sum_b + 0.5_dp*sum_c
          !
         !leftover_e = sum_a*dom_a + 0.5_dp*sum_c*dom_c
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
    ! Merge degenerate is done outside this function for better efficiency
    !CALL merge_degen(nat3, sum_A_out, freq(:,1))
    sum_A_out_modes = sum_A_out
    !
  END FUNCTION sum_A_out_modes_symq
  !
  !
END MODULE variational_tk_symq
