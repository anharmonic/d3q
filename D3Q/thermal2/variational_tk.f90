!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE variational_tk
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
  ! Multiply matrix in diagonal form, like A_out, A_out^-1, A_out^-1/2
  ! with a vector (f, b, etc)
  FUNCTION A_diag_f(A, f, nconf, nat3, nq) RESULT(Af)
    IMPLICIT NONE
    !
    REAL(DP):: Af(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, ix, nu
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nu,it,ix)
    DO iq = 1,nq
!$OMP DO COLLAPSE(3)
    DO nu = 1,nat3
      DO it = 1,nconf
      DO ix = 1,3
        Af(ix,it,nu,iq) = f(ix,it,nu,iq)*A(it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL
    !
  END FUNCTION A_diag_f
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the diagonal matrix A_out
  SUBROUTINE compute_A_out(A_out, input, basis, out_grid, in_grid, S, fc2, fc3)
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
      timer_CALL t_lwphph%start()
      !
      ! Compute the ph-ph and the isotope scattering terms
      A_out_aux = A_out_q(out_grid%xq(:,iq), input%nconf, input%T, sigma_ry, S, in_grid, fc2, fc3, &
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
      ! Compute potentially Casimir linewidth, and get the scattering term as P^be = n(n+1) gamma^be
      IF(input%casimir_scattering) THEN
        timer_CALL t_lwcasi%start()
        lw_casimir = casimir_linewidth_vel( basis%c(:,:,iq), input%casimir_length, input%casimir_dir, S%nat3)
        !
        ! Casimir linewidth is temperature/smearing-independent, sum it to all configurations
        DO it = 1,input%nconf
        DO nu = 1, S%nat3
          A_out(it,nu,iq) = A_out(it,nu,iq) &
                          +lw_casimir(nu) *basis%be(nu,it,iq)*( basis%be(nu,it,iq)+1)
        ENDDO
        ENDDO
        timer_CALL t_lwcasi%stop()
      ENDIF
      !
    ENDDO QPOINT_LOOP
    !
  END SUBROUTINE compute_A_out
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute 3*nat lines of the matrix A_out (those corresponding to the
  ! q-point xq0), including both ph-ph and isotopic scattering processes
  FUNCTION A_out_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, isotopic_disorder)
    USE q_grids,            ONLY : q_basis
    USE fc2_interpolate,    ONLY : forceconst2_grid, bose_phq, set_nu0, &
                                   freq_phq_safe, bose_phq, ip_cart2pat
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : sum_isotope_scattering_modes
    USE input_fc,           ONLY : ph_system_info
    USE mpi_thermal,        ONLY : mpi_bsum
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
    REAL(DP) :: A_out(S%nat3,nconf)
    REAL(DP) :: A_out_isot(S%nat3,nconf)
!     REAL(DP) :: lw(S%nat3,nconf)
    REAL(DP) :: lw_isot(S%nat3)
    REAL(DP) :: P3_isot(S%nat3,S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:), V3Bsq(:,:,:)
    INTEGER :: iq, jq, mu, nu, it, nu0(4), ix
    !
    REAL(DP) :: freq(S%nat3,5), bose(S%nat3,5), xq(3,5)
    REAL(DP),PARAMETER :: epsq = 1.e-12_dp
    !
    ! xq(:,1) -> xq1
    ! xq(:,2) -> xq2
    ! xq(:,3) -> -xq1-xq2
    ! xq(:,4) -> -xq1+xq2
    ! xq(:,5) -> -xq2
    ! And accordingly for U(:,:,i).
    ! Note that U(5) = CONJG(U(2)) and freq(5)=freq(2)


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
    xq(:,1) = xq0
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
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -xq(:,2)-xq(:,1)
      xq(:,4) =  xq(:,2)-xq(:,1)
      xq(:,5) = -xq(:,2) ! => xq4 = -xq5-xq1
      DO jq = 2,4
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
        timer_CALL t_freq%stop()
      !
      ! Interpolate D3(q1,q2,-q1-q2)
        timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
        timer_CALL t_fc3int%stop()
        timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
        timer_CALL t_fc3rot%stop()
        timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
        timer_CALL t_fc3m2%stop()
      !
      ! Interpolate D3(q1,-q2, q2-q1)
      ! For this process, we send q2 -> -q2,
      ! i.e. D2(-q2,q2) -> D2(q2,-q2) = D2(-q2,q2)*
        IF( ALL(ABS(xq(:,2))<epsq) ) THEN
        ! When q2 == 0, just copy over
        V3Bsq = V3sq
        ! xq(:,5) = 0._dp
        ! U(:,:,5) = U(:,:,2)
      ELSE
        U(:,:,5) = CONJG(U(:,:,2))
          timer_CALL t_fc3int%start()
        CALL fc3%interpolate(xq(:,5), xq(:,4), S%nat3, D3)
          timer_CALL t_fc3int%stop()
          timer_CALL t_fc3rot%start()
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,5), U(:,:,4))
          timer_CALL t_fc3rot%stop()
          timer_CALL t_fc3m2%start()
        V3Bsq = REAL( CONJG(D3)*D3 , kind=DP)
          timer_CALL t_fc3m2%stop()
      ENDIF
      !
      CONF_LOOP : &
      DO it = 1,nconf
          timer_CALL t_bose%start()
        DO jq = 1,4
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
          timer_CALL t_bose%stop()
        !
          timer_CALL t_sum%start()
        A_out(:,it) =  A_out(:,it) + grid%w(iq)*sum_A_out_modes( S%nat3, sigma(it), freq, bose, V3sq, V3Bsq, nu0 )
          timer_CALL t_sum%stop()
        !
        IF(isotopic_disorder)THEN
            timer_CALL t_lwisot%start()
          P3_isot = sum_isotope_scattering_modes(S%nat3, S%nat, sigma(it), freq, &
                                              bose, S%ntyp, S%ityp, S%amass_variance, U)

          ! I have to sum the scattering modes on the second index in order to get the
          ! A_out diagonal matrix element
          DO nu = 1, S%nat3
            A_out_isot(:,it) = A_out_isot(:,it) + grid%w(iq)* P3_isot(:,nu)
          ENDDO
            timer_CALL t_lwisot%stop()
        ENDIF
        !
      ENDDO &
      CONF_LOOP
      !
    ENDDO &
    QPOINT_INNER_LOOP
    !
    ! Recollect over MPI processes if necessary
    IF(grid%scattered) THEN
        timer_CALL t_mpicom%start()
      CALL mpi_bsum(S%nat3, nconf, A_out)
        timer_CALL t_mpicom%stop()
    ENDIF
    !
    CALL check_graceful_termination
    !
    A_out_q = A_out + A_out_isot
    !
    DEALLOCATE(U, V3sq, V3Bsq, D3)
    !
  END FUNCTION A_out_q
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the scattering probability required for the matrix A.
  ! The form is a bit more complex than in  PRB 88, 045430 (2013)
  ! because we have rearranged the indexes to only have cohalescence
  ! processes instead of mixed scattering-cohalescence. This trick ensures
  ! that everything is positive definite.
  FUNCTION sum_A_out_modes(nat3, sigma, freq, bose, V3sq, V3Bsq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(nat3,4)
    REAL(DP),INTENT(in) :: bose(nat3,4)
    REAL(DP),INTENT(in) :: V3sq(nat3,nat3,nat3)
    REAL(DP),INTENT(in) :: V3Bsq(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: nu0(4)
    !
    REAL(DP) :: sum_A_out_modes(nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_a, bose_c ! final/initial state populations
    REAL(DP) :: dom_a, dom_c   ! \delta\omega
    REAL(DP) :: ctm_a, ctm_c   !
    REAL(DP) :: norm_a, norm_c
    REAL(DP) :: sum_a, sum_c, sum_ac
    REAL(DP) :: freqm1(nat3,4)
    !
    INTEGER :: i,j,k
    REAL(DP) :: sum_A_out(nat3)
    sum_A_out(:) = 0._dp

    freqm1 = 0._dp
    DO i = 1,nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
      IF(i>=nu0(4)) freqm1(i,4) = 0.5_dp/freq(i,4)
    ENDDO
    !
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP          PRIVATE(i, j, k, bose_a, bose_c, dom_a, dom_c) &
!$OMP          PRIVATE(ctm_a, ctm_c, sum_a, sum_ac, norm_a, norm_c) &
!$OMP          REDUCTION(+: sum_A_out)
!$OMP DO  COLLAPSE(3)
    DO k = 1,nat3
      DO j = 1,nat3
        DO i = 1,nat3
          !
          bose_a = bose(i,1) * bose(j,2) * (bose(k,3)+1)
          bose_c = (bose(i,1)+1) * bose(j,2) * bose(k,4)
          !
          dom_a =  freq(i,1) + freq(j,2) - freq(k,3)
          dom_c = -freq(i,1) + freq(j,2) + freq(k,4)
          !
          ctm_a = bose_a *  f_gauss(dom_a, sigma)
          ctm_c = bose_c *  f_gauss(dom_c, sigma)
          !
          norm_a = pi*freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
          norm_c = pi*freqm1(i,1)*freqm1(j,2)*freqm1(k,4)
          !
          sum_a = norm_a * ctm_a * V3sq(i,j,k)
          sum_c = norm_c * ctm_c * V3Bsq(i,j,k)
          !
          sum_ac = sum_a + 0.5_dp*sum_c
          sum_A_out(i) = sum_A_out(i) + sum_ac
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
    sum_A_out_modes = sum_A_out
    !
  END FUNCTION sum_A_out_modes
  !
  ! Prepares A_out^(-1/2), given A_out (which is diagonal!) in input
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE compute_inv_sqrt_A_out(A, inv_sqrt_A, nconf, nat3, nq)
    IMPLICIT NONE
    !
    REAL(DP):: f(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(out) :: inv_sqrt_A(nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, nu
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nu,it)
    DO iq = 1,nq
!$OMP DO COLLAPSE(2)
    DO nu = 1,nat3
      DO it = 1,nconf
        IF(A(it,nu,iq)>0._dp)THEN
          inv_sqrt_A(it,nu,iq) = 1/DSQRT(A(it,nu,iq))
        ELSE IF (A(it,nu,iq)<0._dp) THEN
          inv_sqrt_A(it,nu,iq) = -1/DSQRT( -A(it,nu,iq) )
          IF(iq/=1.or.(iq==0.and.nu>3)) WRITE(*,*) "Warning: negative A_out", iq, nu, it
        ELSE
          ! This hould be infinity, but it should only happen at Gamma for
          ! acoustic bands, where we can ignore it because everything is zero
          inv_sqrt_A(it,nu,iq) = 0._dp
          IF(iq/=1.or.(iq==1.and.nu>3)) WRITE(*,*) "Suspicious null A_out", iq, nu, it
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL
    !
  END SUBROUTINE compute_inv_sqrt_A_out
  !
  ! Prepares A_out^(-1), given A_out (which is diagonal!) in input
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE compute_inv_A_out(A, inv_A, nconf, nat3, nq)
    IMPLICIT NONE
    !
    REAL(DP):: f(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(out) :: inv_A(nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, nu
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nu,it)
    DO iq = 1,nq
!$OMP DO COLLAPSE(2)
    DO nu = 1,nat3
      DO it = 1,nconf
        IF(A(it,nu,iq)/=0._dp)THEN
          inv_A(it,nu,iq) = 1/A(it,nu,iq)
          ! This should be infinity, but it should only happen at Gamma for
          ! acoustic bands, where we can ignore it because everything is zero
        ELSE
          inv_A(it,nu,iq) = 0._dp
          IF(iq/=1.or.(iq==1.and.nu>3)) WRITE(*,*) "Suspicious null A_out", iq, nu, it
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL
    !
  END SUBROUTINE compute_inv_A_out
  !
  ! \/o\________\\\_________________________________________/^>
  ! Apply \tilde{A} = (1+A_out^{-1/2} A_in A_out^{-1/2}) matrix to a vector f,
  ! A_out^{-1/2} must be already computed and provided in input
  SUBROUTINE tilde_A_times_f(f, Af, inv_sqrt_A_out, input, basis, out_grid, in_grid, S, fc2, fc3)
    USE constants,          ONLY : RY_TO_CMM1
    USE linewidth,          ONLY : linewidth_q
    USE q_grids,            ONLY : q_basis
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)     :: basis
    TYPE(q_grid),INTENT(in)      :: out_grid, in_grid
    REAL(DP),INTENT(in)  ::     f(3,input%nconf,S%nat3,out_grid%nq)
    REAL(DP),INTENT(in)  :: inv_sqrt_A_out(input%nconf,S%nat3,out_grid%nq)
    !
    REAL(DP),INTENT(out) ::    Af(3,input%nconf,S%nat3,out_grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    !
    REAL(DP),ALLOCATABLE :: aux(:,:,:,:)
    INTEGER  :: iq, it, ix, nu
    !
    ALLOCATE(aux(3,input%nconf,S%nat3,out_grid%nq))
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    ! Apply the first 1/sqrt(A_out)
    timer_CALL t_tktld%start()
    DO iq = 1,out_grid%nq
      ! bang!
      DO nu = 1,S%nat3
        DO it = 1,input%nconf
        DO ix = 1,3
          aux(ix,it,nu,iq) = inv_sqrt_A_out(it,nu,iq)*f(ix,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
    ENDDO
    timer_CALL t_tktld%stop()

    DO iq = 1,out_grid%nq
      ! apply A_in
      Af(:,:,:,iq) = A_in_times_f_q(aux, out_grid%xq(:,iq), input%nconf, input%T,&
            sigma_ry,S, basis, in_grid, fc2, fc3, input%isotopic_disorder)
      ! Apply 1+ and the second 1/sqrt(A)
      timer_CALL t_tktld%start()
      DO nu = 1,S%nat3
        DO it = 1,input%nconf
        DO ix = 1,3
          Af(ix,it,nu,iq) = f(ix,it,nu,iq) + inv_sqrt_A_out(it,nu,iq)*Af(ix,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
      timer_CALL t_tktld%stop()
      !
    ENDDO
    !
    DEALLOCATE(aux)
    !
  END SUBROUTINE tilde_A_times_f
  ! \/o\________\\\_________________________________________/^>
  ! Compute 3*nat lines of the matrix A (those corresponding to the
  ! q-point xq0) and apply them to the vector f. In output it returns
  ! 3*nat elements of Af (for each configuration)
  FUNCTION A_in_times_f_q(f, xq0, nconf, T, sigma, S, basis, grid, fc2, fc3, isotopic_disorder)
    USE q_grids,            ONLY : q_basis
    USE fc2_interpolate,    ONLY : forceconst2_grid, bose_phq, set_nu0, &
                                   freq_phq_safe, bose_phq, ip_cart2pat
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : sum_isotope_scattering_modes
    USE input_fc,           ONLY : ph_system_info
    USE mpi_thermal,        ONLY : mpi_bsum
    USE timers
    !USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)     :: basis
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    ! CAREFUL! grid%nq is scattered over MPI, basis%nq is the total number of q vectors
    REAL(DP),INTENT(in) :: f(3,nconf,S%nat3,basis%nq)
    REAL(DP),INTENT(in)  :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    LOGICAL,INTENT(in)  :: isotopic_disorder
    !
    ! FUNCTION RESULT:
    REAL(DP) :: A_in_times_f_q(3,nconf,S%nat3)
    REAL(DP) :: Af_q(3,nconf,S%nat3)
    REAL(DP) :: P3(S%nat3,S%nat3) !aux
    REAL(DP) :: P3_isot(S%nat3,S%nat3) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:), V3Bsq(:,:,:)
    INTEGER :: iq, jq, mu, nu, it, nu0(4), ix
    !
    REAL(DP) :: freq(S%nat3,5), bose(S%nat3,5), xq(3,5)
    REAL(DP),PARAMETER :: epsq = 1.e-12_dp
    !
    ! xq(:,1) -> xq1
    ! xq(:,2) -> xq2
    ! xq(:,3) -> -xq1-xq2
    ! xq(:,4) -> -xq1+xq2
    ! xq(:,5) -> -xq2
    ! And accordingly for U(:,:,i).
    ! Note that U(5) = CONJG(U(2)) and freq(5)=freq(2)


    ALLOCATE(U(S%nat3, S%nat3,5))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(V3Bsq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    Af_q = 0._dp
    P3_isot = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
      timer_CALL t_freq%start()
    xq(:,1) = xq0
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
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -xq(:,2)-xq(:,1)
      xq(:,4) =  xq(:,2)-xq(:,1)
      xq(:,5) = -xq(:,2) ! => xq4 = -xq5-xq1
      DO jq = 2,4
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
        timer_CALL t_freq%stop()
      !
      ! Interpolate D3(q1,q2,-q1-q2)
        timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
        timer_CALL t_fc3int%stop()
        timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
        timer_CALL t_fc3rot%stop()
        timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
        timer_CALL t_fc3m2%stop()
      !
      ! Interpolate D3(q1,-q2, q2-q1)
      ! For this process, we send q2 -> -q2,
      ! i.e. D2(-q2,q2) -> D2(q2,-q2) = D2(-q2,q2)*
        IF( ALL(ABS(xq(:,2))<epsq) ) THEN
        ! When q2 == 0, just copy over
        V3Bsq = V3sq
        ! xq(:,5) = 0._dp
        ! U(:,:,5) = U(:,:,2)
      ELSE
        U(:,:,5) = CONJG(U(:,:,2))
          timer_CALL t_fc3int%start()
        CALL fc3%interpolate(xq(:,5), xq(:,4), S%nat3, D3)
          timer_CALL t_fc3int%stop()
          timer_CALL t_fc3rot%start()
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,5), U(:,:,4))
          timer_CALL t_fc3rot%stop()
          timer_CALL t_fc3m2%start()
        V3Bsq = REAL( CONJG(D3)*D3 , kind=DP)
          timer_CALL t_fc3m2%stop()
      ENDIF
      !
      CONF_LOOP : &
      DO it = 1,nconf
          timer_CALL t_bose%start()
        DO jq = 1,4
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
          timer_CALL t_bose%stop()
        !
        ! P3 is a 3*nat x 3*nat minor of the A matrix, the implicit indexes are
        ! iq0 and iq, the matrix A has dimension (3*nat*nq x 3*nat*nq)
        ! Do not forget the minus sign!!
          timer_CALL t_sum%start()
        P3 =  - sum_A_in_modes( S%nat3, sigma(it), freq, bose, V3sq, V3Bsq, nu0 )
          timer_CALL t_sum%stop()
        !
        IF(isotopic_disorder)THEN
            timer_CALL t_lwisot%start()
          P3_isot = sum_isotope_scattering_modes(S%nat3, S%nat, sigma(it), freq, &
                                              bose, S%ntyp, S%ityp, S%amass_variance, U)
          P3 = P3 + P3_isot
            timer_CALL t_lwisot%stop()
        ENDIF
        !
          timer_CALL t_xain%start()
        ! 3*nat lines of the A matrix are applied now to f to produce 3*nat elements of Af
        DO mu = 1,S%nat3
        DO nu = 1,S%nat3
          DO ix = 1,3
            Af_q(ix,it,nu) = Af_q(ix,it,nu) + P3(nu,mu)*f(ix,it,mu,iq+grid%iq0)*grid%w(iq)
          ENDDO
        ENDDO
        ENDDO
          timer_CALL t_xain%stop()
        !
      ENDDO &
      CONF_LOOP
      !
    ENDDO &
    QPOINT_INNER_LOOP
    !
    ! Recollect over MPI processes if necessary
    IF(grid%scattered) THEN
        timer_CALL t_mpicom%start()
      CALL mpi_bsum(3,nconf,S%nat3, Af_q)
        timer_CALL t_mpicom%stop()
    ENDIF
    !
    CALL check_graceful_termination
    A_in_times_f_q = Af_q
    !
    DEALLOCATE(U, V3sq, V3Bsq, D3)
    !
  END FUNCTION A_in_times_f_q
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the scattering probability required for the matrix A.
  ! The form is a bit more complex than in  PRB 88, 045430 (2013)
  ! because we have rearranged the indexes to only have cohalescence
  ! processes instead of mixed scattering-cohalescence. This trick ensures
  ! that everything is positive definite.
  FUNCTION sum_A_in_modes(nat3, sigma, freq, bose, V3sq, V3Bsq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(nat3,4)
    REAL(DP),INTENT(in) :: bose(nat3,4)
    REAL(DP),INTENT(in) :: V3sq(nat3,nat3,nat3)
    REAL(DP),INTENT(in) :: V3Bsq(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: nu0(4)
    !
    REAL(DP) :: sum_A_in_modes(nat3,nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_a, bose_b, bose_c ! final/initial state populations
    REAL(DP) :: dom_a, dom_b, dom_c   ! \delta\omega
    REAL(DP) :: ctm_a, ctm_b, ctm_c   !
    REAL(DP) :: norm_a, norm_bc
    REAL(DP) :: sum_a, sum_bc, sum_abc
    REAL(DP) :: freqm1(nat3,4)
    !
    INTEGER :: i,j,k
    REAL(DP) :: p(nat3,nat3)
    p(:,:) = 0._dp

    freqm1 = 0._dp
     DO i = 1,nat3
       IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
       IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
       IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
       IF(i>=nu0(4)) freqm1(i,4) = 0.5_dp/freq(i,4)
     ENDDO

    !WHERE(ABS(freq)>eps_freq)
    !  freqm1 = 0.5_dp/freq
    !ELSEWHERE
    !  freqm1 = 0._dp
    !ENDWHERE
    !
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP          PRIVATE(i, j, k, bose_a, bose_b, bose_c, dom_a, dom_b, dom_c) &
!$OMP          PRIVATE(ctm_a, ctm_b, ctm_c, sum_a, sum_bc, sum_abc, norm_a, norm_bc) &
!$OMP          REDUCTION(+: p)
!$OMP DO  COLLAPSE(3)
    DO k = 1,nat3
      DO j = 1,nat3
        DO i = 1,nat3
          !
          bose_a = bose(i,1) * bose(j,2) * (bose(k,3)+1)
          bose_b = bose(i,1) * (bose(j,2)+1) * bose(k,4)
          bose_c = (bose(i,1)+1) * bose(j,2) * bose(k,4)
          !
          dom_a =  freq(i,1) + freq(j,2) - freq(k,3)
          dom_b =  freq(i,1) - freq(j,2) + freq(k,4)
          dom_c = -freq(i,1) + freq(j,2) + freq(k,4)
          !
          ctm_a = bose_a *  f_gauss(dom_a, sigma)
          ctm_b = bose_b *  f_gauss(dom_b, sigma)
          ctm_c = bose_c *  f_gauss(dom_c, sigma)
          !
          norm_a  = pi*freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
          norm_bc = pi*freqm1(i,1)*freqm1(j,2)*freqm1(k,4)
          !
          sum_a  = norm_a  * ctm_a         * V3sq(i,j,k)
          sum_bc = norm_bc * (ctm_b+ctm_c) * V3Bsq(i,j,k)
          !
          sum_abc = - sum_a + sum_bc
          p(i,j) = p(i,j) + sum_abc
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
    sum_A_in_modes = p
    !
  END FUNCTION sum_A_in_modes
  !
  !
  ! Compute thermal conductivity as -2\lambda \over { N T^2 } 1/2 ( f \dot g - f \dot b) )
  ! g = Af-b
  ! f.g = f.Af - b.f
  ! 1/2 (f.g -b.f ) = 1/2 f.Af - b.f
  FUNCTION calc_tk_gf(g, f, b, T, w, Omega, nconf, nat3, nq) RESULT(tk)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE timers
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: g(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: w(nq) ! integration weights (1/nq for regular grids)
    REAL(DP),INTENT(in) :: Omega ! cell volume
    !
    REAL(DP) :: tk(3,3,nconf)
    REAL(DP) :: pref
    INTEGER  :: iq, it, ix, jx, nu
    !
    !
    tk = 0._dp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nu,it,ix,jx)
    DO iq = 1,nq
!$OMP DO COLLAPSE(4)
      DO nu = 1,nat3
        DO it = 1,nconf
          DO jx = 1,3
          DO ix = 1,3
            tk(ix,jx,it) = tk(ix,jx,it) &
                    + w(iq)*f(ix,it,nu,iq) *(g(jx,it,nu,iq)-b(jx,it,nu,iq))
          ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO
    ENDDO
!$OMP END PARALLEL
    DO it = 1,nconf
      pref = -Omega/T(it)**2
      tk(:,:,it) = pref * tk(:,:,it)
    ENDDO
    !
  END FUNCTION calc_tk_gf
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_tk(tk, sigma, T, nconf, name, unit0, iter)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: tk(3,3,nconf)
    REAL(DP),INTENT(in) :: sigma(nconf)
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in) :: nconf
    CHARACTER(len=*),INTENT(in) :: name
    INTEGER,OPTIONAL,INTENT(in) :: unit0, iter
    !
    INTEGER :: it
    IF(present(unit0) .and. .not. present(iter)) CALL errore("print_tk", "wrong args",1)

    ! print to screen
    ioWRITE(stdout,'(2x,a)') name
    IF(present(unit0)) REWIND(unit0)
    DO it = 1,nconf
      ! on screen we only write the diagonal components
      ioWRITE(stdout,'(i6,2f10.4,3(3e17.8,3x))') it, sigma(it), T(it), &
        tk(1,1,it)*RY_TO_WATTMM1KM1, &
        tk(2,2,it)*RY_TO_WATTMM1KM1, &
        tk(3,3,it)*RY_TO_WATTMM1KM1!, &
        !tk(1,2,it)*RY_TO_WATTMM1KM1, &
        !tk(1,3,it)*RY_TO_WATTMM1KM1, &
        !tk(2,3,it)*RY_TO_WATTMM1KM1, &
        !tk(2,1,it)*RY_TO_WATTMM1KM1, &
        !tk(3,1,it)*RY_TO_WATTMM1KM1, &
        !tk(3,2,it)*RY_TO_WATTMM1KM1
      ! on file we only write everything, but in this order: 
      !   Kxx, Kyy, Kzz, Kxy, Kxz, Kyz, Kyx, Kzx, Kzy
      IF(present(unit0)) THEN
      ioWRITE(unit0,'(i6,2f10.4,3(3e17.8,3x))') it, sigma(it), T(it), &
        tk(1,1,it)*RY_TO_WATTMM1KM1, &
        tk(2,2,it)*RY_TO_WATTMM1KM1, &
        tk(3,3,it)*RY_TO_WATTMM1KM1, &
        tk(1,2,it)*RY_TO_WATTMM1KM1, &
        tk(1,3,it)*RY_TO_WATTMM1KM1, &
        tk(2,3,it)*RY_TO_WATTMM1KM1, &
        tk(2,1,it)*RY_TO_WATTMM1KM1, &
        tk(3,1,it)*RY_TO_WATTMM1KM1, &
        tk(3,2,it)*RY_TO_WATTMM1KM1
      ENDIF
    ENDDO
    IF(present(unit0)) FLUSH(unit0)

    !
    IF(present(unit0))THEN
      DO it = 1,nconf
        ioWRITE(unit0+it,'(i6,2f10.4,3(3e17.8,3x))') iter, sigma(it), T(it), &
        tk(1,1,it)*RY_TO_WATTMM1KM1, &
        tk(2,2,it)*RY_TO_WATTMM1KM1, &
        tk(3,3,it)*RY_TO_WATTMM1KM1, &
        tk(1,2,it)*RY_TO_WATTMM1KM1, &
        tk(1,3,it)*RY_TO_WATTMM1KM1, &
        tk(2,3,it)*RY_TO_WATTMM1KM1, &
        tk(2,1,it)*RY_TO_WATTMM1KM1, &
        tk(3,1,it)*RY_TO_WATTMM1KM1, &
        tk(3,2,it)*RY_TO_WATTMM1KM1
        FLUSH(unit0+it)
      ENDDO
    ENDIF
    
    !
  END SUBROUTINE print_tk
  !
  ! \/o\________\\\_________________________________________/^>
  LOGICAL FUNCTION check_gradmod2_tk(g2, name, T, Omega, nconf, nat3, nq, thr)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: g2(3,nconf)
    CHARACTER(len=*),INTENT(in) :: name
    REAL(DP),INTENT(in) :: T(nconf), thr
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: Omega ! cell volume (bohr^3)
    !
    REAL(DP) :: g(3,nconf)
    INTEGER :: it
    !
    DO it = 1,nconf
      g(:,it) = g2(:,it)*RY_TO_WATTMM1KM1!*Omega/nq!/T(it)
    ENDDO
    ioWRITE(stdout,'(2x,a)') name
    check_gradmod2_tk = .true.
    DO it = 1,nconf
      ioWRITE(*,'(12x,i6,3(3e17.8,3x))') it, g(:,it)
      check_gradmod2_tk = check_gradmod2_tk &
                    .and. ALL(g(:,it) < thr)
    ENDDO
    !
  END FUNCTION check_gradmod2_tk
  !
END MODULE variational_tk
