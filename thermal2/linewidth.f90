!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! References:
! [1] Calandra, Lazzeri, Mauri : Physica C 456 (2007) 38-44
! [2] Li. et.al. PhysRevLett.112.175501
#define timer_CALL CALL

MODULE linewidth

#include "mpi_thermal.h"
   USE kinds,            ONLY : DP
   USE mpi_thermal,      ONLY : my_id, num_procs, mpi_bsum, allgather_mat, scatteri_tns, mpi_sum_mat
   USE constants,        ONLY : RY_TO_CMM1, PI
   USE q_grids,          ONLY : q_grid, setup_grid
   USE input_fc,         ONLY : ph_system_info
   USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat_gemm, d3_mixed, ip_cart2pat
   USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
   USE functions,        ONLY : refold_bz
   USE code_input,       ONLY : code_input_type
   USE timers
CONTAINS

   SUBROUTINE check_negative_lw(lw, xq, nat3, nconf, name)
      IMPLICIT NONE
      REAL(DP),INTENT(in) :: lw(nat3, nconf)
      real(DP),INTENT(in) :: xq(3)
      INTEGER,INTENT(in)  :: nat3, nconf
      CHARACTER(*),INTENT(in) :: name
      !
      EXTERNAL errore
      !
      INTEGER :: i,j
      IF(ANY(lw<0._dp))THEN
         DO i = 1,nconf
            DO j = 1,nat3
               IF(lw(j,i) < 0._dp) THEN
                  WRITE(*,*) i, j, lw(j,i)
                  WRITE(*,'(3E13.2)') xq
               ENDIF
            ENDDO
         ENDDO
         CALL errore(name, "negative linewidth", 1)
      ENDIF
      RETURN
   END SUBROUTINE

! <<^V^\\=========================================//-//-//========//O\\//
   SUBROUTINE weights_tetra(grid, xq0, S, fc2, weights_C, weights_X)
      ! USE q_grids,          ONLY : q_grid
      ! USE constants,        ONLY : pi
      ! USE input_fc,         ONLY : ph_system_info
      ! USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      ! USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      ! USE functions,        ONLY : refold_bz
      USE thtetra,          ONLY : tetra_init, tetra_weights_delta_vec
      IMPLICIT NONE
      !
      TYPE(q_grid), INTENT(IN) :: grid
      REAL(DP), INTENT(IN) :: xq0(3)
      TYPE(ph_system_info), INTENT(IN) :: S
      TYPE(forceconst2_grid), INTENT(IN) :: fc2
      REAL(DP), INTENT(OUT) :: weights_C(S%nat3, S%nat3**2, grid%nqtot), weights_X(S%nat3, S%nat3**2, grid%nqtot)
      INTEGER :: iq, ibnd, jbnd, index_double
      REAL(DP),ALLOCATABLE :: freqs_doubled_X(:,:), freqs_doubled_C(:,:)
      REAL(DP), ALLOCATABLE :: gather_doubled_X(:,:), gather_doubled_C(:,:)
      REAL(DP) :: xq(3,3), freq(S%nat3,3)
      !

      ! ALLOCATE(weights_X(nat3, nat3**2, grid%nqtot), weights_C(nat3, nat3**2, grid%nqtot))
      ALLOCATE(freqs_doubled_X(S%nat3**2, grid%nq), freqs_doubled_C(S%nat3**2, grid%nq))
      !
      xq(:,1) = xq0
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1))
      ! this cycle initializes freqs, as a grid of frequencies
      ! this cycle initializes freqs_doubled, which is freq(ibnd,2) +/- freq(jbnd,3)
      timer_CALL t_freqd%start()
      DO iq = 1, grid%nq
         CALL freq_phq_safe(grid%xq(:,iq), S, fc2, freq(:,2))
         ! the third vector (why it's minus?)
         xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
         CALL freq_phq_safe(xq(:,3), S, fc2, freq(:,3))
         DO ibnd = 1, S%nat3
            DO jbnd = 1, S%nat3
               index_double = S%nat3*(ibnd-1) + jbnd
               freqs_doubled_X(index_double,iq) = freq(ibnd,2) + freq(jbnd,3)
               freqs_doubled_C(index_double,iq) = freq(jbnd,3) - freq(ibnd,2)
            ENDDO
         ENDDO
      ENDDO
      IF(grid%scattered) THEN
         ALLOCATE(gather_doubled_X(S%nat3**2, grid%nqtot))
         ALLOCATE(gather_doubled_C(S%nat3**2, grid%nqtot))
         CALL allgather_mat(S%nat3**2, grid%nq, freqs_doubled_X, gather_doubled_X)
         CALL allgather_mat(S%nat3**2, grid%nq, freqs_doubled_C, gather_doubled_C)
      else
         ALLOCATE(gather_doubled_X(S%nat3**2, grid%nq))
         ALLOCATE(gather_doubled_C(S%nat3**2, grid%nq))
         gather_doubled_X = freqs_doubled_X
         gather_doubled_C = freqs_doubled_C
      ENDIF

      timer_CALL t_freqd%stop()
      timer_CALL t_thtetra%start()
      CALL tetra_init(grid%n, S%bg, gather_doubled_X)
      weights_X = tetra_weights_delta_vec(S%nat3, freq(:,1))
      CALL tetra_init(grid%n, S%bg, gather_doubled_C)
      weights_C = tetra_weights_delta_vec(S%nat3, freq(:,1))

      timer_CALL t_thtetra%stop()
      DEALLOCATE(gather_doubled_X, gather_doubled_C)
      DEALLOCATE(freqs_doubled_X, freqs_doubled_C)
   END SUBROUTINE

   PURE function rotate_single_d3(i, j, k, U, D3, invert) result(D3_s2)
      INTEGER, INTENT(IN) :: i, j, k
      LOGICAL, INTENT(IN) :: invert
      COMPLEX(DP), INTENT(IN) :: U(:,:,:), D3(:,:,:)
      !
      INTEGER :: a, b, c, nat3
      COMPLEX(DP) :: D3_s, aux
      REAL(DP) :: D3_s2
      D3_s = 0._dp
      nat3 = SIZE(U, 1)
      !
      DO c = 1, nat3
         DO b = 1, nat3
            if (invert) then
               aux = U(b,i,1)* U(c,k,3)
            else
               aux = U(b,j,2)* U(c,k,3)
            endif
            DO a = 1, nat3
               if (invert) then
                  D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,j,2) )
               else
                  D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,i,1) )
               endif
            END DO
         END DO
      END DO
      D3_s2 = REAL( CONJG(D3_s)* D3_s, kind=DP)
   end function

   FUNCTION linewidth_q(xq1, S, fc2, fc3, grid, input, lw_UN)
      USE functions, ONLY : f_gauss
      USE fc3_interpolate, ONLY : sum_R3, d3_mixed
      USE merge_degenerate,   ONLY : merge_degen

      IMPLICIT NONE
      TYPE(q_grid), INTENT(IN) :: grid
      REAL(DP), INTENT(IN) :: xq1(3)
      TYPE(ph_system_info), INTENT(IN) :: S
      TYPE(forceconst2_grid), INTENT(IN) :: fc2
      CLASS(forceconst3), INTENT(IN) :: fc3
      TYPE(code_input_type), INTENT(IN) :: input
      ! LOGICAL, INTENT(IN), OPTIONAL :: coalescence
      REAL(DP),OPTIONAL,INTENT(out)   :: lw_UN(S%nat3,2,input%nconf) ! Normal/Umklapp contribution
      !
      INTEGER,PARAMETER :: normal=1, umklapp=2
      INTEGER :: iq, jq, j_un, it, nu0(3), i,j,k
      REAL(DP) :: f(3), xq(3,3), freq(S%nat3,3), bose(S%nat3,3)
      REAL(DP) :: freqm1(S%nat3,3), freqtotm1_23, freqtotm1, bose_C, bose_X, dom_C, dom_X, sigma
      COMPLEX(DP) :: aux(S%nat3), linewidth_q(S%nat3,input%nconf)
      complex(DP) :: U(S%nat3,S%nat3,3), D3(S%nat3,S%nat3,S%nat3)
      TYPE(d3_mixed) :: Dqr
      REAL(DP), allocatable :: weights_C(:,:,:), weights_X(:,:,:)
      ! ALLOCATE(U(S%nat3,S%nat3,3), D3(S%nat3,S%nat3,S%nat3))
      REAL(DP) :: V3sq(S%nat3,S%nat3,S%nat3), ctm
      linewidth_q = 0._dp
      IF(present(lw_UN)) lw_UN = 0._dp
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
      timer_CALL t_freq%start()
      xq(:,1) = xq1
      nu0(1) = set_nu0(xq(:,1), S%at)
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      timer_CALL t_freq%stop()

      !> for a two step interpolation to be effective, we need to calculate
      !> fc3%interpolate on xq1 and xq2/xq3. I chose fc3%interpolate(xq1,xq3),
      !> so xq2 and xq1 are inverted then (in ip_cart2pat and in V3sq(j,i,k))
      CALL fc3%sum_R2(xq(:,1), S%nat3, Dqr)

      if (input%delta_approx == 'tetra') then
         allocate(weights_C(S%nat3, S%nat3**2, grid%nqtot))
         allocate(weights_X(S%nat3, S%nat3**2, grid%nqtot))
         CALL weights_tetra(grid, xq1, S, fc2, weights_C, weights_X)
      endif

      DO iq = 1, grid%nq
         !
         timer_CALL t_freq%start()
         ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
         !  IF(PRESENT(coalescence)) THEN
         !  xq(:,2) = grid%xq(:,iq)
         !     xq(:,3) = -(xq(:,2)+xq(:,1))
         !  ELSE
         ! xq(:,2) = [0.5, 0.4, 0.3]
         xq(:,2) = grid%xq(:,iq)
         xq(:,3) = - (xq(:,1) + xq(:,2))
         !  ENDIF
         IF(present(lw_UN)) THEN
            IF(ALL(ABS(refold_bz(xq(:,1), S%bg) &
               +refold_bz(xq(:,2), S%bg) &
               +refold_bz(xq(:,3), S%bg))<1.d-6)) THEN
               j_un = normal
               j_un = umklapp
            ENDIF
         ENDIF

!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
         DO jq = 2,3
            nu0(jq) = set_nu0(xq(:,jq), S%at)
            CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
         ENDDO
!$OMP END PARALLEL DO
         timer_CALL t_freq%stop()
         !
         timer_CALL t_fc3int%start()
         ! old implementation of interpolation:
         !CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
         CALL sum_R3(S, xq(:,3), Dqr, D3)
         timer_CALL t_fc3int%stop()
         timer_CALL t_fc3rot%start()
         !> gemm implementation of rotation
         CALL ip_cart2pat(D3, S%nat3, U(:,:,2), U(:,:,1), U(:,:,3))
         timer_CALL t_fc3rot%stop()
         timer_CALL t_fc3m2%start()
         V3sq = REAL( CONJG(D3)*D3 , kind=DP)
         timer_CALL t_fc3m2%stop()
         DO it = 1,input%nconf
            timer_CALL t_bose%start()
            ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
            DO jq = 1,3
               CALL bose_phq(input%T(it),S%nat3, freq(:,jq), bose(:,jq))
            ENDDO
!$OMP END PARALLEL DO
            timer_CALL t_bose%stop()
            timer_CALL t_sum%start()
            freqm1 = 0._dp
            aux = 0._dp
            DO i = 1,S%nat3
               do j = 1,3
                  IF(i>=nu0(j)) freqm1(i,j) = 0.5_dp/freq(i,j)
               ENDDO
            ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: aux) COLLAPSE(2)
            DO j = 1,S%nat3
               DO k = 1,S%nat3
                  !
                  bose_C = 2* (bose(j,2) - bose(k,3))
                  bose_X = bose(j,2) + bose(k,3) + 1
                  freqtotm1_23= freqm1(j,2) * freqm1(k,3)
                  !
                  DO i = 1,S%nat3
                     !  IF(input%calculation == "cgp" .or. input%calculation == "exact") THEN
                     !     bose_C = bose(i,1) * bose(j,2) * (bose(k,3)+1)
                     !     bose_X = (bose(i,1)+1) * bose(j,2) * bose(k,3) / 2
                     !  ENDIF
                     !
                     !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
                     !
                     freqtotm1 = freqm1(i,1) * freqtotm1_23
                     !IF(freqtot/=0._dp)THEN
                     !
                     sigma = input%sigma(it)/RY_TO_CMM1
                     SELECT CASE(input%delta_approx)
                      CASE ("gauss")
                        dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
                        dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
                        ctm = bose_C * f_gauss(dom_C, sigma) + bose_X * f_gauss(dom_X, sigma)
                      CASE ("tetra")
                        ctm = bose_C * weights_C(i, S%nat3*(j-1) + k, iq+grid%iq0) + &
                              bose_X * weights_X(i, S%nat3*(j-1) + k, iq+grid%iq0)
                        ! CASE ("selfnrg")
                        !   f(1) = freq(i,1)
                        !   f(2) = freq(j,2)
                        !   f(3) = freq(k,3)
                        !   ctm = ctm_selfnrg(sigma, input%T(it), f, bose_C, bose_X)
                        ! CASE ("selfnrg_sp")
                        !   f(1) = freq(i,1)
                        !   f(2) = freq(j,2)
                        !   f(3) = freq(k,3)
                        !   ctm = ctm_selfnrg_spectre(sigma, f, energy, bose_C, bose_X)
                        !  CASE("cgp_C")
                        !   !> the prefactor 2 is to annul the 1/2 in the definition of the linewidth, but maybe
                        !   !> there's another prefactor
                        !   ctm = 2 * bose(i,1) * bose(j,2) * (1.0_dp + bose(k,3)) * weights_C(i,S%nat3*(j-1) + k, iq)
                      CASE DEFAULT
                        CALL errore("sum_rotate_lw", "you should give sigma/tetra_weights", 1)
                     END SELECT
                     !
                     ! IF(TRIM(input%calculation) == "cgp" .or. TRIM(input%calculation) == "exact") THEN
                     !   ctm = ctm * 2 * bose(i,1) * (bose(i,1) + 1)
                     ! ENDIF
                     ! IF(REAL(ctm, DP) < 1) CYCLE
                     ! ci sono tanti negativi che contribuiscono sulla terza cifra, mica da poco
                     !
                     !> rotation of single D3 only where we know that we have scattering
                     timer_CALL t_fc3rot%start()
                     ! V3sq = rotate_single_d3(i,j,k, U, D3, invert=.true.)

                     aux(i) = aux(i) + ctm * V3sq(j,i,k) * freqtotm1
                     ! if (i == 3 .and. j == 4 .and. k == 2 .and. iq == 6) then
                     !    print*, ctm, V3sq(j,i,k), freqtotm1
                     ! endif
                     timer_CALL t_fc3rot%stop()
                  ENDDO
               ENDDO
            ENDDO
!$OMP END PARALLEL DO
            !
            CALL merge_degen(S%nat3, aux, freq(:,1))

            linewidth_q(:,it) = linewidth_q(:,it) + aux * grid%w(iq)
            IF(present(lw_UN)) lw_UN(:,j_un,it) = lw_UN(:,j_un,it) + aux
            timer_CALL t_sum%stop()
         ENDDO
         !
      ENDDO

      timer_CALL t_mpicom%start()
      IF(grid%scattered) CALL mpi_bsum(S%nat3,input%nconf,linewidth_q)
      IF(grid%scattered .and. present(lw_UN)) CALL mpi_bsum(S%nat3,2,input%nconf,lw_UN)
      timer_CALL t_mpicom%stop()

      linewidth_q = linewidth_q * 0.5_dp * pi 
      if (input%delta_approx == 'tetra') linewidth_q = linewidth_q * grid%nqtot
      if (allocated(Dqr%DR3)) deallocate(Dqr%DR3, Dqr%R3)
      if(allocated(weights_C)) deallocate(weights_C, weights_X)


   END FUNCTION

   !> helper function for selfnrg. Not tesed, not used. It would be nice to use 
   !> only linewidth_q and various ctm functions for each case
   FUNCTION ctm_selfnrg_spectre(sigma, freq, ener, bose_C, bose_X)
      USE input_fc,           ONLY : ph_system_info
      USE merge_degenerate,   ONLY : merge_degen
      USE functions,          ONLY : sigma_mgo
      IMPLICIT NONE
      REAL(DP),INTENT(in) :: sigma   ! smearing (regularization) (Ry)
      REAL(DP),INTENT(in) :: freq(3)  ! phonon energies (Ry)
      !
      REAL(DP),INTENT(in) :: ener     ! energies for which to compute the spectral function
      REAL(DP), INTENT(IN) :: bose_C, bose_X
      !
      ! _P -> scattering, _M -> cohalescence
      REAL(DP) :: omega_P,  omega_M   ! \delta\omega
      REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
      COMPLEX(DP) :: ctm_P,ctm_M, reg
      !
      ! Note: using the function result in an OMP reduction causes crash with ifort 14
      COMPLEX(DP) :: ctm_selfnrg_spectre
      !
!     IF(sigma<=0._dp)THEN
!       CALL errore("sum_selfnrg_spectre","spf not implemented in the static limit. "&
!                   //"NEW: To do unshifted spf use 'spf imag'",1)
!     ENDIF
      !
      !
      omega_P  = freq(2)+freq(3)
      omega_P2 = omega_P**2
      !
      omega_M  = freq(2)-freq(3)
      omega_M2 = omega_M**2
      !
      !
      reg = CMPLX(ener, sigma, kind=DP)**2
      !
      ctm_P = 2 * bose_C *omega_P/(omega_P2-reg)
      ctm_M = 2 * bose_X *omega_M/(omega_M2-reg)
      !
      ctm_selfnrg_spectre = ctm_P + ctm_M
      !
   END FUNCTION

!
! \/o\________\\\_________________________________________/^>
! Sum the self energy in a range of frequencies
   FUNCTION ctm_selfnrg(sigma, T, freq, bose_C, bose_X)
      USE input_fc,           ONLY : ph_system_info
      USE functions,          ONLY : df_bose
      USE merge_degenerate,   ONLY : merge_degen
      USE functions,          ONLY : sigma_mgo
      IMPLICIT NONE
      REAL(DP),INTENT(in) :: sigma, T, freq(3), bose_C, bose_X
      REAL(DP) :: omega_X,  omega_C   ! \sigma\omega
      COMPLEX(DP) :: ctm_X, ctm_C, reg, ctm_selfnrg
      !
      ctm_X = 0._dp
      ctm_C = 0._dp
      omega_X  = freq(2)+freq(3)
      omega_C  = freq(2)-freq(3)
      IF(sigma>0._dp)THEN
         reg = CMPLX(freq(1), sigma, kind=DP)**2
         ctm_X = 2 * bose_X *omega_X/(omega_X**2-reg )
         ctm_C = 2 * bose_C *omega_C/(omega_C**2-reg )
         ctm_X = 2 * bose_X *omega_X/(omega_X**2+sigma**2)
         ctm_C = 2 * bose_C *omega_C/(omega_C**2+sigma**2)
         ! In the static limit with sigma=0 case we have to take the
         ! derivative of (n_3-n2)/(w_2-w_3) when w_2 is close to w_3
         IF(omega_X>0._dp)THEN
            ctm_X = 2 * bose_X /omega_X
            ctm_X = 0._dp
         ENDIF
         !
         IF(ABS(omega_C)>1.e-5_dp)THEN
            ctm_C = 2 * bose_C /omega_C
            IF(T>0._dp)THEN
               ctm_C = -2* df_bose(0.5_dp * omega_X, T)
               ctm_C = 0._dp
            ENDIF
         ENDIF
         !
      ENDIF

      ctm_selfnrg = ctm_C + ctm_X
   END FUNCTION
!
!    FUNCTION sum_rotate_lw(S, freq, bose, D3, U, nu0, calc, sigma, T, weights_C, weights_X)
!       USE functions, ONLY : f_gauss => f_gauss
!       USE constants, ONLY : pi, RY_TO_CMM1
!       USE input_fc,           ONLY : ph_system_info
!       USE merge_degenerate,   ONLY : merge_degen
!       IMPLICIT NONE
!       TYPE(ph_system_info),INTENT(in)   :: S
!       REAL(DP),INTENT(in) :: freq(S%nat3,3), bose(S%nat3,3)
!       COMPLEX(DP),INTENT(in) :: D3(S%nat3,S%nat3,S%nat3), U(S%nat3,S%nat3,3)
!       INTEGER,INTENT(in)  :: nu0(3)
!       CHARACTER(10), INTENT(IN) :: calc
!       REAL(DP),INTENT(in), OPTIONAL :: sigma, T
!       REAL(DP), INTENT(IN), OPTIONAL :: weights_C(S%nat3,S%nat3**2), weights_X(S%nat3,S%nat3**2)
!       !
!       COMPLEX(DP) :: sum_rotate_lw(S%nat3)
!       COMPLEX(DP) :: D3_s
!       REAL(DP) :: D3_s2
!       INTEGER :: a,b,c
!       !
!       ! _C -> scattering, _X -> cohalescence
!       REAL(DP) :: bose_C, bose_X ! final/initial state populations
!       REAL(DP) :: dom_C, dom_X   ! \delta\omega
!       COMPLEX(DP) :: ctm   !
!       REAL(DP) :: freqtotm1, freqtotm1_23
!       REAL(DP) :: freqm1(S%nat3,3)
!       COMPLEX(DP) ::  aux
!       !REAL(DP),SAVE :: leftover_e
!       !
!       INTEGER :: i,j,k
!       !
!       freqm1 = 0._dp
!       sum_rotate_lw = 0._dp
!       DO i = 1,S%nat3
!          IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
!          IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
!          IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
!       ENDDO
! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
! !$OMP                     freqtotm1_23,freqtotm1) &
! !$OMP             REDUCTION(+: sum_rotate_lw) COLLAPSE(2)
!       DO k = 1,S%nat3
!          DO j = 1,S%nat3
!             !
!             bose_C = 2* (bose(j,2) - bose(k,3))
!             bose_X = bose(j,2) + bose(k,3) + 1
!             freqtotm1_23= freqm1(j,2) * freqm1(k,3)
!             !
!             DO i = 1,S%nat3
!                !
!                !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
!                !
!                freqtotm1 = freqm1(i,1) * freqtotm1_23
!                !IF(freqtot/=0._dp)THEN
!                !
!                SELECT CASE (calc)
!                 CASE ("gauss")
!                   dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
!                   dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
!                   ctm = (bose_C * f_gauss(dom_C, sigma) + bose_X * f_gauss(dom_X, sigma))
!                 CASE ("tetra")
!                   ctm = bose_C * weights_C(i,S%nat3*(j-1) + k) + bose_X * weights_X(i,S%nat3*(j-1) + k)
!                 CASE ("selfnrg")
!                   ctm = ctm_selfnrg(sigma, T, freq, bose_C, bose_X)
!                 CASE DEFAULT
!                   CALL errore("sum_rotate_lw", "you should give sigma/tetra_weights", 1)
!                END SELECT
!                ! IF(REAL(ctm, DP) < 1) CYCLE
!                ! ci sono tanti negativi che contribuiscono sulla terza cifra, mica da poco
!                !
!                timer_CALL t_fc3rot%start()
!                !> rotation of single D3 only where we know that we have scattering
!                D3_s = 0._dp
!                DO c = 1, S%nat3
!                   DO b = 1, S%nat3
!                      aux = U(b,i,1)* U(c,k,3)
!                      DO a = 1, S%nat3
!                         D3_s = D3_s + D3(a,b,c) * CONJG( aux * U(a,j,2) )
!                      END DO
!                   END DO
!                END DO
!                timer_CALL t_fc3rot%stop()
!                D3_s2 = REAL( CONJG(D3_s)* D3_s, kind=DP)
!                sum_rotate_lw(i) = sum_rotate_lw(i) + ctm * D3_s2 * freqtotm1
!             ENDDO
!          ENDDO
!       ENDDO
! !$OMP END PARALLEL DO
!       !
!       CALL merge_degen(S%nat3, sum_rotate_lw, freq(:,1))
!       !
!    END FUNCTION sum_rotate_lw
! <<^V^\\=========================================//-//-//========//O\\//

   ! <<^V^\\=========================================//-//-//========//O\\//
   ! This function returns the complex self energy from the bubble(?) diagram:
   !  its real part is the order contribution to lineshift
   !  its complex part is the intrinsic linewidth
   FUNCTION selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, freq1, U1)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      !
      REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! FUNCTION RESULT:
      COMPLEX(DP) :: selfnrg_q(S%nat3,nconf)
      COMPLEX(DP) :: se(S%nat3,nconf) ! aux
      !
      COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
      REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
      INTEGER :: iq, jq, it
      INTEGER :: nu0(3)
      !
      REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
      !
      ALLOCATE(U(S%nat3, S%nat3,3))
      ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
      ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
      !
      se = (0._dp, 0._dp)
      !
      timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes at q1
      xq(:,1) = xq0
      nu0(1) = set_nu0(xq(:,1), S%at)
      IF(present(freq1) .and. present(U1)) THEN
         freq(:,1) = freq1
         U(:,:,1)    = U1
      ELSE
         CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      ENDIF
      timer_CALL t_freq%stop()

      DO iq = 1, grid%nq
         !
         timer_CALL t_freq%start()
         xq(:,2) = grid%xq(:,iq)
         xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
         DO jq = 2,3
            nu0(jq) = set_nu0(xq(:,jq), S%at)
            CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
         ENDDO
!$OMP END PARALLEL DO
         timer_CALL t_freq%stop()
         !
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
         DO it = 1,nconf
            ! Compute bose-einstein occupation at q2 and q3
            timer_CALL t_bose%start()
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
            DO jq = 1,3
               CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
            ENDDO
!$OMP END PARALLEL DO
            timer_CALL t_bose%stop()
            !
            timer_CALL t_sum%start()
            se(:,it) = se(:,it) + grid%w(iq)*sum_selfnrg_modes( S, sigma(it), T(it), freq, bose, V3sq, nu0 )
            timer_CALL t_sum%stop()
            !
         ENDDO
         !
      ENDDO
      !
      timer_CALL t_mpicom%start()
      IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,se)
      timer_CALL t_mpicom%stop()
      selfnrg_q = -0.5_dp * se
      !
      DEALLOCATE(U, V3sq)
      !
   END FUNCTION selfnrg_q
   ! \/o\________\\\_________________________________________/^>
   ! Simple spectral function, computed as a superposition of Lorentzian functions
   FUNCTION simple_spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1, shift, input) &
      RESULT(spectralf)
      USE q_grids,            ONLY : q_grid
      USE constants,          ONLY : RY_TO_CMM1
      USE input_fc,           ONLY : ph_system_info
      USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe
      USE fc3_interpolate,    ONLY : forceconst3
      ! USE globals,            ONLY : g_input

      IMPLICIT NONE
      ! ARGUMENTS:
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      type(code_input_type),INTENT(in) :: input
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      LOGICAL,INTENT(in)  :: shift
      !
      REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! FUNCTION RESULT:
      REAL(DP) :: spectralf(ne,S%nat3,nconf)
      !
      COMPLEX(DP) :: selfnrg(S%nat3,nconf)
      !
      INTEGER :: it, i, ie
      REAL(DP) :: gamma(S%nat3,nconf), delta(S%nat3,nconf), omega, denom
      COMPLEX(DP) :: U(S%nat3, S%nat3)
      REAL(DP) :: freq(S%nat3)
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
      IF(present(freq1) .and. present(U1)) THEN
         freq = freq1
         U    = U1
      ELSE
         !CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
         CALL freq_phq_safe(xq0, S, fc2, freq, U)
      ENDIF
      !
      ! use lineshift to compute 3rd order linewidth and lineshift
      IF(shift)THEN
         selfnrg = selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
         gamma(:,:) = -DIMAG(selfnrg(:,:))
         delta(:,:) =   DBLE(selfnrg(:,:))
      ELSE
         gamma = linewidth_q(xq0, S, fc2, fc3, grid, input)
         delta = 0._dp
      ENDIF
      !
      DO it = 1,nconf
         WRITE(*,'(30x,6e12.2,5x,6e12.2)') gamma(:,it)*RY_TO_CMM1,delta(:,it)*RY_TO_CMM1
      ENDDO
      !
      delta = 0._dp
      ! Compute and superpose the Lorentzian functions corresponding to each band
      DO it = 1,nconf
         DO i = 1,S%nat3
            DO ie = 1, ne
               ! Only apply shift (i.e. real part of self nrg) if requested
               omega = freq(i)
               denom =   (ener(ie)**2 -omega**2 -2*omega*delta(i,it))**2 &
                  + 4*omega**2 *gamma(i,it)**2
               IF(ABS(denom)/=0._dp)THEN
                  spectralf(ie,i,it) = 2*omega*gamma(i,it) / denom
               ELSE
                  spectralf(ie,i,it) = 0._dp
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !
   END FUNCTION simple_spectre_q


   ! <<^V^\\=========================================//-//-//========//O\\//
   ! Self energy for all the phonon bands at q in a range of frequencies
   FUNCTION selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
      RESULT(selfnrg_wq)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)           :: grid
      REAL(DP),INTENT(in)               :: sigma(nconf) ! ry
      !
      REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! To interpolate D2 and D3:
      INTEGER :: iq, jq, it
      COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
      REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
      REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
      !
      ! To compute the spectral function from the self energy:
      INTEGER  :: nu0(3)
      COMPLEX(DP),ALLOCATABLE :: selfnrg(:,:,:)
      ! FUNCTION RESULT:
      COMPLEX(DP) :: selfnrg_wq(ne,S%nat3,nconf)
      !
      ALLOCATE(U(S%nat3, S%nat3,3))
      ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
      ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
      ALLOCATE(selfnrg(ne,S%nat3,nconf))
      !
      selfnrg = (0._dp, 0._dp)
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
      timer_CALL t_freq%start()
      xq(:,1) = xq0
      nu0(1) = set_nu0(xq(:,1), S%at)
      IF(present(freq1) .and. present(U1)) THEN
         freq(:,1) = freq1
         U(:,:,1)    = U1
      ELSE
         CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      ENDIF
      timer_CALL t_freq%stop()
      !
      DO iq = 1, grid%nq
         !CALL print_percent_wall(33.333_dp, 300._dp, iq, grid%nq, (iq==1))
         !
         timer_CALL t_freq%start()
         xq(:,2) = grid%xq(:,iq)
         xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
         DO jq = 2,3
            nu0(jq) = set_nu0(xq(:,jq), S%at)
            CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
         ENDDO
!$OMP END PARALLEL DO
         timer_CALL t_freq%stop()
         !
         !
         ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
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
         DO it = 1,nconf
            ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
            timer_CALL t_bose%start()
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
            DO jq = 1,3
               CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
            ENDDO
!$OMP END PARALLEL DO
            timer_CALL t_bose%stop()
            timer_CALL t_sum%start()
            selfnrg(:,:,it) = selfnrg(:,:,it) + grid%w(iq)*sum_selfnrg_spectre( S, sigma(it), freq, bose, V3sq, ne, ener, nu0 )
            timer_CALL t_sum%stop()
            !
         ENDDO
         !
      ENDDO
      !
      timer_CALL t_mpicom%start()
      IF(grid%scattered) CALL mpi_bsum(ne,S%nat3,nconf,selfnrg)
      timer_CALL t_mpicom%stop()
      selfnrg_wq = -0.5_dp * selfnrg
      !
      DEALLOCATE(U, V3sq, D3, selfnrg)
      !
   END FUNCTION selfnrg_omega_q
   !
   ! <<^V^\\=========================================//-//-//========//O\\//
   ! Spectral weight function, computed as in eq. 1 of arXiv:1312.7467v1
   FUNCTION spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1, shift) &
      RESULT(spectralf)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      LOGICAL,INTENT(in)  :: shift ! set to false to drop the real part of the self-energy
      !
      REAL(DP),INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! To compute the spectral function from the self energy:
      INTEGER  :: i, ie, it
      REAL(DP) :: gamma, delta, omega, denom
      COMPLEX(DP) :: selfnrg(ne,S%nat3,nconf)
      ! FUNCTION RESULT:
      REAL(DP)    :: spectralf(ne,S%nat3,nconf)
      !
      ! Once we have the self-energy, the rest is trivial
      selfnrg = selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
      !
      timer_CALL t_mkspf%start()
      DO it = 1,nconf
         DO i = 1,S%nat3
            DO ie = 1, ne
               gamma =  -DIMAG(selfnrg(ie,i,it))
               IF(shift) THEN
                  delta =   DBLE(selfnrg(ie,i,it))
               ELSE
                  delta = 0._dp
               ENDIF
               omega = freq1(i)
               denom =   (ener(ie)**2 -omega**2 -2*omega*delta)**2 &
                  + 4*omega**2 *gamma**2
               IF(ABS(denom)/=0._dp)THEN
                  spectralf(ie,i,it) = 2*omega*gamma / denom
               ELSE
                  spectralf(ie,i,it) = 0._dp
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      timer_CALL t_mkspf%stop()
      !
   END FUNCTION spectre_q

   ! <<^V^\\=========================================//-//-//========//O\\//
   ! Spectral weight function, computed as in eq. 1 of arXiv:1312.7467v1
   ! test, compute this as imaginary part of epsilon(omega)
   FUNCTION spectre2_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
      RESULT(spectralf)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      !LOGICAL,INTENT(in)  :: shift ! set to false to drop the real part of the self-energy
      !
      REAL(DP),INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! To compute the spectral function from the self energy:
      INTEGER  :: i, ie, it
      COMPLEX(DP) :: tepsilon(ne,S%nat3,nconf)
      ! FUNCTION RESULT:
      REAL(DP)    :: spectralf(ne,S%nat3,nconf)
      !
      ! Once we have the self-energy, the rest is trivial
      tepsilon = tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
      !
      timer_CALL t_mkspf%start()
      DO it = 1,nconf
         DO i = 1,S%nat3
            DO ie = 1,ne
               !IF(freq1(i)/=0._dp)THEN
               spectralf(ie,i,it) = DIMAG(tepsilon(ie,i,it))
               !ELSE
               !  spectralf(ie,i,it) = 0._dp
               !ENDIF
            ENDDO
         ENDDO
      ENDDO
      timer_CALL t_mkspf%stop()
      !
   END FUNCTION spectre2_q
   !
   ! <<^V^\\=========================================//-//-//========//O\\//
   !  Complex \tilde{epsilon} in a range of frequencies
   ! Experimental subroutine, assumes the material cubic and takes
   ! epsilon(1,1) as epsilon_infty
   FUNCTION tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
      RESULT(tepsilon)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      USE constants,        ONLY : fpi
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      !
      REAL(DP),INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! To compute the spectral function from the self energy:
      INTEGER  :: i, ie, it
      REAL(DP) :: pref, epsilon_infty, rmass, zeu_i2
      COMPLEX(DP) :: denom
      COMPLEX(DP):: selfnrg(ne,S%nat3,nconf)
      ! FUNCTION RESULT:
      COMPLEX(DP)    :: tepsilon(ne,S%nat3,nconf)
      !
      EXTERNAL errore
      IF(.not.S%lrigid) CALL errore("tepsilon","Cannot compute \tilde{epsilon} without epsilon0", 1)
      ioWRITE(*,*) "BEWARE: tepsilon is only isotropic"
      !
      ! Once we have the self-energy, the rest is trivial
      selfnrg = selfnrg_omega_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
      !
      ! S = 4piZ^2/(volume mass omega0^2)

      rmass = S%amass(1)*S%amass(2)/(S%amass(1)+S%amass(2))
      epsilon_infty = S%epsil(1,1)
      zeu_i2        = S%zeu(1,1,1)**2
      pref =  fpi * zeu_i2 /(S%omega * rmass)
      !pref = 1/freq1(i)**2 ! 4 pi Z^2  / (Volume mu  omega0^2) # mu = reduced mass
!     print*, "rmass", rmass, S%amass(1:2)
!     print*, "epsilon", epsilon_infty
!     print*, "zeu2", zeu_i2
!     print*, "pref", pref
!     print*, "vol", S%omega
!     print*, "denom", rmass*S%omega
!     print*, "fpi * zeu_i2", fpi * zeu_i2
!  denom   1002568.9328649712
!  epsilon   3.1410114605330000
!  fpi * zeu_i2   48.024054760187823
!  pref   4.7901000306236137E-005
!  rmass   8793.6751743770546        22152.652305024902        14582.196429874200
!  vol   114.01023041950073
!  zeu2   3.8216328511998792
!
      timer_CALL t_mkspf%start()
      tepsilon = DCMPLX(0._dp, 0._dp)
      DO it = 1,nconf
         DO i = 1,S%nat3
            IF(freq1(i)==0._dp) CYCLE
            DO ie = 1,ne
               denom = freq1(i)**2 - ener(ie)**2 - 2*freq1(i)*selfnrg(ie,i,it)
               IF(denom==DCMPLX(0._dp,0._dp)) CYCLE
               !tepsilon(ie,i,it) = 1 + pref*freq1(i)**2/denom
               tepsilon(ie,i,it) = epsilon_infty + pref/denom
            ENDDO
         ENDDO
      ENDDO
      timer_CALL t_mkspf%stop()
      !
   END FUNCTION tepsilon_q
   !
   ! <<^V^\\=========================================//-//-//========//O\\//
   ! Infra red reflectivity |sqrt(epsilon)+1/sqrt(epsilon)-1|^2
   ! Experimental! Assumes isoropic material and a bunch of other stuff
   FUNCTION ir_reflectivity_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1) &
      RESULT(reflectivity)
      USE q_grids,          ONLY : q_grid
      USE input_fc,         ONLY : ph_system_info
      USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0
      USE fc3_interpolate,  ONLY : forceconst3, ip_cart2pat
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(in) :: xq0(3)
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
      !
      INTEGER,INTENT(in)  :: ne
      REAL(DP),INTENT(in) :: ener(ne)
      !
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      CLASS(forceconst3),INTENT(in)     :: fc3
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(q_grid),INTENT(in)      :: grid
      REAL(DP),INTENT(in) :: sigma(nconf) ! ry
      !
      REAL(DP),INTENT(in) :: freq1(S%nat3)
      COMPLEX(DP),INTENT(in) :: U1(S%nat3,S%nat3)
      !
      ! To compute the spectral function from the self energy:
      INTEGER  :: i, it
      COMPLEX(DP) :: tepsilon(ne,S%nat3,nconf)
      COMPLEX(DP) :: aux(ne),aux2(ne)
      ! FUNCTION RESULT:
      REAL(DP)    :: reflectivity(ne,S%nat3,nconf)
      !
      ! We use tilde{epsilon}, which itself comes from the self-energy
      tepsilon = tepsilon_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener, freq1, U1)
      !
      timer_CALL t_mkspf%start()
      DO it = 1,nconf
         DO i = 1,S%nat3
            aux  = SQRT(tepsilon(:,i,it))
            aux2 = (aux-1)/(aux+1)
            reflectivity(:,i,it) = DBLE(aux2*CONJG(aux2))
         ENDDO
      ENDDO
      timer_CALL t_mkspf%stop()
      !
   END FUNCTION ir_reflectivity_q
   !
   ! Sum the self energy at the provided ener(ne) input energies
   ! \/o\________\\\_________________________________________/^>
   FUNCTION sum_selfnrg_spectre(S, sigma, freq, bose, V3sq, ne, ener, nu0)
      USE input_fc,           ONLY : ph_system_info
      USE merge_degenerate,   ONLY : merge_degen
      USE functions,          ONLY : sigma_mgo
      IMPLICIT NONE
      TYPE(ph_system_info),INTENT(in)   :: S
      REAL(DP),INTENT(in) :: sigma   ! smearing (regularization) (Ry)
      REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
      REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
      REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
      !
      INTEGER,INTENT(in)  :: ne           ! number of energies...
      REAL(DP),INTENT(in) :: ener(ne)     ! energies for which to compute the spectral function
      INTEGER,INTENT(in)  :: nu0(3)       ! first non-zero phomnon frequency
      !
      ! _P -> scattering, _M -> cohalescence
      REAL(DP) :: bose_P, bose_M      ! final/initial state populations
      REAL(DP) :: factor, sigma_, T=0._dp !freqtotm1
      REAL(DP) :: omega_P,  omega_M   ! \delta\omega
      REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
      COMPLEX(DP) :: ctm_P,ctm_M, reg
      COMPLEX(DP) :: ctm(ne)
      REAL(DP)    :: freqm1(S%nat3,3)  ! phonon energies (Ry)
      !
      INTEGER :: i,j,k, ie
      !
      ! Note: using the function result in an OMP reduction causes crash with ifort 14
      COMPLEX(DP) :: sum_selfnrg_spectre(ne,S%nat3)
      COMPLEX(DP),ALLOCATABLE :: spf(:,:)
      !
!     IF(sigma<=0._dp)THEN
!       CALL errore("sum_selfnrg_spectre","spf not implemented in the static limit. "&
!                   //"NEW: To do unshifted spf use 'spf imag'",1)
!     ENDIF
      !
      ALLOCATE(spf(ne,S%nat3))
      spf = (0._dp, 0._dp)
      !
      DO i = 1,S%nat3
         IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
         IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
         IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
      ENDDO
      !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(ie,i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2, &
!$OMP                     ctm_P,ctm_M,ctm,reg,factor) &
!$OMP             REDUCTION(+: spf) COLLAPSE(2)
      DO k = 1,S%nat3
         DO j = 1,S%nat3
            !
            bose_P   = 1 + bose(j,2) + bose(k,3)
            omega_P  = freq(j,2)+freq(k,3)
            omega_P2 = omega_P**2
            !
            bose_M   = bose(k,3) - bose(j,2)
            omega_M  = freq(j,2)-freq(k,3)
            omega_M2 = omega_M**2
            !
            ! A little optimization: precompute the parts that depends only on the energy
! ! !         DO ie = 1,ne
! ! !           ! regularization:
! ! !           reg = CMPLX(ener(ie), sigma, kind=DP)**2
! ! !           ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
! ! !           ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
! ! !           ctm(ie) = ctm_P + ctm_M
! ! !         ENDDO
            !
            DO i = 1,S%nat3
               !
               factor = V3sq(i,j,k) *freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
               !
               IF(ABS(sigma-666._dp)< 1._dp)THEN
                  sigma_ = sigma_mgo(freq(i,1),T) &
                     +sigma_mgo(freq(j,2),T) &
                     +sigma_mgo(freq(k,3),T)
                  !print*, RY_TO_CMM1*freq(i,1), RY_TO_CMM1*freq(j,2), RY_TO_CMM1*freq(k,3), RY_TO_CMM1*sigma_
               ELSE
                  sigma_ = sigma
               ENDIF

               DO ie = 1, ne
                  ! regularization:
                  reg = CMPLX(ener(ie), sigma_, kind=DP)**2
                  !
                  ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
                  ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
                  ctm(ie) = ctm_P + ctm_M
                  !
                  spf(ie,i) = spf(ie,i) + ctm(ie) * factor !V3sq(i,j,k) * freqtotm1
               ENDDO
! ! !           spf(:,i) = spf(:,i) + ctm(:) * factor !V3sq(i,j,k) * freqtotm1

               !
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !
      CALL merge_degen(ne, S%nat3, spf, freq(:,1))
      sum_selfnrg_spectre = spf
      DEALLOCATE(spf)
      !
   END FUNCTION sum_selfnrg_spectre
   !
   ! \/o\________\\\_________________________________________/^>
   ! Sum the self energy in a range of frequencies
   FUNCTION sum_selfnrg_modes(S, sigma, T, freq, bose, V3sq, nu0)
      USE input_fc,           ONLY : ph_system_info
      USE functions,          ONLY : df_bose
      USE merge_degenerate,   ONLY : merge_degen
      USE functions,          ONLY : sigma_mgo
      IMPLICIT NONE
      TYPE(ph_system_info),INTENT(in)   :: S
      REAL(DP),INTENT(in) :: sigma, T
      REAL(DP),INTENT(in) :: freq(S%nat3,3)
      REAL(DP),INTENT(in) :: bose(S%nat3,3)
      REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
      INTEGER,INTENT(in)  :: nu0(3)
      !
      ! _P -> scattering, _M -> cohalescence
      REAL(DP) :: bose_P, bose_M      ! final/initial state populations
      REAL(DP) :: freqtotm1
      REAL(DP) :: omega_P,  omega_M   ! \sigma\omega
      REAL(DP) :: omega_P2, omega_M2  ! \sigma\omega
      COMPLEX(DP) :: ctm_P, ctm_M, reg
      REAL(DP) :: freqm1(S%nat3,3), sigma_
      !
      INTEGER :: i,j,k
      ! Note: using the function result in an OMP reduction causes crash with ifort 14
      COMPLEX(DP) :: sum_selfnrg_modes(S%nat3)
      COMPLEX(DP) :: se(S%nat3)
      !
      se(:) = (0._dp, 0._dp)
      !
      freqm1 = 0._dp
      DO i = 1,S%nat3
         IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
         IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
         IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
      ENDDO
      !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,&
!$OMP                     omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtotm1) &
!$OMP             REDUCTION(+: se) COLLAPSE(2)
      DO k = 1,S%nat3
         DO j = 1,S%nat3
            !
            bose_P   = 1 + bose(j,2) + bose(k,3)
            bose_M   = bose(k,3) - bose(j,2)
            omega_P  = freq(j,2)+freq(k,3)
            omega_M  = freq(j,2)-freq(k,3)
            !
            IF(sigma>0._dp)THEN
               omega_P2 = omega_P**2
               omega_M2 = omega_M**2
            ELSE IF(sigma<0._dp)THEN ! (static limit)
               ctm_P = 2 * bose_P *omega_P/(omega_P**2+sigma**2)
               ctm_M = 2 * bose_M *omega_M/(omega_M**2+sigma**2)
            ELSE !IF (sigma==0._dp)THEN
               ! In the static limit with sigma=0 case we have to take the
               ! derivative of (n_3-n2)/(w_2-w_3) when w_2 is close to w_3
               IF(omega_P>0._dp)THEN
                  ctm_P = 2 * bose_P /omega_P
               ELSE
                  ctm_P = 0._dp
               ENDIF
               !
               IF(ABS(omega_M)>1.e-5_dp)THEN
                  ctm_M = 2 * bose_M /omega_M
               ELSE
                  IF(T>0._dp)THEN
                     ctm_M = -2* df_bose(0.5_dp * omega_P, T)
                  ELSE
                     ctm_M = 0._dp
                  ENDIF
               ENDIF
               !
            ENDIF
            !
            DO i = 1,S%nat3
               !
               !IF(freq(i,1)<1.d-3 .or. freq(j,2)<1.d-3 .or. freq(k,3)<1.d-3) CYCLE
               ! This comes from the definition of u_qj, Ref. 1. (there is an hidden factor 1/8)
               freqtotm1 = freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
               !
               ! regularization:
               IF(sigma>0._dp)THEN
                  IF(ABS(sigma-666._dp)<1._dp)THEN
                     sigma_ = sigma_mgo(freq(i,1),T) &
                        +sigma_mgo(freq(j,2),T) &
                        +sigma_mgo(freq(k,3),T)
                     !print*, RY_TO_CMM1*freq(i,1), RY_TO_CMM1*freq(j,2), RY_TO_CMM1*freq(k,3), RY_TO_CMM1*sigma_
                  ELSE
                     sigma_ = sigma
                  ENDIF

                  reg = CMPLX(freq(i,1), sigma_, kind=DP)**2
                  ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
                  ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
               ENDIF
               !
               !
               se(i) = se(i) + (ctm_P + ctm_M)*freqtotm1 * V3sq(i,j,k)
               !
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !
      CALL merge_degen(S%nat3, se, freq(:,1))
      sum_selfnrg_modes = se
      !
   END FUNCTION sum_selfnrg_modes
   !
  ! Sum the Imaginary part of the self energy for the phonon modes, uses gaussian function for energy delta
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_linewidth_modes(S, sigma, freq, bose, V3sq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE merge_degenerate,   ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    REAL(DP) :: sum_linewidth_modes(S%nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_C, bose_X ! final/initial state populations 
    REAL(DP) :: dom_C, dom_X   ! \delta\omega
    REAL(DP) :: ctm_C, ctm_X   !
    REAL(DP) :: freqtotm1, freqtotm1_23
    REAL(DP) :: freqm1(S%nat3,3)
    !REAL(DP),SAVE :: leftover_e
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)!, sigma_
    lw(:) = 0._dp
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,&
!$OMP                     freqtotm1_23,freqtotm1) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freqtotm1_23= freqm1(j,2) * freqm1(k,3)
        !
        DO i = 1,S%nat3
          !
          !sigma_= MIN(sigma, 0.5_dp*MAX(MAX(freq(i,1), freq(j,2)), freq(k,3)))
          !
          freqtotm1 = freqm1(i,1) * freqtotm1_23
          !IF(freqtot/=0._dp)THEN
          !
          dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          ctm_C = bose_C * f_gauss(dom_C, sigma)
          !
          dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          ctm_X = bose_X * f_gauss(dom_X, sigma)
          !
          lw(i) = lw(i) - pi * (ctm_C + ctm_X) * V3sq(i,j,k)*freqtotm1
          !
          !leftover_e = pi*freqtotm1 * (ctm_C*dom_C + ctm_X*dom_X) * V3sq(i,j,k)
          !ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(S%nat3, lw, freq(:,1))
    sum_linewidth_modes = lw
    !
  END FUNCTION sum_linewidth_modes

   ! \/o\________\\\_________________________________________/^>
   ! Add the elastic peak of Raman
   SUBROUTINE add_exp_t_factor(nconf, T, ne, nat3, ener, spectralf)
      USE functions,      ONLY : f_bose
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)
      INTEGER,INTENT(in)  :: ne, nat3
      REAL(DP),INTENT(in) :: ener(ne)
      !
      REAL(DP),INTENT(inout) :: spectralf(ne,nat3,nconf)
      !
      REAL(DP) :: factor(ne)
      INTEGER :: it,ie,nu
      !
      DO it = 1,nconf
         factor = (1 + f_bose(ener,T(it))) / ener
         !
         DO nu = 1,nat3
            DO ie = 1,ne
               spectralf(ie,nu,it) = factor(ie)*spectralf(ie,nu,it)
            ENDDO
         ENDDO
         !
      ENDDO
      !
   END SUBROUTINE add_exp_t_factor
   !
   ! \/o\________\\\_________________________________________/^>
   SUBROUTINE gauss_convolution(nconf, T, ne, nat3, ener, spectralf)
      USE functions,      ONLY : f_bose
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)  :: nconf
      REAL(DP),INTENT(in) :: T(nconf)
      INTEGER,INTENT(in)  :: ne, nat3
      REAL(DP),INTENT(in) :: ener(ne)
      !
      REAL(DP),INTENT(inout) :: spectralf(ne,nat3,nconf)
      !
      REAL(DP) :: convol(ne)
      !
      INTEGER :: it,ie,nu
      !
      convol = f_bose(ener,T(it))

      DO it = 1,nconf
         DO nu = 1,nat3
            DO ie = 1,ne
               spectralf(ie,nu,it) = convol(ie)*spectralf(ie,nu,it)
            ENDDO
         ENDDO
      ENDDO
      !
   END SUBROUTINE gauss_convolution
   !
   !
END MODULE linewidth


