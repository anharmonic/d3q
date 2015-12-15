!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! References:
! [1] Calandra, Lazzeri, Mauri : Physica C 456 (2007) 38-44
! [2] arXiv: 1312.7467v1
! [3] R.Biaco et.al. in prep.

MODULE dynbubble

  USE kinds,           ONLY : DP
  USE more_constants,  ONLY : eps_freq
  !USE input_fc,  ONLY : ph_system_info, forceconst2_grid 
  !USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag, ip_cart2pat
  !USE fc3_interpolate, ONLY : forceconst3
  USE timers
#include "mpi_thermal.h"

  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns the complex self energy from the bubble(?) diagram:
  !  its real part is the order contribution to lineshift
  !  its complex part is the intrinsic linewidth
  FUNCTION dynbubble_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info, div_mass_dyn
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq,&
                                 set_nu0, ip_cart2pat, dyn_cart2pat
    USE fc3_interpolate,  ONLY : forceconst3
    USE constants,        ONLY : RY_TO_CMM1
    USE mpi_thermal,      ONLY : mpi_bsum
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
    ! FUNCTION RESULT:
    COMPLEX(DP) :: dynbubble_q(S%nat3,S%nat3,nconf)
    COMPLEX(DP) :: dyn(S%nat3,S%nat3,nconf) ! aux
    COMPLEX(DP) :: dyndio(S%nat3,S%nat3) ! aux
!     REAL(DP) :: sqtrace(S%nat3) ! aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
    INTEGER :: nu0(3)
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    dyn = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes at q1
    ! RAF: inserito un meno
    xq(:,1) = -xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!/nope/!$OMP END PARALLEL DO
      !
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
        !
        dyndio = grid%w(iq)*sum_dynbubble( S%nat3, T(it), sigma(it), freq, bose, D3, nu0 )
!         WRITE(*, *) iq, it
!         WRITE(*, "(3(2f12.4,3x))") -0.5_dp * dyndio/grid%nq
        dyn(:,:,it) = dyn(:,:,it) + dyndio
        !
      ENDDO
      !
    ENDDO
    !
    IF(grid%scattered) CALL mpi_bsum(S%nat3,S%nat3,nconf, dyn)
    
    DO it=1,nconf
        CALL dyn_cart2pat(dyn(:,:,it), S%nat3, U(:,:,1), -1)
    ENDDO
    dynbubble_q = -0.5_dp * dyn
    !
    DEALLOCATE(U)
    !
  END FUNCTION dynbubble_q  
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the the phonon spectral function tracing the self energy
  FUNCTION tracebubble_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener) &
  RESULT (spectralf)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info, div_mass_dyn, write_dyn
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq,&
                                 set_nu0, ip_cart2pat, dyn_cart2pat, fftinterp_mat2
    USE fc3_interpolate,  ONLY : forceconst3
    USE constants,        ONLY : RY_TO_CMM1
    USE mpi_thermal,      ONLY : mpi_bsum
    USE functions,        ONLY : invzmat
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
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: dyn(:,:)!, dyn0(:,:)
    COMPLEX(DP),ALLOCATABLE :: bubble(:,:,:,:), bubbleq(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it, ie, n, m
    INTEGER :: nu0(3)
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3), freqm1(S%nat3)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(dyn(S%nat3,S%nat3)) ! , dyn0(S%nat3,S%nat3))
    ALLOCATE(bubble(ne,S%nat3,S%nat3,nconf), bubbleq(ne,S%nat3,S%nat3))
    !
    dyn = (0._dp, 0._dp)
    !
    !FIXME: I'm doing this interpolation twice in a row
    !CALL fftinterp_mat2(xq0, S%nat3, fc2, dyn0)
    !
    ! Compute eigenvalues, eigenmodes at q1
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    freqm1 = 0._dp
    DO n = 1,S%nat3
      IF(n>nu0(1)) freqm1(n) = 1/freq(n,1)
    ENDDO
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!/nope/!$OMP END PARALLEL DO
      !
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
        !
        bubbleq = -0.5_dp*grid%w(iq)*sum_dynbubble_spectre( S%nat3, T(it), sigma(it), freq, bose, D3, ne, ener, nu0 )
        !
        bubble(:,:,:,it) = bubble(:,:,:,it) + bubbleq(:,:,:)
        !  
      ENDDO
      !
    ENDDO
    !
    IF(grid%scattered) CALL mpi_bsum(S%nat3,S%nat3,ne,nconf, bubble)
    
    !
    ! This is going to take a bit of memory for NE very large (maybe fix in the future)
    DO it=1,nconf
    DO ie = 1,ne
      ! Change the order of the indexes
      FORALL (m=1:S%nat3, n=1:S%nat3)
        dyn(n,m) = bubble(ie,n,m,it)
      END FORALL
!      dyn = 1.*RY_TO_CMM1 !0._dp
      !
      ! Rotate, sum to unperturbed dynmat, invert
      FORALL(n=1:S%nat3) dyn(n,n) = dyn(n,n) &
                          + 0.5_dp * (freq(n,1)**2 - ener(ie)**2)*freqm1(n)
      DO m=1,S%nat3
      DO n=1,S%nat3
        IF(n/=m) dyn(n,m) = 0._dp
      ENDDO
      ENDDO
      !
      ! Is the next line necessary?
      !CALL dyn_cart2pat(dyn, S%nat3, U(:,:,1), -1)
      !
      ! write out the dyn mat to file
      !CALL write_dyn("f", xq0, dyn, S)
      !
      CALL invzmat(S%nat3, dyn)
      !
      ! Change the order of indexes again, take the diagonal (trace will be done later)
      FORALL(n=1:S%nat3) spectralf(ie,n,it) = 2*DIMAG(dyn(n,n))
    ENDDO
    ENDDO
    !
    DEALLOCATE(U)
    !
  END FUNCTION tracebubble_q
  !
  ! \/o\________\\\_________________________________________/^>
  ! Sum the self energy for the phonon modes
  FUNCTION sum_dynbubble(nat3, T, sigma, freq, bose, V3, nu0)
    USE functions,    ONLY : df_bose
    USE constants,    ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nat3
    REAL(DP),INTENT(in) :: T, sigma
    REAL(DP),INTENT(in) :: freq(nat3,1:3)
    REAL(DP),INTENT(in) :: bose(nat3,2:3)
    COMPLEX(DP),INTENT(in) :: V3(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtotm1_23 !nm, freqtotm1_23n, freqtotm1_23
    REAL(DP) :: omega_P,  omega_M, omega_P2, omega_M2
    COMPLEX(DP) :: ctm_P, ctm_M, reg
    !
    INTEGER :: i, j,k, n,m
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_dynbubble(nat3,nat3)
    COMPLEX(DP) :: dyn(nat3,nat3), aux
    REAL(DP)    :: freqm1(nat3,2:3)
    !REAL(DP)    :: sqfreqm1(nat3)
    REAL(DP)    :: sqfreq(nat3)
    !
    ! Prepare some auxiliary
    freqm1 = 0._dp
    !sqfreqm1 = 0._dp
    DO i = 1,nat3
      sqfreq(i) = DSQRT(freq(i,1))
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
    !
    dyn(:,:) = (0._dp, 0._dp)
    !
    DO k = 1,nat3
      DO j = 1,nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        bose_M   = bose(k,3) - bose(j,2)
        omega_P  = freq(j,2)+freq(k,3)
        omega_M  = freq(j,2)-freq(k,3)
        !
        IF(sigma>0._dp)THEN
          ! In this case the factors have to be computed for each omega1, I precompute what I can
          omega_P2 = omega_P**2
          omega_M2 = omega_M**2
        ELSE IF(sigma<0._dp)THEN
          ! In this case I can precompute everything
          ctm_P = 2 * bose_P *omega_P/(omega_P**2+sigma**2)
          ctm_M = 2 * bose_M *omega_M/(omega_M**2+sigma**2)
        ELSE IF (sigma==0._dp)THEN
          ! I can pre-compute everything but it is a bit a mess
          IF(ABS(omega_P)>0._dp)THEN
            ctm_P = 2 * bose_P /omega_P
          ELSE
            ctm_P = 0._dp
          ENDIF
          !
          IF(ABS(omega_M)>1.e-5_dp)THEN
            !bose_M   = bose(k,3) - bose(j,2)
            ctm_M = 2 * bose_M /omega_M
          ELSE
            IF(T>0._dp.and.ABS(omega_P)>0._dp)THEN
              ctm_M = -2* df_bose(0.5_dp * omega_P, T)
              !ctm_M = 0._dp
            ELSE
              ctm_M = 0._dp
            ENDIF
          ENDIF
          !
        ENDIF
        !
        !
        freqtotm1_23 = freqm1(j,2)*freqm1(k,3)
        !
        DO n = 1,nat3
          !
          !freqtotm1_23n = freqtotm1_23 * sqfreqm1(n)
          !
          DO m = 1,nat3
            !
            !freqtotm1_23nm = freqtotm1_23n * sqfreqm1(m)
            !
            IF(sigma>0._dp)THEN
              reg = CMPLX(sqfreq(n)*sqfreq(m), sigma, kind=DP)**2
              ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
              ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
            ENDIF
            ! RAF: scambiato le coniugazioni complesse
            dyn(m,n) = dyn(m,n) + (ctm_P + ctm_M) * freqtotm1_23 &
                                 * V3(m,j,k) * CONJG(V3(n,j,k))
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    sum_dynbubble = dyn
    !
  END FUNCTION sum_dynbubble
  ! \/o\________\\\_________________________________________/^>
  ! Sum the dynamic self energy bubble for a list of energies
  FUNCTION sum_dynbubble_spectre(nat3, T, sigma, freq, bose, V3, ne, ener, nu0)
    USE functions,    ONLY : df_bose
    USE constants,    ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nat3
    REAL(DP),INTENT(in) :: T, sigma
    REAL(DP),INTENT(in) :: freq(nat3,1:3)
    REAL(DP),INTENT(in) :: bose(nat3,2:3)
    COMPLEX(DP),INTENT(in) :: V3(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: ne           ! number of energies...
    REAL(DP),INTENT(in) :: ener(ne)     ! energies for which to compute the spectral function
    INTEGER,INTENT(in)  :: nu0(3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtotm1_23 !nm, freqtotm1_23n, freqtotm1_23
    REAL(DP) :: omega_P,  omega_M, omega_P2, omega_M2
    COMPLEX(DP) :: ctm_P, ctm_M, reg
    !
    INTEGER :: i, j,k, n, m, ie
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_dynbubble_spectre(ne,nat3,nat3)
    COMPLEX(DP) :: dyn(ne,nat3,nat3)
    COMPLEX(DP) :: ctm(ne)
    COMPLEX(DP) :: aux, factor1, factor2
    REAL(DP)    :: freqm1(nat3,2:3)
    REAL(DP)    :: sqfreqm1(nat3)
    !REAL(DP)    :: sqfreq(nat3)
    !
    ! Prepare some auxiliary
    freqm1 = 0._dp
    !
    IF (sigma<0._dp) CALL errore("sum_dynbubble_spectre", "this is dynamic only", 1)
    
    !sqfreqm1 = 0._dp
    DO i = 1,nat3
      !sqfreq(i) = DSQRT(freq(i,1))
      IF(i>=nu0(1)) sqfreqm1(i) = DSQRT(0.5_dp/freq(i,1))
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO
    !
    dyn = (0._dp, 0._dp)
    !
    DO k = 1,nat3
      DO j = 1,nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        bose_M   = bose(k,3) - bose(j,2)
        omega_P  = freq(j,2)+freq(k,3)
        omega_M  = freq(j,2)-freq(k,3)
        !
        DO ie = 1,ne
          ! regularization:
          reg = CMPLX(ener(ie), sigma, kind=DP)**2
          ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
          ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
          ctm(ie) = ctm_P + ctm_M
        ENDDO
        !
        freqtotm1_23 = freqm1(j,2)*freqm1(k,3)
        !
        DO n = 1,nat3
          !
          factor1 = freqtotm1_23 * V3(n,j,k) * sqfreqm1(n)
          DO m = 1,nat3
            !
            ! Do all the energies together 
            ! (I get on output the indexes in the wrong order, but here speed is more critical)
            factor2 = factor1 * CONJG(V3(m,j,k)) * sqfreqm1(m)
            dyn(:,m,n) = dyn(:,m,n) + ctm(:) * factor2
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    sum_dynbubble_spectre = dyn
    !
  END FUNCTION sum_dynbubble_spectre
  !
END MODULE dynbubble
