! References:
! [1] Calandra, Lazzeri, Mauri : Physica C 456 (2007) 38-44
! [2] arXiv: 1312.7467v1

MODULE linewidth

  USE kinds,     ONLY : DP
  USE input_fc,  ONLY : ph_system_info, forceconst2_grid 
  USE interp_fc, ONLY : fftinterp_mat2, mat2_diag, ip_cart2pat
  USE sparse_fc, ONLY : forceconst3

  CONTAINS

! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION linewidth_q(xq0, nconf, T, sigma, S, grid, fc2, fc3) &
  RESULT(lw)
    !
    USE nanoclock
    !
    USE q_grid,         ONLY : q_grid_type
    
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(S%nat3,nconf) ! ry
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(S%nat3,nconf)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !

    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    lw = 0._dp
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      !
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        lw(:,it) = lw(:,it) + sum_modes( S, sigma(:,it), freq, bose, V3sq )
      ENDDO
      !
    ENDDO
    !
    lw = -0.5_dp * lw/grid%nq
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION linewidth_q  
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns the complex self energy from the bubble(?) diagram:
  !  its real part is the order contribution to lineshift
  !  its complex part is the intrinsic linewidth
  FUNCTION selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
    !
    USE q_grid,         ONLY : q_grid_type
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(S%nat3,nconf) ! ry
    !
    ! FUNCTION RESULT:
    COMPLEX(DP) :: selfnrg_q(S%nat3,nconf)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    selfnrg_q = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes at q1
    xq(:,1) = xq0
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      !
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        !
        selfnrg_q(:,it) = selfnrg_q(:,it) + sum_selfnrg_modes( S, sigma(:,it), freq, bose, V3sq )
        !
      ENDDO
      !
    ENDDO
    !
    selfnrg_q = -0.5_dp * selfnrg_q/grid%nq
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION selfnrg_q  
  ! \/o\________\\\_________________________________________/^>
  ! Simple spectral function, computed as a superposition of Lorentzian functions
  FUNCTION simple_spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener) &
  RESULT(spectralf)
    !
    USE nanoclock
    !
    USE q_grid,         ONLY : q_grid_type
    USE constants,      ONLY : RY_TO_CMM1
    
    IMPLICIT NONE
    ! ARGUMENTS:
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(S%nat3,nconf) ! ry
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    ! Local variables
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    COMPLEX(DP) :: selfnrg(S%nat3,nconf)
    !
    INTEGER :: it, i, ie
    REAL(DP) :: gamma, delta, omega, denom
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    REAL(DP) :: freq(S%nat3)
      !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    CALL freq_phq_safe(xq0, S, fc2, freq, U)
    !
    ! use lineshift to compute 3rd order linewidth and lineshift
    selfnrg = selfnrg_q(xq0, nconf, T,sigma, S, grid, fc2, fc3)    
    DO it = 1,nconf
      WRITE(*,'(30x,6e12.2,5x,6e12.2)') -DIMAG(selfnrg(:,it))*RY_TO_CMM1, DBLE(selfnrg(:,it))*RY_TO_CMM1
    ENDDO
    !
    ! Compute and superpose the Lorentzian functions corresponding to each band
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          gamma =  -DIMAG(selfnrg(i,it))
          delta =   DBLE(selfnrg(i,it))
          omega = freq(i)
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
    !
  END FUNCTION simple_spectre_q  
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Full spectral function, computed as in eq. 1 of arXiv:1312.7467v1
  FUNCTION spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener) &
  RESULT(spectralf)
    !
    USE nanoclock
    !
    USE q_grid,         ONLY : q_grid_type
    
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
    TYPE(q_grid_type),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(S%nat3,nconf) ! ry
    !
    ! To interpolate D2 and D3:
    INTEGER :: iq, jq, nu, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie
    REAL(DP) :: gamma, delta, omega, denom
    COMPLEX(DP),ALLOCATABLE :: selfnrg(:,:,:)
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(selfnrg(ne,S%nat3,nconf))
    !
    selfnrg = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      DO it = 1,nconf
        ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!$OMP END PARALLEL DO
        selfnrg(:,:,it) = selfnrg(:,:,it) + sum_selfnrg_spectre( S, sigma(:,it), freq, bose, V3sq, ne, ener )
        !
      ENDDO
      !
    ENDDO
    !
    selfnrg = -0.5_dp * selfnrg/grid%nq
    !
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          gamma =  -DIMAG(selfnrg(ie,i,it))
          delta =   DBLE(selfnrg(ie,i,it))
          omega = freq(i,1)
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
    !
    DEALLOCATE(U, V3sq, D3, selfnrg)
    !
  END FUNCTION spectre_q
  !
  ! Sum the self energy at the provided ener(ne) input energies
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_selfnrg_spectre(S, sigma, freq, bose, V3sq, ne, ener)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma(S%nat3)   ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
    !
    INTEGER,INTENT(in)  :: ne           ! number of energies...
    REAL(DP),INTENT(in) :: ener(ne)     ! energies for which to compute the spectral function
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtot, freqtotm1, eta
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg, num
    !
    INTEGER :: i,j,k, ie
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_selfnrg_spectre(ne,S%nat3)
    COMPLEX(DP),ALLOCATABLE :: spf(:,:)
    !
    ALLOCATE(spf(ne,S%nat3))
    spf = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(ie,i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtot,freqtotm1) &
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
        eta = sigma(k)+sigma(j)
        !
        DO i = 1,S%nat3
          !
          ! This comes from the definition of u_qj, Ref. 1.
          freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
          !
          IF (freqtot/=0._dp) THEN
            freqtotm1 = 1 / freqtot
            !
            DO ie = 1, ne
              ! regularization:
              reg = CMPLX(ener(ie), eta, kind=DP)**2
              !
              ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
              ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
              !
              spf(ie,i) = spf(ie,i) + (ctm_P + ctm_M) * V3sq(i,j,k) * freqtotm1
            ENDDO
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_selfnrg_spectre = spf
    DEALLOCATE(spf)
    !
  END FUNCTION sum_selfnrg_spectre
  !
  ! \/o\________\\\_________________________________________/^>
  ! Sum the self energy for the phonon modes
  FUNCTION sum_selfnrg_modes(S, sigma, freq, bose, V3sq)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma(S%nat3)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtot, freqtotm1, eta
    REAL(DP) :: omega_P,  omega_M   ! \sigma\omega
    REAL(DP) :: omega_P2, omega_M2  ! \sigma\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg
    !
    INTEGER :: i,j,k
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_selfnrg_modes(S%nat3)
    COMPLEX(DP),ALLOCATABLE :: se(:)
    ALLOCATE(se(S%nat3))
    se(:) = (0._dp, 0._dp)
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtot,freqtotm1) &
!$OMP             REDUCTION(+: se) COLLAPSE(2)
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
        eta = sigma(k)+sigma(j)
        !
        DO i = 1,S%nat3
          !
          ! This comes from the definition of u_qj, Ref. 1.
          freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
          !
          IF(freqtot/=0._dp)THEN
            freqtotm1 = 1 / freqtot
            ! regularization:
            reg = CMPLX(freq(i,1), eta, kind=DP)**2
            !
            ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
            ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
            !
            se(i) = se(i) + (ctm_P + ctm_M) * V3sq(i,j,k) * freqtotm1
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_selfnrg_modes = se
    !
  END FUNCTION sum_selfnrg_modes
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_modes(S, sigma, freq, bose, V3sq)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma(S%nat)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    !
    REAL(DP) :: sum_modes(S%nat3)
    !
    ! _X -> scattering, _C -> cohalescence
    REAL(DP) :: bose_X, bose_C ! final/initial state populations 
    REAL(DP) :: dom_X, dom_C   ! \delta\omega
    REAL(DP) :: ctm_X, ctm_C   !
    REAL(DP) :: freqtot, eta
    !
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)
    lw(:) = 0._dp
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,freqtot,bose_X,bose_C,dom_X,dom_C,ctm_X,ctm_C) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_X = bose(j,2) + bose(k,3) + 1
        bose_C = bose(j,2) - bose(k,3)
        eta = sigma(k)+sigma(j)
        !
        DO i = 1,S%nat3
          !
          freqtot = 8*freq(i,1) * freq(j,2) * freq(k,3)
          !
          IF(ABS(freqtot)/=0._dp)THEN
            !
            dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
            ctm_X = bose_C * f_gauss(dom_C, eta)
            !
            dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
            ctm_C = 2*bose_X * f_gauss(dom_X, eta)
            !
            lw(i) = lw(i) + pi/freqtot * (ctm_X + ctm_C) * V3sq(i,j,k)
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_modes = lw
    !
  END FUNCTION sum_modes
  !
  ! Auxiliary subroutines follow:
  ! \/o\________\\\_________________________________________/^>
  ! Interpolate dynamical matrice at q, diagonalize it and compute Bose-Einstein distribution
  SUBROUTINE prepare_phq(xq, T, S, fc2, freq, bose, U)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)  :: xq(3), T
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      REAL(DP),INTENT(out) :: freq(S%nat3), bose(S%nat3)
      COMPLEX(DP),INTENT(out) :: U(S%nat3,S%nat3)
      REAL(DP),PARAMETER :: eps = 1.e-12_dp
      !
      CALL fftinterp_mat2(xq, S, fc2, U)
      CALL mat2_diag(S, U, freq)
      ! Is the following mess really necessary? (3 days later: it is)
      WHERE    (freq >  eps)
        freq = SQRT(freq)
        bose = f_bose(freq, T)
      ELSEWHERE(freq < -eps)
        freq = -SQRT(-freq)
        bose = f_bose(freq, T)
      ELSEWHERE ! i.e. freq=0
        ! this should only happen at Gamma for the 3 acoustic bands
        ! and they normally do not matter for symmetry or velocity = 0
        freq = 0._dp
        bose = 0._dp
      ENDWHERE
      U = CONJG(U)
      !
  END SUBROUTINE prepare_phq
  !
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq_safe(xq, S, fc2, freq, U)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)               :: xq(3)
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      REAL(DP),INTENT(out)              :: freq(S%nat3)
      COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
      REAL(DP),PARAMETER :: epsq = 1.e-6_dp
      REAL(DP) :: cq(3)
      !
      CALL fftinterp_mat2(xq, S, fc2, U)
      CALL mat2_diag(S, U, freq)
      U = CONJG(U)
      
      ! Set patterns and frequency to exactly zero for Gamma (and Gamma+G)
      cq = xq
      CALL cryst_to_cart(1,cq,S%at,-1)
      IF(ALL( ABS(cq-INT(cq))<epsq ) )THEN
        freq(1:3) = 0._dp
!         U(:,1:3) = (0._dp, 0._dp)
      ENDIF
      
      WHERE    (freq >  0._dp)
        freq = SQRT(freq)
      ELSEWHERE(freq < 0._dp)
        freq = -SQRT(-freq)
      ELSEWHERE
        freq = 0._dp
      ENDWHERE
      !
  END SUBROUTINE freq_phq_safe
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq(xq, S, fc2, freq, U)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)               :: xq(3)
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      REAL(DP),INTENT(out)              :: freq(S%nat3)
      COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
      REAL(DP),PARAMETER :: eps = 1.e-16_dp
      !
      CALL fftinterp_mat2(xq, S, fc2, U)
      CALL mat2_diag(S, U, freq)
      U = CONJG(U)
      WHERE    (freq >  eps)
        freq = SQRT(freq)
      ELSEWHERE(freq < -eps)
        freq = -SQRT(-freq)
      ELSEWHERE ! i.e. freq=0
        freq = 0._dp
      ENDWHERE
      !
  END SUBROUTINE freq_phq
  !
  ! Compute Bose-Einstein distribution of freq
  SUBROUTINE bose_phq(T, nat3, freq, bose)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)  :: T
      INTEGER,INTENT(in)   :: nat3
      REAL(DP),INTENT(in)  :: freq(nat3)
      REAL(DP),INTENT(out) :: bose(nat3)
      !
      ! Is the following mess really necessary? (3 days later: it is)
      !WHERE    (freq*T >  eps)
      WHERE    (freq /=  0._dp)
        bose = f_bose(freq, T)
      ELSEWHERE
        bose = 0._dp
      ENDWHERE
      !
  END SUBROUTINE bose_phq
  !
  !
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
    INTEGER :: it,ie,je,nu
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
END MODULE linewidth
