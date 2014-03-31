MODULE linewidth

  USE kinds,    ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid, forceconst3_grid
  USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag, fftinterp_mat3, ip_cart2pat

  CONTAINS

! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION linewidth_q(xq0, T, sigma, S, grid, fc2, fc3) &
  RESULT(lw)
    !
    USE nanoclock
    !
    USE q_grid,         ONLY : q_grid_type
    
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    REAL(DP),INTENT(in) :: T     ! Kelvin
    REAL(DP),INTENT(in) :: sigma ! ry
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu
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
    CALL prepare_phq(xq(:,1), T, S, fc2, freq(:,1), bose(:,1), U(:,:,1))

    !
!     CALL start_nanoclock(lwtot)
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
!       CALL start_nanoclock(i_ph)
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL prepare_phq(xq(:,jq), T, S, fc2, freq(:,jq), bose(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
!       CALL stop_nanoclock(i_ph)
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
!       CALL start_nanoclock(d3time)
      CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
!       CALL stop_nanoclock(d3time)
      !
!       CALL start_nanoclock(c2pat)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
!       CALL stop_nanoclock(c2pat)
      
!       CALL start_nanoclock(tv3sq)
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
!       CALL stop_nanoclock(tv3sq)
      ! ------ end of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      !
!       CALL start_nanoclock(tsum)
      lw = lw + sum_modes( S, sigma, freq, bose, V3sq )
!       CALL stop_nanoclock(tsum)
      !
    ENDDO
    !
    lw = lw/grid%nq
    !
!     CALL stop_nanoclock(lwtot)

!     CALL print_nanoclock(lwtot)
!     CALL print_nanoclock(d3time)
!     !
!     CALL print_nanoclock(i_ph)
!     CALL print_nanoclock(c2pat)
!     CALL print_nanoclock(tv3sq)
!     CALL print_nanoclock(tsum)
!     CALL print_memory()
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION linewidth_q  
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_modes(S, sigma, freq, bose, V3sq)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : RY_TO_CMM1, pi
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma
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
    REAL(DP) :: freqtot
    !
    REAL(DP),PARAMETER :: norm = pi/8
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)
    lw(:) = 0._dp
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,freqtot,bose_X,bose_C,dom_X,dom_C,ctm_X,ctm_C) &
!$OMP             REDUCTION(+: lw) COLLAPSE(3)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        DO i = 1,S%nat3
          !
          freqtot = freq(i,1) * freq(j,2) * freq(k,3)
          !
          IF(ABS(freqtot)>1.d-8)THEN
            bose_X = bose(j,2) + bose(k,3) + 1
            bose_C = bose(j,2) - bose(k,3)
            !
            dom_X =(freq(i,1)+freq(j,2)-freq(k,3))
            dom_C =(freq(i,1)-freq(j,2)-freq(k,3))
            !
            ctm_X = 2 * bose_C * f_gauss(dom_X, sigma)
            ctm_C =     bose_X * f_gauss(dom_C, sigma)
            !
            lw(i) = lw(i) + norm/freqtot * (ctm_X + ctm_C) * V3sq(i,j,k)
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
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns a complex number:
  !  its real part is the "Tadpole" 3rd order contribution to lineshift
  !  its complex part is minus the linewidth
  FUNCTION lineshift_q(xq0, nconf, T, sigma, S, grid, fc2, fc3) &
  RESULT(ls)
    !
    USE q_grid,         ONLY : q_grid_type
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    !
    ! FUNCTION RESULT:
    COMPLEX(DP) :: ls(S%nat3,nconf)
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
    ls = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      !
      CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
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
        !
        ls(:,it) = ls(:,it) + sum_complex_modes( S, sigma(it), freq, bose, V3sq )
        !
      ENDDO
      !
    ENDDO
    !
    ! I think this 1/2 factor comes from 2\gamma = \sum ...
    ls = 0.5_dp * ls/grid%nq
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION lineshift_q  

  ! \/o\________\\\_________________________________________/^>
  ! Simple spectral function, computed as a superposition of Lorentzian functions
  FUNCTION simple_spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener) &
  RESULT(spectralf)
    !
    USE nanoclock
    !
    USE q_grid,         ONLY : q_grid_type
    
    IMPLICIT NONE
    ! ARGUMENTS:
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    ! FUNCTION RESULT:
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    ! Local variables
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    COMPLEX(DP) :: ls(S%nat3,nconf)
    !
    INTEGER :: it, i, ie
    REAL(DP) :: gamma, delta, omega, denom
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    REAL(DP) :: freq(S%nat3)
      !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    CALL freq_phq(xq0, S, fc2, freq, U)
    !
    ! use lineshift to compute 3rd order linewidth and lineshift
    ls = lineshift_q(xq0, nconf, T,sigma, S, grid, fc2, fc3)    
    !
    ! Compute and superpose the Lorentzian functions corresponding to each band
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          gamma = -DIMAG(ls(i,it))
          delta =   DBLE(ls(i,it))
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
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_complex_modes(S, delta, freq, bose, V3sq)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: delta
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M, freqtot      ! final/initial state populations 
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg, num
    !
    INTEGER :: i,j,k
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_complex_modes(S%nat3)
    COMPLEX(DP),ALLOCATABLE :: ls(:)
    ALLOCATE(ls(S%nat3))
    ls(:) = (0._dp, 0._dp)
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtot) &
!$OMP             REDUCTION(+: ls) COLLAPSE(2)
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
        DO i = 1,S%nat3
          !
          freqtot = 4*freq(i,1)*freq(j,2)*freq(k,3)
          !
          IF(ABS(freqtot)>1.d-8)THEN
            ! regularization:
            reg = CMPLX(freq(i,1), delta, kind=DP)**2
            !
            ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
            ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
            !
            ls(i) = ls(i) + (ctm_P + ctm_M) * V3sq(i,j,k) / freqtot
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_complex_modes = - 0.5_dp * ls
    !
  END FUNCTION sum_complex_modes
  !
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
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ener(ne)
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    !
    ! FUNCTION RESULT:
    COMPLEX(DP) :: spf(ne,S%nat3,nconf)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    INTEGER  :: i, ie
    REAL(DP) :: gamma, delta, omega, denom
    REAL(DP) :: spectralf(ne,S%nat3,nconf)
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    spf = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!$OMP END PARALLEL DO
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
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
        spf(:,:,it) = spf(:,:,it) + sum_spectral_function( S, sigma(it), freq, bose, V3sq, ne, ener )
        !
      ENDDO
      !
    ENDDO
    !
    ! I think this 1/2 factor comes from 2\gamma = \sum ...
    spf = 0.5_dp * spf/grid%nq
    !
!     spectralf = -DIMAG(spf)
    DO it = 1,nconf
      DO i = 1,S%nat3
        DO ie = 1, ne
          gamma = -DIMAG(spf(ie,i,it))
          delta =   DBLE(spf(ie,i,it))
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
    DEALLOCATE(U, V3sq, D3)
    !
  END FUNCTION spectre_q  
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_spectral_function(S, sigma, freq, bose, V3sq, ne, ener)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma           ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
    !
    INTEGER,INTENT(in)  :: ne           ! number of energies...
    REAL(DP),INTENT(in) :: ener(ne)     ! energies for which to compute the spectral function
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M, freqtot      ! final/initial state populations 
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg, num
    !
    INTEGER :: i,j,k, ie
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_spectral_function(ne,S%nat3)
    COMPLEX(DP),ALLOCATABLE :: spf(:,:)
    !
    ALLOCATE(spf(ne,S%nat3))
    spf = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(ie,i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtot) &
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
        DO i = 1,S%nat3
          !
          freqtot = 4*freq(i,1)*freq(j,2)*freq(k,3)
          !
          IF(freqtot/=0._dp)THEN
            !
            DO ie = 1, ne
            ! regularization:
              reg = CMPLX(ener(ie), sigma, kind=DP)**2
              !
              ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
              ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
              !
              spf(ie,i) = spf(ie,i) + (ctm_P + ctm_M) * V3sq(i,j,k) / freqtot
            ENDDO
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_spectral_function = - 0.5_dp * spf
    DEALLOCATE(spf)
    !
  END FUNCTION sum_spectral_function
  !
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
  SUBROUTINE freq_phq(xq, S, fc2, freq, U)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)               :: xq(3)
      TYPE(ph_system_info),INTENT(in)   :: S
      TYPE(forceconst2_grid),INTENT(in) :: fc2
      REAL(DP),INTENT(out)              :: freq(S%nat3)
      COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
      REAL(DP),PARAMETER :: eps = 1.e-12_dp
      !
      CALL fftinterp_mat2(xq, S, fc2, U)
      CALL mat2_diag(S, U, freq)
      U = CONJG(U)
      WHERE    (freq >  eps)
        freq = SQRT(freq)
      ELSEWHERE(freq < -eps)
        freq = -SQRT(-freq)
      ELSEWHERE ! i.e. freq=0
        ! this should only happen at Gamma for the 3 acoustic bands
        ! and they normally do not matter for symmetry or velocity = 0
        freq = 0._dp
      ENDWHERE
      !
  END SUBROUTINE freq_phq
  !
  SUBROUTINE bose_phq(T, nat3, freq, bose)
    USE functions,      ONLY : f_bose
    IMPLICIT NONE
      REAL(DP),INTENT(in)  :: T
      INTEGER,INTENT(in)   :: nat3
      REAL(DP),INTENT(in)  :: freq(nat3)
      REAL(DP),INTENT(out) :: bose(nat3)
      REAL(DP),PARAMETER :: eps = 1.e-12_dp
      !
      ! Is the following mess really necessary? (3 days later: it is)
      WHERE    (freq*T >  eps)
        bose = f_bose(freq, T)
      ELSEWHERE
        bose = 0._dp
      ENDWHERE
      !
  END SUBROUTINE bose_phq
  !
END MODULE linewidth
