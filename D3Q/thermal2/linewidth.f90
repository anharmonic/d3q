!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! References:
! [1] Calandra, Lazzeri, Mauri : Physica C 456 (2007) 38-44
! [2] arXiv: 1312.7467v1
#define timer_CALL CALL

MODULE linewidth

  USE kinds,     ONLY : DP
  !USE input_fc,  ONLY : ph_system_info, forceconst2_grid 
  !USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag, ip_cart2pat
  !USE fc3_interpolate, ONLY : forceconst3
  USE timers

  CONTAINS

! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION linewidth_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc3_interpolate,  ONLY : forceconst3
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0, ip_cart2pat
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
    REAL(DP) :: linewidth_q(S%nat3,nconf)
    REAL(DP) :: lw(S%nat3,nconf) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it, nu0(3)
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
      timer_CALL t_freq%start()
    xq(:,1) = xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
      timer_CALL t_freq%stop()
    !
    DO iq = 1, grid%nq
      !
        timer_CALL t_freq%start()
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!/nope/!$OMP END PARALLEL DO
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
          timer_CALL t_bose%start()
        ! Compute bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
          timer_CALL t_bose%stop()
          timer_CALL t_sum%start()
        lw(:,it) = lw(:,it) + sum_linewidth_modes( S, sigma(it), freq, bose, V3sq, nu0 )
          timer_CALL t_sum%stop()
      ENDDO
      !
    ENDDO
    !
    linewidth_q = -0.5_dp * lw/grid%nq
    !
    DEALLOCATE(U, V3sq, D3)
    !
  END FUNCTION linewidth_q  
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns the complex self energy from the bubble(?) diagram:
  !  its real part is the order contribution to lineshift
  !  its complex part is the intrinsic linewidth
  FUNCTION selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0, ip_cart2pat
    USE fc3_interpolate,  ONLY : forceconst3
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
    COMPLEX(DP) :: selfnrg_q(S%nat3,nconf)
    COMPLEX(DP) :: se(S%nat3,nconf) ! aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
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
    ! Compute eigenvalues, eigenmodes at q1
    xq(:,1) = xq0
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
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      DO it = 1,nconf
        ! Compute bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
        !
        se(:,it) = se(:,it) + sum_selfnrg_modes( S, sigma(it), freq, bose, V3sq, nu0 )
        !
      ENDDO
      !
    ENDDO
    !
    selfnrg_q = -0.5_dp * se/grid%nq
    !
    DEALLOCATE(U, V3sq)
    !
  END FUNCTION selfnrg_q  
  ! \/o\________\\\_________________________________________/^>
  ! Simple spectral function, computed as a superposition of Lorentzian functions
  FUNCTION simple_spectre_q(xq0, nconf, T, sigma, S, grid, fc2, fc3, ne, ener) &
  RESULT(spectralf)
    USE q_grids,            ONLY : q_grid
    USE constants,          ONLY : RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe
    USE fc3_interpolate,    ONLY : forceconst3
    
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
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
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
    selfnrg = selfnrg_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)    
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
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq, set_nu0, ip_cart2pat
    USE fc3_interpolate,  ONLY : forceconst3
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
    ! To interpolate D2 and D3:
    INTEGER :: iq, jq, nu, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie
    INTEGER  :: nu0(3)
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
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      DO it = 1,nconf
        ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
        selfnrg(:,:,it) = selfnrg(:,:,it) + sum_selfnrg_spectre( S, sigma(it), freq, bose, V3sq, ne, ener, nu0 )
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
  FUNCTION sum_selfnrg_spectre(S, sigma, freq, bose, V3sq, ne, ener, nu0)
    USE input_fc,           ONLY : ph_system_info
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
    REAL(DP) :: freqtot, freqtotm1
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
    DO k = nu0(3),S%nat3
      DO j = nu0(2),S%nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        omega_P  = freq(j,2)+freq(k,3)
        omega_P2 = omega_P**2
        !
        bose_M   = bose(k,3) - bose(j,2)
        omega_M  = freq(j,2)-freq(k,3)
        omega_M2 = omega_M**2
        !
        DO i = nu0(1),S%nat3
          !
          ! This comes from the definition of u_qj, Ref. 1.
          freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
          !
!           IF (freqtot/=0._dp) THEN
          freqtotm1 = 1 / freqtot
          !
          DO ie = 1, ne
            ! regularization:
            reg = CMPLX(ener(ie), sigma, kind=DP)**2
            !
            ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
            ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
            !
            spf(ie,i) = spf(ie,i) + (ctm_P + ctm_M) * V3sq(i,j,k) * freqtotm1
          ENDDO
!           ENDIF
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
  FUNCTION sum_selfnrg_modes(S, sigma, freq, bose, V3sq, nu0)
    USE input_fc,           ONLY : ph_system_info
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtot
    REAL(DP) :: omega_P,  omega_M   ! \sigma\omega
    REAL(DP) :: omega_P2, omega_M2  ! \sigma\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg
    !
    INTEGER :: i,j,k
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    COMPLEX(DP) :: sum_selfnrg_modes(S%nat3)
    COMPLEX(DP) :: se(S%nat3)
    !
    se(:) = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg,freqtot) &
!$OMP             REDUCTION(+: se) COLLAPSE(2)
    DO k = nu0(3),S%nat3
      DO j = nu0(2),S%nat3
        !
        bose_P   = 1 + bose(j,2) + bose(k,3)
        omega_P  = freq(j,2)+freq(k,3)
        omega_P2 = omega_P**2
        !
        bose_M   = bose(k,3) - bose(j,2)
        omega_M  = freq(j,2)-freq(k,3)
        omega_M2 = omega_M**2
        !
        DO i = nu0(1),S%nat3
          !
          ! This comes from the definition of u_qj, Ref. 1.
          freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
          !
!         IF(freqtot/=0._dp)THEN
          ! regularization:
          reg = CMPLX(freq(i,1), sigma, kind=DP)**2
          !
          ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
          ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
          !
          se(i) = se(i) + (ctm_P + ctm_M)/freqtot * V3sq(i,j,k)
!         ENDIF
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
  ! Sum the Imaginary part of the self energy for the phonon modes, uses gaussian function for energy delta
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_linewidth_modes(S, sigma, freq, bose, V3sq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    USE input_fc,           ONLY : ph_system_info
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
    REAL(DP) :: freqtot, freq23
    !
    REAL(DP),PARAMETER :: norm = pi/8
    !
    INTEGER :: i,j,k
    REAL(DP) :: lw(S%nat3)
    lw(:) = 0._dp
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,freqtot,bose_C,bose_X,dom_C,dom_X,ctm_C,ctm_X,freq23) &
!$OMP             REDUCTION(+: lw) COLLAPSE(2)
    DO k = nu0(3),S%nat3
      DO j = nu0(2),S%nat3
        !
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freq23 = freq(j,2) * freq(k,3)
        !
        DO i = nu0(1),S%nat3
          !
          freqtot = freq(i,1) * freq23
          !IF(freqtot/=0._dp)THEN
          !
          dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          ctm_C = bose_C * f_gauss(dom_C, sigma)
          !
          dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          ctm_X = bose_X * f_gauss(dom_X, sigma)
          !
          lw(i) = lw(i) - norm/freqtot * (ctm_C + ctm_X) * V3sq(i,j,k)
          !ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_linewidth_modes = lw
    !
  END FUNCTION sum_linewidth_modes
  ! Sum the Imaginary part of the self energy for the phonon modes,
  ! NOTE: test function! It computes the full self-energy then it takes the imag part
  !       the interesting part is the conversion between the sigma of a Gaussian
  !       and the width of a Lorentzian
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_linewidth_modes2(S, sigma, freq, bose, V3sq, nu0)
    USE input_fc,           ONLY : ph_system_info
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(S%nat3,3)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3)
    INTEGER,INTENT(in)  :: nu0(3)
    !
    ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
    ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
    ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
    ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
    !  =>  g = 2 sqrt(log(2) c = 0.6 c
    REAL(DP),PARAMETER :: csig =  (2 * DSQRT(DLOG(2._dp)))
    REAL(DP) :: tsigma
    !
    REAL(DP) :: sum_linewidth_modes2(S%nat3)
    COMPLEX(DP) :: se(S%nat3)
    !
    tsigma = csig*sigma
    se = sum_selfnrg_modes(S, tsigma, freq, bose, V3sq, nu0)
    sum_linewidth_modes2  = -DIMAG(se)
    !
  END FUNCTION sum_linewidth_modes2
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
