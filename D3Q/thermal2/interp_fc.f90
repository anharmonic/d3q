!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE interp_fc
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid, forceconst3_grid

  USE nanoclock

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat2(xq, S, fc, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,arg,phase) REDUCTION(+: D)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:, :) = D(:, :) + phase * fc%fc(:, :, i)
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE fftinterp_mat2
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat3(xq2,xq3, S, fc, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst3_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq2(3), xq3(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3, S%nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,arg,phase) REDUCTION(+: D)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:,:,:) = D(:,:,:) + phase * fc%fc(:,:,:, i)
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE fftinterp_mat3
  ! \/o\________\\\_________________________________________/^>
  ! IN PLACE diagonalization of D
  SUBROUTINE mat2_diag(S, D, w2)
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    COMPLEX(DP),INTENT(inout) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(out)      :: w2(S%nat)
    !
    INTEGER  :: nb    ! block size
    INTEGER,save :: nat3=-1, lwork=-1 ! aux. var

    INTEGER :: info      ! flag saying if the exec. of libr. routines was ok

    INTEGER,EXTERNAL ::ILAENV ! function which gives block size
    !
    REAL(DP), ALLOCATABLE    :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    IF ( nat3 /= S%nat3 .or. lwork < 0 ) THEN
      !     check for the block size
      nb = ILAENV( 1, 'ZHETRD', 'U', S%nat3, -1, -1, -1 )
      IF (nb<1) nb = MAX(1,S%nat3)
      IF (nb==1 .or. nb>=S%nat3) then
        lwork=2*S%nat3-1
      ELSE
        lwork = (nb+1)*S%nat3
      ENDIF
      !
      IF(nat3 > 0 ) PRINT*, "WARNING! Redoing ILAENV"
      nat3 = S%nat3
      !
    ENDIF
    !
    ALLOCATE(work (lwork))
    ALLOCATE(rwork (3*S%nat3-2))
    !
    CALL ZHEEV('V','U',S%nat3,D,S%nat3,w2,work,lwork,rwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(rwork)
    DEALLOCATE(work)
    !
  END SUBROUTINE mat2_diag
  ! \/o\________\\\_________________________________________/^>
  ! Returns the cross section of a 3-phonon scattering ;
  ! this is | V^3(q1,q2,q3) |^2 in Fugallo et. al. PRB
  SUBROUTINE scatter_3q(S, fc2, fc3, xq1, xq2, xq3, V3sq)
    USE d3_basis,  ONLY : d3_cart2pat
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in)  :: xq1(3), xq2(3), xq3(3)
    REAL(DP),INTENT(out) :: V3sq(:,:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:)
    REAL(DP) :: xq(3,3)
    INTEGER :: i
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(w2(S%nat3))
    !
    xq(:,1) = xq1 ; xq(:,2) = xq2; xq(:,3) = xq3
    DO i = 1,3
      CALL fftinterp_mat2(xq(:,i), S, fc2, U(:,:,i))
      CALL mat2_diag(S, U(:,:,i), w2)
    ENDDO
    !
    DEALLOCATE(w2)
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
    !
    CALL d3_cart2pat(D3, S%nat, U(:,:,1), U(:,:,2), U(:,:,3))
    V3sq = REAL( CONJG(D3)*D3 , kind=DP)
    !
    DEALLOCATE(U, D3)
    !
  END SUBROUTINE scatter_3q
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE ip_cart2pat(d3in, nat3, u1, u2, u3)
    !   Rotates third derivative of the dynamical basis from cartesian axis
    !   to the basis of the modes. Rotation is not really in place
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! d3 matrix, input: in cartesian basis, output: on the patterns basis
    COMPLEX(DP),INTENT(inout) :: d3in(nat3, nat3, nat3)
    INTEGER,INTENT(in)        :: nat3
    ! patterns (transposed, with respect to what we use in the d3q.x code)
    COMPLEX(DP),INTENT(in)    :: u1(nat3, nat3), u2(nat3, nat3), u3(nat3, nat3) 
    !
    !COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
    COMPLEX(DP)  :: d3tmp(nat3,nat3,nat3)

    INTEGER :: a, b, c, i, j, k
    COMPLEX(DP) :: AUX
    REAL(DP),PARAMETER :: EPS = 1.e-8_dp
    !
    !ALLOCATE(d3tmp(nat3, nat3, nat3))
    d3tmp = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,a,b,c) REDUCTION(+: d3tmp) COLLAPSE(2)
    DO k = 1,nat3
    DO c = 1,nat3
    IF(ABS(u3(c,k))>EPS)THEN
      !
     DO j = 1,nat3
      DO b = 1,nat3
      IF(ABS(u2(b,j))>EPS)THEN
        !
        AUX = u2(b,j) * u3(c,k)
        DO i = 1,nat3
        DO a = 1,nat3
            d3tmp(i, j, k) = d3tmp(i, j, k) &
                            + u1(a,i) * AUX * d3in(a, b, c) 
        ENDDO
        ENDDO
        !
      ENDIF
      ENDDO
      ENDDO
      !
    ENDIF
    ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    d3in  = d3tmp
    !DEALLOCATE(d3tmp)
    !
    RETURN
  END SUBROUTINE ip_cart2pat
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
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,freqtot,bose_X,bose_C,dom_X,dom_C,ctm_X,ctm_C) REDUCTION(+: lw) COLLAPSE(3)
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
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! This function returns a complex number:
  !  its real part is the "Tadpole" 3rd order contribution to lineshift
  !  its complex part is the linewidth
  FUNCTION lineshift_q(xq0, T, sigma, S, grid, fc2, fc3) &
  RESULT(ls)
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
    COMPLEX(DP) :: ls(S%nat3)
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
    ls = (0._dp, 0._dp)
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL prepare_phq(xq(:,1), T, S, fc2, freq(:,1), bose(:,1), U(:,:,1))
    !
!     CALL start_nanoclock(lstot)
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
      ls = ls + sum_complex_modes( S, sigma, freq, bose, V3sq )
!       CALL stop_nanoclock(tsum)
      !
    ENDDO
    !
    ls = ls/grid%nq
    !
!     CALL stop_nanoclock(lstot)

!     CALL print_nanoclock(lstot)
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
  END FUNCTION lineshift_q  
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
    COMPLEX(DP) :: sum_complex_modes(S%nat3)
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    COMPLEX(DP) :: ctm_P, ctm_M, reg, num   !
    !
    INTEGER :: i,j,k
    COMPLEX(DP) :: ls(S%nat3)
    COMPLEX(DP),PARAMETER :: ii = (0._dp, 1._dp)
    ls(:) = (0._dp, 0._dp)
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,ctm_P,ctm_M,reg) &
!$OMP             REDUCTION(+: ls) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        bose_P = 1 + bose(j,2) + bose(k,3)
        omega_P  =(freq(j,2)+freq(k,3))
        omega_P2 = omega_P**2
        !
        bose_M = bose(k,3) - bose(j,2)
        omega_M =(freq(j,2)-freq(k,3))
        omega_M2 = omega_M**2
        !
        DO i = 1,S%nat3
          !
          ! regularization:
          reg = CMPLX(freq(i,1), delta, kind=DP)
          reg = reg**2
          !
!           num = exp(-log(omega_P2-reg))
!           ctm_P = 2 * bose_P *omega_P * num
          ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
!           num = exp(-log(omega_M2-reg))
!           ctm_M = 2 * bose_M *omega_M * num
          ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
          !
          ls(i) = ls(i) + (ctm_P + ctm_M) * V3sq(i,j,k)
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
  SUBROUTINE prepare_phq(xq, T, S, fc2, freq, bose, U)
!     USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag
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
      ! Is the following mess really necessary?
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
END MODULE interp_fc


