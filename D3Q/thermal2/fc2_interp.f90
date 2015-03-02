!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE fc2_interpolate
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid

  USE nanoclock

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
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
      COMPLEX(DP),INTENT(out),OPTIONAL :: U(S%nat3,S%nat3)
      !
      COMPLEX(DP) :: U_(S%nat3,S%nat3)
      REAL(DP),PARAMETER :: eps = 1.e-12_dp
      !
      CALL fftinterp_mat2(xq, S, fc2, U_)
      CALL mat2_diag(S, U_, freq)
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
      IF(present(U)) U = CONJG(U_)
      !
  END SUBROUTINE prepare_phq
  !
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq_safe(xq, S, fc2, freq, U)
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
    USE functions, ONLY : f_bose
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
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE ip_cart2pat(d3in, nat3, u1, u2, u3, d3out)
    !   Rotates third derivative of the dynamical basis from cartesian axis
    !   to the basis of the modes. Rotation is not really in place
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! d3 matrix, input: in cartesian basis, output: on the patterns basis
    COMPLEX(DP),INTENT(inout) :: d3in(nat3, nat3, nat3)
    COMPLEX(DP),OPTIONAL,INTENT(out) :: d3out(nat3, nat3, nat3)
    INTEGER,INTENT(in)        :: nat3
    ! patterns (transposed, with respect to what we use in the d3q.x code)
    COMPLEX(DP),INTENT(in)    :: u1(nat3, nat3), u2(nat3, nat3), u3(nat3, nat3) 
    !
    COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
    !
    INTEGER :: a, b, c, i, j, k
    COMPLEX(DP) :: AUX,AUX2
    REAL(DP),PARAMETER :: EPS = 1.e-8_dp
    !
    ALLOCATE(d3tmp(nat3, nat3, nat3))
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
        AUX2 = 0._dp
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
    IF(present(d3out)) THEN
      d3out = d3tmp
    ELSE
      d3in  = d3tmp
    ENDIF
    DEALLOCATE(d3tmp)
    !
    RETURN
  END SUBROUTINE ip_cart2pat
  !
END MODULE fc2_interpolate


