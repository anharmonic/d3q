!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE interp_fc
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
END MODULE interp_fc


