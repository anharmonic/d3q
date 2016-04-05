!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!

! <<^V^\\=========================================//-//-//========//O\\//
MODULE functions
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  IMPLICIT NONE
  
  REAL(DP),PARAMETER :: one_over_sqrt_2_pi = 1._dp / SQRT( 2*pi)
  REAL(DP),PARAMETER :: one_over_sqrt_pi = 1._dp / SQRT(pi)
  REAL(DP),PARAMETER :: one_over_pi = 1._dp / pi

  CONTAINS
  !  1. NORMALISED GAUSSIAN DISTRIBUTION
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_ngauss(x,s)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,s
    REAL(DP) :: sm1
    sm1 = 1/s
    f_ngauss = (one_over_sqrt_2_pi*sm1) * EXP( -0.5_dp*(sm1*x)**2)
  END FUNCTION
  !
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_ngaussi(x,sm1)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,sm1
    f_ngaussi = (one_over_sqrt_2_pi*sm1) * EXP( -0.5_dp*(sm1*x)**2)
  END FUNCTION
  !
  !  1b. GAUSSIAN DISTRIBUTION
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_gauss(x,s)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,s
    REAL(DP) :: sm1
    sm1 = 1/s
    f_gauss = (one_over_sqrt_pi*sm1) * EXP( -(sm1*x)**2 )
  END FUNCTION
  !
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_gaussi(x,sm1)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,sm1
    f_gaussi = (one_over_sqrt_pi*sm1) * EXP( -(sm1*x)**2 )
  END FUNCTION
  !
  !  2. CAUCHY-LORENTZ DISTRIBUTIONS
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_lorentz(x,g)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,g
    REAL(DP) :: gm1
    gm1 = 1/g
    f_lorentz = one_over_pi*gm1 / (1 + (x*gm1)**2)
  END FUNCTION
  !
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_lorentzi(x,gm1)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,gm1
    f_lorentzi = one_over_pi*gm1 / (1 + (x*gm1)**2)
  END FUNCTION
  !
  !  3. BOSE-EINSTEIN DISTRIBUTIONS
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_bose(x,T)
    USE constants, ONLY : K_BOLTZMANN_RY
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,T
    REAL(DP) :: Tm1
    Tm1 = 1/(T*K_BOLTZMANN_RY)
    f_bose = 1 / (EXP(x*Tm1) - 1)
  END FUNCTION
  !
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_bosei(x,Tm1)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,Tm1
    f_bosei = 1 / (EXP(x*Tm1) - 1)
  END FUNCTION
  !
  !  3b. Derivative of BOSE-EINSTEIN DISTRIBUTIONS
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION df_bose(x,T)
    USE constants, ONLY : K_BOLTZMANN_RY
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,T
    REAL(DP) :: Tm1, expf
    Tm1  = 1/(T*K_BOLTZMANN_RY)
    expf = EXP(x*Tm1)
    df_bose = -Tm1 * expf / (expf-1)**2
  END FUNCTION
  !
  ! Find the periodic copy of xq that's inside the Brillouin zone of lattice bg
  ! Does not keep in account possible degeneracies (i.e. xq on the BZ border)
  FUNCTION refold_bz(xq, bg)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: xq(3), bg(3,3)
    REAL(DP) :: refold_bz(3)
    !
    REAL(DP) :: xp(3), xpmin(3), xpmod, xpmodmin
    INTEGER  :: i,j,k
    INTEGER,PARAMETER :: far = 2
    !
    xpmin = xq
    xpmodmin = SUM(xq**2)
    DO i = -far,far
    DO j = -far,far
    DO k = -far,far
      xp = xq + i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
      xpmod = SUM(xp**2)
      IF(xpmod < xpmodmin)THEN
        xpmin = xp
        xpmodmin = xpmod
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    refold_bz = xpmin
  END FUNCTION
  ! As above, but returns only the vector's length
  FUNCTION refold_bz_mod(xq, bg)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: xq(3), bg(3,3)
    REAL(DP) :: refold_bz_mod
    !
    REAL(DP) :: xp(3), xpmod, xpmodmin
    INTEGER  :: i,j,k
    INTEGER,PARAMETER :: far = 2
    !
    xpmodmin = SUM(xq**2)
    DO i = -far,far
    DO j = -far,far
    DO k = -far,far
      xp = xq + i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
      xpmod = SUM(xp**2)
      IF(xpmod < xpmodmin)THEN
        xpmodmin = xpmod
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    refold_bz_mod = SQRT(xpmodmin)
  END FUNCTION

  SUBROUTINE invzmat (n, a)
    !-----------------------------------------------------------------------
    ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
    ! if the matrix is dimensioned 3x3, it also computes determinant "da"
    ! matrix "a" is unchanged on output - LAPACK
    !
    USE kinds, ONLY : DP
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    COMPLEX(DP),INTENT(inout) :: a(n,n)
    !
    INTEGER :: info, lda, lwork, ipiv(n), nb
    ! info=0: inversion was successful
    ! lda   : leading dimension (the same as n)
    ! ipiv  : work space for pivoting (assumed of length lwork=n)
    COMPLEX(DP),ALLOCATABLE :: work(:) 
    INTEGER,EXTERNAL :: ILAENV
    ! more work space
    !
    lda = n
    !
    nb = ILAENV( 1, 'ZHEEV', 'U', n, -1, -1, -1 )
    lwork=n*nb
    ALLOCATE(work(lwork))
    !
    CALL ZGETRF(n, n, a, lda, ipiv, info)
    CALL errore('invzmat', 'error in DGETRF', ABS(info) )
    CALL ZGETRI(n, a, lda, ipiv, work, lwork, info)
    CALL errore('invzmat', 'error in DGETRI', ABS(info) )
    !
    DEALLOCATE(work)
    !
    RETURN
  END SUBROUTINE invzmat
  
END MODULE functions
! <<^V^\\=========================================//-//-//========//O\\//


