!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Implementaton of quicksort from https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Author: t-nissie@github.com
! Licence GPLv3
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE functions
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  IMPLICIT NONE

  INTERFACE default_if_not_present
    MODULE PROCEDURE default_if_not_present_int
    MODULE PROCEDURE default_if_not_present_logical
  END INTERFACE


  REAL(DP),PARAMETER :: one_over_sqrt_2_pi = 1._dp / DSQRT( 2*pi)
  REAL(DP),PARAMETER :: one_over_sqrt_pi = 1._dp / DSQRT(pi)
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
  ! Pseudo-voigt is a linear combination of a Gaussian and a Lorentian with the
  ! same width, the sigma of the Gussian has to be adapted
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_psvoigt(x,f,eta)
    IMPLICIT NONE
    REAL(DP), PARAMETER :: sigma_to_fwhm = 2*DSQRT(2*DLOG(2._dp))
    REAL(DP), PARAMETER :: gamma_to_fwhm = 2
    REAL(DP),INTENT(in) :: x,f,eta
    REAL(DP) :: sigmam1, gammam1
    sigmam1 = sigma_to_fwhm/f
    gammam1 = gamma_to_fwhm/f
    !
    f_psvoigt = eta * f_lorentzi(x,gammam1) + (1-eta) * f_ngaussi(x,sigmam1)
    !
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
  ! Returns the average length of the thermal oscillation for
  ! and harmonic mode at temperature T and energy x
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_wtoa(x,T)
    USE constants, ONLY : K_BOLTZMANN_RY
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,T
    REAL(DP) :: Tm1
    Tm1 = 1/(T*K_BOLTZMANN_RY)
    f_wtoa = DSQRT(0.5_dp /x /DTANH(0.5_dp * x*Tm1))
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
    refold_bz_mod = DSQRT(xpmodmin)
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

  FUNCTION cross(a,b)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: a(3), b(3)
    REAL(DP) :: cross(3)
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION

  REAL(DP) FUNCTION norm(a)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: a(3)
    norm = DSQRT(SUM(a**2))
  END FUNCTION

  
! Bubble sort is too slow, don't use it
!   SUBROUTINE Bubble_Sort(a)
!   REAL(DP), INTENT(inout), DIMENSION(:) :: a
!     REAL(DP) :: temp
!     INTEGER :: i, j
!     LOGICAL :: swapped
!   
!     DO j = SIZE(a)-1, 1, -1
!       swapped = .FALSE.
!       DO i = 1, j
!         IF (a(i) > a(i+1)) THEN
!          temp = a(i)
!          a(i) = a(i+1)
!          a(i+1) = temp
!          swapped = .TRUE.
!         END IF
!       END DO
!       IF (.NOT. swapped) EXIT
!     END DO
!   END SUBROUTINE  Bubble_Sort
! 
!   SUBROUTINE Bubble_Sort_idx(a, idx)
!   REAL(DP), INTENT(inout), DIMENSION(:) :: a
!   INTEGER, INTENT(out), DIMENSION(:) :: idx
!     REAL(DP) :: temp
!     INTEGER :: i, j, itmp
!     LOGICAL :: swapped
! 
!     IF(size(idx)<size(a)) CALL errore("bubble_idx","not enough room for index",1)
!     FORALL(i=1:size(a)) idx(i) = i
! 
!     DO j = SIZE(a)-1, 1, -1
!       swapped = .FALSE.
!       DO i = 1, j
!         IF (a(i) > a(i+1)) THEN
!          temp = a(i)
!          a(i) = a(i+1)
!          a(i+1) = temp
!          swapped = .TRUE.
!  
!          itmp = idx(i)
!          idx(i) = idx(i+1)
!          idx(i+1) = itmp
!         END IF
!       END DO
!       IF (.NOT. swapped) EXIT
!     END DO
!   END SUBROUTINE  Bubble_Sort_idx

! Implementaton of quicksort from https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Author: t-nissie@github.com
! Licence GPLv3
recursive subroutine quicksort(a, first, last)
  implicit none
  REAL(dp),INTENT(inout) :: a(*)
  INTEGER,INTENT(in) :: first, last
  REAL(dp) :: x, t
  INTEGER :: i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

recursive subroutine quicksort_idx(a, idx, first, last)
  implicit none
  REAL(dp),INTENT(inout) :: a(*)
  INTEGER,INTENT(inout) :: idx(*) 
  INTEGER,INTENT(in) :: first, last
  REAL(DP) :: x, t
  INTEGER :: i, j, n
  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     n = idx(i);  idx(i) = idx(j);  idx(j) = n
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort_idx(a, idx, first, i-1)
  if (j+1 < last)  call quicksort_idx(a, idx, j+1, last)
end subroutine quicksort_idx

  INTEGER FUNCTION default_if_not_present_int(deft, arg) &
          RESULT(default_if_not_present)
    INTEGER,INTENT(in) :: deft
    INTEGER,OPTIONAL,INTENT(in) :: arg
    IF(present(arg)) THEN
      default_if_not_present = arg
    ELSE
      default_if_not_present = deft
    ENDIF
    
  END FUNCTION
  LOGICAL FUNCTION default_if_not_present_logical(deft, arg) &
          RESULT(default_if_not_present)
    LOGICAL,INTENT(in) :: deft
    LOGICAL,OPTIONAL,INTENT(in) :: arg
    IF(present(arg)) THEN
      default_if_not_present = arg
    ELSE
      default_if_not_present = deft
    ENDIF
    
  END FUNCTION
  
  REAL(DP) FUNCTION sigma_mgo(w, T)
    !USE constants, ONLY : RY_TO_CMM1 ==> 109737.31570111268
    IMPLICIT NONE
    REAL(DP), INTENT(in) :: w, T
!     REAL(DP) :: A = (42.1109/RY_TO_CMM1)**3
!     REAL(DP) :: B = 778.684/RY_TO_CMM1
!    sigma_mgo = A/(B-w)**2
      REAL(DP) :: a = -0.481937_dp 
      REAL(DP) :: b = 0.0115197_dp 
      sigma_mgo = a*w**2 + b*w
  END FUNCTION

!----------------------------------------------------------------------------
  SUBROUTINE cdiag_serial( n, h, ldh, e, v )
    !----------------------------------------------------------------------------
    !! Calculates all the eigenvalues and eigenvectors of a complex
    !! hermitean matrix H. On output, the matrix is unchanged.
    !
    USE kinds,            ONLY : DP
    USE mp_bands,         ONLY : nbgrp, me_bgrp, root_bgrp, intra_bgrp_comm
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    ! ... on INPUT
    !
    INTEGER :: n
    !! Dimension of the matrix to be diagonalized
    INTEGER :: ldh
    !! Leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP) :: h(ldh,n)
    !! Matrix to be diagonalized
    !
    ! ... on OUTPUT
    !
    REAL(DP)    :: e(n)
    !! eigenvalues
    COMPLEX(DP) :: v(ldh,n)
    !! eigenvectors (column-wise)
    !
    ! ... local variables for LAPACK 
    !
    INTEGER                  :: lwork, nb, info
    REAL(DP),    ALLOCATABLE :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    INTEGER, EXTERNAL :: ILAENV
    ! ILAENV returns optimal block size "nb"
    ! ... check for the block size
    !
    nb = ILAENV( 1, 'ZHETRD', 'U', n, - 1, - 1, - 1 )
    IF ( nb < 1 .OR. nb >= n ) THEN
       lwork = 2*n
    ELSE
       lwork = ( nb + 1 )*n
    END IF
    !
    v = h
    !
    ALLOCATE( work( lwork ) )    
    ALLOCATE( rwork( 3 * n - 2 ) )    
    !
    CALL ZHEEV( 'V', 'U', n, v, ldh, e, work, lwork, rwork, info )
    !
    CALL errore( 'cdiagh', 'diagonalization (ZHEEV) failed', ABS( info ) )
    !
    ! ... deallocate workspace
    !
    DEALLOCATE( rwork )
    DEALLOCATE( work )
    !
    RETURN
    !
  END SUBROUTINE cdiag_serial

END MODULE functions
! <<^V^\\=========================================//-//-//========//O\\//


