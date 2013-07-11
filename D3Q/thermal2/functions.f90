!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE functions
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  
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
    f_ngauss = (one_over_sqrt_2_pi*sm1) * EXP( -0.5_dp*sm1*(x**2) ) 
  END FUNCTION
  !
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  FUNCTION f_ngaussi(x,sm1)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: x,sm1
    f_ngaussi = (one_over_sqrt_2_pi*sm1) * EXP( -0.5_dp*sm1*(x**2) ) 
  END FUNCTION
  REAL(DP) ELEMENTAL & ! <`\.......''..','
  !
  !  1b. GAUSSIAN DISTRIBUTION
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
END MODULE functions
! <<^V^\\=========================================//-//-//========//O\\//


