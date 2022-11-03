!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!

MODULE more_constants
  USE kinds, ONLY : DP
#include "mpi_thermal.h"
  REAL(DP),PARAMETER :: RY_TO_JOULE =  0.5* 4.35974394e-18
  REAL(DP),PARAMETER :: RY_TO_SECOND = 2* 2.418884326505e-17
  REAL(DP),PARAMETER :: RY_TO_METER = 5.2917721092e-11
  CHARACTER(len=3),PARAMETER :: INVALID = '///'
  REAL(DP),PARAMETER :: DHUGE = -1.E+30_dp
  REAL(DP),PARAMETER :: MASS_DALTON_TO_RY = 0.5_dp*1822.88839_dp
  REAL (DP), PARAMETER :: RY_TO_WATT = RY_TO_JOULE / RY_TO_SECOND
  REAL(DP),PARAMETER :: RY_TO_WATTMM1KM1 = RY_TO_WATT / RY_TO_METER
  ! added
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_DP ! m
  REAL(DP), PARAMETER :: HARTREE_SI       = 4.35974394E-18_DP   ! J
  REAL(DP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_DP   ! J
  REAL(DP), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_DP   ! J s, this is the Planck constant h, not \hbar
  REAL(DP), PARAMETER :: ryau_sec = H_PLANCK_SI/RYDBERG_SI
  REAL(DP), PARAMETER :: ryvel_si = BOHR_RADIUS_SI/ryau_sec 

  !
  REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
  REAL(DP),PARAMETER :: eps_freq = 0._dp !1.e-8_dp
  COMPLEX(DP), PARAMETER :: complex_i = cmplx(0.0,1.0)
  CONTAINS
  !
!   CHARACTER(len=256) &
!   FUNCTION sigma_file_name(prefix, nq1, nq2, nq3, cT, csigma)
!     IMPLICIT NONE
!     CHARACTER(len=256),INTENT(in) :: prefix
!     INTEGER,INTENT(in) :: nq1, nq2, nq3
!     CHARACTER(len=6),INTENT(in) ::  cT, csigma
!     !
!     CHARACTER (LEN=6), EXTERNAL :: int_to_char
!     !
!     sigma_file_name= TRIM(prefix)//&
!                 "."//TRIM(int_to_char(nq1))// &
!                 "."//TRIM(int_to_char(nq2))// &
!                 "."//TRIM(int_to_char(nq3))// &
!                 "."//TRIM(cT)// &
!                 "."//TRIM(csigma)//".out"
!     !
!   END FUNCTION sigma_file_name
!   !
  ! Write a number from a list using as many digits after the dot as the longest in the list.
  CHARACTER(len=6) &
  FUNCTION write_conf(it,nt,T) RESULT(str)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: it, nt
    REAL(DP),INTENT(in) :: T(nt)
    
    INTEGER :: max_digit_left, max_digit_right, jt, ls
    REAL(DP) :: Tx, Tfrac
    CHARACTER(len=64):: fmt1='', fmt2=''
    
    max_digit_right=0 
    max_digit_left=0 
   
    DO jt = 1,nt      
      max_digit_left = MAX( max_digit_left , CEILING(LOG10(T(jt))) )
      IF(T(jt)>0) THEN
        Tfrac = T(jt)-INT(T(jt))
        IF(Tfrac>0._dp) THEN
          Tx = 1._dp/Tfrac
          max_digit_right = MAX( max_digit_right , CEILING(LOG10(Tx)) )
        !rint*, ">>", jt, T(jt), Tx,  CEILING(LOG10(Tx))
        ENDIF
      ENDIF
    ENDDO
    
    str=""
    WRITE(fmt1,'(i6)') max_digit_right
    fmt2 = "(1f6."//TRIM(ADJUSTL(fmt1))//")"
    !print*, fmt1, fmt2, max_digit_left, max_digit_right
   
    IF(max_digit_right>0)THEN 
      WRITE(str,fmt2) T(it)
    ELSE
      WRITE(str,'(i6)') INT(T(it))
    ENDIF
    str=TRIM(ADJUSTL(str))
    
  END FUNCTION write_conf
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_citations_linewidth
    ioWRITE(stdout,*)
    ioWRITE(stdout,'(5x,a)') &
        " ",&
        "For 2n+1 calculations of force constants please cite:",&
        " 1. Lorenzo Paulatto, Francesco Mauri, and Michele Lazzeri",&
        "    Phys. Rev. B 87, 214303 (2013)"
    ioWRITE(stdout,'(5x,a)') &
        " ",&
        "For thermal transport calculations please cite:",&
        " 2. Giorgia Fugallo, Michele Lazzeri, Lorenzo Paulatto, and Francesco Mauri",&
        "    Phys. Rev. B 88, 045430 (2013)", &
        " 3. A. Cepellotti, G. Fugallo, L. Paulatto, M. Lazzeri, F. Mauri, N. Marzari,",&
        "    Nature communications 6 (2015)",&
        " 4. G. Fugallo, A. Cepellotti, L. Paulatto, M. Lazzeri, N. Marzari, F. Mauri,",&
        "    Nano letters 14 (11), 6109-6114 (2014)",&  
        " 5. M. Simoncelli, N. Marzari, F. Mauri, Unified theory of thermal transport in crystals and glasses,",&
        "    Nature Physics 15 (8), 809-813 (2019)",&    
        " 6. M. Simoncelli, N. Marzari, F. Mauri,",&
        "    Wigner formulation of thermal transport in solids, Phys. Rev. X 12, 041011 (2022)"
    ioWRITE(stdout,'(5x,a)') &
        " ",&
        "For spectral function calculations also cite:",&
        " 7. Lorenzo Paulatto, Ion Errea, Matteo Calandra, and Francesco Mauri,",&
        "    Phys. Rev. B 91, 054304 (2015)"
    ioWRITE(stdout,*)

  END SUBROUTINE print_citations_linewidth
  !
END MODULE more_constants
!
