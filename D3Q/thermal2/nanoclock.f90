!
! Copyright Lorenzo Paulatto 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
!
MODULE nanoclock
  !
  USE kinds, ONLY : DP
  USE iso_c_binding, ONLY : c_double
  !
!   TYPE(nanoclock) :: timer
!   timer%name = "MAIN"
!   CALL start_nanoclock(timer)
!   CALL stop_nanoclock(timer)
!   CALL print_nanoclock(timer)
  !
  TYPE nanotimer
    REAL(kind=c_double) :: t0 = -1._c_double
    REAL(kind=c_double) :: tot = 0._c_double
    CHARACTER(len=32)   :: name = "unknown"
  END TYPE nanotimer
  !
  INTERFACE
    FUNCTION c_nanosec() BIND(C,name="c_nanosec")
      USE iso_c_binding, ONLY : c_double
      IMPLICIT NONE
      REAL(kind=c_double) :: c_nanosec
    END FUNCTION
  END INTERFACE
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  !
  FUNCTION f_nanosec()
    IMPLICIT NONE
    REAL(DP) :: f_nanosec
    f_nanosec = REAL(c_nanosec(), kind=DP)
    RETURN
  END FUNCTION f_nanosec
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE start_nanoclock(timer)
    IMPLICIT NONE
    TYPE(nanotimer) :: timer
    IF( .not. timer%t0 < 0) CALL errore("start_nanoclock", &
                                        TRIM(timer%name)//" clock is already running", 1)
    timer%t0 = c_nanosec()
  END SUBROUTINE start_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE stop_nanoclock(timer)
    IMPLICIT NONE
    TYPE(nanotimer) :: timer
    IF( timer%t0 < 0) CALL errore("stop_nanoclock", &
                                  TRIM(timer%name)//" clock is not running", 1)
    timer%tot = timer%tot + (c_nanosec()-timer%t0)
    timer%t0 = -1._c_double
  END SUBROUTINE stop_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE print_nanoclock(timer)
    IMPLICIT NONE
    TYPE(nanotimer) :: timer
    !
    IF(timer%t0>0) THEN
      ! print a running clock
      PRINT*, TRIM(timer%name), timer%tot + (c_nanosec()-timer%t0), "(r)"
    ELSE
      PRINT*, TRIM(timer%name), timer%tot
    ENDIF
    
  END SUBROUTINE print_nanoclock
  ! \/o\________\\\_________________________________________/^>
  
END MODULE

