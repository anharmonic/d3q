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
    CHARACTER(len=16)   :: name = "unknown"
    !
    REAL(kind=c_double) :: t0 = -1._c_double
    REAL(kind=c_double) :: tot = 0._c_double
    INTEGER :: calls = 0
  END TYPE nanotimer
  !
  TYPE(nanotimer) :: c2pat = nanotimer("c2pat")
  TYPE(nanotimer) :: tv3sq = nanotimer("v3sq")
  TYPE(nanotimer) :: i_ph  = nanotimer("i_ph")
  TYPE(nanotimer) :: lwtot = nanotimer("LW_2tot")
  TYPE(nanotimer) :: tsum = nanotimer("sum_bands")
  TYPE(nanotimer) :: d3time= nanotimer("D3")
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
    TYPE(nanotimer),INTENT(inout) :: timer
    IF( .not. timer%t0 < 0) CALL errore("start_nanoclock", &
                                        TRIM(timer%name)//" clock is already running", 1)
    timer%t0 = c_nanosec()
    timer%calls = timer%calls +1
  END SUBROUTINE start_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE stop_nanoclock(timer)
    IMPLICIT NONE
    TYPE(nanotimer),INTENT(inout) :: timer
    IF( timer%t0 < 0) CALL errore("stop_nanoclock", &
                                  TRIM(timer%name)//" clock is not running", 1)
    timer%tot = timer%tot + (c_nanosec()-timer%t0)
    timer%t0 = -1._c_double
  END SUBROUTINE stop_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE print_line()
    IMPLICIT NONE
    WRITE(*,'(2x," * ",16("*")," * ",12("*"),"***** * ", 12("*")," *")') 
  END SUBROUTINE print_line
  SUBROUTINE print_header()
    IMPLICIT NONE
    WRITE(*,'(2x," * ",a16," * ",5x,a12," * ", a6," *")') &
    "TIMER NAME", " TIME (s) ", "CALLS"
  END SUBROUTINE print_header
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_nanoclock(timer)
    IMPLICIT NONE
    TYPE(nanotimer),INTENT(in),OPTIONAL :: timer
    !
    IF(.not. present(timer)) THEN
      CALL print_header()
      CALL print_line()
      RETURN
    ENDIF
    
    IF(timer%t0>0) THEN
      ! print a running clock
      WRITE(*,'(2x," * ",a16," * ",f12.6,"s (r) * ", i12," *")') &
      TRIM(timer%name), timer%tot + (c_nanosec()-timer%t0), "(r)", timer%calls
    ELSE
      WRITE(*,'(2x," * ",a16," * ",f12.6,"s     * ", i12," *")') &
      TRIM(timer%name), timer%tot, timer%calls
    ENDIF
    
  END SUBROUTINE print_nanoclock
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_memory()
    USE iso_c_binding,  ONLY : c_int
    IMPLICIT NONE
    INTEGER(kind=c_int) :: kb
    CHARACTER(len=2)    :: unit
    !
    CALL memstat(kb)
    unit = "kB"
    IF(kb>10000) THEN
      kb = kb/1000
      unit = "MB"
    ENDIF
    WRITE(*,'(2x," * ",6x,a16," : ",i8,a2,16x," *")') "Memory used ", kb, unit
    !
  END SUBROUTINE print_memory

  
END MODULE

