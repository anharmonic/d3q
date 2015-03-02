!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
!
! The nanotimer type describes a nanonsecond-grained timer, it can be started,
! stopped and retrieved via start_nanosecond and stop_nanosecond. Item %tot contains
! the total elapsed time, %calls the number of calls and %t0 the starting time
! (for running clocks)
! f_nanosec returns a system nanoclocs timer, relative to the first call of 
! the c_nanotimer functions, it is not used here (we always use the C routine directly).
! 
MODULE nanoclock
  !
  USE kinds,            ONLY : DP
  USE iso_c_binding,    ONLY : c_double
  IMPLICIT NONE
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
    
    CONTAINS
      procedure :: read  => read_nanoclock
      procedure :: print => print_nanoclock
      procedure :: start => start_nanoclock
      procedure :: stop  => stop_nanoclock
  END TYPE nanotimer
  !
  ! You can define timers like this:
!   TYPE(nanotimer) :: my_timer = nanotimer("Timer description")
!   ...
!   CALL my_timer%start()
!   ...
!   CALL my_timer%stop()
!   CALL my_timer%print()
  !
  ! c_nanosec returns the nanoseconds relative to the first time it was called
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
    CLASS(nanotimer),INTENT(inout) :: timer
    IF( .not. timer%t0 < 0) CALL errore("start_nanoclock", &
                                        TRIM(timer%name)//" clock is already running", 1)
    timer%t0 = c_nanosec()
    timer%calls = timer%calls +1
  END SUBROUTINE start_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE stop_nanoclock(timer)
    IMPLICIT NONE
    CLASS(nanotimer),INTENT(inout) :: timer
    IF( timer%t0 < 0) CALL errore("stop_nanoclock", &
                                  TRIM(timer%name)//" clock is not running", 1)
    timer%tot = timer%tot + (c_nanosec()-timer%t0)
    timer%t0 = -1._c_double
  END SUBROUTINE stop_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  REAL(DP) FUNCTION read_nanoclock(timer)
    IMPLICIT NONE
    CLASS(nanotimer),INTENT(in) :: timer
    !
    IF(timer%t0>0) THEN
      read_nanoclock = REAL(timer%tot + (c_nanosec()-timer%t0), kind=DP)
    ELSE
      read_nanoclock = REAL(timer%tot, kind=DP)
    ENDIF
    
  END FUNCTION read_nanoclock
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE print_nanoclock(timer)
    IMPLICIT NONE
    CLASS(nanotimer),INTENT(in) :: timer
    !
    IF(timer%calls<=0) RETURN
    !
    IF(timer%t0>0) THEN
      ! print a running clock
      WRITE(*,'(2x," * ",a16," * ",f15.6," s     * ",f15.6," ms/call * ", f8.3, " % wtime * ", i12," calls * RUNNING!")') &
      TRIM(timer%name), timer%tot+ (c_nanosec()-timer%t0), (1000*timer%tot)/timer%calls, &
      100*(timer%tot+ (c_nanosec()-timer%t0))/c_nanosec(), timer%calls
    ELSE
      WRITE(*,'(2x," * ",a16," * ",f15.6," s     * ",f15.6," ms/call * ", f8.3, " % wtime * ", i12," calls *")') &
      TRIM(timer%name), timer%tot, (1000*timer%tot)/timer%calls, &
      100*timer%tot/c_nanosec(), timer%calls
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
    WRITE(*,'(2x," * ",6x,a16," : ",i8,a2,27x," *")') "Memory used ", kb, unit
    !
  END SUBROUTINE print_memory
  ! \/o\________\\\_________________________________________/^>
  
END MODULE

