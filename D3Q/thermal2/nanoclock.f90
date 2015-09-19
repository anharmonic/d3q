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
  USE mpi_thermal,      ONLY : ionode
#include "para_io.h"  
  IMPLICIT NONE
  !
!   TYPE(nanoclock) :: timer
!   timer%name = "MAIN"
!   CALL start_nanoclock(timer)
!   CALL stop_nanoclock(timer)
!   CALL print_nanoclock(timer)
  !
  TYPE nanotimer
    CHARACTER(len=24)   :: name = "unknown"
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
  ! Init function is optional, it is only needed if you want to start the WALL
  ! count before any specific clock
  SUBROUTINE init_nanoclock()
    IMPLICIT NONE
    REAL(DP) :: dummy
    dummy = REAL(c_nanosec(), kind=DP)
    RETURN
  END SUBROUTINE init_nanoclock
  ! Returns the number of nanoseconds passed since the first call to c_nanosec
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
    USE mpi_thermal, ONLY : num_procs, mpi_ipl_sum
    IMPLICIT NONE
    CLASS(nanotimer),INTENT(in) :: timer
    REAL(DP) :: tot
    INTEGER  :: calls
    !
    IF(timer%calls<=0) RETURN
    !
    calls = timer%calls
    CALL mpi_ipl_sum(calls)
    tot = REAL(timer%tot,kind=DP)
    ! CAREFUL! clock could be running on some processor but not on all of them!
    IF(timer%t0>0) tot = tot+ (c_nanosec()-timer%t0) 
    CALL mpi_ipl_sum(tot)
    !
    ioWRITE(*,'(2x," * ",a24," * ",f15.6," * ",f15.6," * ",f15.3," * ",f15.6," * ", f8.3, " * ", i12," *")') &
    TRIM(timer%name), 1000*tot/num_procs, (1000*tot)/(calls*num_procs), 1000*tot, (1000*tot)/calls, &
    100*tot/c_nanosec()/num_procs, calls
    
  END SUBROUTINE print_nanoclock
  
  FUNCTION get_wall()
    REAL(kind=DP) :: get_wall
    get_wall= REAL(c_nanosec(), kind=DP)
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_timers_header()
    IMPLICIT NONE
    !
    ioWRITE(*,'(2x," * ",24x," * ",12x," ms * ",7x," ms/call * ",8x," ms*cpu * ",3x," ms*cpu/call * ", " % wtime * ",6x," calls *")')
    !
  END SUBROUTINE print_timers_header
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

