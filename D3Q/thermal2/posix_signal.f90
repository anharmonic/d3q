!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE posix_signal
! This module is a Fortran 2003 interface to the customize_signals.c C file
! Compatible with Intel/PGI/Gcc(>=4.3) compilers
#include "mpi_thermal.h"
  ! This module is compiled only if the following preprocessing option
  ! is enabled
#ifdef __TERMINATE_GRACEFULLY
  USE iso_c_binding

  IMPLICIT NONE

  INTEGER,VOLATILE::signal_trapped=0
  LOGICAL,VOLATILE::sigint_trapped=.false.
  INTEGER(kind=c_int),PARAMETER :: SIGINT = 2_c_int

  ! Interface to the c function is clib/signal_handler.c
  INTERFACE
    FUNCTION init_TERMINATE_GRACEFULLY(new_handler) BIND(c, name = "init_TERMINATE_GRACEFULLY")
      USE iso_c_binding
      TYPE(C_FUNPTR),VALUE,INTENT(IN):: new_handler
      INTEGER(C_INT)::init_TERMINATE_GRACEFULLY
    END FUNCTION init_TERMINATE_GRACEFULLY
  END INTERFACE
#endif

  CONTAINS

  SUBROUTINE set_TERMINATE_GRACEFULLY()
    USE iso_c_binding
    TYPE(C_FUNPTR),TARGET::ptr
!    INTERFACE
!       SUBROUTINE routine(signal) bind(C)
!         USE iso_c_binding
!         INTEGER(C_INT),VALUE, INTENT(IN)::signal
!       END SUBROUTINE routine
!    END INTERFACE

#ifdef __TERMINATE_GRACEFULLY
    ptr = C_FUNLOC(set_graceful_termination)

    IF (init_TERMINATE_GRACEFULLY(ptr) .NE. 0) THEN
       CALL errore("set_TERMINATE_GRACEFULLY", "The association of signals INT or TERM failed!", 1)
    ENDIF
#endif

  END SUBROUTINE set_TERMINATE_GRACEFULLY
  !
  ! Sets the signal_trapped flag on all nodes/processors
  ! Only the master will use the signal, though
#ifdef __TERMINATE_GRACEFULLY
  SUBROUTINE set_graceful_termination(signum) BIND(c)
    USE iso_c_binding
    USE mpi_thermal, ONLY : abort_mpi
    INTEGER(C_INT),VALUE,INTENT(IN):: signum
    ! Double CTRL-C will stop immediately;
    ! Some implementation of MPI send SIGTERM to every process when
    ! SIGINT (aka CTRL-C) is received which causes the processes
    ! to receive both SIGINT AND SIGTERM. We track SIGINT separately and
    ! only stop immediately if we receive SIGINT twice
    IF(sigint_trapped.and.signum==SIGINT) THEN
      WRITE(stderr, '(5x,a,i6)') "**** SIGINT ALREADY TRAPPED ONCE: terminating immediately!!", signum
      CALL ptrace()
      CALL abort_mpi(signum)
      STOP 1
    ELSE
      WRITE(stderr, '(5x,a,i6)') "**** Trapped signal: trying to terminate gracefully", signum
      signal_trapped = INT(signum)
      !
      ! if the signal was SIGNIT (aka CTRL-C) we track it separately
      IF(signum==SIGINT) THEN
        WRITE(stderr, '(5x,a)') "**** press CTRL-C again to terminate immediately (no restart possible!)"
        sigint_trapped = .true.
      ENDIF
      !
    ENDIF
    !
  END SUBROUTINE set_graceful_termination
#endif
  !
  ! This subroutine must be called by all the mpi processes
  SUBROUTINE check_graceful_termination()
    USE mpi_thermal, ONLY : mpi_bsum, stop_mpi, mpi_any, abort_mpi
    USE timers,      ONLY : print_all_timers, check_time_limit
    IMPLICIT NONE
    LOGICAL :: time_is_out
    !
#ifdef __TERMINATE_GRACEFULLY
    CALL mpi_bsum(signal_trapped)
    ! time limit may be a few second out of sync over mpi, let's
    ! be sure that everyone stops
    time_is_out = check_time_limit()
!    print*, "time_is_out", time_is_out
    CALL mpi_any(time_is_out)
!    print*, "time_is_out 2", time_is_out
    !
    IF(signal_trapped>0 .or. time_is_out) THEN
      CALL print_all_timers()
      !CALL stop_mpi()
      CALL abort_mpi(signal_trapped)
      STOP 0
    ENDIF
#endif
  END SUBROUTINE check_graceful_termination


END MODULE posix_signal


