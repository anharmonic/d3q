!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE timers
  USE nanoclock, ONLY : nanotimer, get_wall, print_timers_header, init_nanoclock
  USE kinds,     ONLY : DP
#include "mpi_thermal.h"

  INTEGER :: limit_seconds = -1

  TYPE(nanotimer) :: t_freq =   nanotimer("ph interp & diag"), &
                     t_bose =   nanotimer("bose distrib"), &
                     t_sum  =   nanotimer("sum modes"), &
                     t_fc3int = nanotimer("fc3 interpolate"), &
                     t_fc3m2  = nanotimer("fc3 modulus sq"), &
                     t_fc3rot = nanotimer("fc3 rotate"), &
                     t_lwphph = nanotimer("lw ph-ph"), &
                     t_lwisot = nanotimer("lw isotopes"), &
                     t_lwcasi = nanotimer("lw casimir"), &
                     t_velcty = nanotimer("ph group velocity"), &
                     t_mpicom = nanotimer("mpi communication"), &
                     t_readdt = nanotimer("read fc data")
!    *      lw isotopes *       48.455432 s     *   272.221529 ms/call *    7.365 % wtime *          178 calls *
!    *       lw casimir *        0.092465 s     *     0.519465 ms/call *    0.014 % wtime *          178 calls *
!    *         lw ph-ph *      608.567777 s     *  3418.920093 ms/call *   92.498 % wtime *          178 calls *
!   * Contributions to ph-ph linewidth time:
!    *      freq & bose *      205.748674 s     *     0.070546 ms/call *   31.272 % wtime *      2916530 calls *
!    *        sum modes *       17.402060 s     *     0.007956 ms/call *    2.645 % wtime *      2187264 calls *
!    *  fc3 interpolate *      361.900451 s     *     0.496374 ms/call *   55.006 % wtime *       729088 calls *
!    *   fc3 modulus sq *        0.967292 s     *     0.001327 ms/call *    0.147 % wtime *       729088 calls *
!    *       fc3 rotate *       17.824916 s     *     0.024448 ms/call *    2.709 % wtime *       729088 calls *

  TYPE(nanotimer) :: t_tksma   = nanotimer("sma thermalk"), &
                     t_tksum   = nanotimer("sum of tk terms"), &
                     t_tkaout  = nanotimer("prepare A_out and prec"), &
                     t_tkain   = nanotimer("prepare and apply A_in"), &
                     t_tkcg    = nanotimer("computer CG step"), &
                     t_xain    = nanotimer("A_in * f (A_in given)"), &
                     t_lwinout = nanotimer("lw input/output")

  CONTAINS
  SUBROUTINE print_all_timers()
    IMPLICIT NONE

    CALL print_timers_header()
    CALL t_tksma%print()
    CALL t_tksum%print()
    CALL t_tkaout%print()
    CALL t_tkain%print()
    CALL t_tkcg%print()
    CALL t_freq%print()
    CALL t_bose%print()
    CALL t_sum%print()
    CALL t_fc3int%print()
    CALL t_fc3m2 %print()
    CALL t_fc3rot%print()
    CALL t_lwphph%print()
    CALL t_lwisot%print()
    CALL t_lwcasi%print()
    CALL t_velcty%print()
    CALL t_mpicom%print()
    CALL t_readdt%print()
    CALL t_lwinout%print()
    CALL t_xain%print()
  END SUBROUTINE
  !
  SUBROUTINE set_time_limit(max_seconds, max_time)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: max_seconds
    REAL(DP),INTENT(in) :: max_time
    REAL(DP) :: aux
    ! let' start time tracking in case it has not already
    CALL init_nanoclock()
    ! max_time is in hh.mmss format, we do not check if minute>59 or seconds>59
    IF(max_time>0) THEN
      limit_seconds = NINT(100*(max_time*100 - INT(max_time*100)))
      aux = DBLE(INT(max_time*100))/100._dp
      limit_seconds = limit_seconds &
                     +NINT(60*(100*aux-INT(aux)*100))&
                     +3600*NINT(max_time)

      !print*, "lll", max_time, 100*(max_time*100 - INT(max_time*100)), &
      !0*(100*aux-INT(aux)*100), 3600*INT(max_time)
    ELSE IF(max_seconds>0)THEN
      limit_seconds = max_seconds
    ENDIF
    IF(limit_seconds>0)THEN
      ioWRITE(stdout, '(2x,a,i12,a)') "Code will try to stop after", limit_seconds, "s"
    ENDIF
  END SUBROUTINE
  !
  LOGICAL FUNCTION check_time_limit()
    IMPLICIT NONE
    !print*, "check", limit_seconds, INT(get_wall())
    check_time_limit = (INT(get_wall())>limit_seconds)&
                     .and. (limit_seconds>0)
  END FUNCTION


END MODULE
