!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
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
                     t_fc3dint= nanotimer("fc3 double intp"), &
                     t_fc3m2  = nanotimer("fc3 modulus sq"), &
                     t_fc3rot = nanotimer("fc3 rotate"), &
                     t_lwphph = nanotimer("lw ph-ph"), &
                     t_lwisot = nanotimer("lw isotopes"), &
                     t_lwcasi = nanotimer("lw casimir"), &
                     t_velcty = nanotimer("ph group velocity"), &
                     t_mpicom = nanotimer("mpi communication"), &
                     t_readdt = nanotimer("read fc data"), &
                     t_merged = nanotimer("merge degenerate")

  TYPE(nanotimer) :: t_tksma   = nanotimer("sma thermalk"), &
                     t_tksum   = nanotimer("sum of tk terms"), &
                     t_tkaout  = nanotimer("prepare A_out"), &
                     t_tkain   = nanotimer("prepare and apply A_in"), &
                     t_tkcg    = nanotimer("compute CG step"), &
                     t_tkprec  = nanotimer("precond./initialize"), &
                     t_tktld   = nanotimer("1/sqrt(A_out)"), &
                     t_xain    = nanotimer("A_in * f (A_in given)"), &
                     t_lwchk   = nanotimer("check A_out>0"), &
                     t_lwinout = nanotimer("lw input/output"), &
                     t_restart = nanotimer("restart data i/o")
  TYPE(nanotimer) :: t_asr3a   = nanotimer("asr3 iteration"), &
                     t_asr3s   = nanotimer("asr3 symmetrize"), &
                     t_asr3io  = nanotimer("asr3 in/output"), &
                     t_asr3idx = nanotimer("asr3 index")
                     
  TYPE(nanotimer) :: t_spf          = nanotimer("spectral function"), &
                     t_qresolved    = nanotimer("q-resolved spf"), &
                     t_qresolved_io = nanotimer("q-resolved i/o & comm"), &
                     t_qsummed      = nanotimer("q-resolved spf"), &
                     t_qsummed_io   = nanotimer("q-resolved i/o & comm")
                     
  TYPE(nanotimer) :: t_optimize     = nanotimer("optimize grid")
                    

  CONTAINS
  SUBROUTINE print_all_timers()
    IMPLICIT NONE

    ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
    CALL print_timers_header()
    CALL t_tksma%print()
    CALL t_tksum%print()
    CALL t_tkaout%print()
    CALL t_tkain%print()
    CALL t_tkcg%print()
    ioWRITE(stdout,'("   * lw contribs ")')
    CALL t_lwphph%print()
    CALL t_lwisot%print()
    CALL t_lwcasi%print()
    ioWRITE(stdout,'("   * fine grain terms ")')
    CALL t_freq%print()
    CALL t_bose%print()
    CALL t_sum%print()
    CALL t_fc3int%print()
    CALL t_fc3m2 %print()
    CALL t_fc3rot%print()
    CALL t_velcty%print()
    CALL t_mpicom%print()
    CALL t_readdt%print()
    CALL t_lwinout%print()
    CALL t_restart%print()
    CALL t_tkcg%print()
    CALL t_tkprec%print()
    CALL t_tktld%print()
    CALL t_xain%print()
    CALL t_merged%print()
    ! ASR timers:
    CALL t_asr3a%print()
    CALL t_asr3s%print()
    CALL t_asr3io%print()
    CALL t_asr3idx%print()
    ! q-resolved spectra function
    CALL t_qresolved%print()
    CALL t_qresolved_io%print()
    CALL t_qsummed%print()
    CALL t_qsummed_io%print()
    !
    CALL t_optimize%print()
    
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
