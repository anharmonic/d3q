!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE timers
  USE nanoclock, ONLY : nanotimer, get_wall, init_nanoclock
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
                     t_merged = nanotimer("merge degenerate"), &
                     t_mkspf  = nanotimer("spectral function"), &
                     t_iodata = nanotimer("read data")

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
    CALL t_iodata%print()
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
    CALL t_mkspf%print()
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
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_percent_wall(gran_pc, gran_sec, i, n, reset)
    USE mpi_thermal, ONLY : mpi_any, abort_mpi
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: gran_pc, gran_sec
    INTEGER,INTENT(in) :: i, n
    LOGICAL,INTENT(in) :: reset
    !
    REAL(DP),SAVE :: last_print = 0._dp
    REAL(DP),SAVE :: iter_start = 0._dp
    REAL(DP) :: wall, fact, pc, iter_time, iter_end_time
    LOGICAL ::print_now, time_is_out
    !
    fact = 100._dp / gran_pc
    !ioWRITE(*,*) fact*DBLE(i)/n,INT(DBLE(fact*(i-.5_dp))/n), fact/n
    print_now = fact*DBLE(i)/n-INT(DBLE(fact*(i-.5_dp))/n)  <  1.49_dp*fact/n

    wall = get_wall()
    print_now = print_now .or. (wall-last_print)>gran_sec
!    print_now = print_now .and. (wall-last_print)>gran_sec/10._dp

    print_now = print_now .or. (i==n)
    print_now = print_now .or. (i==2) ! give an estimate as soon as possible

    IF(reset) iter_start = wall
    !
    last_print = wall
    pc = 100*DBLE(i-1)/(n-1)
    iter_time = (wall-iter_start)
    !
    IF(pc>0._dp) THEN
      iter_end_time = 100*iter_time/pc
      IF(print_now.and.ionode) WRITE(stdout,'(f12.1,"% | STEP TIME:",f12.1,"s | STEP END: ",f12.1,"s &
                    &| WALL:",f12.1,"s ")') &
                    pc, iter_time,iter_end_time, wall
    ELSE
      iter_end_time = 0
      IF(print_now .and. ionode) WRITE(stdout,'(f12.1,"% | STEP TIME:",f12.1,"s | STEP END: ",9x,"-.-s &
                    &| WALL:",f12.1,"s ")') &
                    pc, iter_time, wall
    ENDIF
    !
    time_is_out = check_time_limit(iter_end_time+iter_start)
    ! the following line causes the code to hang when the number of q-point per CPU is not equal
    ! and print_percent_wall is called inside a q-point loop (i.e. like in the SPF subroutine)
    !CALL mpi_any(time_is_out)
    IF(time_is_out)THEN
      ioWRITE(stdout,'(a,2/,a,2/,a)') "********","Cannot finish before time limit: giving up to save resources","********"
      CALL print_all_timers()
      CALL abort_mpi(255)
    ENDIF
    !
    FLUSH(stdout)
  END SUBROUTINE
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_timers_header()
    IMPLICIT NONE
    !
    ioWRITE(*,'(2x," * ",24x," * ",12x," ms * ",7x," ms/call * ",8x," ms*cpu * "'&
    //',3x," ms*cpu/call * ", " % wtime * ",6x," calls *")')
    !
  END SUBROUTINE print_timers_header
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_memory()
    USE clib_wrappers,           ONLY : memstat
    IMPLICIT NONE
    INTEGER :: kb
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
  !
  ! \/o\________\\\_________________________________________/^>
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
      ioWRITE(stdout, '(2x,a,i12,a)') "Code will try to stop if it cannot finish in", limit_seconds, "s"
    ENDIF
  END SUBROUTINE
  !
  ! \/o\________\\\_________________________________________/^>
  !
  LOGICAL FUNCTION check_time_limit(predicted_final_time)
    IMPLICIT NONE
    REAL(DP),OPTIONAL,INTENT(in) :: predicted_final_time
    REAL(DP) :: time
    IF(present(predicted_final_time)) THEN
       time = predicted_final_time
    ELSE
       time = get_wall()
    ENDIF
    !print*, "check", limit_seconds, INT(get_wall())
    check_time_limit = (INT(time)>limit_seconds)&
                     .and. (limit_seconds>0)
  END FUNCTION

 
END MODULE
