!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto - Nov 2011
MODULE d3_restart

  ! Number of processors in every possible parallelisation level (several not supported by D3)
  ! restart is only tested with exactly the same parallel configuration (changing nproc_pool MAY work)
  USE mp_world,  ONLY : cur_nproc => nproc
!  USE mp_global, ONLY : cur_nproc_pool => nproc_pool !, &
!                         cur_nproc_image => nproc_image, & ! not in D3 (mpi)
!                         cur_nproc_ortho => nproc_ortho, & ! not in D3 (mpi)
!                         cur_nogrp => nogrp, &             ! not in D3 (openmp)
!                         cur_npgrp => npgrp                ! not in D3 (openmp)

  LOGICAL :: done_dwfc(-3:3,-3:3) = .false. ! set in generate_dwfc2 (solve_linter_d3q.f90), after computing d_i2 wfc_i1
  LOGICAL :: done_pdvp(-3:3,-3:3) = .false. ! as above, for <psi|d_i2 V|psi k+q_i1>
  LOGICAL :: done_lmetq0          = .false. ! as above, after computing E_fermi shift
  LOGICAL :: done_nscf            = .false. ! set in run_nscf_d3 after computing ground-state wavefunctions

  CHARACTER(len=8),PARAMETER :: d3_rfile_suffix = '.d3rfile'

CONTAINS
!
SUBROUTINE d3_from_scratch()
  ! Set up the restart flags to initial data, i.e. "nothing has been done"
  IMPLICIT NONE
  done_dwfc(-3:3,-3:3) = .false.
  done_pdvp(-3:3,-3:3) = .false.
  done_lmetq0          = .false.
  done_nscf            = .false.
END SUBROUTINE d3_from_scratch
!
SUBROUTINE d3_check_restart(what)
  ! Write, Read OR Cleanup the restart file (tmp_dir_d3/prefix.rfile)
  ! all the data stored in the file is present in this module.  
  !
  ! 'read': after reading also broadcast to all nodes, if the data cannot be read it
  !        will call d3_from_scratch in order to initialize it.
  !
  ! 'write': after writing it will check if the limit time has been reached and, if it
  !         is the case, stop the code.
  !
  ! 'clean': delete the restart file
  !
  USE mp,         ONLY : mp_bcast
  USE io_files,   ONLY : prefix
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE d3_iofiles, ONLY : tmp_dir_d3
  USE mp_world,   ONLY : world_comm
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(in) :: what 
  !
  CHARACTER(len=256) :: rfile
  INTEGER :: ios, runit
  CHARACTER(len=12),PARAMETER :: sub='restart_dwfc'
  INTEGER :: rst_nproc, rst_nproc_pool
  INTEGER,EXTERNAL :: find_free_unit
  !
  ! REMEMBER TO UPDATE mp_bcast AND d3_check_restart IF YOU CHANGE THE NAMELIST:
  NAMELIST / d3_restart_info / &
      done_dwfc, done_pdvp, done_lmetq0, done_nscf, &
      rst_nproc, rst_nproc_pool
  !
  !WRITE(stdout, "(7x,a)") "REMARK: restart not tested, skipping"
  !RETURN ! <------------ FIXME
  !
  ! FIXME: use seqopn_d3 for this:
  IF(ionode)THEN
    rfile=TRIM(tmp_dir_d3)//"/"//TRIM(prefix)//TRIM(d3_rfile_suffix)
    runit=find_free_unit()
  ENDIF
  !
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRITE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  IF (TRIM(what)=='write') THEN
    !
    IF (ionode) THEN
      rst_nproc      = cur_nproc
!      rst_nproc_pool = cur_nproc_pool
      !
      OPEN(UNIT  = runit, ACCESS= 'sequential',  FILE  = TRIM(rfile), &
            FORM  ='formatted', status='unknown', iostat=ios)
      IF (ios==0) WRITE(runit,d3_restart_info, iostat=ios)
      IF (ios==0) CLOSE(runit, iostat=ios)

    ENDIF
    CALL mp_bcast(ios, ionode_id, world_comm)
    IF(ios/=0) CALL errore(sub,'cannot write restart file "'//TRIM(rfile)//'"',1)
    !
    ! this is a nice spot to stop: check time!
    CALL d3_check_time()
    !
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> READ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ELSE IF (TRIM(what)=='read') THEN
    !
    ! Initialize the variables:
    CALL d3_from_scratch()
    !
    WRITE(stdout, "(7x,a)") "REMARK: restart not tested, skipping"
    RETURN
    !
    IF (ionode) THEN
      OPEN(UNIT  = runit, ACCESS= 'sequential',  FILE  = TRIM(rfile), &
            FORM  ='formatted', status='old', iostat=ios)
      IF (ios==0) READ(runit,d3_restart_info,iostat=ios)
      CLOSE(runit)
    ENDIF
    !
    CALL mp_bcast(ios, ionode_id, world_comm)
    IF (ios/=0) THEN
      WRITE(stdout,'(5x,a)') "REMARK: recover data not found or invalid: restarting from scratch."
      ! Reinitialize, in case some garbage was read
      CALL d3_from_scratch()
    ELSE
      WRITE(stdout,'(5x,a)') "REMARK: recover data has been read succesfully."
      CALL mp_bcast(done_nscf,   ionode_id, world_comm)
      CALL mp_bcast(done_dwfc,   ionode_id, world_comm)
      CALL mp_bcast(done_pdvp,   ionode_id, world_comm)
      CALL mp_bcast(done_lmetq0, ionode_id, world_comm)
      ! check the number of processors
      CALL mp_bcast(rst_nproc,      ionode_id, world_comm)
      CALL mp_bcast(rst_nproc_pool, ionode_id, world_comm)
!      IF(rst_nproc/=cur_nproc .or. rst_nproc_pool/=cur_nproc_pool) &
        CALL errore(sub, 'Trying to restart with a different number of processor or pools.', 2)
    ENDIF
    !
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CLEAN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ELSE IF (TRIM(what)=='clean') THEN
    IF (ionode) THEN
      OPEN(UNIT  = runit, ACCESS= 'sequential',  FILE  = TRIM(rfile), &
           FORM  ='formatted', status='unknown', iostat=ios)
      IF (ios==0) CLOSE(runit, iostat=ios, status='DELETE')
    ENDIF
    CALL mp_bcast(ios, ionode_id, world_comm)
    IF(ios/=0) CALL errore(sub,'cannot delete restart file "'//TRIM(rfile)//'"',3)  
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ..... >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ELSE
    CALL errore(sub, "can only do 'read' or 'write' or 'clean'.",9)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE d3_check_restart
!
SUBROUTINE d3_check_time
  USE kinds, ONLY : DP
  USE d3_control,         ONLY : code, max_seconds
  USE d3_reset_module,    ONLY : d3_reset
  USE stop_d3_module,     ONLY : stop_d3
  USE io_global,          ONLY : stdout
  !
  IMPLICIT NONE
  REAL(DP) :: t1
  REAL(DP),EXTERNAL :: get_clock
  ! Time limit not set: just continue
  IF(max_seconds<0) RETURN
  !
  t1 = get_clock(code)
  ! Time limit not reached yet: continue
  IF(t1<max_seconds) RETURN
  !
  ! Time limit reached: STOP!
  write( stdout, '(/,5x,"===================================================")')
  write( stdout,   '(5x,"=   Wall time:",i6," > limit:",i6," => STOPPING!  =")') INT(t1), max_seconds
  write( stdout, '(  5x,"===================================================",/)')
  !
  CALL d3_reset(print_clock=.false., cleanup=.false.)
  CALL print_clock_d3()
  CALL stop_d3()
  FLUSH(stdout)
  !
  RETURN
  !
END SUBROUTINE d3_check_time
!
END MODULE d3_restart
