!
! Copyright (C) 2001-2008 Quantum-ESPRSSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the ionode_id directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3_readin_module
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_readin()
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq. T
  !    input is read from unit 5. A namelist is used on the machine which
  !    allows it. A second routine read_file reads the variables saved
  !    on the data file by the self-consistent program.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : ntyp => nsp, amass
  USE uspp,             ONLY : okvan
  USE pwcom,            ONLY : lsda
  USE run_info,         ONLY : title
  USE control_flags,    ONLY : iverbosity
  USE phcom,            ONLY : lgamma
  USE d3com,            ONLY : ethr_ph, istop, d3_mode, fild3dyn, max_seconds
  USE d3_control,       ONLY : print_star, print_perm, print_trev,&
                               restart, safe_io, d3dir
  USE noncollin_module, ONLY : noncolin
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : ionode, stdout
  USE constants,        ONLY : eps8
  USE parameters,       ONLY : npk
  USE control_ph,       ONLY : tmp_dir_ph
  USE save_ph,          ONLY : tmp_dir_save
  USE nscf_d3,          ONLY : run_nscf_d3
  USE d3_iofiles,       ONLY : fild1rho, fild2rho, fild3rho, fildrho_dir
  USE d3_grid,          ONLY : d3_grid_init, d3_single_point_init, d3_grid_slice
  USE d3_kgrid,         ONLY : nk1=>d3_nk1, nk2=>d3_nk2, nk3=>d3_nk3, &
                               k1=>d3_k1, k2=>d3_k2, k3=>d3_k3, &
                               degauss=>d3_degauss
  USE d3_debug,         ONLY : read_d3_debug, bcast_d3_debug
  !
  IMPLICIT NONE
  !
  REAL(DP):: xq1(3),xq2(3),xq3(3)
  INTEGER :: nq1, nq2, nq3
  !
  INTEGER :: ios, ipol, it
  INTEGER :: first, last, step, offset
  ! counters
  CHARACTER(len=256) :: outdir, fildrho, mode
  CHARACTER(len=9),PARAMETER :: sub='d3_readin'
  REAL(DP),PARAMETER :: Gamma(3) = (/ 0._dp, 0._dp, 0._dp /)
  REAL(DP) :: max_time
  !
  CHARACTER(len=256),EXTERNAL :: trimcheck
  LOGICAL,EXTERNAL :: eqvect

  NAMELIST / inputph / ethr_ph, amass, prefix, & ! controls, pw.x prefix
       outdir, fildrho_dir, d3dir, & ! directories (of $prefix.save, of fildrho files, for scratch)
       fild3dyn, fildrho, fild1rho, fild2rho, fild3rho, & ! output/input
       istop, mode, iverbosity, &  ! legacy, to be removed or fixed
       first, last, step, offset, & ! partial grid definition
       safe_io, restart, max_time, max_seconds, & ! restart controls (sort of)
       print_star, print_perm, print_trev, &  ! fildrho files to write
       nk1, nk2, nk3, k1, k2, k3, degauss
  ! convergence threshold
  ! atomic masses
  ! write control
  ! directory for temporary files
  ! the punch file produced by pwscf
  ! file with the dynamical matrix (output)
  ! file with deltarho at q
  ! file with deltarho at p (gamma in the old code)
  ! file with deltarho at -p-q (new)
  ! list of the q=0 modes to be computed
  ! .true.==> writes some auxiliary
  ! .true.==> this is a recover run
  ! to stop the program at a given point
  ! variables used for testing purposes
  CALL start_clock('d3_readin')
  !
  ONLY_IONODE : &
  IF ( ionode ) THEN
     !
     CALL input_from_file()
     !    Read the first line of the input file
     READ (5, '(a)', iostat = ios) title
     IF(ios/=0) CALL errore (sub, 'reading title ', ABS (ios) )
     !
     !   set default values for variables in namelist
     !
     ethr_ph = 1.d-5
     iverbosity = 0
     CALL get_env( 'ESPRESSO_TMPDIR', outdir )
     IF ( TRIM( outdir ) == ' ' ) outdir = './'
     CALL get_env( 'ESPRESSO_D3DIR', d3dir )
     IF ( TRIM( d3dir ) == ' ' ) d3dir=outdir
     CALL get_env( 'ESPRESSO_FILDRHO_DIR', fildrho_dir )
     IF ( TRIM( fildrho_dir ) == ' ' ) fildrho_dir=outdir
     outdir= trimcheck(outdir)
     d3dir = trimcheck(d3dir)
     fildrho_dir = trimcheck(fildrho_dir)

     
     prefix = 'pwscf'
     mode   = 'single triplet'
     fild3dyn = 'd3dyn'
     fildrho = ' '
     fild1rho = ' '
     fild2rho = ' '
     fild3rho = ' '
     ! dimensions of regular MP grid: 
     nq1=0
     nq2=0
     nq3=0
     ! definition of triplets to compute, for manual parallelisation:
     first=-1
     last=-1
     step=1
     offset=0
     !
     restart = .false.    ! true  : scan for fild3dyn files and do not recompute (grid calc only)
     safe_io = .false.    ! true  : close and reopen dpsi and psidHpsi files to guarantee partial restart
     print_star = .true.  ! false : do not print the star of the q triplet
     print_perm = .false. ! true  : print also inequivalent permutations of q1, q2 and q3
     print_trev = .true.  ! false : do not use time-reversal to print D3 of -q1,-q2,-q3 when they're not in the star
     istop = 0
     !
     max_seconds = -1     ! stop after this Wall time
     max_time    = -1._dp ! as above, in the hh.mm format, 
                          ! note: only one can be specified
     nk1=-1; nk2=-1; nk3=-1 ! do not change the k-point grid
     k1=-1;  k2=-1;  k3=-1
     degauss=-1._dp
     !
     !     reading the namelist inputph
     !
     WRITE(stdout, '(5x,a)') "Waiting for input..."
     CALL flush_unit( stdout )

     READ (5, inputph, iostat = ios)
     IF(ios/=0) CALL errore (sub, 'reading inputph namelist', ABS (ios) )
     !
     !    reads the q-point
     !
     d3_mode = TRIM(mode)
     IF ( mode(1:10) == "gamma-only" ) THEN
        WRITE(stdout, '(5x,a)') "Doing Gamma point calculation"
        xq1 = 0._dp
        xq2 = 0._dp
        xq3 = 0._dp
        lgamma = .true.
     ELSE IF ( mode(1:7) == "gamma-q" ) THEN
        WRITE(stdout, '(5x,a)') "Doing 0,-q,q calculation"
        READ (5, *, iostat = ios) (xq3 (ipol), ipol = 1, 3)
        IF(ios/=0) CALL errore (sub, 'reading xq1', ABS (ios) )
        xq2 = -xq3
        xq1 = 0._dp
        lgamma = SUM(xq1**2) < eps8 .and. SUM(xq2**2) < eps8 .and. SUM(xq3**2) < eps8
     ELSE IF ( mode(1:6) == "single" ) THEN
        WRITE(stdout, '(5x,a)') "Doing single-triplet calculation"
        READ (5, *, iostat = ios) (xq1 (ipol), ipol = 1, 3)
        IF(ios/=0) CALL errore (sub, 'reading xq1', ABS (ios) )
        READ (5, *, iostat = ios) (xq2 (ipol), ipol = 1, 3)
        IF(ios/=0) CALL errore (sub, 'reading xq2', ABS (ios) )
        READ (5, *, iostat = ios) (xq3 (ipol), ipol = 1, 3)
        IF(ios/=0) CALL errore (sub, 'reading xq3', ABS (ios) )
        lgamma = SUM(xq1**2) < eps8 .and. SUM(xq2**2) < eps8 .and. SUM(xq3**2) < eps8
     ELSE IF ( mode(1:7) == "partial" ) THEN
        WRITE(stdout, '(5x,a)') "Doing partial-grid dispersion calculation"
        READ (5, *, iostat = ios) (xq1 (ipol), ipol = 1, 3)
        IF(ios/=0) CALL errore (sub, 'reading xq', ABS (ios) )
        READ (5, *, iostat = ios) nq1, nq2, nq3
        IF(ios/=0) CALL errore (sub, 'reading nq1, nq2, nq3', ABS (ios) )
     ELSE IF ( mode(1:4) == "full" ) THEN
        WRITE(stdout, '(5x,a)') "Doing grid dispersion calculation"
        READ (5, *, iostat = ios) nq1, nq2, nq3
        IF(ios/=0) CALL errore (sub, 'reading nq1, nq2, nq3', ABS (ios) )
     ELSE
      CALL errore(sub, "Possible calculations: 'gamma-only', 'gamma-q', 'single triplet', "//&
                       "'partial grid', 'full grid'", 1)
     ENDIF
     !
     ! Read debug flags
     CALL read_d3_debug(5)
     !
     ! This is the temporary directory used by pw.x, it MUST be left unchanged for read_file to work:
     tmp_dir = trimcheck(outdir)
     !
     !     Check all namelist variables
     !
     IF (ethr_ph.LE.0.d0) CALL errore (' d3_readin', ' Wrong ethr_ph ', 1)
     IF (iverbosity.NE.0.AND.iverbosity.NE.1) &
          CALL errore (sub, 'Wrong iverbosity', 1)

     IF (fildrho /= ' ') THEN !CALL errore (sub, 'WARNING! This is the new code: use fild{1,2,3}rho instead of fildrho', 1)
        IF (fildrho(1:5) /= 'auto:' .and. .not. lgamma) &
           CALL errore(sub, 'fildrho must start with "auto:" (except for gamma-only calculation)',1)
        IF (fild1rho/=' ' .or. fild2rho/=' ' .or. fild3rho/=' ') &
          CALL errore(sub, 'you can input either fildrho or fildXrho (X=1,2,3) NOT both!', 5)
        fild1rho = fildrho
        fild2rho = fildrho
        fild3rho = fildrho
     ELSE
        IF (fild1rho==' ') CALL errore (sub, 'Wrong fild1rho '//TRIM(fild1rho), 2)
        IF (fild2rho==' ') CALL errore (sub, 'Wrong fild2rho '//TRIM(fild2rho), 3)
        IF (fild3rho==' ') CALL errore (sub, 'Wrong fild3rho '//TRIM(fild3rho), 4)
     ENDIF
     !
     IF (max_time>0._dp) THEN
       IF(max_seconds>0) CALL errore(sub, 'You can specify either max_seconds or max_time, not both', 8)
       max_seconds = NINT(60*100*(max_time-INT(max_time)))
       IF(max_seconds>=3600) CALL errore(sub, 'max_time must be in the hh.mm format',9)
       max_seconds = max_seconds + 3600*INT(max_time)
     ENDIF
     IF(max_seconds>0) WRITE(stdout, '(5x,a,i10,a)') "Max time:",max_seconds,"s"
     !
  END IF &
  ONLY_IONODE
  !
  ! >>>>>>>>>> BROADCAST <<<<<<<<<<
  CALL bcast_d3_input()
  CALL bcast_d3_debug()
  !
  !IF( degauss /= 0._dp .and. ( eqvect(xq1,Gamma,Gamma, 1.d-5) &
  !                         .or.eqvect(xq2,Gamma,Gamma, 1.d-5) &
  !                         .or.eqvect(xq3,Gamma,Gamma, 1.d-5) ) )&
  !    CALL errore(sub, 'One q is equivalent to Gamma, but is not Gamma, and system is metallic:'//&
  !                             ' Efermi shift not implemented!',1)
  !
  ! *** GLOBAL VARIABLES (is this still necessary??):
  ! ***  USE save_ph,     ONLY : tmp_dir_save  <--- prefix used inside pw.x
  ! ***  USE control_ph,  ONLY : tmp_dir_ph    <--- prefix used inside ph.x 
  !                             (messy: except for lgamma .AND. .not. ldisp, or for -nimage>1)
  ! ***  USE io_files,    ONLY : tmp_dir       <--- prefix actually used by diropn by default
  tmp_dir_save=tmp_dir
  tmp_dir_ph = TRIM(tmp_dir)//'_ph0'
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !   NOTE: here we read the simple PW files, not the ones modified by the
  !   phonon code!
  !
  CALL read_file()
  !
!   kplusq( 1)%xq =  xq1
!   kplusq( 2)%xq =  xq2
!   kplusq( 3)%xq =  xq3
  !
  ! Set up the points on which D3 shall be computed:
  IF ( mode(1:7) == "partial" ) THEN
    CALL d3_grid_init(nq1,nq2,nq3, xq1)
  ELSE IF ( mode(1:4) == "full" ) THEN
    CALL d3_grid_init(nq1,nq2,nq3)
  ELSE
    CALL d3_single_point_init(xq1, xq2, xq3)
  ENDIF
  !
  CALL d3_grid_slice(first, last, step, offset)
  !
  IF (lsda)     CALL errore (sub, 'lsda not implemented', 1)
  IF (okvan)    CALL errore (sub, 'US not implemented', 1)
  IF (noncolin) CALL errore (sub, 'd3 is not working in the noncolinear case', 1)

  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
!   CALL allocate_part()

  DO it = 1, ntyp
     IF (amass (it) <= 0.d0) CALL errore (sub, 'Wrong masses', it)
  ENDDO
  !
  CALL stop_clock('d3_readin')
  !
  RETURN
  !
  ! I do not like contained subroutines, but here it is very practical..
  CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_d3_input
    !-----------------------------------------------------------------------
    !
    !     In this routine the first processor sends input data to all
    !     the other processors
    !
#ifdef __MPI
    USE mp,             ONLY : mp_bcast
    USE io_global,      ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    CALL mp_bcast(mode,    ionode_id)
    CALL mp_bcast(tmp_dir, ionode_id)
    CALL mp_bcast(d3dir,   ionode_id)
    CALL mp_bcast(prefix,  ionode_id)
    !
    CALL mp_bcast(xq1, ionode_id)
    CALL mp_bcast(xq2, ionode_id)
    CALL mp_bcast(xq3, ionode_id)
    !
    CALL mp_bcast(nq1, ionode_id)
    CALL mp_bcast(nq2, ionode_id)
    CALL mp_bcast(nq3, ionode_id)
    !
    CALL mp_bcast(fildrho,  ionode_id)
    CALL mp_bcast(fild1rho, ionode_id)
    CALL mp_bcast(fild2rho, ionode_id)
    CALL mp_bcast(fild3rho, ionode_id)
    !
    CALL mp_bcast(first,   ionode_id)
    CALL mp_bcast(last,    ionode_id)
    CALL mp_bcast(step,    ionode_id)
    CALL mp_bcast(offset,  ionode_id)
    CALL mp_bcast(safe_io, ionode_id)
    !
    CALL mp_bcast(fild3dyn,   ionode_id)
    CALL mp_bcast(ethr_ph,    ionode_id)
    CALL mp_bcast(amass,      ionode_id)
    CALL mp_bcast(iverbosity, ionode_id)
    CALL mp_bcast(lgamma,     ionode_id)
    !
    CALL mp_bcast(restart,    ionode_id)
    CALL mp_bcast(max_seconds,ionode_id)
    CALL mp_bcast(print_star, ionode_id)
    CALL mp_bcast(print_perm, ionode_id)
    !
    CALL mp_bcast(nk1, ionode_id)
    CALL mp_bcast(nk2, ionode_id)
    CALL mp_bcast(nk3, ionode_id)
    CALL mp_bcast(k1, ionode_id)
    CALL mp_bcast(k2, ionode_id)
    CALL mp_bcast(k3, ionode_id)
    CALL mp_bcast(degauss, ionode_id)
    !
#endif
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE bcast_d3_input
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
END SUBROUTINE d3_readin
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_readin_module
!-----------------------------------------------------------------------

