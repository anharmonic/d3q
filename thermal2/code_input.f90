!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Code contributions from Raffaello Bianco
!
! <<^V^\\=========================================//-//-//========//O\\//
! This module is a big mess that read input and check if there are no conflicts
MODULE code_input
  USE kinds,    ONLY : DP
  USE mpi_thermal, ONLY : ionode
  USE timers
#include "mpi_thermal.h"
  !
  REAL(DP) :: default_sigma = 10._dp
  
  ! NOTE: energies in the LWINPUT_TYPE structure are assumed to be in CM^-1
  !       in the rest of the code energies are in Rydberg!
  TYPE code_input_type
    !
    CHARACTER(len=16) :: calculation ! lw=linewidth, spf=spectral function
    CHARACTER(len=16) :: mode        ! "full" or "simple" spectral function 
    CHARACTER(len=256) :: outdir
    CHARACTER(len=256) :: prefix     ! put this in front of file names
    !
    LOGICAL :: exp_t_factor ! add elastic peak (sort of not working)
    CHARACTER(len=7) :: sort_freq ! only applies to "full" calculation: 
                                 ! sort w+w_shift when saving to file. Instead of w,lw,ls 
                                 ! you have w,lw,ls+w with the last two blocks sorted differently 
                                 ! than the first one to avoid unesthetical jumps in band plots
    !
    CHARACTER(len=256) :: file_mat3
    CHARACTER(len=256) :: file_mat2
    CHARACTER(len=256) :: file_mat2_final
    CHARACTER(len=8)   :: asr2
    !
    INTEGER            :: skip_q
    INTEGER            :: nconf
    REAL(DP),ALLOCATABLE :: T(:), sigma(:)
    ! for spectral function:
    INTEGER :: ne
    REAL(DP) :: e0, de
    REAL(DP) :: sigma_e 
    ! for final state:
    INTEGER  :: nu_initial
    REAL(DP) :: e_initial
    REAL(DP) :: q_initial(3)
    LOGICAL  :: q_resolved, q_summed
    REAL(DP) :: sigmaq
    ! intrinsic ph-ph scattering, this is true by default:
    LOGICAL  :: intrinsic_scattering
    ! for isotope contribution to lw
    LOGICAL  :: isotopic_disorder
    ! for border scattering, grain size effects
    LOGICAL  :: mfp_cutoff
    LOGICAL  :: casimir_scattering
    ! NOTE: Casimir length must also include the structure factor (usually 0.5)
    REAL(DP) :: sample_length
    REAL(DP) :: sample_dir(3)
    !
    INTEGER :: nk(3), nk_in(3)
    REAL(DP) :: xk0(3), xk0_in(3)
    !
    ! only for tk:
    REAL(DP) :: thr_tk
    INTEGER  :: niter_max
    !
    CHARACTER(len=6) :: grid_type
    CHARACTER(len=6) :: grid_type_in
    ! for dynbubble and r2q:
    LOGICAL :: print_dynmat
    LOGICAL :: print_velocity
    LOGICAL :: print_neutron_cs
    !
    LOGICAL :: store_lw ! for tk-sma, save the linewidth to file
    !
    LOGICAL :: optimize_grid
    REAL(DP) :: optimize_grid_thr
    !
    LOGICAL :: restart
  END TYPE code_input_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(code, input, qpts, S, fc2, fc3)
    !USE io_global,      ONLY : stdout
    USE q_grids,              ONLY : q_grid, setup_path, setup_grid, setup_plane_grid
    USE constants,            ONLY : RY_TO_CMM1, BOHR_RADIUS_SI
    USE more_constants,       ONLY : INVALID, DHUGE, MASS_DALTON_TO_RY
    USE clib_wrappers,             ONLY : f_mkdir_safe
    USE fc3_interpolate,      ONLY : forceconst3
    USE nist_isotopes_db,     ONLY : compute_gs
    USE input_fc,             ONLY : div_mass_fc2, forceconst2_grid, ph_system_info
    USE mpi_thermal,          ONLY : ionode, mpi_broadcast
    USE cmdline_param_module, ONLY : cmdline_to_namelist
    USE timers
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*),INTENT(in)    :: code
    TYPE(code_input_type),INTENT(out) :: input
    TYPE(q_grid),INTENT(out)  :: qpts
    TYPE(forceconst2_grid),INTENT(out) :: fc2
    CLASS(forceconst3),POINTER,OPTIONAL,INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(out)   :: S
    !
    ! Input variable, and defaul values:
    CHARACTER(len=16)  :: calculation = "" ! "spf"
    CHARACTER(len=256) :: file_mat3  = INVALID ! no default
    CHARACTER(len=256) :: file_mat2  = INVALID ! no default
    CHARACTER(len=256) :: file_mat2_final  = INVALID ! default = file_mat2
    CHARACTER(len=256) :: prefix     = INVALID ! default: calculation.mode
    !
    CHARACTER(len=256) :: outdir = './'              ! where to write output files
    CHARACTER(len=8)   :: asr2 = "no"                ! apply sum rule to phonon force constants
    INTEGER            :: nconf = -1                 ! number of smearing/temperature couples
    INTEGER            :: nq = -1                    ! number of q-point to read, only for lw 
    INTEGER            :: skip_q = 0                 ! skip this many points when computing a BZ path
    INTEGER            :: nk(3) = (/-1, -1, -1/)     ! integration grid for lw, db and tk, (the outer one for tk_sma)
    REAL(DP)           :: xk0(3) = 0._dp             ! grid shift as fraction of half grid step
    REAL(DP)           :: xk0_in(3) = (/ DHUGE, DHUGE, DHUGE/)          ! grid shift as fraction of half grid step
    INTEGER            :: nk_in(3) = (/-1, -1, -1/)  ! inner integration grid, only for tk_sma
    LOGICAL            :: exp_t_factor = .false.     ! add elastic peak of raman, only in spectre calculation
    CHARACTER(len=7)   :: sort_freq = "default"      ! how to sort frequencies (default, overlap, shifted)
    CHARACTER(len=6)   :: grid_type="simple"         ! "simple" uniform qpoints grid, or "bz" symmetric BZ grid
    CHARACTER(len=6)   :: grid_type_in=INVALID       ! "simple" uniform qpoints grid, or "bz" symmetric BZ grid
    LOGICAL            :: print_dynmat = .false.     ! print the dynamical matrix for each q (only r2q and dynbubble code)
    LOGICAL            :: print_velocity   = .true.    ! print the phonon group velocity for each q (only r2q code)
    LOGICAL            :: print_neutron_cs = .false.    ! print the phonon cross-section for inelastic neutron scattering (only r2q code)
    LOGICAL            :: store_lw = .false.         ! for tk-sma calculation, store the lw for the grid point to file
    !
    ! The following variables are used for spectre and final state calculations
    INTEGER  :: ne = -1                 ! number of energies on which to sample the spectral decomposition
    REAL(DP) :: de = 1._dp, e0 = 0._dp  ! energy step and minimum
    REAL(DP) :: sigma_e  = -1._dp       ! smearing used for delta of e in plots, like jdos, phdos and final state decomposition
    INTEGER  :: nu_initial = 0          ! initial mode for final state decomposition, set to
                                        ! zero to get the sum of all modes at e_initial
    REAL(DP) :: e_initial = -1._dp      ! initial energy for final state decomposition
                                        ! may default to omega_{nu_initial} if not not set
    REAL(DP) :: q_initial(3) = 0._dp    ! initial q for final state decomp
    LOGICAL  :: q_resolved = .false.    ! save the final state function of q and e as well
                                        ! in different files
    LOGICAL  :: q_summed = .false.      ! save the final state of q as well
    REAL(DP) :: sigmaq = 0.1_dp         ! reciprocal space smearing for final q decomposition
    !
    ! Compute intrinsic ph-ph scattering, setting this to false is only useful
    ! to re-compute different casimir or isotopic parameters and then recompute
    ! SMA using tools/recompute_sma.m
    LOGICAL  :: intrinsic_scattering = .true.
    ! Use isotopic disorder (only for tk calculations)
    LOGICAL  :: isotopic_disorder = .false.
    REAL(DP),ALLOCATABLE :: isotopes_mass(:), isotopes_conc(:), auxm(:), auxs(:)
    INTEGER :: n_isotopes, atomic_N
    !
    ! Border scattering by cutting off phonons with mfp>sample_length
    LOGICAL  :: mfp_cutoff = .false.
    ! Border scattering via Casimir model (only tk calculations)
    LOGICAL  :: casimir_scattering = .false.
    ! Applies to both casimit_scattering and mfp_cutoff:
    REAL(DP) :: sample_length_au = -1._dp ! length in bohr
    REAL(DP) :: sample_length_mu = -1._dp ! length in micrometres
    REAL(DP) :: sample_length_mm = -1._dp ! length in millimitres
    REAL(DP) :: sample_dir(3) = 0._dp
    !
    REAL(DP) :: volume_factor = 1._dp ! volume scaling, e.g. for 2D slabs
    !
    INTEGER  :: max_seconds = -1
    REAL(DP) :: max_time    = -1._dp
    !
    REAL(DP) :: thr_tk = 1.d-2
    INTEGER  :: niter_max = 1000
    !
    LOGICAL :: restart = .false.
    !
    LOGICAL :: optimize_grid = .false.
    REAL(DP) :: optimize_grid_thr = 1.d-2
    !
    ! Local variables use to read the list or grid of q-points required by lw
    REAL(DP) :: xq(3), xq0(3), e1(3), e2(3)
    INTEGER  :: ios, ios2, i, j, ij, naux, nq1, nq2, nq3, nT_aux, nsigma_aux
    REAL(DP),ALLOCATABLE :: T_aux(:), sigma_aux(:)
    !
    CHARACTER(len=1024) :: line, word
    CHARACTER(len=16)   :: word2, word3
    CHARACTER(len=512)  :: input_file
    CHARACTER(LEN=256), EXTERNAL :: TRIMCHECK
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    INTEGER             :: input_unit, aux_unit, PASS
    INTEGER,EXTERNAL :: find_free_unit
    !
    LOGICAL :: qpoints_ok=.false.,  &! true after reading QPOINTS
               configs_ok=.false.,  &! true after reading CONFIGS
               isotopes_ok=.false., &! true after reading ISOTOPES
               do_grid=.false.       ! is true, construct a regular grid of q-points
    !
    NAMELIST  / lwinput / &
      calculation, outdir, prefix, &
      file_mat2, file_mat3, asr2, &
      nconf, skip_q, nq, nk, grid_type, xk0, &
      optimize_grid, optimize_grid_thr, &
      ne, de, e0, sigma_e, &
      nu_initial, e_initial, q_initial, q_resolved, q_summed, sigmaq,&
      exp_t_factor, sort_freq, &
      isotopic_disorder, &
      casimir_scattering,  &
      sample_length_au, sample_length_mu, sample_length_mm, sample_dir,&
      max_seconds, max_time

    NAMELIST  / tkinput / &
      calculation, outdir, prefix, &
      file_mat2, file_mat3, asr2, &
      thr_tk, niter_max, &
      nconf, nk, nk_in, grid_type, grid_type_in, xk0, xk0_in, &
      optimize_grid, optimize_grid_thr, &
      intrinsic_scattering, &
      isotopic_disorder, store_lw, &
      casimir_scattering, mfp_cutoff, sample_dir, &
      sample_length_au, sample_length_mu, sample_length_mm, &
      volume_factor, &
      max_seconds, max_time, restart

    NAMELIST  / dbinput / &
      calculation, outdir, prefix, &
      file_mat2, file_mat3, file_mat2_final, asr2, &
      nconf, nk, nq, xk0, grid_type, print_dynmat, &
      ne, de, e0, sigma_e, &
      max_seconds, max_time

    NAMELIST  / r2qinput / &
      calculation, outdir, prefix, &
      file_mat2, asr2, nq, sort_freq, &
      ne, de, e0, sigma_e, &
      nk, nconf, grid_type, xk0, &
      isotopic_disorder, &
      casimir_scattering,  &
      sample_length_au, sample_length_mu, sample_length_mm, sample_dir,&
      print_dynmat, print_velocity, print_neutron_cs
      !
      
    timer_CALL t_iodata%start()
    
    IONODE_READS_INPUT_FILE : &
    IF(ionode)THEN
      input_file="input."//TRIM(code)
      CALL parse_command_line(input_file)
      !
      PASSES : DO PASS = 1,2
        IF(PASS==1) THEN
          IF(TRIM(input_file)=="-")THEN
            ioWRITE(stdout,'(2x,3a)') "Warning! Reading standard input will probably not work with MPI"
            input_unit = 5
          ELSE
            ioWRITE(stdout,'(2x,3a)') "Reading input file '", TRIM(input_file), "'"
            input_unit = find_free_unit()
            OPEN(unit=input_unit, file=input_file, status="OLD", action="READ")
          ENDIF
          aux_unit = input_unit
        ELSE IF (PASS==2) THEN
          aux_unit = find_free_unit()
!           IF(ionode) THEN
            WRITE(stdout,'(2x,3a)') "merging with command line arguments"
            OPEN(unit=aux_unit, file=TRIM(input_file)//".tmp~", status="UNKNOWN", action="READWRITE")
            !OPEN(unit=aux_unit, status="SCRATCH", action="READWRITE")
            CALL cmdline_to_namelist(TRIM(code)//"input", aux_unit)
!             CALL mpi_wbarrier()
            REWIND(aux_unit)
!           ELSE
!             CALL mpi_wbarrier()
!             OPEN(unit=aux_unit, file="."//TRIM(input_file)//".tmp", status="UNKNOWN", action="READ")
!           ENDIF
        ENDIF
        !
        IF(code=="LW")THEN
          IF(PASS==1) calculation="lw"
          READ(aux_unit, lwinput)
          IF(PASS==2.and.ionode) WRITE(*, lwinput)
        ELSE IF (code=="TK")THEN
          IF(PASS==1) calculation="sma"
          do_grid = .true.
          READ(aux_unit, tkinput)
          IF(PASS==2.and.ionode) WRITE(*, tkinput)
        ELSE IF (code=="DB")THEN
          IF(PASS==1) calculation="db"
          READ(aux_unit, dbinput)
          IF(PASS==2.and.ionode) WRITE(*, dbinput)
        ELSE IF (code=="R2Q")THEN
          IF(PASS==1) calculation="freq"
          READ(aux_unit, r2qinput)
          IF(PASS==2.and.ionode) WRITE(*, r2qinput)
          IF(calculation=="freq" .and. PASS==1)THEN
            ! Do not read, nconf and configs in the R2Q-freq case
            IF(ANY(nk<=0)) nk=1
            IF(nconf<=0)THEN
              nconf=1
              configs_ok = .true.
            ENDIF
          ELSEIF(calculation=="rms" .or. calculation=="dos" &
                 .and. PASS==1)THEN
            ! Do not read, QPOINTS in the R2Q-rms and jdos case
            IF(nq<=0) nq = 1
            qpoints_ok = .true.
          ENDIF
        ELSE
          CALL errore('READ_INPUT', 'Wrong code', 1)
        ENDIF
        !
!         IF(PASS==2.and..not.ionode) CLOSE(aux_unit,STATUS="KEEP")
!         CALL mpi_wbarrier()
!         IF(PASS==2.and.ionode) CLOSE(aux_unit,STATUS="DELETE")
        IF(PASS==2) CLOSE(aux_unit,STATUS="DELETE")
      ENDDO PASSES
      !
    ENDIF &
    IONODE_READS_INPUT_FILE 
    !
    CALL broadcast_namelist_variables()
    CALL mpi_broadcast(do_grid)
    CALL mpi_broadcast(configs_ok)
    CALL mpi_broadcast(qpoints_ok)
    !
    IF(TRIM(file_mat2) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    IF(TRIM(file_mat2_final) == INVALID ) file_mat2_final = file_mat2
    IF(TRIM(file_mat3) == INVALID .and. present(fc3)) &
        CALL errore('READ_INPUT', 'Missing file_mat3', 1)
    IF(ANY(nk<0)) CALL errore('READ_INPUT', 'Missing nk', 1)    
    !IF(nconf<0)   CALL errore('READ_INPUT', 'Missing nconf', 1)    

    CALL set_time_limit(max_seconds, max_time)
    
    input%file_mat2 = file_mat2
    input%file_mat2_final = file_mat2_final
    input%file_mat3 = file_mat3
    input%outdir    = TRIMCHECK(outdir)
    input%asr2      = asr2
    input%skip_q   = skip_q
    input%nconf     = nconf
    input%nk        = nk
    input%grid_type = grid_type
    input%grid_type_in = grid_type_in
    input%exp_t_factor = exp_t_factor
    input%sort_freq    = sort_freq
    input%print_dynmat     = print_dynmat
    input%print_velocity   = print_velocity
    input%print_neutron_cs = print_neutron_cs
    input%store_lw         = store_lw 
    !
    input%intrinsic_scattering = intrinsic_scattering
    input%isotopic_disorder    = isotopic_disorder
    input%casimir_scattering   = casimir_scattering
    input%mfp_cutoff           = mfp_cutoff
    input%sample_dir           = sample_dir
    IF(casimir_scattering .and. mfp_cutoff) &
      CALL errore("code_input","don't use both casimir_scattering and mfp_cutoff",1)
    !
    input%thr_tk = thr_tk
    input%niter_max = niter_max
    input%restart = restart
    !
    input%optimize_grid = optimize_grid
    input%optimize_grid_thr = optimize_grid_thr
    !
    IF(ANY(nk_in<0)) nk_in = nk
    input%nk_in = nk_in
    IF(xk0_in(1)==DHUGE) xk0_in(1) = xk0(1)
    IF(xk0_in(2)==DHUGE) xk0_in(2) = xk0(2)
    IF(xk0_in(3)==DHUGE) xk0_in(3) = xk0(3)
    input%xk0_in    = 0.5_dp*xk0_in/input%nk_in
    IF(TRIM(input%grid_type_in)==INVALID) input%grid_type_in = grid_type
    input%xk0    = 0.5_dp*xk0/input%nk
    !
    ios = f_mkdir_safe(input%outdir)
    IF(ios>0) CALL errore('READ_INPUT', 'cannot create directory: "'//TRIM(input%outdir)//'"',1)
    !
    ! read data before reading the q-point, because we need the unit
    ! cell to go to/from crystal coords etc
    CALL READ_DATA(input, s, fc2, fc3)
    !
    IF(volume_factor/=1._dp)THEN
      s%Omega = s%Omega*volume_factor
      ioWRITE(*,'(2x,"Unit-cell volume rescaled by:",f12.6," new volume: ",f12.6," bohr^3")') &
      volume_factor, s%Omega
    ENDIF
    !
    READ(calculation,*,iostat=ios) input%calculation, input%mode
    IF(ios/=0) THEN
      input%calculation = calculation
      input%mode = "full"
    ENDIF
    ioWRITE(*,*) "calculation: ", input%calculation, input%mode
    !
!     IF(nq<0.and.TRIM(input%calculation)/="grid".and.code=="LW") &
!         CALL errore('READ_INPUT', 'Missing nq', 1)    
    CALL cryst_to_cart(1,input%xk0,S%bg, +1)
    xk0 = input%xk0 ! put it back there to do "outer grid" later (HACK)
    CALL cryst_to_cart(1,input%xk0_in,S%bg, +1)
    !
    IF(TRIM(prefix)==INVALID)THEN
      input%prefix = TRIM(input%calculation)//"_"//TRIM(input%mode)
    ELSE
      input%prefix = prefix
    ENDIF
    !
    IF(TRIM(input%calculation) == 'spf' .and. ne < 0) &
      CALL errore('READ_INPUT', 'Missing ne for spf calculation', 1)    

    input%ne = ne
    input%de = de
    input%e0 = e0
    IF(sigma_e<0._dp)THEN
      input%sigma_e = 5*de
    ELSE
      input%sigma_e = sigma_e
    ENDIF
    !
    IF(TRIM(input%calculation) == 'final' .and. (e_initial < 0 .and. nu_initial==0)) &
      CALL errore('READ_INPUT', 'Missing e_initial for final state calculation', 1)    
    input%nu_initial = nu_initial
    input%e_initial  = e_initial
    input%q_initial  = q_initial
    ! if we also want the q of the final state:
    input%q_resolved = q_resolved
    input%q_summed = q_summed
    input%sigmaq  = sigmaq
    IF(calculation == 'cgp' .and. (grid_type == 'bz' .or. grid_type == 'ws'))&
      CALL errore('READ_INPUT', "CGP calculation isn't properly implemented with bz grids", 1)
    !
    IF(calculation == 'cgp' .and. (grid_type == 'random'))&
      CALL errore('READ_INPUT', "CGP thermal conductivity is not variational with"//&
      "random grids: convergence is not guaranteed!!!", 1)
    
    IF(input%casimir_scattering .or. input%mfp_cutoff)THEN
      IF(COUNT((/sample_length_au,sample_length_mu,sample_length_mm/)>=0._dp)>1) &
        CALL errore('READ_INPUT', "You cannot specify more than one: sample_length_{au,mu,mm}",1)
      IF(sample_length_au<0._dp .and. sample_length_mu<0._dp .and. sample_length_mm<0._dp) &
        CALL errore('READ_INPUT', "You must specify one of: sample_length_{au,mu,mm}",2)
      IF(sample_length_au>0._dp) input%sample_length = sample_length_au
      IF(sample_length_mu>0._dp) input%sample_length = sample_length_mu/(BOHR_RADIUS_SI*1.d+6)
      IF(sample_length_mm>0._dp) input%sample_length = sample_length_mm/(BOHR_RADIUS_SI*1.d+3)
      ioWRITE(*,'(5x,a,1f12.2)') "Sample length (bohr)", input%sample_length
    ELSE
      input%sample_length = 0._dp
    ENDIF
    !
    IF(ionode) READ(input_unit,'(a1024)', iostat=ios) line
    CALL mpi_broadcast(line)
    CALL mpi_broadcast(ios)
    !
    READ_CARDS : &
    DO WHILE (ios==0)
      !
      !print*, ionode, ">>", TRIM(line)
      READ(line,*,iostat=ios2) word
      IF(ios2/=0) THEN
        IF(ionode) READ(input_unit,'(a1024)', iostat=ios) line
        CALL mpi_broadcast(line)
        CALL mpi_broadcast(ios)
        !
        CYCLE READ_CARDS 
      ENDIF
      !
      CARD_FOUND_CASE : &
      SELECT CASE (TRIM(word))
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE ("QPOINTS")
        IF(qpoints_ok) THEN
          !CALL errore("READ_INPUT", "Won't reads q-points twice", 1)
          ioWRITE(*,*) "WARNING! Ignoring QPOINTS (duplicated, or not required)"
          CYCLE READ_CARDS
        ENDIF
        qpoints_ok = .true.
        !
        IF(code=="TK") THEN
          ioWRITE(*,*) "Ignoring QPOINTS (code TK)"
          CYCLE READ_CARDS
        ENDIF
        !
        ioWRITE(*,*) "Reading QPOINTS"
        !
        qpts%basis = 'cartesian'
        READ(line,*,iostat=ios) word2, word3, xk0
        IF(ios/=0) THEN
          READ(line,*,iostat=ios) word2, word3
          xk0 = 0._dp
        ENDIF
        IF(ios==0) qpts%basis = TRIM(word3(1:9))
        !
        IF(TRIM(qpts%basis) == "grid" .or. TRIM(qpts%basis) == "bz" &
           .or. TRIM(qpts%basis)=="ws" .or.  TRIM(qpts%basis)=="randws" &
           .or. TRIM(qpts%basis)=="bxsf" .or. TRIM(qpts%basis)=="xsf" &
           .or. TRIM(qpts%basis)=="random" ) THEN
          !
          grid_type = qpts%basis
          qpts%basis = "cartesian"
          do_grid = .true.
          !
          IF(ionode) READ(input_unit,*,iostat=ios) nq1, nq2, nq3
          CALL mpi_broadcast(ios)
          CALL mpi_broadcast(nq1)
          CALL mpi_broadcast(nq2)
          CALL mpi_broadcast(nq3)
          xq0 = 0.5_dp * xk0 / (/ nq1, nq2, nq3 /)
          CALL cryst_to_cart(1,xq0,S%bg, +1)
          IF(ios/=0) CALL errore("READ_INPUT", "Reading QPOINTS nq1, nq2, nq3 for grid calculation", 1)
          line=''
          !this is done later
          !CALL setup_grid(input%grid_type, S%bg, nq1,nq2,nq3, qpts)
          CYCLE READ_CARDS
        ELSE IF(TRIM(qpts%basis) == "plane") THEN
          !
          grid_type = qpts%basis
          qpts%basis = "cartesian"
          do_grid = .false.
          !
          IF(ionode) READ(input_unit,*,iostat=ios) nq1, nq2
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT", "plane 0", 1)
          CALL mpi_broadcast(nq1)
          CALL mpi_broadcast(nq2)
          
          IF(ionode) READ(input_unit,*,iostat=ios) e1
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT", "plane 1", 1)
          CALL mpi_broadcast(3, e1)
          
          IF(ionode) READ(input_unit,*,iostat=ios) e2
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT", "plane 1", 1)
          CALL mpi_broadcast(3, e2)

          IF(ionode) READ(input_unit,*,iostat=ios) xq0
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT", "plane 1", 1)
          CALL mpi_broadcast(3, xq0)
          
          line=''
          !this is done later
          CALL setup_plane_grid(grid_type, S%bg, nq1,nq2, e1, e2, xq0, qpts)
          CYCLE READ_CARDS
        ENDIF
        !
        ! Read the number of q-points, if not specified in the namelist
        IF(nq<0) THEN 
          READ(input_unit,*, iostat=ios) nq
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT","Expecting number of q-points.", 1)
          CALL mpi_broadcast(nq)
        ENDIF
        !
        ! Read the q-points, one by one, add to path
        QPOINT_LOOP : & ! ..............................................................
        DO i = 1, nq
          IF(ionode) READ(input_unit,'(a1024)', iostat=ios) line
          CALL mpi_broadcast(ios)
          IF(ios/=0) CALL errore("READ_INPUT","Expecting q point: input error.", 1)
          CALL mpi_broadcast(line)
          !
          ! Try to read point and number of points
          READ(line,*, iostat=ios2) xq(1), xq(2), xq(3), naux
          IF(ios2==0) THEN
            IF(TRIM(qpts%basis) == "crystal")  CALL cryst_to_cart(1,xq,S%bg, +1)
            CALL setup_path(xq, naux, qpts, S%at)
            CYCLE QPOINT_LOOP
          ENDIF
          !
          ! Try to read just the point 
          READ(line,*, iostat=ios2) xq(1), xq(2), xq(3)
          IF(ios2==0) THEN
            IF(TRIM(qpts%basis) == "crystal")  CALL cryst_to_cart(1,xq,S%bg, +1)
            CALL setup_path(xq, 1, qpts, S%at)
            CYCLE QPOINT_LOOP
          ENDIF
          !
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting q point, got: '"//TRIM(line)//"'.", 1)
          EXIT QPOINT_LOOP
          !
        ENDDO &
        QPOINT_LOOP ! .................................................................
        
        ioWRITE(*,"(2x,a,i4,a,i6,a)") "Read", nq," lines, set-up ",qpts%nq,&
                                 " q-points, "//TRIM(qpts%basis)//" coordinates"
        !
        IF(TRIM(qpts%basis) == "crystal")  THEN
          !CALL cryst_to_cart(qpts%nq,qpts%xq,S%bg, +1)
          qpts%basis = "cartesian"
          ioWRITE(*,"(4x,a)") "q-points converted to cartesian coordinates (2pi/alat)"
        ENDIF
        ioWRITE(*,*)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE ("CONFIGS")
        IF(configs_ok) CALL errore("READ_INPUT", "Won't reads configs twice", 1)
        configs_ok = .true.
        !
        READ(line,*,iostat=ios) word2, word3
        IF(ios/=0) word3 = "list"
        !
        IF(TRIM(word3)=="matrix")THEN
          !
          nT_aux=-1
          nsigma_aux = -1
          READ_CONFIGS_MATRIX : &
          DO
            IF(ionode) READ(input_unit,"(a1024)", iostat=ios2) line
            CALL mpi_broadcast(ios2)
            IF(ios2/= 0) CALL errore("READ_INPUT","Error configs matrix.", 1)
            CALL mpi_broadcast(line)
            !
            READ(line,*, iostat=ios2)  word2, naux
            IF(ios2/=0)THEN
              READ(line,*, iostat=ios2)  word2
              IF(ios2/= 0) CALL errore("READ_INPUT","Error configs matrix line.", 2)
              naux = 1
            ENDIF
            !
            IF(TRIM(word2)=="T")THEN
              nT_aux = naux
              IF(nT_aux < 0) CALL errore("READ_INPUT","Bad number of temperatures.", 1)
              ALLOCATE(T_aux(nT_aux))
              IF(ionode) READ(input_unit,*, iostat=ios2) T_aux
              CALL mpi_broadcast(ios2)
              IF(ios/= 0) CALL errore("READ_INPUT","Error reading temperatures.", 1)
              CALL mpi_broadcast(nT_aux, T_aux)
            ELSE IF (TRIM(word2)=="sigma") THEN
              nsigma_aux = naux
              IF(nsigma_aux < 0) CALL errore("READ_INPUT","Bad number of sigmas.", 1)
              ALLOCATE(sigma_aux(nsigma_aux))
              IF(ionode) READ(input_unit,*, iostat=ios2) sigma_aux
              CALL mpi_broadcast(ios2)
              IF(ios/= 0) CALL errore("READ_INPUT","Error reading sigmas.", 1)
              CALL mpi_broadcast(nsigma_aux, sigma_aux)
            ELSE IF (TRIM(word2)=="" .or. word2(1:1) =='!')THEN
              CYCLE READ_CONFIGS_MATRIX
            ELSE
              CALL errore("READ_INPUT","Unexpected line in CONFIGS matrix: '"//TRIM(word2)//"'", 1)
            ENDIF
            !
            IF(nT_aux>0 .and. nsigma_aux>0) EXIT READ_CONFIGS_MATRIX
          ENDDO &
          READ_CONFIGS_MATRIX
          !
          nconf = nT_aux * nsigma_aux
          ioWRITE(stdout,'(2x,4(a,i4))') "Building matrix of ", nconf, " configurations, from",&
                                          nT_aux," temperature and ", nsigma_aux," smearings."
          input%nconf = nconf
          ALLOCATE(input%sigma(nconf), input%T(nconf))
          ij = 0
          DO i = 1,nsigma_aux
            DO j  = 1,nT_aux
              ij = ij+1
              input%T(ij) = T_aux(j)
              input%sigma(ij)     = sigma_aux(i)
            ENDDO
          ENDDO
          DEALLOCATE(T_aux, sigma_aux)
          !
        ELSE IF (TRIM(word3)=="list")THEN
          !
          IF(nconf<0) THEN
            IF(ionode) READ(input_unit,*, iostat=ios2) nconf
            CALL mpi_broadcast(ios2)
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting number of configs.", 1)
            CALL mpi_broadcast(nconf)
            IF(nconf<1) CALL errore("READ_INPUT","Wrong number of configs.", 1)
            input%nconf = nconf
          ENDIF
          !
          ALLOCATE(input%sigma(nconf), input%T(nconf))
          !
          ioWRITE(stdout,*) "Reading CONFIGS", nconf
          DO i = 1,nconf
            IF(ionode) READ(input_unit,'(a1024)', iostat=ios2) line
            CALL mpi_broadcast(ios2)
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration: input error.", 1)
            CALL mpi_broadcast(line)
            !
            ! try to read sigma and temperature
            READ(line,*,iostat=ios2) input%sigma(i), input%T(i)
            ! If it fails, read just temperature
            IF(ios2/=0) THEN
              READ(line,*,iostat=ios2) input%T(i)
              ! If this fails, complain
              IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration, got: '"//TRIM(line)//"'.", 1)
              ! reuse previous value of sigma if we read jus ttemperature
              IF(i>1) THEN
                input%sigma(i) = input%sigma(i-1)
              ELSE
                CALL errore("READ_INPUT","I need at least one value of sigma", 1)
              ENDIF
            ENDIF
          ENDDO
          !
        ELSE
          CALL errore("READ_INPUT","CONFIGS can be 'list' (default) or 'matrix'.", 1)
        ENDIF
        !
        ioWRITE(*,'(2x,a,/,100(8f9.1,/))') "Temperatures:", input%T
        ioWRITE(*,'(2x,a,/,100(8f9.3,/))') "Smearings:   ", input%sigma
        ioWRITE(*,*)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE ("ISOTOPES")
        IF(isotopes_ok) CALL errore("READ_INPUT", "Won't reads isotopes twice", 1)
        isotopes_ok = .true.
        !
        ioWRITE(*,*) "Reading ISOTOPES", S%ntyp
        IF(.not.isotopic_disorder) THEN
          ioWRITE(*,*) "WARNING! you did not set isotopic_disorder to true!"
        ENDIF
        ALLOCATE(auxs(S%ntyp), auxm(S%ntyp))
        !
        ISOTOPE_TYPE_LOOP : &
        DO i = 1,S%ntyp
          IF(ionode) READ(input_unit,'(a1024)', iostat=ios2) line
          CALL mpi_broadcast(ios2)
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope: input error.", 1)
          CALL mpi_broadcast(line)
          !
          ! Try to read isotope name and method
          READ(line,*,iostat=ios2) word2, word3
          !
          ! If this fails, complain
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope, got: '"//TRIM(line)//"'.", 1)
          IF(word2 /= S%atm(i)) THEN
            ioWRITE(*,"(a,i3,2a)") "WARNING: isotope name from input does not match FC file", i, word2, S%atm(i)
          ENDIF
          !
          !
          IF (word3=="natural") THEN
            CALL compute_gs(auxm(i), auxs(i), word2, 0, 0)
            
          ELSE IF (word3=="N") THEN
            READ(line,*,iostat=ios2) word2, word3, atomic_N
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope atomic number: input error.", 2)
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, 0)
            
          ELSE IF (word3=="M") THEN
            n_isotopes = 1
            ALLOCATE(isotopes_mass(n_isotopes), isotopes_conc(n_isotopes))
            isotopes_conc = 1._dp
            READ(line,*,iostat=ios2) word2, word3, isotopes_mass(1)
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope atomic mass: input error.", 3)
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, n_isotopes, isotopes_mass, isotopes_conc)
            DEALLOCATE(isotopes_mass, isotopes_conc)
            
          ELSE IF (word3=="isotopes") THEN
            READ(line,*,iostat=ios2) word2, word3, n_isotopes
            ALLOCATE(isotopes_mass(n_isotopes), isotopes_conc(n_isotopes))
            DO j = 1,n_isotopes
              IF(ionode) READ(input_unit,*, iostat=ios2) isotopes_mass(j), isotopes_conc(j)
              CALL mpi_broadcast(ios2)
              IF(ios2/=0) CALL errore("READ_INPUT","Reading isotope line: input error.", 4)
            ENDDO
            CALL mpi_broadcast(n_isotopes, isotopes_mass)
            CALL mpi_broadcast(n_isotopes, isotopes_conc)
            !
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, n_isotopes, isotopes_mass, isotopes_conc)
            DEALLOCATE(isotopes_mass, isotopes_conc)
          ELSE IF (word3=="manual") THEN
            READ(line,*,iostat=ios2) word2, word3, auxm(i), auxs(i)
              IF(ios2/=0) CALL errore("READ_INPUT","Reading isotope line: input error.", 5)
          ELSE
            CALL errore("READ_INPUT","I do not understand this isotope choice '"//TRIM(word2)//"'.", 1)
          ENDIF
        ENDDO &
        ISOTOPE_TYPE_LOOP 
        !
        ! Replace atomic masses with those from input:
        S%amass(1:S%ntyp) = auxm(1:S%ntyp) * MASS_DALTON_TO_RY
        S%amass_variance(1:S%ntyp) = auxs(1:S%ntyp)
        DEALLOCATE(auxs, auxm)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE DEFAULT
        IF(TRIM(line) /= '') THEN
          ioWRITE(*,*) "Skip:", TRIM(line)
        ENDIF
      END SELECT CARD_FOUND_CASE 
      !
      IF(ionode) READ(input_unit,'(a1024)', iostat=ios) line
      CALL mpi_broadcast(ios)
      CALL mpi_broadcast(line)
      word = ''
      !
    ENDDO &
    READ_CARDS
    !
    ! Setup a grid of q-points for which to compute the linewidth
    ! if no "Q_POINTS" section has been found, use the same grid
    ! as we used for integrating the phonon scattering
    ! (i.e. inner and outer grids are the same)
    IF(do_grid) THEN
      IF(.not.qpoints_ok .or. code=="TK") THEN
        nq1=nk(1); nq2=nk(2); nq3=nk(3)
        grid_type = input%grid_type
        xq0       = input%xk0
        qpoints_ok = .true.
      ENDIF
      !
      ioWRITE(*,*) "--> Setting up outer grid"
      CALL setup_grid(grid_type, S%bg, nq1, nq2, nq3, qpts, scatter=.false., xq0=xq0)
      input%prefix = TRIM(input%prefix)//&
                "."//TRIM(int_to_char(nq1))// &
                "x"//TRIM(int_to_char(nq2))// &
                "x"//TRIM(int_to_char(nq3))
      IF(LEN(TRIM(grid_type))==2) input%prefix=TRIM(input%prefix)//"@"//TRIM(grid_type)
      !DO i = 1,qpts%nq
      !  WRITE(*,'(i3,3f12.6,f16.3)') i, qpts%xq(:,i), qpts%w(i)
      !ENDDO
    ENDIF
    !
    ! Set natural isotope concentration for every species, if not read from input
    ! note that, if we're not doing isotope scattering, the values from input are taken.
    IF(input%isotopic_disorder.and..not.isotopes_ok)THEN
      ioWRITE(*,*) "Setting default isotopes"
      DO j = 1,S%ntyp
        CALL compute_gs(S%amass(j), S%amass_variance(j), S%atm(j), 0, 0)
      ENDDO
      S%amass = S%amass * MASS_DALTON_TO_RY
    ENDIF
    ! Finally, divide the FCs by the sqrt of these masses
    CALL div_mass_fc2(S, fc2)
    IF(present(fc3)) CALL fc3%div_mass(S)
    !
    IF(.not.qpoints_ok) CALL errore("READ_INPUT", "QPOINTS card not found &
                                    &(it may have been eaten by the preceeding card)", 1)
    IF(.not.configs_ok) CALL errore("READ_INPUT", "CONFIGS card not found &
                                    &(it may have been eaten by the preceeding card)", 1)
    !
    IF(input_unit/=5 .and. input_unit/=0) CLOSE(input_unit)
    !
!     WRITE(*,*) "=========== being input ==========="
!     WRITE(*,*) input
!     WRITE(*,*) "===========  end input  ==========="
    !
    timer_CALL t_iodata%stop()
    !
    CONTAINS 
    !
    SUBROUTINE broadcast_namelist_variables()
      USE mpi_thermal, ONLY : mpi_broadcast
      IMPLICIT NONE
      CALL mpi_broadcast(asr2)
      CALL mpi_broadcast(calculation)
      CALL mpi_broadcast(casimir_scattering)
      CALL mpi_broadcast(de)
      CALL mpi_broadcast(e0)
      CALL mpi_broadcast(sigma_e)
      CALL mpi_broadcast(nu_initial)
      CALL mpi_broadcast(e_initial)
      CALL mpi_broadcast(exp_t_factor)
      CALL mpi_broadcast(file_mat2)
      CALL mpi_broadcast(file_mat2_final)
      CALL mpi_broadcast(file_mat3)
      CALL mpi_broadcast(grid_type)
      CALL mpi_broadcast(grid_type_in)
      CALL mpi_broadcast(intrinsic_scattering)
      CALL mpi_broadcast(isotopic_disorder)
      CALL mpi_broadcast(max_seconds)
      CALL mpi_broadcast(max_time)
      CALL mpi_broadcast(mfp_cutoff)
      CALL mpi_broadcast(nconf)
      CALL mpi_broadcast(ne)
      CALL mpi_broadcast(niter_max)
      CALL mpi_broadcast(3,nk)
      CALL mpi_broadcast(3,nk_in)
      CALL mpi_broadcast(nq)
      CALL mpi_broadcast(outdir)
      CALL mpi_broadcast(prefix)
      CALL mpi_broadcast(print_dynmat)
      CALL mpi_broadcast(print_velocity)
      CALL mpi_broadcast(print_neutron_cs)
      CALL mpi_broadcast(3,q_initial)
      CALL mpi_broadcast(q_resolved)
      CALL mpi_broadcast(q_summed)
      CALL mpi_broadcast(restart)
      CALL mpi_broadcast(3,sample_dir)
      CALL mpi_broadcast(sample_length_au)
      CALL mpi_broadcast(sample_length_mm)
      CALL mpi_broadcast(sample_length_mu)
      CALL mpi_broadcast(sigmaq)
      CALL mpi_broadcast(sort_freq)
      CALL mpi_broadcast(skip_q)
      CALL mpi_broadcast(store_lw)
      CALL mpi_broadcast(thr_tk)
      CALL mpi_broadcast(volume_factor)
      CALL mpi_broadcast(3,xk0)
      CALL mpi_broadcast(3,xk0_in)
      CALL mpi_broadcast(optimize_grid)
      CALL mpi_broadcast(optimize_grid_thr)
    END SUBROUTINE
  END SUBROUTINE READ_INPUT
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_DATA(input, s, fc2, fc3)
    USE input_fc,           ONLY : same_system, read_fc2, aux_system, &
                                   forceconst2_grid, ph_system_info
    USE asr2_module,        ONLY : impose_asr2
    !USE io_global,          ONLY : stdout
    USE fc3_interpolate,    ONLY : read_fc3, forceconst3
    USE clib_wrappers,           ONLY : memstat
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)        :: input
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    CLASS(forceconst3),POINTER,OPTIONAL,INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(inout)   :: s
    TYPE(ph_system_info) :: s3
    !
    !INTEGER(kind=c_int) :: kb
    INTEGER :: kb
    !
      timer_CALL t_readdt%start()
    CALL read_fc2(input%file_mat2, S,  fc2)
    !print*, S
    IF(present(fc3)) THEN
      fc3 => read_fc3(input%file_mat3, S3)
      !
      IF(.not.same_system(S, S3)) THEN
        ioWRITE(stdout,*) "WARNING! FC2 and FC3 systems DO NOT MATCH !!!"
      ENDIF
    ENDIF
    !
    CALL aux_system(S)
    !
    CALL memstat(kb)
    ioWRITE(stdout,*) "Reading : done."
    ioWRITE(stdout,*) "Memory used : ", kb/1000, "Mb"
    !
    CALL impose_asr2(input%asr2, S%nat, fc2, S%zeu)
    ! NOTE: we now divide by the mass in READ_INPUT, as the masses
    !       read from input may be different (i.e. when including isotope scattering)
      timer_CALL t_readdt%stop()
    !
  END SUBROUTINE READ_DATA
  !
  SUBROUTINE parse_command_line(input_file)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(inout) :: input_file ! left unchanged if "-in XXXX" not found
    INTEGER :: count, i, ierr
    CHARACTER(len=512) :: argv
    count = command_argument_count()
    i = 0
    DO
      i = i +1
      IF(i>count) EXIT 
      CALL get_command_argument(i,argv)
      IF(argv=='-in')THEN
        i = i +1
        CALL get_command_argument(i,argv,status=ierr)
        IF(ierr>0) CALL errore('parse_command_line','cannot read input file name',1)
        IF(ierr<0) CALL errore('parse_command_line','input file name too long (max 512 chars)',1)
        input_file = argv
      ENDIF
    ENDDO
  END SUBROUTINE

END MODULE code_input
! <<^V^\\=========================================//-//-//========//O\\//
