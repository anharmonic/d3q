 ! This module contains subroutines to read QE output 
 ! read_md() -> read atomic positions, force and energy from the output of pw.x 
!               (single calculation, md run or even a cioncatenation of output files)
 !of molecular dynamics
! 10th Feb., 2021
!
!-------------------------------------------------------------------------
MODULE read_md_module
#include "mpi_thermal.h"
  !
 CONTAINS
  !
 SUBROUTINE new_sc(S,fc2,ta_sc)
  !-----------------------------------------------------------------------
  ! Make supercell
  !
  USE Kinds, ONLY    : DP
  USE input_fc, ONLY : read_fc2, forceconst2_grid, ph_system_info, &  
                       read_system, aux_system  ! use only: .. to know precisely what comes from where, a good practice
  !
  IMPLICIT NONE
  !
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc2
  INTEGER          :: i, j, k, jj, kk, nu, mu, nat_tot
  REAL(DP),ALLOCATABLE   :: ta_sc(:,:)
  !
  nat_tot = S%nat * fc2%nq(1) * fc2%nq(2) * fc2%nq(3)            
  ALLOCATE(ta_sc(3,nat_tot))
  OPEN(unit=112, file='new.txt', status='unknown')  
  k = 0   ! initialize
  DO j = 1, fc2%n_R
  DO i = 1, S%nat             
  k = k + 1     ! accummulate
  ta_sc(:,k) = S%tau(:,i)+fc2%xR(:,j)         
    ioWRITE(112,*) ta_sc(:,k)       
  ENDDO
  ENDDO
  CLOSE(112)
  !
  END SUBROUTINE
  !
  ! Generate supercell that correspond to force constants FC
  SUBROUTINE fc_to_supercell(S, fc2, aa, bb, nat_tot, tau_sc)
    USE input_fc,  ONLY : forceconst2_grid, ph_system_info
    USE kinds, ONLY : DP
    IMPLICIT NONE
    TYPE(ph_system_info)   :: S   ! system information 
    TYPE(forceconst2_grid) :: fc2 ! force constants describing the harmonic system  ALLOCATE(aa(3,3))
    REAL(DP),INTENT(out) :: aa(3,3), bb(3,3)
    INTEGER,INTENT(out)  :: nat_tot
    REAL(DP),INTENT(out),ALLOCATABLE :: tau_sc(:,:)
    INTEGER :: i,j,k, uni
    aa(:,1) = S%at(:,1)*fc2%nq(1)
    aa(:,2) = S%at(:,2)*fc2%nq(2)
    aa(:,3) = S%at(:,3)*fc2%nq(3)
    bb(:,1) = S%bg(:,1)/fc2%nq(1) ! vectors in reciprocal lattice of the sc 
    bb(:,2) = S%bg(:,2)/fc2%nq(2)
    bb(:,3) = S%bg(:,3)/fc2%nq(3)

    nat_tot = S%nat * fc2%nq(1) * fc2%nq(2) * fc2%nq(3)
    ! generate equilibrium positions of the atoms in bohr units

    OPEN(newunit=uni,file="sc.dat",action="WRITE",form="formatted")
    ALLOCATE(tau_sc(3,nat_tot))
    WRITE(uni,'("CELL_PARAMETERS bohr")') 
    WRITE(uni,'(3f14.8)') aa*S%alat
    WRITE(uni,'("ATOMIC_POSITIONS bohr")') 
    k = 0   ! initialize
    DO j = 1, fc2%n_R           !
      DO i = 1, S%nat           ! loop over cart. coord in file fc2
      k = k + 1                 ! 
      tau_sc(:,k) = (S%tau(:,i)+fc2%xR(:,j))*S%alat     !     
      WRITE(uni,'(a3,3f14.8)') S%atm(S%ityp(i)), tau_sc(:,k)
      ENDDO
    ENDDO
    CLOSE(uni)

  END SUBROUTINE 
  !

  SUBROUTINE read_pioud(file_tau, file_for, file_toten, toten0, nat_tot, alat, first_step, n_skip, n_steps,&
                    tau_md, force_md, toten_md, tau0, u_disp)
    !----------------------------------------------------------------------
    ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
    ! at finite temperature md.out. 
    !
    USE kinds,     ONLY : DP
    USE mpi_thermal, ONLY : my_id, num_procs, mpi_bsum
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*),INTENT(in) :: file_tau, file_for, file_toten
    REAL(DP),INTENT(in)  :: toten0 ! energy of unperturbed system, for reference
    INTEGER,INTENT(in)   :: nat_tot ! number of atoms to read
    REAL(DP),INTENT(in)  :: alat ! alat  in bohr and cell size in units of alat
    INTEGER,INTENT(in)   :: first_step, n_skip ! first step to read, number of steps to skip between two read
    INTEGER,INTENT(inout) :: n_steps         ! input: maximum number of steps to read, 
                                            ! output : number of steps actually read
    REAL(DP),ALLOCATABLE,INTENT(out)   :: tau_md(:,:,:), toten_md(:), force_md(:,:,:)
    REAL(DP),INTENT(in),OPTIONAL :: tau0(3,nat_tot) ! equilibrium atomic positions, only used to compute displacements
    REAL(DP),INTENT(out),ALLOCATABLE,OPTIONAL :: u_disp(:,:,:) ! displacements


  END SUBROUTINE

  SUBROUTINE read_md(md_file, toten0, nat_tot, alat, aa, first_step, n_skip, n_steps,&
                    tau_md, force_md, toten_md, tau0, u_disp)
  !----------------------------------------------------------------------
  ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
  ! at finite temperature md.out. 
  !
  USE kinds,     ONLY : DP
  ! USE input_fc,  ONLY : read_fc2, forceconst2_grid, ph_system_info, & 
  !                       read_system, aux_system ! use only: .. to know precisely what comes from where, a good practice
  USE mpi_thermal, ONLY : my_id, num_procs, mpi_bsum
  USE constants, ONLY : BOHR_RADIUS_ANGS
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*),INTENT(in) :: md_file
  REAL(DP),INTENT(in)  :: toten0 ! energy of unperturbed system, for reference
  INTEGER,INTENT(in)   :: nat_tot ! number of atoms to read
  REAL(DP),INTENT(in)  :: alat, aa(3,3) ! alat  in bohr and cell size in units of alat
  INTEGER,INTENT(in)   :: first_step, n_skip ! first step to read, number of steps to skip between two read
  INTEGER,INTENT(inout) :: n_steps         ! input: maximum number of steps to read, 
                                           ! output : number of steps actually read
  REAL(DP),ALLOCATABLE,INTENT(out)   :: tau_md(:,:,:), toten_md(:), force_md(:,:,:)
  REAL(DP),INTENT(in),OPTIONAL :: tau0(3,nat_tot) ! equilibrium atomic positions, only used to compute displacements
  REAL(DP),INTENT(out),ALLOCATABLE,OPTIONAL :: u_disp(:,:,:) ! displacements
  !
  REAL(DP),ALLOCATABLE :: vel_md(:,:,:) ! the velocity in bohr/time
  INTEGER :: i, j, k, nu, mu, beta, n_steps0, &
                      jat, kat, iat, ios, uni, uni_f, i_step, k_step
  REAL(DP)		:: dt, avg_x, avg_y, avg_z, avgsq_x, avgsq_y, avgsq_z
  CHARACTER(len=1024)    :: line
  CHARACTER(len=8)   :: dummy, dummyc, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
  LOGICAL,EXTERNAL   :: matches 
  LOGICAL      :: look_for_forces, look_for_toten_md, actually_read_the_forces, &
                  is_crystal_coords, is_alat_units, is_bohr_units, is_angstrom_units
  REAL(DP),PARAMETER :: ANGS_TO_BOHR = 1/BOHR_RADIUS_ANGS
  INTEGER :: first_step_=-1, n_skip_=-1, n_steps_=-1
  !
  !n_steps = 100
  IF(.not.ALLOCATED(force_md)) ALLOCATE(force_md(3,nat_tot,n_steps))
  ALLOCATE(tau_md(3,nat_tot,n_steps))
  ALLOCATE(toten_md(n_steps))

  OPEN(newunit=uni_f, file=TRIM(md_file)//".extract", action="READ",form="formatted")
  READ(uni_f, *, iostat=ios) dummy, first_step_, n_skip_, n_steps_
  IF(ios==0 .and. dummy == "EXTRACT:" .and. &
          first_step_==first_step .and. n_skip_==n_skip .and. n_steps_==n_steps) THEN
    ioWRITE(*,*) "NOTICE: Reading from "//TRIM(md_file)//".extract"
    uni = uni_f
    OPEN(newunit=uni_f, file="/dev/null", action="WRITE",form="formatted")
  ELSE
    CLOSE(uni_f)
    OPEN(newunit=uni,file=md_file,action="READ",form="formatted")
    OPEN(newunit=uni_f, file=TRIM(md_file)//".extract", action="WRITE",form="formatted")
  ENDIF
  WRITE(uni_f,*) "EXTRACT:", first_step, n_skip, n_steps
  !
  ! total number of steps to read
  n_steps0 = n_steps
  ! number of steps to read on this CPU
  n_steps = n_steps/num_procs
  IF( my_id<(n_steps0-n_steps*num_procs) ) n_steps = n_steps+1
  i_step = n_steps
  CALL mpi_bsum(i_step)
  IF(i_step/=n_steps0) CALL errore("md_read","parallel steps distribution not ok",1)
  !
  i_step = 0
  k_step = 0
  look_for_forces=.false.
  look_for_toten_md =.false.
  actually_read_the_forces = .false.
  ! 
  READ_LOOP : &
  DO
   READ(uni,"(a1024)", iostat=ios) line ! reads a line to variable "line"
    IF(ios /= 0 ) EXIT READ_LOOP      ! stop reading if there is an error
! READ atomic positions from initial configuration before md run 
! pw.x always prints initial positions in alat units    
    IF(matches("positions (alat units)",line))THEN
    !IF(     matches("positions (cryst. coord.)",line) &
    !   .or. matches("positions (alat units)",line) )THEN
      
      !is_crystal_coords = matches("positions (cryst. coord.)", line)

      k_step = k_step+1
      IF(k_step >= first_step+n_skip*my_id .and. MODULO(k_step-first_step, n_skip*num_procs) == 0)THEN
         actually_read_the_forces = .true.
         WRITE(uni_f,*) TRIM(line)
      ELSE
         actually_read_the_forces = .false.
         WRITE(uni_f,*) TRIM(line), " skip"
      ENDIF
      !
      IF(actually_read_the_forces)THEN
        IF(i_step >= n_steps) THEN
          !ioWRITE(*,*) "nstep reached... exiting"
          EXIT READ_LOOP
        ENDIF        
        i_step = i_step+1
        !print*, "Reading step", k_step, " as ", i_step
        WRITE(*,*) "Reading initial coords step", k_step, " as ", i_step, " on cpu", my_id
        look_for_toten_md = .true.
        look_for_forces = .true.
        DO k=1, nat_tot
          READ(uni,"(a1024)", iostat=ios) line
          WRITE(uni_f,*) TRIM(line)
          READ(line, *) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, tau_md(:,k, i_step)
        ENDDO

        !IF(is_crystal_coords)) CALL cryst_to_cart(fc2%n_R*S%nat,tau_md(:,:, i_step), aa, +1)
        tau_md(:,:, i_step) = tau_md(:,:, i_step)*alat
        !
       ENDIF
! READ atomic positions after md steps
    ELSEIF(matches("ATOMIC_POSITIONS",line) &
          .and. .not. (look_for_forces .or. look_for_toten_md))THEN
      
      ! Find out the unit of measure/coordinate system
      is_crystal_coords = matches("crystal",line)
      is_alat_units = matches("alat",line)
      is_bohr_units = matches("bohr",line)
      is_angstrom_units = matches("angstrom",line)
      IF(COUNT((/is_crystal_coords, is_alat_units, is_bohr_units, is_angstrom_units/))/=1) &
         CALL errore("read_md","unknown atomic positions format",1)

      IF(look_for_forces .or. look_for_toten_md) CALL errore("read_md","i found coordinates twice",1)

      k_step = k_step+1
      !IF(k_step >= first_step .and. MODULO(k_step, n_skip) == 0)THEN
      IF(k_step >= first_step+n_skip*my_id .and. MODULO(k_step-first_step-n_skip*my_id, n_skip*num_procs) == 0)THEN
         actually_read_the_forces = .true.
         WRITE(uni_f,*) TRIM(line)
      ELSE
         actually_read_the_forces = .false.
         WRITE(uni_f,*) TRIM(line), " skip"
      ENDIF
          
      IF(actually_read_the_forces)THEN
        IF(i_step >= n_steps) THEN
          !print*, "nstep reached... exiting"
          EXIT READ_LOOP
        ENDIF

        i_step = i_step + 1
        WRITE(*,*) "Reading step", k_step, " as ", i_step, " on cpu ", my_id

        look_for_toten_md = .true.  ! this...
        look_for_forces = .true. ! ..and this moved after the IF
        DO k=1, nat_tot
          READ(uni,"(a1024)", iostat=ios) line
          WRITE(uni_f,*) TRIM(line), actually_read_the_forces
          READ(line, *) dummyc, tau_md(:,k, i_step)
        ENDDO
        IF(is_crystal_coords) CALL cryst_to_cart(nat_tot,tau_md(:,:, i_step), aa, +1)
        IF(is_crystal_coords.or.is_alat_units) tau_md(:,:, i_step) = tau_md(:,:, i_step)*alat
        IF(is_angstrom_units)  tau_md(:,:, i_step) = tau_md(:,:, i_step)*ANGS_TO_BOHR
      ENDIF
! READ forces on atoms after md steps
    ELSE IF(matches("Forces acting on atoms",line) .and. look_for_forces) THEN 
      WRITE(uni_f,*) TRIM(line), actually_read_the_forces
      WRITE(uni_f,*)
      look_for_forces = .false.
      READ(uni,*)
      DO k=1, nat_tot
        READ(uni,"(a1024)", iostat=ios) line
        WRITE(uni_f,*) TRIM(line), actually_read_the_forces
        READ(line,*) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, force_md(:, k, i_step)
      ENDDO
! READ total energy after md steps
    ELSE IF(matches("!    total energy ",line) .and. look_for_toten_md) THEN
      look_for_toten_md = .false.
      WRITE(uni_f,*) TRIM(line), actually_read_the_forces
      READ(line,*)  dummy1, dummy2, dummy3, dummy4, toten_md(i_step)
    ENDIF
  ENDDO READ_LOOP

  CLOSE(uni)
  CLOSE(uni_f)

  !i_step = n_steps
  n_steps = i_step
  CALL mpi_bsum(i_step)
  ioWRITE(stdout,'(2x,a,i8)') "Total number of steps read among all CPUS", i_step

  IF(i_step<n_steps0)THEN
    ioWRITE(*,'(2x,a,i8,a,i8,a)') "NOTICE: Looking for ", n_steps0, "  steps, only", i_step, "  found."
    !n_steps = i_step
  ENDIF

  IF(look_for_toten_md .or. look_for_forces)THEN
    ! error
    CALL errore('compute_frc', 'did not find matching frc or nrg', 1)
  ENDIF
 !--------------------------------------------------------------------------------------
 ! compute F_aimd
  OPEN(441, file="F_aimd.dat", status="unknown")
  DO i_step = 1, n_steps
    DO k=1, nat_tot
      ioWRITE(441,'(3f14.9)') force_md(:,k,i_step)     
    ENDDO
  ioWRITE(441,*)
  ENDDO
  CLOSE(441)
 !
  !--------------------------------------------------------------------------------------
  OPEN(442, file="e_aimd.dat", status="unknown")
  DO i_step = 1, n_steps
    toten_md(i_step) = toten_md(i_step) - toten0
    ioWRITE(442,'(E14.6)') toten_md(i_step)
  !ioWRITE(442,*)
  ENDDO
  CLOSE(442)

  IF(present(u_disp) .and. present(tau0)) THEN
    IF(.not.ALLOCATED(u_disp)) ALLOCATE(u_disp(3,nat_tot,n_steps))
    !-------------------------------------------------------------------------------------
    !  Compute Displacement from Molecular Dynamics 
    !            Feb. 4th - 12th: 
    !PRINT*, "---------------------------------------------------------------------"
      i_step = 0
    OPEN(2244,file="new_disp-md.dat", status="unknown")
    DO i_step = 1, n_steps
      DO k=1, nat_tot
        ! ...cartesian components of displacement for each atom in supercell
        u_disp(1,k,i_step) = tau_md(1, k, i_step) - tau0(1, k) 
        u_disp(2,k,i_step) = tau_md(2, k, i_step) - tau0(2, k) 
        u_disp(3,k,i_step) = tau_md(3, k, i_step) - tau0(3, k)
        ioWRITE(2244,'(3(3f10.5,10x))') u_disp(:,k,i_step),tau_md(:, k, i_step),tau0(:, k)
        ! MSD = 1./n_steps*(DABS(u())**2
      ENDDO
      ioWRITE(2244,*)
      u_disp(:,:,i_step) = u_disp(:,:,i_step)
    ENDDO
    CLOSE(2244)
  ENDIF
  !
 !-------------------------------------------------------------------------------------
 !  Compute velocity from position of atoms tau_md()
 !            June 5th, 2021 
   i_step = 0
   dt = 0.4838
   !IF(present(vel_md))THEN
   ALLOCATE(vel_md(3,nat_tot,n_steps))
   OPEN(2240,file="velocity.dat", status="unknown")
   DO i_step = 1, n_steps
    DO k=1, nat_tot
      ! ...cartesian components of displacement for each atom in supercell
      IF(i_step==1)THEN
	      vel_md(:,k,i_step) = (tau_md(:, k, i_step) - tau_md(:, k, i_step+1))/dt
	    ELSE IF(i_step==n_steps)THEN
	      vel_md(:,k,i_step) = (tau_md(:, k, i_step) - tau_md(:, k, i_step-1))/dt
	    ELSE
        vel_md(:,k,i_step) = (tau_md(:, k, i_step-1) - tau_md(:, k, i_step+1))/(2*dt)
        !v(:,k,i_step) = v(:,k,i_step)
      ENDIF
      ioWRITE(2240,'(3f14.8,10x)') vel_md(:, k, i_step)

      ENDDO
      ioWRITE(2240,*)
    ENDDO
  CLOSE(2240)
  !
  !
  !----------------------------------------------------
  ! Compute standard deviation of atomic positions 13th July
  !
   i_step = 0
   OPEN(2245,file="stdev.m", status="unknown")
   DO i_step = 1, n_steps
    DO k=1, nat_tot
      tau_md(:, k, i_step) = tau_md(:, k, i_step) !*1.0_dp
    ENDDO
      !ioWRITE(2245,*)
      !
  ENDDO
  ! 
  !ioWRITE(2245,'(3f10.5,10x)') tau_md(:,1,n_steps)
  !
  avg_x = SUM(tau_md(1,1,:))/n_steps
  avg_y = SUM(tau_md(2,1,:))/n_steps
  avg_z = SUM(tau_md(3,1,:))/n_steps
  !
  avgsq_x = SUM(tau_md(1,1,:)**2)/n_steps
  avgsq_y = SUM(tau_md(2,1,:)**2)/n_steps
  avgsq_z = SUM(tau_md(3,1,:)**2)/n_steps
  
  !ioWRITE(2245,'(3f10.5,10x)') tau_md(:, 1, :) !/S%alat
  ioWRITE(2245,*) "===========STANDARD DEVIATION==========="
  ioWRITE(2245,'(3f10.5,10x)') SQRT(avgsq_x -avg_x**2 ), SQRT(avgsq_y -avg_y**2 ), SQRT(avgsq_z -avg_z**2 )
      
   CLOSE(2245)
   !
  !
 END SUBROUTINE
  !------------------------------------------------------------------------ 
  !
END MODULE
