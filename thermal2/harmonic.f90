 ! This module contains subroutines and functions for compute_force.90 code
! 10th Feb., 2021
!
!-------------------------------------------------------------------------
MODULE harmonic_module
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
  INTEGER          :: i, j, k, jj, kk, nu, mu, nat_sc
  REAL(DP),ALLOCATABLE   :: ta_sc(:,:)
  !
  nat_sc = S%nat * fc2%nq(1) * fc2%nq(2) * fc2%nq(3)            
  ALLOCATE(ta_sc(3,nat_sc))
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
 SUBROUTINE read_md(first_step, n_skip, n_steps,S,fc2,u_disp,force_sc)
  !----------------------------------------------------------------------
  ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
  ! at finite temperature md.out. 
  !
  USE kinds,     ONLY : DP
  USE input_fc,  ONLY : read_fc2, forceconst2_grid, ph_system_info, & 
                        read_system, aux_system ! use only: .. to know precisely what comes from where, a good practice
  USE mpi_thermal, ONLY : my_id, num_procs, mpi_bsum
  !
  IMPLICIT NONE
  !
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc2
  INTEGER,INTENT(in) :: first_step, n_skip
  INTEGER,INTENT(inout) :: n_steps
  INTEGER          :: i, j, k, nu, mu, beta, nat_sc, n_steps0, &
                      jat, kat, iat, ios, uni, i_step, k_step
  INTEGER, ALLOCATABLE   :: idx_R_map(:,:) 
  REAL(DP),ALLOCATABLE   :: tau_scmd(:,:,:), tot_ene(:), force_sc(:,:,:), &
          u_disp(:,:,:), tau_sc(:,:), aa(:,:), bb(:,:), &
          h_force(:,:,:)
  CHARACTER(len=1024)    :: line
  CHARACTER(len=8)   :: dummy, dummyc, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
  LOGICAL,EXTERNAL   :: matches 
  LOGICAL      :: look_for_forces, look_for_toten, actually_read_the_forces
  !
  !n_steps = 100
  ALLOCATE(force_sc(3,S%nat*fc2%n_R,n_steps))
  ALLOCATE(tau_scmd(3,S%nat*fc2%n_R,n_steps))
  ALLOCATE(u_disp(3,S%nat*fc2%n_R,n_steps))
  ALLOCATE(tot_ene(n_steps))

  ALLOCATE(aa(3,3))
  aa(:,1) = S%at(:,1)*fc2%nq(1)
  aa(:,2) = S%at(:,2)*fc2%nq(2)
  aa(:,3) = S%at(:,3)*fc2%nq(3)
  ALLOCATE(bb(3,3))
  bb(:,1) = S%bg(:,1)/fc2%nq(1) ! vectors in reciprocal lattice of the sc 
  bb(:,2) = S%bg(:,2)/fc2%nq(2)
  bb(:,3) = S%bg(:,3)/fc2%nq(3)

  nat_sc = S%nat * fc2%nq(1) * fc2%nq(2) * fc2%nq(3)
  ! generate equilibrium positions of the atoms in bohr units

  uni=1123
  OPEN(unit=uni,file="sc.dat",action="WRITE",form="formatted")
  ALLOCATE(tau_sc(3,nat_sc))
  WRITE(uni,'("CELL_PARAMETERS bohr")') 
  WRITE(uni,'(3f14.8)') aa*S%alat
  WRITE(uni,'("ATOMIC_POSITIONS bohr")') 
  k = 0   ! initialize
  DO j = 1, fc2%n_R           !
    DO i = 1, S%nat           ! loop over cart. coord in file fc2
    k = k + 1     ! accummulate
    tau_sc(:,k) = (S%tau(:,i)+fc2%xR(:,j))*S%alat     !     
    WRITE(uni,'(a3,3f14.8)') S%atm(S%ityp(i)), tau_sc(:,k)
    ENDDO
  ENDDO
  CLOSE(uni)


  uni = 1122
  OPEN(uni,file="md.out",action="READ",form="formatted")
  i_step = 0
  k_step = 0
  look_for_forces=.false.
  look_for_toten =.false.
  actually_read_the_forces = .false.

  
  READ_LOOP : &
  DO
   READ(uni,"(a1024)", iostat=ios) line ! reads a line to variable "line"
    IF(ios /= 0 ) EXIT READ_LOOP      ! stop reading if there is an error
! READ atomic poistions from initail configuration before md run 
    !IF(matches("positions (alat units)",line))THEN
    IF(matches("positions (cryst. coord.)",line))THEN
      k_step = k_step+1
      IF(k_step >= first_step+n_skip*my_id .and. MODULO(k_step-first_step, n_skip*num_procs) == 0)THEN
         actually_read_the_forces = .true.
      ELSE
         actually_read_the_forces = .false.
      ENDIF
      !
      IF(actually_read_the_forces)THEN
        IF(i_step >= n_steps) THEN
          print*, "nstep reached... exiting"
          EXIT READ_LOOP
        ENDIF        
        i_step = i_step+1
        !print*, "Reading step", k_step, " as ", i_step
        WRITE(*,*) "Reading initial coords step", k_step, " as ", i_step, " on cpu", my_id
        look_for_toten = .true.
        look_for_forces = .true.
        k=0
        DO j = 1, fc2%n_R
        DO iat = 1, S%nat
        k = k + 1
          READ(uni,"(a1024)", iostat=ios) line
          !print*, TRIM(line) ! remove spaces
          READ(line, *) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, tau_scmd(:,k, i_step)
        ENDDO
        ENDDO
        CALL cryst_to_cart(fc2%n_R*S%nat,tau_scmd(:,:, i_step), aa, +1)
        tau_scmd(:,:, i_step) = tau_scmd(:,:, i_step)*S%alat
       ENDIF
! READ atomic positions after md steps
    ELSEIF(matches("ATOMIC_POSITIONS",line) &
          .and. .not. (look_for_forces .or. look_for_toten))THEN
      
      IF(.not.matches("crystal",line)) CALL errore("read_md","can only read crystal coords",1)
      IF(look_for_forces .or. look_for_toten) CALL errore("read_md","i found coordinates twice",1)

      k_step = k_step+1
      !IF(k_step >= first_step .and. MODULO(k_step, n_skip) == 0)THEN
      IF(k_step >= first_step+n_skip*my_id .and. MODULO(k_step-first_step, n_skip*num_procs) == 0)THEN
         actually_read_the_forces = .true.
      ELSE
         actually_read_the_forces = .false.
      ENDIF
          
      IF(actually_read_the_forces)THEN
        IF(i_step >= n_steps) THEN
          print*, "nstep reached... exiting"
          EXIT READ_LOOP
        ENDIF

        i_step = i_step + 1
        WRITE(*,*) "Reading step", k_step, " as ", i_step, " on cpu ", my_id

        look_for_toten = .true.  ! this...
        look_for_forces = .true. ! ..and this moved after the IF
        k=0 
        DO j = 1, fc2%n_R
        DO iat = 1, S%nat
          k = k + 1
          READ(uni,"(a1024)", iostat=ios) line
          !print*, TRIM(line)
          READ(line, *) dummyc, tau_scmd(:,k, i_step)
        ENDDO
        ENDDO
        CALL cryst_to_cart(fc2%n_R*S%nat,tau_scmd(:,:, i_step), aa, +1)
        tau_scmd(:,:, i_step) = tau_scmd(:,:, i_step)*S%alat
      ENDIF
! READ forces on atoms after md steps
    ELSE IF(matches("Forces acting on atoms",line) .and. look_for_forces) THEN 
      look_for_forces = .false.
      READ(uni,*)
      k=0
      DO j = 1, fc2%n_R
      DO iat = 1, S%nat
        k=k+1
        READ(uni,"(a1024)", iostat=ios) line
        !print*, TRIM(line)
        READ(line,*) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, force_sc(:, k, i_step)
      ENDDO
      ENDDO
! READ total energy after md steps
    ELSE IF(matches("!    total energy ",line) .and. look_for_toten) THEN
      look_for_toten = .false.
      !print*, TRIM(line)
      READ(line,*)  dummy1, dummy2, dummy3, dummy4, tot_ene(i_step)
    ENDIF
  ENDDO READ_LOOP

  !i_step = n_steps
  n_steps = i_step
  CALL mpi_bsum(i_step)
  ioWRITE(stdout,'(2x,a,i8)') "Total number of steps read among all CPUS", i_step

  IF(i_step<n_steps0*num_procs)THEN
    ioWRITE(*,'(2x,a,i8,a,i3)') "Looking for ", n_steps, "I only found", i_step, " on cpu", my_id
    n_steps = i_step
  ENDIF

  IF(look_for_toten .or. look_for_forces)THEN
    ! error
    CALL errore('compute_frc', 'did not find matching frc or nrg', 1)
  ENDIF
 !--------------------------------------------------------------------------------------
 ! compute F_aimd
  OPEN(441, file="F_aimd.dat", status="unknown")
  DO i_step = 1, n_steps
     k=0 
      DO j = 1, fc2%n_R ! loop over all atoms in the super cell
      DO jat = 1, S%nat
        k=k+1
        ioWRITE(441,'(3f14.9)') force_sc(:,k,i_step)     
      ENDDO
      ENDDO
  ioWRITE(441,*)
  ENDDO
  CLOSE(441)
 !-------------------------------------------------------------------------------------
 !  Compute Displacement from Molecular Dynamics 
 !            Feb. 4th - 12th: 
 !PRINT*, "---------------------------------------------------------------------"
   i_step = 0
   OPEN(2244,file="new_disp-md.dat", status="unknown")
   DO i_step = 1, n_steps
      k=0
      DO j = 1, fc2%n_R
      DO iat = 1, S%nat
        ! ...cartesian components of displacement for each atom in supercell
        k=k+1
        u_disp(1,k,i_step) = tau_scmd(1, k, i_step) - tau_sc(1, k) 
        u_disp(2,k,i_step) = tau_scmd(2, k, i_step) - tau_sc(2, k) 
        u_disp(3,k,i_step) = tau_scmd(3, k, i_step) - tau_sc(3, k)
        ioWRITE(2244,'(3(3f10.5,10x))') u_disp(:,k,i_step),tau_scmd(:, k, i_step),tau_sc(:, k)
        ! MSD = 1./n_steps*(DABS(u())**2
      ENDDO
      ENDDO
      !CALL cryst_to_cart(fc2%n_R*S%nat, u_disp(:,:,i_step), aa, +1)
      ioWRITE(2244,*)
      !u_disp(:,:,i_step) = u_disp(:,:,i_step)/S%alat ! bohr units
      !u_disp(:,:,i_step) = u_disp(:,:,i_step)
    ENDDO
  !DEALLOCATE(tot_ene)
  !DEALLOCATE(tau_scmd)
  !DEALLOCATE(u_disp)
  !DEALLOCATE(aa)
  !DEALLOCATE(force_sc)
  CLOSE(2244)
  CLOSE(uni)
  !
  !
 END SUBROUTINE
  !------------------------------------------------------------------------
  !
  !
 SUBROUTINE harmonic_force(n_steps,S,fc2,u_disp,h_force)
  !----------------------------------------------------------------------
  ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
  ! at finite temperature md.out. 
  !
  USE Kinds, ONLY     : DP
  USE input_fc, ONLY  : read_fc2, forceconst2_grid, ph_system_info, & 
                          read_system, aux_system ! use only: .. to know precisely what comes from where, a good practice
  !USE harmonic_subroutines
  !
  IMPLICIT NONE
  !
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc2
  INTEGER          :: i, j, k, jj, kk, nu, mu, beta, nat_sc, &
          i_step, n_steps, jat, kat, cR(3)
  INTEGER, ALLOCATABLE,SAVE   :: idx_R_map(:,:) 
  REAL(DP),ALLOCATABLE   :: h_force(:,:,:), u_disp(:,:,:) !, force_sc(:,:,:)

  IF(.NOT.ALLOCATED(h_force)) ALLOCATE(h_force(3,S%nat*fc2%n_R,n_steps))

  !n_steps = 100
  IF(.not.ALLOCATED(idx_R_map)) THEN
    !
    ALLOCATE(idx_R_map(fc2%n_R,fc2%n_R))
    !IF(.not.ALLOCATED(f_harm)) ALLOCATE(f_harm(....))

    ! Build map (i,j)->i such that 
    !    (R_j-R_k) = (R_i + RR),
    ! for R_i, R_j and R_k unit-cell vectors in the super-cell and RR a super-lattice vector
    idx_R_map = -1
    !fc2%i_0 
    DO j = 1, fc2%n_R
    DO k = 1, fc2%n_R
      cR = fc2%yR(:,j) - fc2%yR(:,k) ! R~ = R-R'
      cR(1) = MODULO( cR(1), fc2%nq(1))
      cR(2) = MODULO( cR(2), fc2%nq(2))
      cR(3) = MODULO( cR(3), fc2%nq(3))   ! R^ is R~ but taken inside the grid
      DO i = 1, fc2%n_R
        IF(ALL(cR==fc2%yR(:,i))) THEN
          ! is the index that we are looking for
          idx_R_map(j,k) = i
          EXIT ! <- stop the loop and exit
        ENDIF
      ENDDO
      !WRITE(666,'(3(3i4,5x))') fc2%yR(:,j), fc2%yR(:,k), fc2%yR(:,idx_R_map(j,k))
    IF(idx_R_map(j,k)==-1) CALL errore("harm_force", "could not find some R,R'", 1)
    ENDDO
    ENDDO
    !
  ENDIF

  ! compute force
  ! OPEN(331, file="F_harm.dat", status="unknown")
   h_force(:,:,:) = 0._dp
   DO i_step = 1, n_steps
     jj=0 
     DO j = 1, fc2%n_R  ! loop over all atoms in the super cell
     DO jat = 1, S%nat
        jj=jj+1      
        nu = (jat-1)*3
        kk=0
        DO k = 1, fc2%n_R   !loop over index of vector of the cell in the supercell
        DO kat = 1, S%nat ! ... over all atoms in the unit cell
        kk=kk+1
        DO beta = 1,3 ! ... over cartesian axes x,y,z for each atom
          mu = (kat-1)*3 + beta
          i = idx_R_map(j,k)
          h_force(:,jj, i_step) = h_force(:,jj, i_step) &
                                - fc2%FC(nu+1:nu+3,mu,i) * u_disp(:,kk,i_step)
        ENDDO
        ENDDO
        ENDDO
        !WRITE(331,'(3f14.9)') h_force(:,jj,i_step) !, force_sc(:,jj,i_step) !/force(:,jj,i_step)
     END DO
     END DO
     !WRITE(331,*)
    ENDDO
 !  CLOSE(331)
 !
 END SUBROUTINE 
 !
END MODULE
