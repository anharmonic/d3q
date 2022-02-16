!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE tdph_module
  USE kinds, ONLY : DP
#include "mpi_thermal.h"
  !
  INTEGER :: nfar=0

  TYPE tdph_input_type
    !
    CHARACTER(len=256) :: ai = 'md'
    CHARACTER(len=256) :: fmd = 'md.out'
    CHARACTER(len=256) :: ftau, fforce, ftoten
    CHARACTER(len=256) :: file_mat2 = 'mat2R.periodic'
    CHARACTER(len=8) :: fit_type = "force"
    CHARACTER(len=8) :: minimization = "acrs"
    CHARACTER(len=9) :: basis = "mu"
    INTEGER :: nfirst, nskip, nmax, nprint
    REAL(DP) :: e0, thr, T, randomization
    !
  END TYPE tdph_input_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT_TDPH(input)
    USE code_input,           ONLY : parse_command_line
    USE cmdline_param_module, ONLY : cmdline_to_namelist
    !
    IMPLICIT NONE
    !
    TYPE(tdph_input_type),INTENT(out) :: input
    CHARACTER(len=256)  :: input_file
    INTEGER :: ios
    !
    ! Input variable, and defaul values:
    CHARACTER(len=256) :: ai  = 'md'
    CHARACTER(len=256) :: fmd = 'md.out'
    CHARACTER(len=256) :: ftau   = 'positions.dat'
    CHARACTER(len=256) :: fforce = 'forces.dat'
    CHARACTER(len=256) :: ftoten = 'pioud.dat'
    CHARACTER(len=256) :: file_mat2 = 'mat2R.periodic'
    CHARACTER(len=8) :: fit_type = "force"
    CHARACTER(len=8) :: minimization = "acrs"
    CHARACTER(len=9) :: basis = "mu"
    INTEGER :: nfirst=1000, nskip=100, nmax=-1, nprint=1000, nread=-1

    INTEGER :: input_unit, aux_unit, err1, err2
    !CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER,EXTERNAL :: find_free_unit
    REAL(DP) :: e0 = 0._dp, thr = 1.d-8, T=-1._dp, randomization=0._dp
    !
    NAMELIST  / tdphinput / &
        fmd, fforce, ftau, ftoten, &
        ai, file_mat2, fit_type, &
        nfirst, nskip, nmax, nprint, nread, &
        e0, minimization, thr, T, basis, &
        randomization

    ioWRITE(*,*) "Waiting for input"
    !
    IF(ionode)THEN
      input_file="input.TDPH"
      CALL parse_command_line(input_file)
      IF(TRIM(input_file)=="-")THEN
        ioWRITE(stdout,'(2x,3a)') "Warning! Reading standard input will probably not work with MPI"
        input_unit = 5
        err1=0
      ELSE
        ioWRITE(stdout,'(2x,3a)') "Reading input file '", TRIM(input_file), "'"
        !input_unit = find_free_unit()
        OPEN(newunit=input_unit, file=input_file, status="OLD", action="READ",iostat=err1)
      ENDIF
      !
      IF(err1==0)THEN
        aux_unit = find_free_unit()
        READ(input_unit, tdphinput, iostat=err1)
        ioWRITE(stdout,'(2x,3a)') "merging with command line arguments"
      ELSE
        !ioWRITE(stdout,'(2x,3a)') "no input file, trying with command line arguments"
        CALL errore("tdph", "bad input file", 1)
      ENDIF
      OPEN(unit=aux_unit, file=TRIM(input_file)//".tmp~", status="UNKNOWN", action="READWRITE")
      CALL cmdline_to_namelist("tdphinput", aux_unit)
      REWIND(aux_unit)
      READ(aux_unit, tdphinput, iostat=err2)
      IF(err2==0)CLOSE(aux_unit, status="DELETE")
      
      IF(err1/=0 .and. err2/=0 .and. ionode)  WRITE(stdout,'(2x,3a)') "Warning: no input file or parameters!"

      ioWRITE(stdout, tdphinput)
    ENDIF
    !
    CALL bcast_namelist_variables()
    !
    !IF(ANY(<1))  CALL errore("READ_INPUT_TDPH","Invalid nk",1)

    IF(e0==0._dp .and. fit_type=='energy') &
        CALL errore("tdph","need zero energy to fit energy difference", 1)

    IF(nmax>0 .and. nread>0) CALL errore("tdph", "You cannot specify both nread and nmax",1)
    IF(nmax<0 .and. nread>0) nmax = nfirst+nskip*nread
    IF(nmax<0 .and. nread<0) nmax = 50000 ! default value when nothing specified
    IF(ANY((/nfirst,nmax,nskip/)<1)) &
        CALL errore("tdph","wrong parameters", 1)
    !
    input%ai            = ai
    input%fmd           = fmd
    input%ftau          = ftau
    input%fforce        = fforce
    input%ftoten        = ftoten
    input%file_mat2     = file_mat2
    input%fit_type      = fit_type
    input%minimization  = minimization
    input%nfirst        = nfirst
    input%nskip         = nskip
    input%nmax          = nmax
    input%nprint        = nprint
    input%e0            = e0
    input%thr           = thr
    input%T             = T
    input%basis         = basis
    input%randomization = randomization
    !
    CONTAINS 
    SUBROUTINE bcast_namelist_variables()
        USE mpi_thermal, ONLY : mpi_broadcast
        IMPLICIT NONE
        CALL mpi_broadcast(ai)
        CALL mpi_broadcast(fmd)
        CALL mpi_broadcast(ftau)
        CALL mpi_broadcast(fforce)
        CALL mpi_broadcast(ftoten)
        CALL mpi_broadcast(file_mat2)
        CALL mpi_broadcast(fit_type)
        CALL mpi_broadcast(nfirst)
        CALL mpi_broadcast(nskip)
        CALL mpi_broadcast(nmax)
        CALL mpi_broadcast(nprint)
        CALL mpi_broadcast(nread)
        CALL mpi_broadcast(e0)
        CALL mpi_broadcast(minimization)
        CALL mpi_broadcast(thr)
        CALL mpi_broadcast(T)
        CALL mpi_broadcast(basis)
        CALL mpi_broadcast(randomization)
    END SUBROUTINE
  END SUBROUTINE READ_INPUT_TDPH
  !
  SUBROUTINE set_qe_global_geometry(Si)
    USE input_fc,           ONLY : ph_system_info
    USE cell_base,          ONLY : at, bg, celldm, ibrav, omega
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
    USE noncollin_module,   ONLY : m_loc, nspin_mag
    !
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in) :: Si
    !
    ! Quantum-ESPRESSO symmetry subroutines use the global variables
    ! we copy the system data from structure S
    ntyp   = Si%ntyp
    nat    = Si%nat
    IF(allocated(tau)) DEALLOCATE(tau)
    ALLOCATE(tau(3,nat))
    IF(allocated(ityp)) DEALLOCATE(ityp)
    ALLOCATE(ityp(nat))
    celldm = Si%celldm
    at     = Si%at
    bg     = Si%bg
    omega  = Si%omega
    atm(1:ntyp)    = Si%atm(1:ntyp)
    amass(1:ntyp)  = Si%amass(1:ntyp)
    tau(:,1:nat)   = Si%tau(:,1:nat)
    ityp(1:nat)    = Si%ityp(1:nat)
  
  END SUBROUTINE

END MODULE

!
!----------------------------------------------------------------------------
PROGRAM tdph
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : amu_ry, K_BOLTZMANN_SI, K_BOLTZMANN_RY, RYTOEV !ry_to_kelvin
  USE parameters,         ONLY : ntypx
  !USE mp,                 ONLY : mp_bcast
  !USE mp_global,          ONLY : mp_startup, mp_global_end
  !USE mp_world,           ONLY : world_comm
  !USE environment,        ONLY : environment_start, environment_end
  ! symmetry
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, &
                                  irt, copy_sym, nrot, inverse_s, t_rev
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq, invsymq
  USE control_lr,         ONLY : lgamma
  USE decompose_d2,       ONLY : smallg_q_fullmq, find_d2_symm_base, sym_and_star_q, dotprodmat, &
                                 make_qstar_d2, allocate_sym_and_star_q, tr_star_q, dynmat_basis
  USE cmdline_param_module
  USE input_fc,           ONLY : forceconst2_grid, ph_system_info, read_system, aux_system, read_fc2, &
                                 div_mass_fc2, multiply_mass_dyn, write_fc2
  USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE asr2_module,        ONLY : impose_asr2
  USE quter_module,       ONLY : quter
  USE mpi_thermal,        ONLY : start_mpi, stop_mpi, num_procs, mpi_broadcast
  USE tdph_module
  ! harmonic
  USE read_md_module,     ONLY : read_md, read_pioud, fc_to_supercell
  USE harmonic_module,    ONLY : harmonic_force_md
  ! minipack
  USE lmdif_module,       ONLY : lmdif0
  USE lmdif_p_module,     ONLY : lmdif_p0
  USE timers
  USE random_numbers, ONLY : randy
  !
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: CODE="TDPH"
  CHARACTER(len=256) :: fildyn, filout
  INTEGER :: ierr, nargs
  !
  INTEGER       :: nq1, nq2, nq3, nqmax, nq_wedge, nqq, nq_done

  REAL(DP),ALLOCATABLE      :: x_q(:,:), w_q(:)
  ! for harmonic_module
  REAL(DP),ALLOCATABLE      :: u(:,:,:), force_harm(:,:,:), force_md(:,:,:), tau_md(:,:,:), tau_sc(:,:), &
                               force_diff(:), diff_tot(:), toten_md(:), h_energy(:)
  ! for lmdf1
  INTEGER                   :: n_steps, n_steps_tot, j_steps, nat_sc, ulog
  INTEGER :: mdata, mdata_tot !, mfcn
  !
  REAL(DP) :: xq(3), syq(3,48), at_sc(3,3), bg_sc(3,3), force_ratio(3), aux
  LOGICAL :: sym(48), lrigid_save, skip_equivalence, time_reversal
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:,:,:,:), d2(:,:), w2(:,:), &
                             star_wdyn(:,:,:,:, :), star_dyn(:,:,:)
  REAL(DP),ALLOCATABLE :: decomposition(:), xqmax(:,:), ph_coefficients(:), &
                          ph_coefficients0(:), metric(:)
  INTEGER :: i,j, icar,jcar, na,nb, iq, nph, iph, iswitch, first_step, n_skip
  INTEGER,ALLOCATABLE :: rank(:)
  TYPE(ph_system_info) :: Si
  TYPE(forceconst2_grid) :: fc, fcout
  TYPE(dynmat_basis),ALLOCATABLE :: dmb(:)
  TYPE(sym_and_star_q),ALLOCATABLE :: symq(:)
  TYPE(tdph_input_type) :: input
  !
  CALL start_mpi()
  CALL remove_stack_limit()
  !
  ! Read namelist tdphinput
  CALL READ_INPUT_TDPH(input)

  CALL read_fc2(input%file_mat2, Si, fc)
  lrigid_save = Si%lrigid
  !Si%lrigid = .false.
  CALL impose_asr2("simple", Si%nat, fc, Si%zeu)
  CALL aux_system(Si)
  CALL div_mass_fc2(Si, fc)
  !
  CALL set_qe_global_geometry(Si)
  !
  ! ######################### symmetry setup #########################
  ! Symmetry setup uses global variable, which is not ideal, but
  ! not a huge problem since we only need to do this part once
  ! at the beginning.
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
  CALL set_sym_bl ( )
  ioWRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  IF(.not.allocated(m_loc))  THEN
    ALLOCATE(m_loc(3,Si%nat))
    m_loc = 0._dp
  ENDIF
  
  CALL find_sym ( Si%nat, Si%tau, Si%ityp, .false., m_loc )
  ioWRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !
  ! Find the reduced grid of q-points:
  skip_equivalence = .FALSE.
  time_reversal    = .TRUE.
  nq1 = fc%nq(1)
  nq2 = fc%nq(2)
  nq3 = fc%nq(3)
  nqmax = nq1*nq2*nq3
  ALLOCATE(x_q(3,nqmax), w_q(nqmax))
  CALL kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, Si%bg, nqmax,&
                    0,0,0, nq1,nq2,nq3, nq_wedge, x_q, w_q )
  !
  ioWRITE(stdout, *) "Generated ", nq_wedge, "points"
  
  ALLOCATE(rtau( 3, 48, Si%nat), d2(3*Si%nat,3*Si%nat))

  ! Variable to hold the dyn matrix and q-points of the entire grid
  nq_done = 0
  ALLOCATE(star_wdyn(3,3,Si%nat,Si%nat, nqmax))
  ALLOCATE(xqmax(3,nqmax))

  ! For every q-point in the irreducible wedge, we find its symmetry
  ! and the basis of the space of symmetry-constrained dynfactoramical matrices
  ! Again, this part uses gloabl variables which is a annoying, but once
  ! we have the set of basis matrices, we don't have to touch it
  ! anymore.
  ALLOCATE(dmb(nq_wedge))
  ALLOCATE(rank(nq_wedge))
  ALLOCATE(symq(nq_wedge))

  CALL t_init%start()

  Q_POINTS_LOOP : &
  DO iq = 1, nq_wedge
    ioWRITE(stdout, *) "____[[[[[[[", iq, "]]]]]]]]____"
    ioWRITE(stdout, '(i6, 3f12.4)') iq, x_q(:,iq)
    !
    ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
    ! part 1: call smallg_q and the copy_sym, 
    xq = x_q(:,iq)
    minus_q = .true.
  
    sym = .false.
    sym(1:nsym) = .true.
    CALL smallg_q_fullmq(xq, 0, Si%at, Si%bg, nsym, s, sym, minus_q)
    nsymq = copy_sym(nsym, sym)
    ! recompute the inverses as the order of sym.ops. has changed
    CALL inverse_s ( ) 
  
    ! part 2: this computes gi, gimq
    call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!    WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
!    IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G"
    !
    ! finally this does some of the above again and also computes rtau...
    CALL sgam_lr(Si%at, Si%bg, nsym, s, irt, Si%tau, rtau, Si%nat)
    !
    ! Now, I copy all the symmetry definitions to a derived type
    ! which I'm going to use EXCLUSIVELY from here on
    CALL allocate_sym_and_star_q(Si%nat, symq(iq))
    symq(iq)%xq = xq
    symq(iq)%nrot = nrot
    symq(iq)%nsym  = nsym
    symq(iq)%nsymq = nsymq
    symq(iq)%minus_q = minus_q
    symq(iq)%irotmq = irotmq
    symq(iq)%s = s
    symq(iq)%invs = invs
    symq(iq)%rtau = rtau
    symq(iq)%irt= irt

    !integer :: nrot, nsym, nsymq, irotmq
    !integer :: nq, nq_tr, isq (48), imq
    !real(DP) :: sxq (3, 48)
    !
    ! the next subroutine uses symmetry from global variables to find the basis of crystal-symmetric
    ! matrices at this q point
    CALL fftinterp_mat2(xq, Si, fc, d2)
    d2 = multiply_mass_dyn(Si,d2)

    IF(ionode)THEN
      CALL find_d2_symm_base(xq, rank(iq), dmb(iq)%basis, &
         Si%nat, Si%at, Si%bg, symq(iq)%nsymq, symq(iq)%minus_q, &
         symq(iq)%irotmq, symq(iq)%rtau, symq(iq)%irt, symq(iq)%s, symq(iq)%invs, &
         d2, input%basis)
    ENDIF
    ! 
    CALL t_comm%start()
    call mpi_broadcast(rank(iq))
    IF(.not.ionode) ALLOCATE(dmb(iq)%basis(Si%nat3, Si%nat3, rank(iq)))
    CALL mpi_broadcast(Si%nat3, Si%nat3, rank(iq), dmb(iq)%basis)
    CALL t_comm%stop()
      !
    ! Calculate the list of points making up the star of q and of -q
    CALL tr_star_q(symq(iq)%xq, Si%at, Si%bg, symq(iq)%nsym, symq(iq)%s, symq(iq)%invs, &
                   symq(iq)%nq_star, symq(iq)%nq_trstar, symq(iq)%sxq, &
                   symq(iq)%isq, symq(iq)%imq, .false. )

    ioWRITE(stdout, '(5x,a,2i5)') "Found star of q and -q", symq(iq)%nq_star, symq(iq)%nq_trstar
    syq = symq(iq)%sxq
    call cryst_to_cart(symq(iq)%nq_trstar, syq, Si%at, -1)
    DO i = 1, symq(iq)%nq_trstar
       syq(1,i) = MODULO(syq(1,i), 1._dp)
       syq(2,i) = MODULO(syq(2,i), 1._dp)
       syq(3,i) = MODULO(syq(3,i), 1._dp)
       syq(:,i) = syq(:,i) * (/nq1, nq2, nq3/)
       ioWRITE(stdout,'(i4,3i3,l2)') i, NINT(syq(:,i)), (i>symq(iq)%nq_star)
    ENDDO

  ENDDO Q_POINTS_LOOP
  !
  ! Number of degrees of freedom for the entire grid:
  nph = SUM(rank)
  ioWRITE(stdout, '("=====================")')
  ioWRITE(stdout, '(5x,a,2i5)') "TOTAL number of degrees of freedom", nph  
  
  ! Allocate a vector to hold the decomposed phonons over the entire grid
  ! I need single vector in order to do minimization, otherwise a derived
  ! type would be more handy
  ALLOCATE(ph_coefficients(nph), ph_coefficients0(nph))
  iph = 0
  Q_POINTS_LOOP2 : &
  DO iq = 1, nq_wedge
    xq = symq(iq)%xq
    !IF(iq==1) xq=xq+1.d-6
    !
    ! Interpolate the system dynamical matrix at this q
    CALL fftinterp_mat2(xq, Si, fc, d2)
    ! Remove the mass factor, I cannot remove it before because the effective
    ! charges/long range interaction code assumes it is there
    d2 = multiply_mass_dyn(Si,d2)
    !
    ! Decompose the dynamical matrix over the symmetric basis at this q-point
    ioWRITE(stdout,'(2x,a)') "== DECOMPOSITION =="
    DO i = 1,rank(iq)
      iph = iph +1
      ph_coefficients(iph) = dotprodmat(3*Si%nat,d2, dmb(iq)%basis(:,:,i))
      ioWRITE(stdout,"(i3,1f12.6)") i, ph_coefficients(iph)
    ENDDO
    !
  ENDDO Q_POINTS_LOOP2
  !
  ph_coefficients0 = ph_coefficients
!-----------------------------------------------------------------------
  ! Variables that can be adjusted according to need ...
  !
  n_steps = input%nmax ! total molecular dynamics steps TO READ
  first_step = input%nfirst ! start reading from this step
  n_skip = input%nskip !        ! number of steps to skip

  CALL fc_to_supercell(Si, fc, at_sc, bg_sc, nat_sc, tau_sc)
  CALL t_init%stop()
!###################  end of initialization ####################################################

  ALLOCATE(tau_md(3,nat_sc,n_steps))
  ALLOCATE(force_md(3,nat_sc,n_steps))
  ALLOCATE(toten_md(n_steps))
  CALL t_read%start()
  IF(input%ai=="md") THEN
  CALL read_md(input%fmd, input%e0, nat_sc, Si%alat, at_sc, first_step, n_skip, n_steps, &
               n_steps_tot, tau_md, force_md, toten_md, tau_sc, u)
  ELSE IF(input%ai=="pioud")THEN
    CALL read_pioud(input%ftau, input%fforce, input%ftoten, input%e0, nat_sc, Si%alat, &
                    first_step, n_skip, n_steps, n_steps_tot, tau_md, force_md, toten_md, tau_sc, u)
  ELSE
    CALL errore("tdph","unknown input format", 1)
  ENDIF
  CALL t_read%stop()

!###################  end of data input ####################################################

  mdata = 3*nat_sc*n_steps
  mdata_tot = 3*nat_sc*n_steps_tot
  ALLOCATE(force_diff(3*nat_sc))
  ALLOCATE(diff_tot(mdata_tot))
  ALLOCATE(force_harm(3,nat_sc,n_steps))
  nfar = 0
  
  IF(input%randomization/=0._dp)THEN
    IF(ionode)THEN
      !
      aux = 2*SUM(ABS(ph_coefficients))/nph*ABS(input%randomization)
      !
      IF(input%randomization>0._dp)THEN
        WRITE(*,"(2x,a,f12.6)") "Adding random noise to initial phonons up to ", aux
        DO i = 1, nph
          ph_coefficients(i) = ph_coefficients(i) + aux*(randy()-.5_dp)
        ENDDO
      ELSE IF(input%randomization<0._dp)THEN
        WRITE(*,"(2x,a,f12.6)") "Randomizing initial phonons up to ", aux
        DO i = 1, nph
          ph_coefficients(i) = aux*(randy()-.5_dp)
        ENDDO
      ENDIF
      !
    ENDIF
    ! be sure that every CPU has the same coefficients
    CALL mpi_broadcast(nph, ph_coefficients)
  ENDIF

  ! FIXME: Compute and save to file the forces before minimization, should be removed
  iswitch = 0
  CALL chi_lmdif(mdata_tot, nph, ph_coefficients, diff_tot, iswitch)
  OPEN(118,file="h_enr.dat0",status="unknown") 
  DO i = 1, n_steps
  !WRITE(117,*) "i_step = ",i
        WRITE(118,'(E16.8)') h_energy(i) !
  END DO
  CLOSE(118)
  !  -- end of FIXME


  !ph_coefficients(3) =   ph_coefficients(3) *0.0_dp
  ulog=1119
  OPEN(unit=ulog, file="tdph.log", form="formatted", status="unknown")

  CALL t_minim%start()
  SELECT CASE (input%minimization)
  CASE("acrs")
    CALL ACRS0(nph,ph_coefficients, force_diff, chisq_acrs)
  CASE("lmdif")
    CALL lmdif_p0(chi_lmdif, mdata_tot, nph, ph_coefficients, diff_tot, input%thr, iswitch)
  CASE("none")
    ! Do nothing
  CASE DEFAULT
    CALL errore("tdph","unknown minimization engine requested",1)
  END SELECT
  CALL t_minim%stop()
  !
!###################  end of minimization ####################################################
! Now write to file the file matrices in "centered" and "periodic" form

  nfar = 0
  iswitch = 0
  CALL chi_lmdif(mdata_tot, nph, ph_coefficients, diff_tot, iswitch)
  Si%lrigid = lrigid_save
  CALL write_fc2("matOUT.periodic", Si, fcout)

  ! Force ratio
  OPEN(116,file="force_ratio.dat",status="unknown") 
  !
  DO i = 1, n_steps
   WRITE(116,*) "i_step = ",i
   DO j = 1, nat_sc
     WHERE(force_md(:,j,i)/=0._dp) 
       force_ratio = force_harm(:,j,i)/force_md(:,j,i)
     ELSEWHERE
       force_ratio = 0._dp
     END WHERE
 
     ioWRITE(116,'(1f14.6,5x,3(3f14.6,5x))') DSQRT(norm2(force_ratio)),&
                             force_ratio, force_harm(:,j,i), force_md(:,j,i)
     END DO
  END DO
  CLOSE(116)
  ! 
  ! Harmonic energy
  ! 
  OPEN(117,file="h_enr.dat",status="unknown") 
  DO i = 1, n_steps
  !WRITE(117,*) "i_step = ",i
    ioWRITE(117,'(i5,3(E13.6, 3x))') i, h_energy(i), toten_md(i), &
                                         EXP(-toten_md(i)/(K_BOLTZMANN_RY*300.0_DP))
  END DO
  CLOSE(117)
   !
  nfar = 2
  CALL chi_lmdif(mdata_tot, nph, ph_coefficients, diff_tot, iswitch)
  CALL write_fc2("matOUT.centered", Si, fcout)

  CLOSE(ulog)
  ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
  CALL print_timers_header()
  CALL t_read%print()
  CALL t_init%print()
  CALL t_minim%print()
  ioWRITE(stdout,'("   * inside minimization ")')
  CALL t_chi2%print()
  ioWRITE(stdout,'("   * inside chi2 ")')
  CALL t_comm%print()
  CALL t_force%print()
  CALL t_recom%print()

  CALL stop_mpi()

  !write fcout to file (final result)
 ! ---- the program is over ---- !
 !
 CONTAINS

! Penalty function in the form required by LMDIF : return an array of size
! at least as large as the number of degrees of freedom containing the penalty
! (NOT SQUARED) for each dimension. 
!-----------------------------------------------------------------------
 SUBROUTINE chi_lmdif(mdata_tot, nph, ph_coef, diff_tot, iswitch)
  !-----------------------------------------------------------------------
  ! Calculates the square difference, fdiff2, btw harmonic and ab-initio
  ! forces for n_steps molecur dyanmics simulation 
  !
  USE tdph_module,  ONLY : nfar
  USE mpi_thermal,  ONLY : mpi_bsum, allgather_vec, my_id
  USE decompose_d2, ONLY : recompose_fc
  IMPLICIT NONE
  INTEGER,INTENT(in)    :: mdata_tot, nph
  REAL(DP),INTENT(in)   :: ph_coef(nph)
  REAL(DP),INTENT(out)  :: diff_tot(mdata_tot)
  REAL(DP)  :: diff(mdata)
  INTEGER,INTENT(inout) :: iswitch

  INTEGER :: nq_done, iph, iq, i, j, k
  INTEGER,SAVE :: iter = 0
  CHARACTER (LEN=6),  EXTERNAL :: int_to_char
  REAL(DP) :: chi2, e0
  REAL(DP),SAVE :: last_save = 0._dp

  CALL t_chi2%start()

  CALL recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_coef,&
                    nq1, nq2, nq3, nqmax, nfar, fcout)

  IF(nfar.ne.0) RETURN
  !
  ! READ atomic positions, forces, etc, and compute harmonic force and energy
  CALL harmonic_force_md(n_steps, nat_sc, Si,fcout,u,force_harm,h_energy)
    !
    SELECT CASE(input%fit_type)
    CASE('force', 'forces')
      diff = RESHAPE( (force_harm(:,:,:) - force_md(:,:,:)), (/ mdata /) )
    CASE('energy')
      ! not implemented with MPI_GATHER
      STOP 100
      !DO i = 1, n_steps
      !diff = diff + (h_energy(i)-toten_md(i))**2
      !ENDDO
    CASE('thforce')
      ! not implemented with MPI_GATHER
      STOP 100
      !diff = diff + RESHAPE( (ABS(force_harm(:,:,i) - force_md(:,:,i))**2 *EXP(-toten_md(i)) ), (/ mfc /) )
    CASE DEFAULT
      CALL errore("tdph", 'unknown chi2 method', 1)
    END SELECT
  
  CALL t_comm%start() 
  CALL allgather_vec(mdata, diff, diff_tot)
  CALL t_comm%stop() 

  IF(iswitch==1)THEN
    iter = iter+1
    chi2 = SQRT(SUM(diff_tot**2))
    ioWRITE(*,'(i10,e12.2)') iter, chi2

    !chi2 = (SUM( (force_harm(1:3,1:nat_sc,1:n_steps) - force_md(1:3,1:nat_sc,1:n_steps))**2 ))
    !CALL mpi_bsum(chi2)
    !chi2=SQRT(chi2)

    ioWRITE(ulog, "(i10,f14.6,9999f12.6)") iter, chi2, ph_coef
    ! Every input%nprint steps have a look
    IF(MODULO(iter,input%nprint)==0) &
      CALL write_fc2("matOUT.iter_"//TRIM(int_to_char(iter)), Si, fcout)
  ENDIF

  CALL t_chi2%stop()
  !
  END SUBROUTINE chi_lmdif
  !

  ! Penalty function as required by ACRS: return a single positive real number to minimize.
  !-----------------------------------------------------------------------
  SUBROUTINE chisq_acrs(ph_coef, nph, fdiff2)
    !-----------------------------------------------------------------------
    ! Calculates the square difference, fdiff2, btw harmonic and ab-initio
    ! forces for n_steps molecur dyanmics simulation 
    !
    USE tdph_module,  ONLY : nfar
    USE mpi_thermal,  ONLY : mpi_bsum
    USE decompose_d2, ONLY : recompose_fc
    IMPLICIT NONE
    INTEGER,INTENT(in)    :: nph
    REAL(DP),INTENT(in)   :: ph_coef(nph)
    REAL(DP),INTENT(out)  :: fdiff2
  
    INTEGER :: nq_done, iph, iq, i, j, k
    INTEGER,SAVE :: iter = 0
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    REAL(DP) :: chi2

    CALL t_minim%start()
  
    CALL recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_coef,&
                      nq1, nq2, nq3, nqmax, nfar, fcout)
          !
    IF(nfar.ne.0) RETURN
    !
    CALL harmonic_force_md(n_steps, nat_sc, Si,fcout,u,force_harm,h_energy)
    !
    fdiff2 = 0._dp
  
    DO i = 1, n_steps
      SELECT CASE(input%fit_type)
      CASE('force', 'forces')
        fdiff2 = fdiff2 + SUM((force_harm(:,:,i) - force_md(:,:,i))**2 )
      CASE('energy')
        fdiff2 = fdiff2 + (h_energy(i)-toten_md(i))**2
      CASE('thforce')
        fdiff2 = fdiff2 + SUM( (force_harm(:,:,i) - force_md(:,:,i))**2 )*EXP(-toten_md(i))
      CASE DEFAULT
        CALL errore("tdph", 'unknown chi2 method', 1)
      END SELECT
    ENDDO
    
    CALL  t_comm%start()
    CALL mpi_bsum(fdiff2) 
    CALL  t_comm%stop()

    iter = iter+1
    chi2 = SQRT(fdiff2)/n_steps_tot
    ioWRITE(*,'(i10,e12.2)') iter, chi2
    ioWRITE(ulog, "(i10,' ',f14.6,'  ',9999f12.6)") iter, chi2, ph_coef
    ! Every 1000 steps have a look
    IF(input%nprint>0 .and. MODULO(iter,input%nprint)==0)&
       CALL write_fc2("matOUT.iter_"//TRIM(int_to_char(iter)), Si, fcout)
    !ENDIF
  
    CALL t_minim%stop()
    !
    END SUBROUTINE chisq_acrs
  
  !----------------------------------------------------------------------------
 END PROGRAM tdph
!------------------------------------------------------------------------------
