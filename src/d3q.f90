!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program d3q
  !-----------------------------------------------------------------------
  !
  USE pwcom,         ONLY : lgauss
  USE phcom,         ONLY : all_done
  USE d3com,         ONLY : code, ethr_ph, restart, fild3dyn
  USE ions_base,     ONLY : nat
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : gamma_only
  USE mp_global,     ONLY : mp_startup
  USE mp_world,      ONLY : mpime
  USE environment,   ONLY : environment_start
  USE d3_basis,      ONLY : patq
  USE kplus3q,       ONLY : kplusq
  USE linter_d3q,    ONLY : generate_dwfc2
  USE uspp,          ONLY : nlcc_any
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE cell_base,     ONLY : at

  ! D3 subroutines:
  USE d3_iofiles,       ONLY : openfild3, openfile_drho, d3_add_rho_core,&
                               setup_d3_iofiles
  USE d3_exc_module
  USE d3_exc_gc_module
  USE dpsi1dv2dpsi3_module
  USE dpsi1dpsi2dv3_module
  USE dq1rhodq23v_module
  USE rhodq123v_module
  USE d3ionq_module
  USE d3_nlcc_module,   ONLY : d3_nlcc_0, d3_nlcc_123
  USE d3_valence_module
  USE d3matrix_module
  USE dpsi1dv2psi_module
  !
  USE d3_shuffle,       ONLY : nperms, d3perms, d3_shuffle_equiv, d3_check_permutations
  !USE davcio_debug

  USE nscf_d3,            ONLY : run_nscf_d3
  USE d3_grid,            ONLY : d3_triplets, n_triplets, n_triplets_todo, &
                                 i_triplet_first, i_triplet_last, i_triplet_offset, i_triplet_step
  USE d3_readin_module,   ONLY : d3_readin
  USE d3_setup_module,    ONLY : d3_setup_q_independent, d3_setup
  USE d3_init_module,     ONLY : d3_init
  USE d3_reset_module,    ONLY : d3_reset
  USE allocate_d3_module, ONLY : allocate_d3
  USE stop_d3_module,     ONLY : stop_d3
  USE d3matrix_io2,       ONLY : read_d3dyn_xml2
  USE d3_open,            ONLY : listu_d3
  USE d3_restart,         ONLY : d3_check_restart, d3_check_time, d3_from_scratch
  USE d3_debug
  USE mp_world,           ONLY : world_comm
  !USE funct,              ONLY : dft_is_gradient
  USE xc_lib,             ONLY : xclib_dft_is
  !USE dft_par_mod,  ONLY: isgradient
  USE gc_d3,            ONLY : setup_d3gc
!
  USE mp_pools,        ONLY : intra_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm, inter_bgrp_comm
!  USE mp_diag,         ONLY : mp_start_diag
  USE command_line_options,  ONLY : input_file_, ndiag_

  IMPLICIT NONE
  include 'laxlib.fh'

  TYPE d3matrix_with_permutations
    COMPLEX(DP),ALLOCATABLE :: dyn(:,:,:)
  END TYPE
  TYPE(d3matrix_with_permutations) :: d3(nperms)

  REAL(DP) :: t0=0._dp, t1=0._dp
  REAL(DP),EXTERNAL :: get_clock
  REAL(DP) :: xq1(3), xq2(3), xq3(3)

  COMPLEX(DP),ALLOCATABLE :: d3tmp(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: &
      d3dyn(:,:,:),        &
      d3dyn_dpdvdp(:,:,:), &
      d3dyn_dpdpdv(:,:,:), &
      d3dyn_drd2v(:,:,:),  &
      d3dyn_rd3v(:,:,:),   &
      d3dyn_ion(:,:,:),    &
      d3dyn_exc(:,:,:),    &
      d3dyn_nlcc(:,:,:),   &
      d3dyn_smear(:,:,:)
  INTEGER :: ip=0, jp=0, i_triplet, n_triplets_done
  LOGICAL :: found
  !
  !
  gamma_only = .false.
  all_done=.false.
  !
  ! Initialize MPI, clocks, print initial messages
  !
!#ifdef __MPI
!  CALL mp_startup ( )
!  CALL mp_startup ( start_images=.false. )
!  CALL mp_start_diag ( ndiag_, world_comm, intra_bgrp_comm, &
!       do_distr_diag_inside_bgrp_ = .true. )
!  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
!       inter_bgrp_comm )
!  CALL mp_startup ( start_images=.true. )
!  CALL laxlib_start ( ndiag_, world_comm, intra_bgrp_comm, &
!       do_distr_diag_inside_bgrp_ = .true. )
!  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
!       inter_bgrp_comm )
!#endif
  CALL environment_start ( code )
  !
  CALL start_clock('D3TOTEN')
  !
#ifndef __XLF
  WRITE(stdout, '(4(5x,a,/) )') &
      " ___  ____   __    ____  ____  ____  _  _  ___    ____  ___  ",&
      "/ __)(_  _) /__\  (  _ \(_  _)(_  _)( \( )/ __)  (  _ \(__ ) ",&
      "\__ \  )(  /(__)\  )   /  )(   _)(_  )  (( (_-.   )(_) )(_ \ ",&
      "(___/ (__)(__)(__)(_)\_) (__) (____)(_)\_)\___/  (____/(___/ "
#endif
  !
  ! Initialization routines
  !
  CALL d3_readin()
    ! It reads standard input and pw.x save file, plus it runs a few sanity checks
    !         CALL get_env( 'ESPRESSO_TMPDIR', outdir )
    !         CALL bcast_d3_input()
    !         CALL read_file()
    !         CALL d3_grid_init .or. d3_single_point_init

  CALL d3_setup_q_independent()
  CALL setup_d3gc()

    !           CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid)
    !           CALL dmxc_spin (rhoup, rhodw, dmuxc (ir, 1, 1), &
  !  _  _  _       ___  ____  ____  ____     ____  ____  ___  ____  _  _      _  _  _
  ! ( \( \( \     / __)(  _ \(_  _)(  _ \   (  _ \( ___)/ __)(_  _)( \( )    / )/ )/ )
  !  \ \\ \\ \   ( (_-. )   / _)(_  )(_) )   ) _ < )__)( (_-. _)(_  )  (    / // // /
  !   \_)\_)\_)   \___/(_)\_)(____)(____/   (____/(____)\___/(____)(_)\_)  (_/(_/(_/
  !
  n_triplets_done = 0
  GRID_CALCULATION_QPOINT_LOOP : &
  DO i_triplet = i_triplet_first+i_triplet_offset, i_triplet_last, i_triplet_step
    n_triplets_done = n_triplets_done+1
    xq1 = d3_triplets(i_triplet)%xq1
    xq2 = d3_triplets(i_triplet)%xq2
    xq3 = d3_triplets(i_triplet)%xq3
    !
    IF(restart)THEN
      CALL read_d3dyn_xml2(fild3dyn, xq1,xq2,xq3, at=at, found=found, seek=.true.)
      IF(found) THEN
        write( stdout, '(/,5x,"===================================================")')
        write( stdout,   '(5x,"=       triplet ",i4," already done (skipping)      =")') i_triplet
        write( stdout,   '(5x,"===================================================",/)')
        CYCLE GRID_CALCULATION_QPOINT_LOOP
      ENDIF
      !
    ENDIF
    !
    write( stdout, '(/,5x,"===================================================")')
    write( stdout,   '(5x,"=             Starting D3 calculation             =")')
    write( stdout,   '(5x,"===================================================")')
    IF (n_triplets>1) THEN
    write( stdout,   '(5x,"=     triplet ",i4," ( ",i4," of ",i4,", ",i3,"% done)     =")') &
                    i_triplet, n_triplets_done, n_triplets_todo, &
                    NINT( 100*(REAL(n_triplets_done-1)/REAL(n_triplets_todo) ))
    write( stdout,   '(5x,"===================================================")')
    ENDIF
    write( stdout,   '(5x,"=    q1 = ( ",3f11.6," )   =")') xq1
    write( stdout,   '(5x,"=    q2 = ( ",3f11.6," )   =")') xq2
    write( stdout,   '(5x,"=    q3 = ( ",3f11.6," )   =")') xq3
    write( stdout,   '(5x,"===================================================")')
    FLUSH( stdout )
    !
    CALL setup_d3_iofiles(xq1, xq2, xq3)
    !
    IF(restart)THEN
      CALL d3_check_restart('read')
    ELSE
      CALL d3_from_scratch()
    ENDIF
    FLUSH( stdout )
    !
    CALL d3_setup(xq1, xq2, xq3)
    FLUSH( stdout )
    ! It reads the fild?rho files and sets up symmetry
    !  This subroutine prepares several variables which are needed in the
      !           CALL find_sym ( nat, tau, ityp, nr1, nr2, nr3, .FALSE., &
      !           CALL smallg_q (kplusq(iq)%xq, modenum, at, bg, nsym, s, ftau, symq(iq)%sym, symq(iq)%minus_q)
      !           CALL sgam_ph (at, bg, symq(2)%nsymq, s, irt, tau, rtau, nat, sym)
      !           CALL allocate_d3_pattern(nat, patq(iq))
      !           CALL allocate_d3_symmetry(nat, -1, symq(iq))
      !           CALL io_pattern(nat,fild1rho,symq(iq)%nirr,symq(iq)%npert,patq(iq)%u,-1)
      !           CALL allocate_d3_symmetry(nat, symq(iq)%npertx, symq(iq))
      !           CALL set_sym_irr (nat, at, bg, kplusq(iq)%xq, s, invs, nsym, rtau, &
      !           CALL allocate_d3_pattern(nat, patq(iq))
      !           CALL set_irr_nosym (nat, at, bg, kplusq(iq)%xq, s, invs, nsym, rtau, &
    !
    write( stdout, '(/,5x,"=============== run_nscf_d3 start",i3," ==============",/)') 
    CALL run_nscf_d3(.true.)
!     CALL run_nscf_d3(.false.)
    write( stdout, '(/,5x,"=============== run_nscf_d3 done",i3," ===============",/)') 
    FLUSH( stdout )
    ! inside run_nscf_d3 subroutine setup_nscf_d3 is CALLed which increases the number
    ! of kpoints in order to generate k, k+/-q1, k+/-q2, k+/-q3
    ! wfcs are then recalculated for all these points
      !           CALL clean_pw( .FALSE. )
      !           CALL close_files()
      !           CALL restart_from_file()
      !           CALL setup_nscf_d3 (q1,q2,q3)
      !                   CALL sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
      !                   CALL irreducible_BZ (nrot, s, nsymq, minus_q, at, bg, npk, nkstot, xk, wk, &
      !                   CALL set_kplus3q( xk, wk, q1,q2,q3, nkstot, npk, nat )
      !               CALL init_run()
      !               CALL non_scf()
      !               CALL punch( 'all' )
      !               CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
      !               CALL close_files()
    write( stdout, '(/,5x,"=============== allocate_d3 start",i3," ==============",/)') 
    CALL allocate_d3() !(~allocate_phq)
    write( stdout, '(/,5x,"=============== allocate_d3 done",i3," ===============",/)') 
    !CALL d3_summary()
    !
    write( stdout, '(/,5x,"============== openfile_drho start",i3," =============",/)') 
    CALL openfile_drho()
    write( stdout, '(/,5x,"============== openfile_drho done",i3," ==============",/)') 
    FLUSH( stdout )
      !           CALL drho_change_q
      !                 CALL drho_add_phase
    write( stdout, '(/,5x,"================ openfild3 start",i3," ===============",/)') 
    CALL openfild3()
    FLUSH( stdout )
    write( stdout, '(/,5x,"================ openfild3 done",i3," ================",/)') 
    !
    CALL listu_d3(-1)
    FLUSH( stdout )
    !
    write( stdout, '(/,5x,"================= d3_init start",i3," ================",/)') 
    CALL d3_init()
    write( stdout, '(/,5x,"================= d3_init done",i3," =================",/)') 
    FLUSH( stdout )
      !           IF ( nlcc_any ) CALL set_drhoc(kplusq(iq)%xq, d3c(iq)%drc)
      !           CALL setlocq_coul( kplusq(iq)%xq, upf(nt)%zp, tpiba2, ngm, g, omega,d3v(iq)%loc(1,nt) )
      !           CALL setlocq( kplusq(iq)%xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
      !           CALL write_igkq_d3(ik)
    write( stdout, '(/,5x,"========== d3_check_permutations start",i3," =========",/)') 
    CALL d3_check_permutations()
    write( stdout, '(/,5x,"========== d3_check_permutations done",i3," ==========",/)') 
    FLUSH( stdout )
    !
    CALL print_clock(code)
    CALL mp_barrier(world_comm)
    !
    ! Temporary space to save stuff on file!
    DO ip = 1, nperms
      ALLOCATE(d3(ip)%dyn(3*nat, 3*nat, 3*nat))
    ENDDO
    !
    ALLOCATE(d3tmp(3*nat,3*nat,3*nat))
    !
    ALLOCATE(        d3dyn(3*nat,3*nat,3*nat) )
    ALLOCATE( d3dyn_dpdvdp(3*nat,3*nat,3*nat) )
    ALLOCATE( d3dyn_dpdpdv(3*nat,3*nat,3*nat) )
    ALLOCATE(  d3dyn_drd2v(3*nat,3*nat,3*nat) )
    ALLOCATE(   d3dyn_rd3v(3*nat,3*nat,3*nat) )
    ALLOCATE(    d3dyn_ion(3*nat,3*nat,3*nat) )
    ALLOCATE(  d3dyn_smear(3*nat,3*nat,3*nat) )
    ALLOCATE(    d3dyn_exc(3*nat,3*nat,3*nat) )
    ALLOCATE(   d3dyn_nlcc(3*nat,3*nat,3*nat) )
    d3dyn(:,:,:) = (0.d0, 0.d0)
    !
    !___________________________________________________________________________________
    !
    ! d3_add_rho_core(+1) adds the variation of the core charge
    ! to the charge derivative and writes it to a different file.
    ! The variation of the charge, modified this way is used by
    ! the routines d3_exc and d3dyn_cc.
    !
    write( stdout, '(/,5x,"================== add core start",i3," ==============",/)') 
    WRITE( stdout, '(/,5x,"Adding derivative of core charge")')
    !printed = .false.
    IF(dbg_add_core) THEN
       CALL d3_add_rho_core( +1._dp )
    ELSE
       WRITE(stdout,'(5x,"WARNING: Adding null core")')
       CALL d3_add_rho_core( 0._dp ) ! DEBUG: add a null core, multiplied by zero
    ENDIF
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"d3_add_rho_core   cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== add core done",i3," ===============",/)') 
    FLUSH( stdout )

    !___________________________________________________________________________________
    !
    !printed = .false.
    t0 = get_clock (code)
    write( stdout, '(/,5x,"===================================================")')
    write( stdout, '  (5x,"= Nscf calculation of the perturbed wavefunctions =")')
    write( stdout,   '(5x,"===================================================")')
    IF(dbg_do_dwfc) &
    CALL generate_dwfc2()
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    write( stdout,'(//,5x,"generate_dwfc   cpu time:",f9.2," sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== nscf dpsi done =================",/)')
    !
    !___________________________________________________________________________________
    !
    ! write on files terms of the type: <dpsi| dH | psi>, that
    ! will be used for the metallic case
    !
    write( stdout, '(/,5x,"================== precomp start",i3," ===============",/)') 
    !printed = .false.
    WRITE( stdout, '(/,5x,"Pre-computing < Pc dpsi_(k+X)/du(-X)| dH/du(Y) | psi_k-Y >")')
    IF(dbg_do_dpdvp) &
      CALL gen_dpsi1dv2psi()
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"gen_dpsi1dv2psi     cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== precomp done",i3," ================",/)') 
    FLUSH( stdout )
    !___________________________________________________________________________________
    !
    ! All the ingredients are ready here, time to compute the actual 3rd order terms
    !___________________________________________________________________________________
    !
    ! calculate the terms < dpsi| dH | dpsi >
    !
    DBG_dpdvdp : IF(dbg_do_dpdvdp) THEN 
    write( stdout, '(/,5x,"================== dpdvdp start",i3," ================",/)') 
    WRITE( stdout, '(/,5x,"Calculating the matrix elements <dpsi |dH |dpsi>")')
    d3tmp = (0._dp, 0._dp)
    DO ip = 1,nperms
      !printed = .false.
      d3(ip)%dyn = (0._dp, 0._dp)
      IF (d3perms(ip)%todo) THEN
        CALL dpsi1dv2dpsi3(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
      ELSE
        jp = d3perms(ip)%shuffle_from
        CALL d3_shuffle_equiv(nat, d3perms(jp)%i,d3perms(jp)%j,d3perms(jp)%k, &
                              d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, &
                              d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
      ENDIF
      d3tmp = d3tmp + d3(ip)%dyn
      CALL dbgwrite_d3dyn(2*d3(ip)%dyn, 'dpdvdp.'//d3perms(ip)%name, 1)
    ENDDO
    CALL dbgwrite_d3dyn(d3tmp, 'dpdvdp', 1)
    d3dyn_dpdvdp = d3tmp
    !
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    !
    WRITE( stdout, '(5x,"dpsi1dv2dpsi3 cpu time:",f9.2, &
          &   " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== dpdvdp done",i3," =================",/)') 
    FLUSH( stdout )
    ENDIF DBG_dpdvdp
    !___________________________________________________________________________________
    !
    ! calculate the term < dpsi| dpsi > < psi | dH | psi>
    DBG_dpdpdv : IF(dbg_do_dpdpdv) THEN 
    write( stdout, '(/,5x,"================== dpdpdv start",i3," ================",/)') 
    WRITE( stdout, '(/,5x,"Calculating the matrix elements <dpsi|dpsi>< psi|dH|psi> ")')
    d3tmp = (0._dp, 0._dp)
    DO ip = 1,nperms
      !printed = .false.
      d3(ip)%dyn = (0._dp, 0._dp)
      IF (d3perms(ip)%todo) THEN
        CALL dpsi1dpsi2dv3(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
      ELSE
        jp = d3perms(ip)%shuffle_from
        CALL d3_shuffle_equiv(nat, d3perms(jp)%i,d3perms(jp)%j,d3perms(jp)%k, &
                              d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, &
                              d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
      ENDIF
      d3tmp = d3tmp + d3(ip)%dyn
      CALL dbgwrite_d3dyn(d3(ip)%dyn, 'dpdpdv.'//d3perms(ip)%name, 1)
    ENDDO
    CALL dbgwrite_d3dyn(d3tmp, 'dpdpdv', 1)
    d3dyn_dpdpdv = d3tmp
    !
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    !
    WRITE( stdout, '(5x,"dpsi1dpsi2dv3 cpu time:",f9.2, &
          &   " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== dpdpdv done",i3," =================",/)') 
    FLUSH( stdout )
    ENDIF DBG_dpdpdv
    !___________________________________________________________________________________
    !
    ! calculate the term   drho * d2V
    !
    DBG_drhod2v : IF(dbg_do_drhod2v) THEN
    write( stdout, '(/,5x,"================== dpd2v start",i3," =================",/)') 
    WRITE( stdout, '(/,5x,"Calculating the matrix elements <psi |d^2 v |dpsi>")')
    !
    d3tmp = (0._dp, 0._dp)
    DO ip = 1,nperms
      !printed = .false.
      d3(ip)%dyn = (0._dp, 0._dp)
      !
      IF (d3perms(ip)%todo .and. d3perms(ip)%todo_first) THEN
        CALL dq1rhodq23v(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
      ELSE
        IF ( .not. d3perms(ip)%todo_first) THEN
          !print*, "doing first"
          jp         = d3perms(ip)%first_from
          d3(ip)%dyn = d3(jp)%dyn
        ELSE
          !print*, "doing shuffle"
          jp         = d3perms(ip)%shuffle_from
          CALL d3_shuffle_equiv(nat, d3perms(jp)%i, d3perms(jp)%j, d3perms(jp)%k, &
                                d3perms(ip)%i, d3perms(ip)%j, d3perms(ip)%k, &
                                d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
        ENDIF
      ENDIF
      d3tmp = d3tmp + d3(ip)%dyn
      !WRITE(10077,'(//,a)') 'drd2v.'//d3perms(ip)%name
      !WRITE(10077,'(2i3,2x,3i3,2x,3i3,l2)') ip,jp, &
      !d3perms(jp)%i, d3perms(jp)%j, d3perms(jp)%k, &
      !d3perms(ip)%i, d3perms(ip)%j, d3perms(ip)%k, d3perms(ip)%shuffle_conjg
      !WRITE(10077,'(3(2f12.6,3x))') d3(jp)%dyn
      !WRITE(10077,'(/,a)')  
      !WRITE(10077,'(3(2f12.6,3x))') d3(ip)%dyn
      CALL dbgwrite_d3dyn(2*d3(ip)%dyn, 'drd2v.'//d3perms(ip)%name, 1)
    ENDDO
    !stop 1
    !
    CALL dbgwrite_d3dyn(d3tmp, 'drd2v', 1)
    d3dyn_drd2v = d3tmp
    !
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"drhod2v       cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== dpd2v done",i3," ==================",/)') 
    FLUSH( stdout )
    ENDIF DBG_drhod2v
    !___________________________________________________________________________________
    !
    ! It calculates the term   rho * d3V
    !
    DBG_rhod3v : IF(dbg_do_rhod3v) THEN
    write( stdout, '(/,5x,"=================== rho d3v start =================",/)')
    WRITE( stdout, '(/,5x,"Calculating the matrix elements <psi |d^3v |psi>")')
      d3dyn_rd3v = (0._dp, 0._dp)
      !printed = .false.
    CALL rhodq123v(d3dyn_rd3v)
      CALL dbgwrite_d3dyn (d3dyn_rd3v, 'rd3v', 1)
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"d3vrho        cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"=================== rho d3v done ==================",/)')
    FLUSH( stdout )
    ENDIF DBG_rhod3v
    !___________________________________________________________________________________
    !
    ! It calculates the contribution due to ionic term
    !
    DBG_ion : IF(dbg_do_ion) THEN
    write( stdout, '(/,5x,"================== ewald start",i3," =================",/)') 
    WRITE( stdout, '(/,5x,"Calculating the Ewald contribution")')
    d3dyn_ion = (0._dp, 0._dp)
    CALL d3ionq( kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq, &
                 patq(1)%u, patq(2)%u, patq(3)%u, ethr_ph, d3dyn_ion )
    CALL dbgwrite_d3dyn (d3dyn_ion, 'd3ionq', 1)
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"d3ionq        cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================== ewald done",i3," ==================",/)') 
    FLUSH( stdout )
    ENDIF DBG_ion
    !___________________________________________________________________________________
    !
    ! In the metallic case, calculate some additional terms
    !
    d3dyn_smear = (0._dp, 0._dp)
    !
    ADD_SMEARING_CONTRIBUTION : &
    IF (lgauss .and. dbg_do_smearing) THEN
      d3dyn_smear = 0._dp
      !
      DBG_smr1 : IF(dbg_do_smr_ijk) THEN
      write( stdout, '(/,5x,"================== valence start",i3," ==============",/)') 
      WRITE( stdout, '(/,5x,"Calculating the valence contribution")')
       ! first term
      d3tmp = (0._dp, 0._dp)
      DO ip = 1,nperms
        !printed = .false.
        d3(ip)%dyn = (0._dp, 0._dp)
        IF (d3perms(ip)%todo) THEN
          CALL d3_valence_ijk(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
        ELSE
          jp = d3perms(ip)%shuffle_from
          CALL d3_shuffle_equiv(nat, d3perms(jp)%i,d3perms(jp)%j,d3perms(jp)%k, &
                                d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, &
                                d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
        ENDIF
        d3tmp = d3tmp + d3(ip)%dyn
        CALL dbgwrite_d3dyn(6*d3(ip)%dyn, 'smr_ijk.'//d3perms(ip)%name, 1)
      ENDDO
      CALL dbgwrite_d3dyn(d3tmp, 'smr_ijk', 1)
      d3dyn_smear = d3tmp
      ENDIF DBG_smr1
      !_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
      ! second term
      ! only do perms that are not equivalent (tricky to do properly):
      DBG_smr2 : IF(dbg_do_smr_ij) THEN
      d3tmp = (0._dp, 0._dp)
      DO ip = 1,nperms
        !printed = .false.
        d3(ip)%dyn = (0._dp, 0._dp)
        !
        IF (d3perms(ip)%todo .and. d3perms(ip)%todo_first) THEN
          CALL d3_valence_ij(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
        ELSE
          IF ( .not. d3perms(ip)%todo_first) THEN
            jp         = d3perms(ip)%first_from
            d3(ip)%dyn = d3(jp)%dyn
          ELSE
            jp         = d3perms(ip)%shuffle_from
            CALL d3_shuffle_equiv(nat, d3perms(jp)%i, d3perms(jp)%j, d3perms(jp)%k, &
                                  d3perms(ip)%i, d3perms(ip)%j, d3perms(ip)%k, &
                                  d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
          ENDIF
        ENDIF
        d3tmp = d3tmp + d3(ip)%dyn
        CALL dbgwrite_d3dyn(d3(ip)%dyn, 'smr_ij.'//d3perms(ip)%name, 1)
      ENDDO
      CALL dbgwrite_d3dyn(d3tmp, 'smr_ij', 1)
      d3dyn_smear = d3dyn_smear + d3tmp
      ENDIF DBG_smr2
      !_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
      ! terms 3 and 4, only for q1=q2=q3=Gamma
      IF(dbg_do_smr_g) THEN
      d3tmp = (0._dp,0._dp)
      CALL d3_valence_gamma(d3tmp)
      CALL dbgwrite_d3dyn (d3tmp, 'smr_g', 1)
      d3dyn_smear = d3dyn_smear+d3tmp
      ENDIF
      !
      !
      CALL dbgwrite_d3dyn (d3dyn_smear, 'smr', 1)
      t1 = get_clock (code) - t0
      t0 = get_clock (code)
      WRITE( stdout, '(5x,"d3_valence    cpu time:",f9.2, &
          &         " sec    Total time:",f12.2," sec")') t1, t0
      write( stdout, '(/,5x,"================== valence done",i3," ===============",/)') 
      FLUSH( stdout )
    ENDIF ADD_SMEARING_CONTRIBUTION
    !___________________________________________________________________________________
    !
    ! d3_add_rho_core(+1) adds the variation of the core charge
    ! to the charge derivative and writes it to a different file.
    ! The variation of the charge, modified this way is used by
    ! the routines d3_exc and d3dyn_cc.
    !
!    write( stdout, '(/,5x,"================== add core start",i3," ==============",/)') 
!    WRITE( stdout, '(/,5x,"Adding derivative of core charge")')
!    !printed = .false.
!    IF(dbg_add_core) THEN
!       CALL d3_add_rho_core( +1._dp )
!    ELSE
!       WRITE(stdout,'(5x,"WARNING: Adding null core")')
!       CALL d3_add_rho_core( 0._dp ) ! DEBUG: add a null core, multiplied by zero
!    ENDIF
!    t1 = get_clock (code) - t0
!    t0 = get_clock (code)
!    WRITE( stdout, '(5x,"d3_add_rho_core   cpu time:",f9.2, &
!        &         " sec    Total time:",f12.2," sec")') t1, t0
!    write( stdout, '(/,5x,"================== add core done",i3," ===============",/)') 
!    FLUSH( stdout )
    !___________________________________________________________________________________
    !
    ! calculate d3Ei * drho * drho * drho, where drho is the variation
    ! of the charge and d3Ei is the third derivative of the
    ! Kohn-Sham-Energy term depending on the charge density.
    !
    DBG_exc : IF(dbg_do_exc) THEN
    write( stdout, '(/,5x,"================ exc contrib start",i3," =============",/)') 
    d3dyn_exc = (0._dp, 0._dp)
    !printed = .false.
    IF ( xclib_dft_is('gradient') .and. dbg_exc_do_gga) THEN
      WRITE( stdout, '(/,5x,"Calculating the exchange-correlation contribution with GGA")')
      d3tmp = (0._dp, 0._dp)
      DO ip = 1,nperms
        !printed = .false.
        d3(ip)%dyn = (0._dp, 0._dp)
        IF (d3perms(ip)%todo) THEN
          CALL d3_exc_gc(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
        ELSE
          jp = d3perms(ip)%shuffle_from
          CALL d3_shuffle_equiv(nat, d3perms(jp)%i,d3perms(jp)%j,d3perms(jp)%k, &
                                d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, &
                                d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
        ENDIF
        d3tmp = d3tmp + d3(ip)%dyn
        CALL dbgwrite_d3dyn(d3(ip)%dyn, 'exc_GGA.'//d3perms(ip)%name, 1)
      ENDDO
      !
      ! DEBUG DEBUG DEBUG!
      d3dyn_exc = d3tmp/6._dp
!       CALL dbgwrite_d3dyn (d3dyn_exc, 'exc_GGA', 1)
      !
!       d3tmp = (0._dp, 0._dp)
!       CALL d3_exc(d3tmp)
!       CALL dbgwrite_d3dyn (d3tmp, 'exc_LDA', 1)
!       !
!       d3dyn_exc = d3dyn_exc + d3tmp
      !
    ELSE
      WRITE( stdout, '(/,5x,"Calculating the exchange-correlation contribution")')
      d3dyn_exc = 0._dp
      CALL d3_exc(d3dyn_exc)
      CALL dbgwrite_d3dyn (d3dyn_exc, 'exc_LDA', 1)
    ENDIF
    CALL dbgwrite_d3dyn (d3dyn_exc, 'exc', 1)
    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"d3_exc        cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"================ exc contrib done",i3," =============",/)') 
    ENDIF DBG_exc
    FLUSH( stdout )
    !___________________________________________________________________________________
    !
    ! calculate additional terms due to non_linear-core-corrections
    !
    !
    d3dyn_nlcc = (0._dp, 0._dp)
    !
    ADD_NLCC_CORRECTION : &
    IF (nlcc_any .and. dbg_do_nlcc) THEN
      write( stdout, '(/,5x,"=============== nlcc contrib start",i3," ============",/)') 
      WRITE( stdout, '(/,5x,"Calculating the core-correction contribution")')
      !
      ! One diagonal term (correction to rho d3v):
      !printed = .false.
      IF(dbg_do_nlcc_0) THEN
      CALL d3_nlcc_0(d3dyn_nlcc)
      CALL dbgwrite_d3dyn(d3dyn_nlcc, 'nlcc_0', 1)
      ENDIF
      !_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
      !
      DBG_nlcc_123 : IF(dbg_do_nlcc_123) THEN
      d3tmp = (0._dp, 0._dp)
      DO ip = 1,nperms
        !printed = .false.
        d3(ip)%dyn = (0._dp, 0._dp)
        IF (d3perms(ip)%todo) THEN
          CALL d3_nlcc_123(d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, d3(ip)%dyn)
        ELSE
          jp = d3perms(ip)%shuffle_from
          CALL d3_shuffle_equiv(nat, d3perms(jp)%i,d3perms(jp)%j,d3perms(jp)%k, &
                                d3perms(ip)%i,d3perms(ip)%j,d3perms(ip)%k, &
                                d3perms(ip)%shuffle_conjg, d3(jp)%dyn, d3(ip)%dyn)
        ENDIF
        d3tmp = d3tmp + d3(ip)%dyn
        CALL dbgwrite_d3dyn(2*d3(ip)%dyn, 'nlcc_123.'//d3perms(ip)%name, 1)
      ENDDO
      CALL dbgwrite_d3dyn(d3tmp, 'nlcc_123', 1)
      d3dyn_nlcc = d3dyn_nlcc + d3tmp
      !
      CALL dbgwrite_d3dyn(d3dyn_nlcc, 'nlcc', 1)
      !
      t1 = get_clock (code) - t0
      t0 = get_clock (code)
      WRITE( stdout, '(5x,"d3_nlcc_*      cpu time:",f9.2, &
          &         " sec    Total time:",f12.2," sec")') t1, t0
      FLUSH( stdout )
      ENDIF DBG_nlcc_123
      !
      write( stdout, '(/,5x,"=============== nlcc contrib done",i3," =============",/)') 
    ENDIF ADD_NLCC_CORRECTION
    !___________________________________________________________________________________
    !
    ! Sum all contributions
    !
    write( stdout, '(/,5x,"=============== computing D3 start",i3," ============",/)') 
    !
    d3dyn = d3dyn_dpdvdp & ! <dpsi|dV|dpsi>
           +d3dyn_dpdpdv &
           +d3dyn_drd2v  & ! drho d2V
           +d3dyn_rd3v   & ! rho d3V
           +d3dyn_ion    & ! Ewald sum
           +d3dyn_smear  & ! metals only, 3 terms
           +d3dyn_exc    & ! d3E_xc/drho d1rho d2rho d3rho
           +d3dyn_nlcc     ! correction to exc and drho d2v terms
    CALL dbgwrite_d3dyn(d3dyn, 'dyn', 1)
    !
    ! Finally, symmetrize and compute the star of q points
    !
    WRITE( stdout, '(/,5x,"Symmetrizing and writing the tensor to disc")')
    WRITE( stdout, '(/,5x,"CALLing d3matrix")')
    CALL d3matrix(d3dyn,fild3dyn)

    t1 = get_clock (code) - t0
    t0 = get_clock (code)
    WRITE( stdout, '(5x,"d3matrix      cpu time:",f9.2, &
        &         " sec    Total time:",f12.2," sec")') t1, t0
    write( stdout, '(/,5x,"=============== computing D3 done",i3," =============",/)') 
    FLUSH( stdout )
    !
    write( stdout, '(/,5x,"================== cleanup start",i3," ==============",/)') 
    ! Deallocate
    DO ip = 1, nperms
      DEALLOCATE(d3(ip)%dyn)
    ENDDO
    DEALLOCATE(d3tmp)
    DEALLOCATE(d3dyn, d3dyn_dpdvdp, d3dyn_dpdpdv, d3dyn_drd2v, d3dyn_rd3v, &
              d3dyn_ion, d3dyn_smear, d3dyn_exc, d3dyn_nlcc )
    !
    CALL d3_check_time()
    !
    CALL d3_reset(print_clock=(n_triplets_done<n_triplets_todo), cleanup=.true.)
    !
    write( stdout, '(/,5x,"================== cleanup done",i3," ===============",/)') 
    FLUSH( stdout )
    !
  ENDDO &
  GRID_CALCULATION_QPOINT_LOOP
  !
  !  _  _  _       ___  ____  ____  ____     ____  _  _  ____       _  _  _
  ! ( \( \( \     / __)(  _ \(_  _)(  _ \   ( ___)( \( )(  _ \     / )/ )/ )
  !  \ \\ \\ \   ( (_-. )   / _)(_  )(_) )   )__)  )  (  )(_) )   / // // /
  !   \_)\_)\_)   \___/(_)\_)(____)(____/   (____)(_)\_)(____/   (_/(_/(_/
  !
  CALL stop_clock('D3TOTEN')
  !
  CALL laxlib_end()
  CALL stop_d3()
  !
#ifndef __XLF
  WRITE(stdout, '(4(5x,a,/))') &
        " ____  ___    ____  _____  _  _  ____",  &
        "(  _ \(__ )  (  _ \(  _  )( \( )( ___)", &
        " )(_) )(_ \   )(_) ))(_)(  )  (  )__)",  &
        "(____/(___/  (____/(_____)(_)\_)(____)"
#endif
  !
END PROGRAM d3q

