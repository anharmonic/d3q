!
! Copyright (C) 2001-2011 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3_setup_module
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_setup_q_independent()
  !-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE lsda_mod,         ONLY : lsda, nspin
  USE gvecs,            ONLY : doublegrid
  USE eqv,              ONLY : dmuxc
  USE scf,              ONLY : rho, rho_core, v, vltot, vrs, kedtau
!  USE funct,            ONLY : dmxc, dmxc_spin
  USE fft_base,         ONLY : dfftp
  !USE gc_d3,            ONLY : setup_d3gc
  !
  IMPLICIT NONE
  INTEGER :: ir
  REAL (DP) :: rhotot, rhoup, rhodw
  !
  ! 1) Computes the total local potential (external+scf) on the smoot grid
  !
  CALL start_clock('d3_setup_noq')

  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Computes the derivative of the xc potential
  !    it only depends on the ground-state density, hence it is equal for 
  !    q1, q2 and q3.
  !
  ALLOCATE(dmuxc(dfftp%nnr, nspin_mag, nspin_mag))
  CALL setup_dmuxc()
  !dmuxc (:,:,:) = 0.d0
  !IF (lsda) THEN
  !   DO ir = 1, dfftp%nnr
  !      rhoup = rho%of_r(ir, 1) + 0.5d0 * rho_core (ir)
  !      rhodw = rho%of_r(ir, 2) + 0.5d0 * rho_core (ir)
  !      CALL dmxc_spin (rhoup, rhodw, dmuxc (ir, 1, 1), &
  !           dmuxc(ir, 2, 1), dmuxc(ir, 1, 2), dmuxc(ir, 2, 2) )
  !   ENDDO
  !ELSE
  !   DO ir = 1, dfftp%nnr
  !      rhotot = rho%of_r(ir, 1) + rho_core(ir)
  !      IF (rhotot > 1.d-30)   dmuxc(ir, 1, 1) =  dmxc(rhotot)
  !      IF (rhotot < - 1.d-30) dmuxc(ir, 1, 1) = -dmxc(-rhotot)
  !   ENDDO
  !ENDIF
  !
  !CALL setup_d3gc()
  !
  CALL stop_clock('d3_setup_noq')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_setup(xq1, xq2, xq3)
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  d3toten program:
  !  1) computes the total local potential (external+scf) on the smoot
  !     grid to be used in h_psi and similia
  !  2) computes dmuxc 3.1) with GC if needed
  !  3) for metals sets the occupated bands
  !  4) computes alpha_pv
  !  5.1) computes the variables needed to pass to the pattern representat
  !       of the small group of q
  !     u      the patterns
  !     t      the matrices of the small group of q on the pattern basis
  !     tmq    the matrix of the symmetry which sends q -> -q + G
  !     gi     the G associated to each symmetry operation
  !     gimq   the G of the q -> -q+G symmetry
  !     irgq   the small group indices
  !     nsymq  the order of the small group of q
  !     irotmq the index of the q->-q+G symmetry
  !     nirr   the number of irreducible representation
  !     npert  the dimension of each irreducible representation
  !     nmodes the number of modes
  !     minus_q true if there is a symmetry sending q -> -q+G
  !  5.2) computes the variables needed to pass to the pattern representat
  !       of the group of the crystal
  !     ug0     the patterns
  !     tg0     the matrices of the group on the pattern basis
  !     nsymg0  the order of the group of the crystal
  !     nirrg0  the number of irreducible representation
  !     npertg0 the dimension of each irreducible representation
  !  6) set the variables needed to deal with nlcc
  !  7) set the variables needed to distribute one loop between pools
  !  8) set the variables needed to calculate only selected q=0 modes
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout, ionode_id, ionode
  USE lr_symm_base,     ONLY : rtau, irgq, minus_q, irotmq, nsymq
  USE modes,            ONLY : nmodes, npertx
  USE d3com,            ONLY : npert_i, npert_f
  !
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,        ONLY : at, bg
  USE fft_base,         ONLY : dfftp
  USE symm_base,        ONLY : nsym, s, irt, invs, inverse_s, &
                               s_axis_to_cart, find_sym, copy_sym, sname
  USE uspp_param,       ONLY : upf
  USE control_flags,    ONLY : modenum
  USE constants,        ONLY : degspin, pi
  USE uspp,             ONLY : nlcc_any
  USE d3_basis,         ONLY : patq, allocate_d3_pattern
  USE d3_symmetry,      ONLY : symq, allocate_d3_symmetry, sym_gamma, &
                               d3_set_sym_irr, minus_3q
  USE kplus3q,          ONLY : q_special_cases, kplusq
  USE io_files,         ONLY : tmp_dir
  USE control_ph,       ONLY : tmp_dir_ph
  !
  USE d3_iofiles,       ONLY : fildrho_q, tmp_dir_d3, fildrho_dir
  USE mp,               ONLY : mp_bcast, mp_barrier
  USE mp_world,         ONLY : world_comm
  USE noncollin_module, ONLY : m_loc
  USE extfield,         ONLY : gate
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! minimum band energy
  ! maximum band energy
  ! working array

  INTEGER :: iq, ii, ipert
  ! counters
  LOGICAL :: sym(48), sym2(48), magnetic_sym
  ! the symmetry operations
  REAL (DP) :: mdum(3) = 0._dp
  CHARACTER(len=256)   :: fildrho_tmp
  CHARACTER(len=8),PARAMETER :: sub = 'd3_setup'
  !
  CALL start_clock (sub)
  WRITE(stdout, "(5x,a)") "================== setup started =================="
  !
  ! 5) set all the variables needed to use the pattern representation
  !    kplusq is set here (except the k+q grids!):
  CALL q_special_cases(xq1, xq2, xq3, at)
  !
  ! Find symmetries of the crystal (i.e. bravais lattic + ions base) and reorder them
  ! to be the first symm_base:nrot ones. This changes several global variables in module symm_base
  modenum = 0
  magnetic_sym = .false.
  ALLOCATE(m_loc(3,nat))
  m_loc = 0._dp
  CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc, gate )
  sym(:)      = .false.
  sym(1:nsym) = .true.
  !
  ! For each q we find it's small group and set up accordingly the variables in its symq type
  !
  WRITE(stdout,'(7x,a)') "Symmetry operations:"
  DO iq = 1,3
    symq(iq)%sym = sym
!    call smallg_q(kplusq(iq)%xq, modenum, at, bg, nsym, s, ftau, symq(iq)%sym, symq(iq)%minus_q)
    call smallg_q(kplusq(iq)%xq, modenum, at, bg, nsym, s, symq(iq)%sym, symq(iq)%minus_q)
!    CALL smallg_q (xq, 0, at, bg, nsym, s, sym, minus_q)
    ! FIXME! minus_q is currently not available in D3 code    
    symq(iq)%minus_q = .false. 
    ! The phonon mechanism for minus_q doesn ot work in D3 because, i.e. in the case of G,q,-q,
    ! G will get irotmq = 1, while q and -q will get the same irotmq: I do not know what to do
    ! in this case. A specific subroutine that looks for a triple irotmq is required..
    
    WRITE(stdout,'(9x,a,i1,":",i3,2x,48l1)') "q", iq, symq(iq)%nsymq, symq(iq)%sym
  ENDDO
  ! 
  ! Find the intersection of the small groups of the 3 q vectors
  sym = symq(1)%sym.and.symq(2)%sym.and.symq(3)%sym
  WRITE(stdout,'(7x,a)') "Symmetry common to all q vectors:"
  WRITE(stdout,'(9x,4x,48l1)') sym
  !
  ! Reorder the symmetry operations in the intersection of the small groups:
  ! DANGER DANGER DANGER!!! copy_sym fools around with global variables!
  nsymq  = copy_sym( nsym, sym ) ! <-- count sym.ops. and reorder them:
                                 ! 1..nsym  = symmetries of the crystal
                                 ! 1..nsymq = small group of q (note: nsym >= nsymq)
  ! set irgq (the index of symmetry operation of the small group of q) to identity
  ! as they have been reordered in copy_sym...
  irgq = 0
  DO ii = 1,nsymq
    irgq(ii) = ii
  ENDDO
  WRITE(stdout,'(9x,a,i2,a)') "-->", nsymq, " symmetry operation(s) have been found"
  !
  WRITE(stdout,*) "REMARK: q -> -q symmetry yet unsupported in d3q.x"
  minus_q = symq(1)%minus_q .and. symq(2)%minus_q .and. symq(3)%minus_q
  WRITE(stdout,'(9x,a,3l,a)') "individual q --> -q operations: ", symq(1)%minus_q, symq(2)%minus_q, symq(3)%minus_q
  !
  ! Look if there is a single sym.op. that links ALL THREE q's to their opposites
  CALL minus_3q(kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq, at, bg, nsym, s, minus_q, irotmq)
  IF(.not. minus_q) THEN
!     ! set individual minus_q to .false. if global minus_q is false
!     symq(1)%minus_q = .false.
!     symq(2)%minus_q = .false.
!     symq(3)%minus_q = .false.
  ELSE
    WRITE(stdout,'(9x,3a)') "a global q --> -q  symm.op. has been found (",TRIM(sname(irotmq)),")"
 ENDIF
  !
  CALL inverse_s( ) ! <-- computes invs from symm_base
  CALL s_axis_to_cart ( )
  !
  !  the first nsymq matrices are symmetries of the small group of q
  !
  ! 5.1) Finds the variables needeed for the pattern representation
  !      of the small group of q
  !
  ALLOCATE(rtau( 3, 48, nat))
  sym2(1:nsym)=.true.  !<-- BAD HACK! fix sgam_ph instead so that it does something sensible
  CALL sgam_ph(at, bg, nsym, s, irt, tau, rtau, nat, sym2)

  nmodes = 3 * nat
  !
  ALLOCATE(patq(-3:3))
  !
  IF (modenum .ne. 0) THEN
      CALL errore(sub, 'TO BE FIXED!', 1)
  ELSE
!      ANY_CRYSTAL_SYMMETRY : &
!      IF (nsym > 1) THEN
        ! read from file the patterns used in the 3 phonon calculations
        DO iq = 1,3
          !
          CALL allocate_d3_pattern(nat, patq(iq))
          CALL allocate_d3_symmetry(nat, -1, symq(iq))
          !
          IF(kplusq(iq)%ldrho) THEN
            fildrho_tmp = fildrho_q(kplusq(iq)%drho_from)%name
          ELSE IF(kplusq(iq)%ldrho_cc) THEN
            fildrho_tmp = fildrho_q(kplusq(iq)%drho_cc_from)%name
          ELSE
            CALL errore(sub, 'This should never happen!', 2)
          ENDIF
          !
          ! FIXME: workaround for filename mess - needed to find where the patterns are
          tmp_dir = tmp_dir_ph
          !
          WRITE(stdout, '(7x,a,i1)') "Reading patterns and q-point for q",iq
          IF(ionode) &
            CALL io_pattern(nat, fildrho_tmp, symq(iq)%nirr, symq(iq)%npert,&
                            patq(iq)%u, kplusq(iq)%xq_drho, fildrho_dir,-1)
          CALL mp_bcast(symq(iq)%nirr,      ionode_id, world_comm)
          CALL mp_bcast(symq(iq)%npert,     ionode_id, world_comm)
          CALL mp_bcast(kplusq(iq)%xq_drho, ionode_id, world_comm)
          CALL mp_bcast(patq(iq)%u,         ionode_id, world_comm)
          !
          tmp_dir = tmp_dir_d3
          !
          IF ( kplusq(iq)%ldrho_is_mine .and. &
              SUM(ABS(kplusq(iq)%xq_drho-kplusq(iq)%xq))>1.e-6_dp) THEN
              WRITE(stdout, '(9x,a,i1,a,3f8.4,a,3f8.4,a)') "--> drho of q",iq,&
                            " (",kplusq(iq)%xq,") was computed by ph.x for (",kplusq(iq)%xq_drho,")"
          ENDIF
          !
          IF(kplusq(iq)%ldrho_cc) THEN
            WRITE(stdout, '(9x,a,i1)') "--> taking c.c. of patterns for q",iq
            patq(iq)%u = CONJG(patq(iq)%u)
          ENDIF
          !
          ! the maximum number of perturbation amongs all irreducible
          ! representations for each q vector
          !
          symq(iq)%npertx = MAXVAL(symq(iq)%npert(1:symq(iq)%nirr))
#ifdef D3_EXTREME_VERBOSITY
          WRITE(stdout, '(/,7x,a,i1)') "* patterns of q",iq
          DO ipert = 1,3*nat
            WRITE(stdout, '(9x,"* pert. ",i3)') ipert
            WRITE(stdout, '((9x,3(3x,2f10.5)))') patq(iq)%u(:,ipert)
          ENDDO
#endif
        ENDDO
        !
        ! global npertx is the maximum of them all
        npertx = MAXVAL((/symq(1)%npertx, symq(2)%npertx, symq(3)%npertx/))

        DO iq = 1,3
          CALL allocate_d3_symmetry(nat, symq(iq)%npertx, symq(iq))
          !
          CALL d3_set_sym_irr(nat, at, bg, kplusq(iq)%xq, s, invs, nsym, rtau, &
              irt, symq(iq)%irgq, symq(iq)%nsymq, symq(iq)%minus_q, symq(iq)%irotmq, &
              symq(iq)%t, symq(iq)%tmq, symq(iq)%npertx, patq(iq)%u, symq(iq)%npert, &
              symq(iq)%nirr, symq(iq)%gi, symq(iq)%gimq)
        ENDDO
!      ELSE ANY_CRYSTAL_SYMMETRY
!         WRITE(stdout, '(7x,a)') "Not using symmetry! Patterns are all 1-atom cartesian displacements."
!         npertx=1
!         DO iq = 1,3
!           !
!           CALL allocate_d3_pattern(nat, patq(iq))
!           CALL allocate_d3_symmetry(nat, -1, symq(iq))
!           CALL allocate_d3_symmetry(nat,  1, symq(iq))
!           !
!           symq(iq)%npertx=1
!           CALL set_irr_nosym (nat, at, bg, kplusq(iq)%xq, s, invs, nsym, rtau, &
!               irt, symq(iq)%irgq, symq(iq)%nsymq, symq(iq)%minus_q, symq(iq)%irotmq, &
!               symq(iq)%t, symq(iq)%tmq, npertx, patq(iq)%u, symq(iq)%npert, &
!               symq(iq)%nirr, symq(iq)%gi, symq(iq)%gimq, iverbosity)
! 
!         ENDDO
!         npertx = 1
!      ENDIF &
!      ANY_CRYSTAL_SYMMETRY
     !
     sym_gamma => null()
     DO iq = 1,3
        ! Set up pattern for -q points
        CALL allocate_d3_pattern(nat, patq(-iq))
        patq(-iq)%u = CONJG(patq(iq)%u)
        !
        ! If we are using Gamma, set up its symmetry
        IF(SUM(kplusq(iq)%xq(:)**2) < 1.e-5_dp .and. .not. ASSOCIATED(sym_gamma) ) THEN
            WRITE(stdout, '(7x,a,i1)') "Symmetry information for Gamma associated with q",iq
            sym_gamma => symq(iq)
        ENDIF
     ENDDO
  ENDIF
  !
  WRITE(stdout, '(/,5x,a)') "Symmetry analysis done"
  !
  ! Set non linear core correction stuff
  !
  nlcc_any = ANY( upf(1:ntyp)%nlcc )
  !
  !IF (nlcc_any) ALLOCATE (drc( ngm, ntyp))
  !
  ! 7) Sets up variables needed to distribute one loop between pools
  !
  ! FIXME: the code is not ready to be used with npert_i and npert_f, grep all the places where they are 
  !        used and if needed add a sum inter_pool
  !        
! #ifdef __MPI
!   CALL block_distribute( 3*nat, my_pool_id, npool, npert_i, npert_f, nresto )
! #else
  npert_i = 1
  npert_f = 3 * nat
! #endif
  !
  !
  WRITE(stdout, "(5x,a,//)") "=================== setup done ===================="
  CALL stop_clock (sub)

  CALL mp_barrier(world_comm)

  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_setup
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
END MODULE d3_setup_module
!-----------------------------------------------------------------------
