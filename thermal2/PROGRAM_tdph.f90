!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE tdph_module
  USE kinds, ONLY : DP
  !
  TYPE dynmat_basis
     COMPLEX(DP),ALLOCATABLE :: basis(:,:,:)
  END TYPE

  CONTAINS

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
  USE constants,          ONLY : amu_ry
  USE parameters,         ONLY : ntypx
  USE mp,                 ONLY : mp_bcast
  USE mp_global,          ONLY : mp_startup, mp_global_end
  USE mp_world,           ONLY : world_comm
  USE io_global,          ONLY : ionode_id, ionode, stdout
  USE environment,        ONLY : environment_start, environment_end
  ! symmetry
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, &
                                  irt, copy_sym, nrot, inverse_s, t_rev
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq, invsymq
  USE control_lr,         ONLY : lgamma
  USE decompose_d2
  USE cmdline_param_module
  USE input_fc,           ONLY : forceconst2_grid, ph_system_info, read_system, aux_system, read_fc2, &
                                 div_mass_fc2, multiply_mass_dyn, write_fc2
  USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE asr2_module,        ONLY : impose_asr2
  USE quter_module,       ONLY : quter
  USE tdph_module


  !
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: CODE="MKWEDGE"
  CHARACTER(len=256) :: fildyn, filout
  INTEGER :: ierr, nargs
  !
  INTEGER       :: nq1, nq2, nq3, nqmax, nq_wedge, nqq
  INTEGER       :: nq_done
  REAL(DP),ALLOCATABLE      :: x_q(:,:), w_q(:)
  REAL(DP) :: xq(3), syq(3,48)
  !
  LOGICAL :: sym(48), lrigid, skip_equivalence, time_reversal
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:,:,:,:), d2(:,:), w2(:,:), &
                             star_wdyn(:,:,:,:, :), star_dyn(:,:,:)
  REAL(DP),ALLOCATABLE :: decomposition(:), xqmax(:,:)
  INTEGER :: i,j, icar,jcar, na,nb, iq, ndf
  INTEGER,ALLOCATABLE :: rank(:)
  TYPE(ph_system_info) :: Si
  TYPE(forceconst2_grid) :: fc, fcout
  TYPE(dynmat_basis),ALLOCATABLE :: dmb(:)
  TYPE(sym_and_star_q),ALLOCATABLE :: symq(:)
  !
  CALL mp_startup()
  CALL environment_start(CODE)
  !
  ! setup output
  !fileout = cmdline_param_char("o", TRIM(fildyn)//".out")
  !
  ! Read system info from a mat2R file
  !OPEN(unit=999,file="mat2R",action='READ',status='OLD')
  !CALL read_system(999, Si)
  !CALL aux_system(Si)
  !CLOSE(999)

  CALL read_fc2("mat2R", Si, fc)
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
  WRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  IF(.not.allocated(m_loc))  THEN
    ALLOCATE(m_loc(3,Si%nat))
    m_loc = 0._dp
  ENDIF
  
  CALL find_sym ( Si%nat, Si%tau, Si%ityp, .false., m_loc )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !
  ! Find the reduced grid of q-points:
  skip_equivalence = .FALSE.
  time_reversal    = .TRUE.
  nq1 = fc%nq(1)
  nq2 = fc%nq(2)
  nq3 = fc%nq(3)
  nqmax = nq1*nq2*nq3
  ALLOCATE(x_q(3,nqmax), w_q(nqmax))
  call kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, Si%bg, nqmax,&
                    0,0,0, nq1,nq2,nq3, nq_wedge, x_q, w_q )
  !
  WRITE(stdout, *) "Generated ", nq_wedge, "points"
  
  ALLOCATE(rtau( 3, 48, Si%nat), d2(3*Si%nat,3*Si%nat))

  ! Variable to hold the dyn matrix and q-points of the entire grid
  nq_done = 0
  ALLOCATE(star_wdyn(3,3,Si%nat,Si%nat, nqmax))
  ALLOCATE(xqmax(3,nqmax))

  ! For every q-point in the irreducible wedge, we find its symmetry
  ! and the basis of the space of symmetry-constrained dynamical matrices
  ! Again, this part uses gloabl variables which is a annoying, but once
  ! we have the set of basis matrices, we don't have to touch it
  ! anymore.
  ALLOCATE(dmb(nq_wedge))
  ALLOCATE(rank(nq_wedge))
  ALLOCATE(symq(nq_wedge))

  Q_POINTS_LOOP : &
  DO iq = 1, nq_wedge
    WRITE(stdout, *) "____[[[[[[[", iq, "]]]]]]]]____"
    WRITE(stdout, '(i6, 3f12.4)') iq, x_q(:,iq)
    !
    ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
    ! part 1: call smallg_q and the copy_sym, 
    xq = x_q(:,iq)
    minus_q = .true.
  
    sym = .false.
    sym(1:nsym) = .true.
    CALL smallg_q(xq, 0, Si%at, Si%bg, nsym, s, sym, minus_q)
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
    ! the next subroutine uses symmetry from global variables to find he basis of crystal-symmetric
    ! matrices at this q point
    CALL find_d2_symm_base(xq, rank(iq), dmb(iq)%basis, &
       Si%nat, Si%at, Si%bg, symq(iq)%nsymq, symq(iq)%minus_q, &
       symq(iq)%irotmq, symq(iq)%rtau, symq(iq)%irt, symq(iq)%s, symq(iq)%invs )
    !
    !
    ! Calculate the list of points making up the star of q and of -q
    CALL tr_star_q(symq(iq)%xq, Si%at, Si%bg, symq(iq)%nsym, symq(iq)%s, symq(iq)%invs, &
                   symq(iq)%nq_star, symq(iq)%nq_trstar, symq(iq)%sxq, &
                   symq(iq)%isq, symq(iq)%imq, .false. )

    WRITE(stdout, '(5x,a,2i5)') "Found star of q and -q", symq(iq)%nq_star, symq(iq)%nq_trstar
    syq = symq(iq)%sxq
    call cryst_to_cart(symq(iq)%nq_trstar, syq, Si%at, -1)
    DO i = 1, symq(iq)%nq_trstar
       syq(1,i) = MODULO(syq(1,i), 1._dp)
       syq(2,i) = MODULO(syq(2,i), 1._dp)
       syq(3,i) = MODULO(syq(3,i), 1._dp)
       syq(:,i) = syq(:,i) * (/nq1, nq2, nq3/)
       WRITE(stdout,'(i4,3i3,l2)') i, NINT(syq(:,i)), (i>symq(iq)%nq_star)
    ENDDO

  ENDDO Q_POINTS_LOOP

  !
  ! Number of degrees of freedom for the entire grid:
  ndf = SUM(rank)
  WRITE(stdout, '("=================")')
  WRITE(stdout, '(5x,a,2i5)') "TOTAL number of degrees of freedom", ndf  

  nq_done = 0
  Q_POINTS_LOOP2 : &
  DO iq = 1, nq_wedge
    xq = symq(iq)%xq
    !
    ! Interpolate the system dynamical matrix at this q
    CALL fftinterp_mat2(xq, Si, fc, d2)
    ! Remove the mass factor, I cannot remove it before becaus the effective
    ! charges/long range interaction code assumes it is there
    d2 = multiply_mass_dyn(Si,d2)
    !
    IF(nq_done+symq(iq)%nq_trstar> nqmax) CALL errore("tdph","too many q-points",1)
    !
    ! Rotate the dynamical matrices to generate D(q) for every q in the star
    ALLOCATE(star_dyn(3*Si%nat,3*Si%nat, symq(iq)%nq_trstar))
    CALL make_qstar_d2 (d2, Si%at, Si%bg, Si%nat, symq(iq)%nsym, symq(iq)%s, &
                        symq(iq)%invs, symq(iq)%irt, symq(iq)%rtau, &
                        symq(iq)%nq_star, symq(iq)%sxq, symq(iq)%isq, &
                        symq(iq)%imq, symq(iq)%nq_trstar, star_dyn, &
                        star_wdyn(:,:,:,:,nq_done+1:nq_done+symq(iq)%nq_trstar))

    ! rebuild the full list of q vectors in the grid by concatenating all the stars
    xqmax(:,nq_done+1:nq_done+symq(iq)%nq_trstar) &
        = symq(iq)%sxq(:,1:symq(iq)%nq_trstar)

    nq_done = nq_done + symq(iq)%nq_trstar
    !
    ! Just a simple check
!    DO i = 1, nq_trstar
!       WRITE(stdout,'(i4,3f12.4,l2, 1f15.9)') i, sxq(:,i), (i>nq_star), &
!                       dotprodmat(3*nat,star_dyn(:,:,i), star_dyn(:,:,i))
!    ENDDO
    DEALLOCATE(star_dyn)
  
    !CALL compact_dyn(nat, d2, phi)
    print*, "== DECOMPOSITION =="
    DO i = 1,rank(iq)
      WRITE(stdout,"(i3,1f12.6)") i, dotprodmat(3*Si%nat,d2, dmb(iq)%basis(:,:,i))
    ENDDO
    !
    !d2 = 0
    ! print*, "== RECOMPOSITION =="
    !DO i = 1,rank
    !  d2 = d2 + decomposition(i)*basis(:,:,i)
    !ENDDO

  ENDDO Q_POINTS_LOOP2
  !
  CALL quter(nq1, nq2, nq3, Si%nat,Si%tau,Si%at,Si%bg, star_wdyn, xqmax, fcout, 2)
  CALL write_fc2("matOUZ", Si, fcout)


  ! print*, d2

  !DEALLOCATE(phi, d2, w2)
  !DEALLOCATE(rtau, tau, ityp)
  !DEALLOCATE(irt) ! from symm_base
  !----------------------------------------------------------------------------
 END PROGRAM tdph
!----------------------------------------------------------------------------
!
