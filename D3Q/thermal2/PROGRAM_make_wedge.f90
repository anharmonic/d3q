!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! A small utility that reads the first q from a dynamical matrix file (either xml or plain text),
! recomputes the system symmetry (starting from the lattice) and generates the star of q.
!
! Useful for debugging and for producing the star of the wannier-phonon code output.
!
! Syntax:
!   q2qstar.x filein [fileout]
!
! fileout default: rot_filein (old format) or rot_filein.xml (new format)
!
!----------------------------------------------------------------------------
PROGRAM make_wedge
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
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, irt, copy_sym, nrot, inverse_s, t_rev
  ! for reading the dyn.mat.
  USE cell_base,          ONLY : at, bg, celldm, ibrav, omega
  USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
  ! as above, unused here
  USE control_ph,         ONLY : xmldyn
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  !
  USE dynmat,             ONLY : w2
  !
  ! for non-xml file only:
  USE dynamicalq,         ONLY : dq_phiq => phiq, dq_tau => tau, dq_ityp => ityp, zeu 
  ! fox xml files only
  USE io_dyn_mat,         ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                                 read_dyn_mat, read_dyn_mat_tail, &
                                 write_dyn_mat_header
  ! small group symmetry
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq, invsymq
  USE control_lr,         ONLY : lgamma
  USE decompose_d2
  USE cmdline_param_module
  USE input_fc,           ONLY : forceconst2_grid, ph_system_info, read_system, aux_system, read_fc2, &
                                 div_mass_fc2, multiply_mass_dyn
  USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE asr2_module,        ONLY : impose_asr2


  !
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: CODE="MKWEDGE"
  CHARACTER(len=256) :: fildyn, filout
  INTEGER :: ierr, nargs
  !
  INTEGER       :: nq1, nq2, nq3, nqmax, nqs, isq (48), imq, nqq,nq_out
  REAL(DP),ALLOCATABLE      :: x_q(:,:), w_q(:)
  REAL(DP) :: xq(3), sxq(3,48)
  !
  LOGICAL :: sym(48), lrigid, skip_equivalence, time_reversal
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:,:,:,:), d2(:,:), basis(:,:,:), star_dyn(:,:,:)
  REAL(DP),ALLOCATABLE :: decomposition(:)
  INTEGER :: i,j, icar,jcar, na,nb, rank, iq
  TYPE(ph_system_info) :: Sinfo
  TYPE(forceconst2_grid) :: fc
  !
  CALL mp_startup()
  CALL environment_start(CODE)
  !
  ! setup output
  !fileout = cmdline_param_char("o", TRIM(fildyn)//".out")
  !
  ! Read system info from a mat2R file
  !OPEN(unit=999,file="mat2R",action='READ',status='OLD')
  !CALL read_system(999, Sinfo)
  !CALL aux_system(Sinfo)
  !CLOSE(999)

  CALL read_fc2("mat2R", Sinfo, fc)
  CALL impose_asr2("simple", Sinfo%nat, fc, Sinfo%zeu)
  CALL aux_system(Sinfo)
  CALL div_mass_fc2(Sinfo, fc)
  !
  ! Quantum-ESPRESSO symmetry subroutines use the global variables
  ! we copy the system data from structure S
  ntyp   = Sinfo%ntyp
  nat    = Sinfo%nat
  ALLOCATE(tau(3,nat))
  ALLOCATE(ityp(nat))
  celldm = Sinfo%celldm
  at     = Sinfo%at
  bg     = Sinfo%bg
  omega  = Sinfo%omega
  atm(1:ntyp)    = Sinfo%atm(1:ntyp)
  amass(1:ntyp)  = Sinfo%amass(1:ntyp)
  tau(:,1:nat)   = Sinfo%tau(:,1:nat)
  ityp(1:nat)    = Sinfo%ityp(1:nat)
  !CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  !at = at / celldm(1)  !  bring at in units of alat
  !CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
  !CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  ! ######################### symmetry setup #########################
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
  CALL set_sym_bl ( )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  IF(.not.allocated(m_loc))  THEN
    ALLOCATE(m_loc(3,nat))
    m_loc = 0._dp
  ENDIF
  
  CALL find_sym ( nat, tau, ityp, .false., m_loc )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !
  ! Find the reduced grid of q-points:
  skip_equivalence = .FALSE.
  time_reversal    = .TRUE.
  nq1 = 8
  nq2 = 8
  nq3 = 8
  nqmax = nq1*nq2*nq3
  ALLOCATE(x_q(3,nqmax), w_q(nqmax))
  call kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, bg, nqmax,&
                    0,0,0, nq1,nq2,nq3, nqs, x_q, w_q )
  !
  WRITE(stdout, *) "Generated ", nqs, "points"

  ALLOCATE(rtau( 3, 48, nat), d2(3*nat,3*nat))
  ALLOCATE(w2(3*nat))
  Q_POINTS_LOOP : &
  DO iq = 1, nqs
    WRITE(stdout, *) "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
    WRITE(stdout, '(i6, 3f12.4)') iq, x_q(:,iq)

  
    ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
    ! part 1: call smallg_q and the copy_sym, 
    xq = x_q(:,iq)
    minus_q = .true.
  
    sym = .false.
    sym(1:nsym) = .true.
    CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q)
    nsymq = copy_sym(nsym, sym)
    ! recompute the inverses as the order of sym.ops. has changed
    CALL inverse_s ( ) 
  
  
  
    ! part 2: this computes gi, gimq
    call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
    WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
    IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G:"
    !
    ! finally this does some of the above again and also computes rtau...
    CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
    !
    ! the next subroutine uses symmetry from global variables to find he basis of crystal-symmetric
    ! matrices at this q point
    CALL find_d2_symm_base(xq, rank, basis)
    !
    ! Calculate the list of points making up the star of q
    CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, .true. )
    IF(imq/=0) THEN
       nq_out = nqs
       ALLOCATE(star_dyn(3*nat,3*nat, nqs))
    ELSE ! In this cas -q is not in the star of q, but I can recover it by time-reversal D(-q) = D(q)^*
       nq_out = 2*nqs
       IF(nq_out>48) CALL errore("make_wedge","unexpected imq=0 and nqs>48/2",1)
       ALLOCATE(star_dyn(3*nat,3*nat, 2*nqs))
       sxq(:,nqs+1:2*nqs) = -sxq(:,1:nqs)
    ENDIF

    ! Interpolate the system dynamical matrix at this q
    CALL fftinterp_mat2(xq, Sinfo, fc, d2)
    ! Remove the mass factor, I cannot remove it before because al the effective
    ! charges code assumes it is there
    d2 = multiply_mass_dyn(Sinfo,d2)
    !
    star_dyn = 6666._dp
!      ! Rotate the dynamical matrices making up the star of q it on file (i actually want to have it in a variable!)
    CALL make_qstar_d2 (d2, at, bg, nat, nsym, s, invs, irt, rtau, &
                         nqs, sxq, isq, imq, nq_out, star_dyn)
    DO i = 1, nq_out
       WRITE(stdout,'(i4,3f12.4,l2, 1f15.9)') i, sxq(:,i), (i>nqs), dotprodmat(3*nat,star_dyn(:,:,i), star_dyn(:,:,i))
    ENDDO
    DEALLOCATE(star_dyn)
!    ENDDO
  
    !CALL compact_dyn(nat, d2, phi)
    !  print*, d2
    print*, "== DECOMPOSITION =="
    DO i = 1,rank
      print*, dotprodmat(3*nat,d2, basis(:,:,i))
    ENDDO
    !
    !d2 = 0
    !DO i = 1,rank
    !  d2 = d2 + decomposition(i)*basis(:,:,i)
    !ENDDO

  ENDDO Q_POINTS_LOOP
  ! print*, d2

  !DEALLOCATE(phi, d2, w2)
  !DEALLOCATE(rtau, tau, ityp)
  !DEALLOCATE(irt) ! from symm_base
  !----------------------------------------------------------------------------
 END PROGRAM make_wedge
!----------------------------------------------------------------------------
!
