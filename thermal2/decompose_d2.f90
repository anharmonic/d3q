MODULE decompose_d2 
   USE kinds, ONLY : dp
#include "mpi_thermal.h"

   TYPE sym_and_star_q
      real(DP) :: xq (3)
      ! xq: q vector
      integer :: s (3, 3, 48), invs(48)
      ! invs: list of inverse operation indices
      integer :: nrot, nsym, nsymq, irotmq
      ! nrot  symmetry operations of the lattice
      ! nsym  symmetry operations of the crystal
      ! nsymq symmetry operations of the q-point
      ! index of the rotation that send q -> -q
      logical :: minus_q
      !
      integer :: nq_star, nq_trstar, isq (48), imq
      ! nq  : degeneracy of the star of q
      ! nq_tr  : degeneracy of the star of q and -q
      ! isq : index of q in the star for a given sym
      ! imq : index of -q in the star (0 if not present)
   
      real(DP) :: sxq (3, 48)
      ! list of vectors in the star
      !
      real(DP),allocatable :: rtau(:,:,:) !3,48,nat
      ! position of the rotated atoms
      integer,allocatable :: irt(:,:)
      ! the rotated of each atom
   END TYPE
   !
   TYPE dynmat_basis
      COMPLEX(DP),ALLOCATABLE :: basis(:,:,:)
   END TYPE


CONTAINS

SUBROUTINE allocate_sym_and_star_q(nat, symq)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: nat
   TYPE(sym_and_star_q),INTENT(inout) :: symq
   ALLOCATE(symq%rtau(3,48,nat))
   ALLOCATE(symq%irt(48,nat))
END SUBROUTINE


subroutine generate_simple_base(ndim, mtx, nx)
   implicit none
   integer,intent(in) :: ndim
   integer,intent(out) :: nx
   complex(DP),intent(out) :: mtx(ndim, ndim, ndim**2)
   integer :: i,j

   nx = 0  
   !first matrices with a single 1 along the diagonal
   mtx = 0._dp
   DO i = 1, ndim
     nx=nx+1
     mtx(i,i, nx) = 1._dp
   ENDDO
   !second, matrices with a single non-zero off-diagonal real term
   DO i = 1,ndim
   DO j = i+1, ndim
     nx=nx+1
     mtx(i,j,nx) = dsqrt(0.5_dp)
     mtx(j,i,nx) = dsqrt(0.5_dp)
   ENDDO
   ENDDO
   !third, matrices with a single non-zero off-diagonal imaginary term
   DO i = 1,ndim
   DO j = i+1, ndim
     nx=nx+1
     mtx(i,j,nx) =  CMPLX(0._dp, dsqrt(0.5_dp), kind=dp)
     mtx(j,i,nx) = -CMPLX(0._dp, dsqrt(0.5_dp), kind=dp)
   ENDDO
   ENDDO
   IF(nx>ndim**2) CALL errore("gen_sbase", "too many matrices?", 1)

   ioWRITE(stdout,'(2x,a,i8)') "Initial basis size:", nx

end subroutine

subroutine generate_mu_base(ndim, mtx, u0, nx)
   USE functions,          ONLY : cdiag_serial
   implicit none
   integer,intent(in)      :: ndim
   complex(dp),INTENT(in)  :: u0(ndim,ndim)
   integer,intent(out)     :: nx
   complex(DP),intent(out) :: mtx(ndim, ndim, ndim**2)
   !
   integer :: i,j
   real(dp)  :: e(ndim)
   complex(dp)  :: u(ndim,ndim), ev(ndim,ndim), v1(ndim,1), v2(1,ndim)

   ! diagonalize u0
   u=u0
   call cdiag_serial (ndim, u, ndim, e, ev)

   nx = 0  
   !first matrices with a single 1 along the diagonal
   mtx = 0._dp
   DO i = 1, ndim
   DO j = 1, ndim
     IF(ABS(e(i))>1.d-6 .and. ABS(e(j))>1.d-6 )THEN
      nx=nx+1
      v1(:,1) = ev(:,i)
      v2(1,:) = CONJG(ev(j,:))
      mtx(:,:, nx) = matmul(v1,v2)
     ENDIF
   ENDDO
   ENDDO

   IF(nx>ndim**2) CALL errore("gen_sbase", "too many matrices?", 1)

   ioWRITE(stdout,'(2x,a,i8)') "Initial basis size:", nx

end subroutine


! Generate a base for the space of dynamical matrices that respect the
! crystal symmetry
!---------------------------------------------------------------------
subroutine find_d2_symm_base(xq, rank, basis, nat, at, bg, &
                             nsymq, minus_q, irotmq, rtau, irt, s, invs, u0 )
!---------------------------------------------------------------------

  USE kinds,     ONLY : DP
  implicit none
!
!   first the dummy variables
!
  real(DP), INTENT(IN) :: xq (3)
  integer,INTENT(out) :: rank
  complex(DP),ALLOCATABLE,INTENT(out) :: basis(:,:,:)
  integer,INTENT(in)  :: nat
  real(DP),INTENT(in) :: at(3,3), bg(3,3)
  integer,INTENT(in)  :: nsymq
  logical,INTENT(in)  :: minus_q
  integer,INTENT(in)  :: irotmq
  real(DP),INTENT(in) :: rtau(3,48,nat)
  integer,INTENT(in)  :: irt(48,nat), s(3,3,48), invs(48)
  complex(dp),INTENT(in),optional :: u0(3*nat,3*nat)

! input: the q point

!  INTEGER, INTENT(OUT) :: npert(3*nat), nirr
  REAL(DP)   :: eigen(3*nat)
  complex(DP):: u(3*nat, 3*nat)
  
  integer :: i, j, k, nx, jx, na, nb

  complex(DP) :: wdyn (3, 3, nat, nat), phi (3 * nat, 3 * nat)
  complex(DP) :: mtx(3*nat, 3*nat, 9*nat**2)
  real(DP)    :: normtx
  print*, 10000040

   ! build an initial trivial basis for the hermitean matrices space
   IF(present(u0))THEN
      call generate_mu_base(3*nat, mtx, u0, nx)
   ELSE
      call generate_simple_base(3*nat, mtx, nx)
   ENDIF
   print*, 10000050

   ! symmetrize each matrix
   DO i = 1,nx
     CALL scompact_dyn(nat, mtx(:,:,i), wdyn)
     do na = 1, nat
        do nb = 1, nat
           call trntnsc( wdyn(:,:,na,nb), at, bg, -1 )
        enddo
     enddo
     call symdynph_gq_new (xq,wdyn,s,invs,rtau,irt,nsymq,nat,irotmq,minus_q)
     do na = 1, nat
        do nb = 1, nat
           call trntnsc( wdyn(:,:,na,nb), at, bg, +1 )
        enddo
     enddo
     CALL compact_dyn(nat, mtx(:,:,i), wdyn)
     !WRITE(*,'(6(2f7.3,3x),5x)') mtx(:,:,i)
   ENDDO
   print*, 10000051

   ! Some matrices can be zero at this point, we throw them away
   jx = 0
   DO i = 1,nx
     normtx = dotprodmat(3*nat, mtx(:,:, i),mtx(:,:, i))
     IF(normtx>1.d-12)THEN
      jx = jx+1
      mtx(:,:,jx) = mtx(:,:, i)/ DSQRT(normtx)
      !print*, i, jx, normtx
     ENDIF
   ENDDO
   ioWRITE(stdout,'(2x,a,i8)') "Number of non-zero symmetric matrices: ", jx
   nx = jx

   ! Graham-Schmidt
   jx = 0
   DO i = 1, nx
     phi = mtx(:,:,i)
     DO j = 1, jx !i-1 ! I only orthogonalize w.r.t the non-zero matrices that I have found
         normtx = dotprodmat(3*nat, mtx(:,:,j), mtx(:,:,j))
         mtx(:,:,i) = mtx(:,:,i) - mtx(:,:,j) * dotprodmat(3*nat, mtx(:,:,i), mtx(:,:,j))
!        mtx(:,:,i) = mtx(:,:,i) - mtx(:,:,j) * dotprodmat(3*nat, phi, mtx(:,:,j))
     ENDDO
     normtx = dotprodmat(3*nat, mtx(:,:,i), mtx(:,:,i))
     IF(normtx > 1.d-12)THEN
       jx = jx+1
       mtx(:,:,jx) = mtx(:,:,i)/DSQRT(normtx)
       IF(jx<i) mtx(:,:,i) = 0._dp
       !print*, jx, i, normtx
     ENDIF
   ENDDO
   ioWRITE(stdout,'(2x,a,i8)') "Number of orthonormal matrices:", jx
   nx = jx

  rank = nx
  ALLOCATE(basis(3*nat,3*nat,rank))
  basis(:,:,1:rank) = mtx(:,:,1:rank)

#ifdef __DEBUG
  phi = 0._dp
   ioWRITE(*,*) "++++++++++"
   DO i = 1, nx
     IF(ANY(ABS(mtx(:,:,i))>1.d-8))THEN
       phi  = phi+mtx(:,:,i)
       CALL scompact_dyn(nat, mtx(:,:,i), wdyn)
       ioWRITE(*,*) "+", i
       DO j = 1,nat
       DO k = 1,nat
            ioWRITE(*,*) j,k
            ioWRITE(*,'(3(2f7.3,3x),5x)') wdyn(:,:,j,k)
       ENDDO
       ENDDO
       call cdiagh (3 * nat, phi, 3 * nat, eigen, u)
       DO k = 1,3*nat
         IF(abs(eigen(k))>1.d-8) THEN
            ioWRITE(*,'("a", f10.6, 3x, 6(2f7.3,3x))') eigen(k), u(:,k)
         ENDIF
       ENDDO
     ENDIF
   ENDDO

   ioWRITE(*,*) "+", "Decomposition done"
#endif


!   WRITE(*,*) "total:" 
!   phi = phi / DSQRT(DBLE(nx))
!   CALL scompact_dyn(nat, phi, wdyn)
!   DO j = 1,nat
!   DO k = 1,nat
!        WRITE(*,*) j,k
!        WRITE(*,'(3(2f7.3,3x),5x)') wdyn(:,:,j,k)
!   ENDDO
!   ENDDO
!   call cdiagh (3 * nat, phi, 3 * nat, eigen, u)
!   DO k = 1,3*nat
!     !IF(abs(eigen(k))>1.d-8) THEN
!        WRITE(*,'("a", f10.6, 3x, 6(2f7.3,3x))') eigen(k), u(:,k)
!     !ENDIF
!   ENDDO

   return
end subroutine find_d2_symm_base

function dotprodmat(n,a,b) result(r)
  use kinds, only : dp
  implicit none
  integer,intent(in) :: n
  complex(dp),intent(in) :: a(n,n), b(n,n)
  !complex(dp) :: c(n,n)
  real(dp) :: r
  complex(dp) :: z
  integer :: i,j
  !c = 0._dp
  !call zgemm('N','C',n,n,n,1.0d0,a,n,b,n,0.0d0,c,n)
  !c = matmul(a,conjg(transpose(b)))
  !c = matmul(conjg(transpose(b)),a)
  z=0._dp
  !do i = 1,n
  ! z = z+c(i,i)
  !enddo
  do i = 1, n
  do j = 1, n
    z = z + CONJG(b(j,i))*a(j,i)
  enddo
  enddo
  IF(ABS(imag(z))>1.d-8) CALL errore("dotprodmat", "unexpected imaginary part", 1)
  r = dble(z)
end function

!
! Given D(q), comutes the D matrics in the star of q
!-----------------------------------------------------------------------
subroutine make_qstar_d2 (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
     nq, sxq, isq, imq, nq_trstar, star_dyn, star_wdyn)
  !-----------------------------------------------------------------------
  ! Generates the dynamical matrices for the star of q a
  ! If there is a symmetry operation such that q -> -q +G then imposes on
  ! dynamical matrix those conditions related to time reversal symmetry.
  ! It will return nq_trstar dynamical matrices, with
  !   nq_trstar = nq      if -q is in the star of q
  !   nq_trstar = 2*nq    otherwise
  !
  USE kinds, only : DP
  USE io_dyn_mat, only : write_dyn_mat
  USE control_ph, only : xmldyn
  implicit none
  ! input variables
  integer,INTENT(in) :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       nq, isq (48), imq, nq_trstar
  ! number of atoms in the unit cell
  ! number of symmetry operations
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom
  ! degeneracy of the star of q
  ! symmetry op. giving the rotated q
  ! index of -q in the star (0 if non present)
  ! unit number
  complex(DP),INTENT(in) :: dyn (3 * nat, 3 * nat)
  ! the input dynamical matrix. if imq.ne.0 the
  complex(DP),OPTIONAL,INTENT(out) :: star_dyn (3 * nat, 3 * nat, nq_trstar)
  complex(DP),OPTIONAL,INTENT(out) :: star_wdyn (3,3,nat,nat, nq_trstar)
  ! output matrices 

  real(DP),INTENT(in) :: at (3, 3), bg (3, 3), rtau (3, 48, nat), sxq (3, 48)
  ! direct lattice vectors
  ! reciprocal lattice vectors
  ! for each atom and rotation gives the R vector involved
  ! list of q in the star
  !
  !  local variables
  integer :: na, nb, iq, nsq, isym, icar, jcar, i, j, counter
  ! counters
  ! nsq: number of sym.op. giving each q in the list

  complex(DP) :: phi (3, 3, nat, nat), phi2 (3, 3, nat, nat)
  ! work space
  counter=0
  !
  ! Sets number of symmetry operations giving each q in the list
  !
  nsq = nsym / nq
  if (nsq * nq /= nsym) then
     print*, nsq, nq, nsym
     call errore ('q2star_ph', 'wrong degeneracy', 1)
  endif
  !
  ! Writes dyn.mat. dyn(3*nat,3*nat) on the 4-index array phi(3,3,nat,nat)
  !
  CALL scompact_dyn(nat, dyn, phi)
  !
  ! Go to crystal coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, - 1)
     enddo
  enddo
  !
  ! If -q is in the list impose first of all the conditions coming from
  ! time reversal symmetry
  !
  if (imq /= 0) then
     phi2 (:,:,:,:) = (0.d0, 0.d0)
     isym = 1
     do while (isq (isym) /= imq)
        isym = isym + 1
     enddo
     call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
          rtau, sxq (1, imq) )
     do na = 1, nat
        do nb = 1, nat
           do i = 1, 3
              do j = 1, 3
                 phi (i, j, na, nb) = 0.5d0 * (phi (i, j, na, nb) + &
                                        CONJG(phi2(i, j, na, nb) ) )
              enddo
           enddo
        enddo
     enddo
     !phi2 (:,:,:,:) = phi (:,:,:,:)
     !!
     !! Back to cartesian coordinates
     !!
     !do na = 1, nat
     !   do nb = 1, nat
     !      call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
     !   enddo
     !enddo
     !
     ! Saves 4-index array phi2(3,3,nat,nat) on the dyn.mat. dyn(3*nat,3*nat)
     !
     !CALL compact_dyn(nat, dyn, phi2)
  endif
  !
  ! For each q of the star rotates phi with the appropriate sym.op. -> phi
  !
  do iq = 1, nq
     phi2 (:,:,:,:) = (0.d0, 0.d0)
     do isym = 1, nsym
        if (isq (isym) == iq) then
           call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
                rtau, sxq (1, iq) )
        endif
     enddo
     phi2 (:,:,:,:) = phi2 (:,:,:,:) / DBLE (nsq)
     !
     ! Back to cartesian coordinates
     !
     do na = 1, nat
        do nb = 1, nat
           call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
        enddo
     enddo
     !
     ! Compact and store for output
     IF(present(star_dyn))  CALL compact_dyn(nat, star_dyn(:,:,iq), phi2)
     IF(present(star_wdyn)) star_wdyn(:,:,:,:, iq)= phi2
     !write(*,'(a,i5,999(" (",2f12.7,") "))') "iq done", iq, star_dyn(:,:,iq)
     !
     ! TODOO : also produce the star of -q when it is not in the star of q
     if (imq == 0) then
        !
        ! if -q is not in the star recovers its matrix by time reversal
        !
        do na = 1, nat
           do nb = 1, nat
              do i = 1, 3
                 do j = 1, 3
                    phi2 (i, j, na, nb) = CONJG(phi2 (i, j, na, nb) )
                 enddo
              enddo
           enddo
        enddo
        !
        ! and writes it 
        !print*, "iq done", iq+nq
        IF(present(star_dyn)) CALL compact_dyn(nat, star_dyn(:,:,iq+nq), phi2)
        IF(present(star_wdyn)) star_wdyn(:,:,:,:, iq)= phi2
     endif
  enddo
  !
  return
end subroutine make_qstar_d2

! Compute the star of q (by calling star_q) and if time reversal is not
! a symmetriic of the small group of q, all apply it to obtain the star
! of -q
!-----------------------------------------------------------------------
subroutine tr_star_q (xq, at, bg, nsym, s, invs, nq, nq_tr, sxq, isq, imq, verbosity)
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE kinds, only : DP
  implicit none
  integer, intent(in) :: nsym, s (3, 3, 48), invs(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  real(DP), intent(in) :: xq (3), at (3, 3), bg (3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  integer, intent(out) :: nq, nq_tr, isq (48), imq
  ! nq  : degeneracy of the star of q
  ! nq  : degeneracy of the star of q and -q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  real(DP), intent(out) :: sxq (3, 48)
  ! list of vectors in the star of q
  logical, intent(in) :: verbosity
  ! if true prints several messages.

  CALL star_q(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, verbosity )
 
  IF(imq==0) THEN
     nq_tr = 2*nq
     IF(nq>48) CALL errore("make_wedge","unexpected imq=0 and nq_star>48/2",1)
     sxq(:,nq+1:2*nq) = -sxq(:,1:nq)
  ELSE
     nq_tr = nq
  ENDIF

end subroutine


! The original smallg_q, only looks for minus_q in a very specific case, 
! this subroutine always does.
!-----------------------------------------------------------------------
SUBROUTINE smallg_q_fullmq (xq, modenum, at, bg, nrot, s, sym, minus_q)
   !-----------------------------------------------------------------------
   !
   ! This routine selects, among the symmetry matrices of the point group
   ! of a crystal, the symmetry operations which leave q unchanged.
   ! Furthermore it checks if one of the above matrices send q --> -q+G.
   ! In this case minus_q is set true.
   !
   !  input-output variables
   !
   USE kinds, ONLY : DP
   USE symm_base, ONLY : t_rev
   
   implicit none
 
   real(DP), parameter :: accep = 1.e-5_dp
 
   real(DP), intent(in) :: bg (3, 3), at (3, 3), xq (3)
   ! input: the reciprocal lattice vectors
   ! input: the direct lattice vectors
   ! input: the q point of the crystal
 
   integer, intent(in) :: s (3, 3, 48), nrot, modenum
   ! input: the symmetry matrices
   ! input: number of symmetry operations
   ! input: main switch of the program, used for
   !        q<>0 to restrict the small group of q
   !        to operation such that Sq=q (exactly,
   !        without G vectors) when iswitch = -3.
   logical, intent(inout) :: sym (48), minus_q
   ! input-output: .true. if symm. op. S q = q + G
   ! output: .true. if there is an op. sym.: S q = - q + G
   !
   !  local variables
   !
 
   real(DP) :: aq (3), raq (3), zero (3)
   ! q vector in crystal basis
   ! the rotated of the q vector
   ! the zero vector
 
   integer :: irot, ipol, jpol
   ! counter on symmetry op.
   ! counter on polarizations
   ! counter on polarizations
 
   logical :: eqvect
   ! logical function, check if two vectors are equa
   !
   ! return immediately (with minus_q=.true.) if xq=(0,0,0)
   !
   minus_q = .true.
   if ( (xq (1) == 0.d0) .and. (xq (2) == 0.d0) .and. (xq (3) == 0.d0) ) &
        return
   !
   !   Set to zero some variables
   !
   minus_q = .false.
   zero(:) = 0.d0
   !
   !   Transform xq to the crystal basis
   !
   aq = xq
   call cryst_to_cart (1, aq, at, - 1)
   !
   !   Test all symmetries to see if this operation send Sq in q+G or in -q+G
   !
   do irot = 1, nrot
      if (.not.sym (irot) ) goto 100
      raq(:) = 0.d0
      do ipol = 1, 3
         do jpol = 1, 3
            raq(ipol) = raq(ipol) + DBLE( s(ipol,jpol,irot) ) * aq( jpol)
         enddo
      enddo
      IF (t_rev(irot)==1) raq=-raq
      sym (irot) = eqvect (raq, aq, zero, accep)
      !
      !  if "iswitch.le.-3" (modenum.ne.0) S must be such that Sq=q exactly !
      !
      if (modenum.ne.0 .and. sym(irot) ) then
         do ipol = 1, 3
            sym(irot) = sym(irot) .and. (abs(raq(ipol)-aq(ipol)) < 1.0d-5)
         enddo
      endif
      if (.not.minus_q) then
 !     if (sym(irot).and..not.minus_q) then
         raq = - raq
         minus_q = eqvect (raq, aq, zero, accep)
      endif
 100  continue
   enddo
   !
   !  if "iswitch.le.-3" (modenum.ne.0) time reversal symmetry is not included !
   !
   if (modenum.ne.0) minus_q = .false.
   !
   return
   !
 END SUBROUTINE smallg_q_fullmq
 
 SUBROUTINE recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_coef, nq1, nq2, nq3, nqmax, nfar, fcout)
   USE kinds, ONLY : DP
   USE input_fc, ONLY : read_fc2, forceconst2_grid, ph_system_info
   USE quter_module,       ONLY : quter

   IMPLICIT NONE
   TYPE(ph_system_info),INTENT(in)   :: Si
   TYPE(sym_and_star_q),INTENT(in) :: symq(nq_wedge)
   TYPE(forceconst2_grid),INTENT(inout) :: fcout
   TYPE(dynmat_basis),INTENT(in) :: dmb(nq_wedge)
   INTEGER,INTENT(in) :: rank(nq_wedge)
   REAL(DP),INTENT(in) :: ph_coef(nph)

   INTEGER, INTENT(in) :: nq1, nq2, nq3, nqmax, nq_wedge, nfar, nph
   !
   INTEGER :: i, iph, iq, nq_done
   COMPLEX(DP) :: d2(Si%nat3, Si%nat3)
   !
   COMPLEX(DP),ALLOCATABLE :: star_wdyn(:,:,:,:,:), star_dyn(:,:,:)
   REAL(DP),ALLOCATABLE :: xqmax(:,:)

   ALLOCATE(star_wdyn(3,3,Si%nat,Si%nat, nqmax))
   ALLOCATE(xqmax(3,nqmax))
   
   ! Reconstruct the dynamical matrix from the coefficients
   nq_done = 0
   iph = 0
   Q_POINTS_LOOP3 : &
   DO iq = 1, nq_wedge
      d2 = 0._dp
      DO i = 1,rank(iq)
         iph = iph+1
         d2 = d2+ ph_coef(iph)*dmb(iq)%basis(:,:,i)
      ENDDO
      ! WRITE(999,'(i3,3f12.6)') iq,symq(iq)%xq
      ! WRITE(999,'(3(2f12.6,4x))') d2

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

      DEALLOCATE(star_dyn)

   ENDDO Q_POINTS_LOOP3

   IF(iph.ne.nph) CALL errore("minimize", "wrong iph", 1)
   !
   CALL quter(nq1, nq2, nq3, Si%nat,Si%tau,Si%at,Si%bg, star_wdyn, xqmax, fcout, nfar)
   !
   DEALLOCATE(xqmax, star_wdyn)

   IF(nfar.ne.0) RETURN
   END SUBROUTINE

END MODULE

