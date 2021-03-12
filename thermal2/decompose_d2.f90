
MODULE decompose_d2 
   USE kinds, ONLY : dp

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

CONTAINS

SUBROUTINE allocate_sym_and_star_q(nat, symq)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: nat
   TYPE(sym_and_star_q),INTENT(inout) :: symq
   ALLOCATE(symq%rtau(3,48,nat))
   ALLOCATE(symq%irt(48,nat))
END SUBROUTINE

! Generate a base for the space of dynamical matrices that respect the
! crystal symmetry
!---------------------------------------------------------------------
subroutine find_d2_symm_base(xq, rank, basis, nat, at, bg, &
                             nsymq, minus_q, irotmq, rtau, irt, s, invs )
!---------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  !USE ions_base, ONLY : nat !, tau, ntyp => nsp, ityp, amass
  !USE cell_base, ONLY : at, bg
  !USE symm_base, ONLY : s, invs, irt
  !USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, rtau

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

! input: the q point

!  INTEGER, INTENT(OUT) :: npert(3*nat), nirr
  REAL(DP)   :: eigen(3*nat)
  complex(DP):: u(3*nat, 3*nat)
  
  integer :: i, j, k, nx, jx, na, nb

  complex(DP) :: wdyn (3, 3, nat, nat), phi (3 * nat, 3 * nat)
  complex(DP) :: mtx(3*nat, 3*nat, 9*nat**2)
  real(DP)    :: normtx
!
!call write_matrix('random matrix',wdyn,nat)
!
! symmetrize the random matrix with the little group of q
!
   nx = 0  
   !first matrices with a single 1 along the diagonal
   mtx = 0._dp
   DO i = 1, 3*nat
     nx=nx+1
     mtx(i,i, nx) = 1._dp
   ENDDO
   !second, matrices with a single non-zero off-diagonal real term
   DO i = 1,3*nat
   DO j = i+1, 3*nat
!   DO j = 1, 3*nat
     nx=nx+1
     mtx(i,j,nx) = dsqrt(0.5_dp)
     mtx(j,i,nx) = dsqrt(0.5_dp)
   ENDDO
   ENDDO
   !third, matrices with a single non-zero off-diagonal imaginary term
   DO i = 1,3*nat
   DO j = i+1, 3*nat
     nx=nx+1
     mtx(i,j,nx) =  CMPLX(0._dp, dsqrt(0.5_dp), kind=dp)
     mtx(j,i,nx) = -CMPLX(0._dp, dsqrt(0.5_dp), kind=dp)
   ENDDO
   ENDDO
   WRITE(*,*) "Initial basis size:", nx
  
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
   WRITE(*,*) "Number of non-zero symmetric matrices::", jx
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
   PRINT*, "Number of orthonormal matrices:", jx
   nx = jx

  rank = nx
  ALLOCATE(basis(3*nat,3*nat,rank))
  basis(:,:,1:rank) = mtx(:,:,1:rank)

#ifdef __DEBUG
  phi = 0._dp
   WRITE(*,*) "++++++++++"
   DO i = 1, nx
     IF(ANY(ABS(mtx(:,:,i))>1.d-8))THEN
       phi  = phi+mtx(:,:,i)
       CALL scompact_dyn(nat, mtx(:,:,i), wdyn)
       WRITE(*,*) "+", i
       DO j = 1,nat
       DO k = 1,nat
            WRITE(*,*) j,k
            WRITE(*,'(3(2f7.3,3x),5x)') wdyn(:,:,j,k)
       ENDDO
       ENDDO
       call cdiagh (3 * nat, phi, 3 * nat, eigen, u)
       DO k = 1,3*nat
         IF(abs(eigen(k))>1.d-8) THEN
            WRITE(*,'("a", f10.6, 3x, 6(2f7.3,3x))') eigen(k), u(:,k)
         ENDIF
       ENDDO
     ENDIF
   ENDDO

   WRITE(*,*) "+", "Decomposition done"
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
        print*, "iq done", iq+nq
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
  USE io_global,  ONLY : stdout
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

END MODULE

