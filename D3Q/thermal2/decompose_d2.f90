
MODULE decompose_d2 


CONTAINS
! Generate a base for the space of dynamical matrices that respect the
! crystal symmetry
!---------------------------------------------------------------------
subroutine find_d2_symm_base(xq, rank, basis)
!---------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  USE ions_base, ONLY : nat, tau, ntyp => nsp, ityp, amass
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, sr, invs, nsym, irt, t_rev
!  USE modes,     ONLY : num_rap_mode, name_rap_mode
!  USE noncollin_module, ONLY : noncolin, nspin_mag
!  USE spin_orb,  ONLY : domag
!  USE constants, ONLY: tpi
!  USE control_ph, ONLY : search_sym
!  USE control_flags, ONLY : iverbosity
!  USE random_numbers, ONLY : randy
!  USE rap_point_group, ONLY : name_rap

!  use mp, only: mp_bcast
!  use io_global, only : ionode_id
!  use mp_images, only : intra_image_comm

  USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, gi, gimq, rtau
  USE control_lr,   ONLY : lgamma

  implicit none
!
!   first the dummy variables
!
  real(DP), INTENT(IN) :: xq (3)
  integer,INTENT(out) :: rank
  complex(DP),ALLOCATABLE,INTENT(out) :: basis(:,:,:)
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
   WRITE(*,*) "NX::", nx
  
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
      print*, i, jx, normtx
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
       print*, jx, i, normtx
     ENDIF
   ENDDO
   PRINT*, "Number of orthonormal matrices:", jx
   nx = jx

  rank = nx
  ALLOCATE(basis(3*nat,3*nat,rank))
  basis(:,:,1:rank) = mtx(:,:,1:rank)

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
!   WRITE(*,*) "-", 1
       ENDDO
!   WRITE(*,*) "-", 2
       call cdiagh (3 * nat, phi, 3 * nat, eigen, u)
!   WRITE(*,*) "-", 3
       DO k = 1,3*nat
         !IF(abs(eigen(k))>1.d-8) THEN
            WRITE(*,'("a", f10.6, 3x, 6(2f7.3,3x))') eigen(k), u(:,k)
         !ENDIF
       ENDDO
     ENDIF
   ENDDO

   WRITE(*,*) "+", "Decomposition done"

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

END MODULE
