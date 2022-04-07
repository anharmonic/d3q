!
! Written by Lorenzo Paulatto (2013-2022) IMPMC @ UPMC / CNRS UMR7590
!  and Ibrahim G. garba
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE decompose_zstar
   USE kinds, ONLY : dp
#include "mpi_thermal.h"

CONTAINS

   ! Copy of symtensor that does not use global variables
   !--------------------------------------------------------------------------
SUBROUTINE symtensor_zstar( nat, tens, at, bg, nsym, s, irt )
   !-----------------------------------------------------------------------
   !! Symmetrize a function \(f(i,j,na)\) (e.g. the effective charges in 
   !! cartesian axis), where \(i,j\) are the cartesian components and \(na\)
   !! is the atom index.
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: nat
   !! number of atoms
   REAL(DP), INTENT(INOUT) :: tens(3,3,nat)
   REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
   !! tensor function to symmetrize
   INTEGER :: nsym
   integer,INTENT(in)  :: irt(48,nat), s(3,3,48)
   !
   ! ... local variables
   !
   INTEGER :: na, isym, nar, i,j,k,l
   REAL(DP), ALLOCATABLE :: work (:,:,:)
   !
   IF (nsym == 1) RETURN
   !
   ! bring tensor to crystal axis
   !
   DO na=1,nat
      CALL cart_to_crys ( tens (:,:,na), at )
   END DO
   !
   ! symmetrize in crystal axis
   !
   ALLOCATE (work(3,3,nat))
   work (:,:,:) = 0.0_dp
   DO na = 1, nat
      DO isym = 1, nsym
         nar = irt (isym, na)
         DO i = 1, 3
            DO j = 1, 3
               DO k = 1, 3
                  DO l = 1, 3
                     work (i,j,na) = work (i,j,na) + &
                        s (i,k,isym) * s (j,l,isym) * tens (k,l,nar)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   tens (:,:,:) = work (:,:,:) / DBLE(nsym)
   DEALLOCATE (work)
   !
   ! bring tensor back to cartesian axis
   !
   DO na=1,nat
      CALL crys_to_cart ( tens (:,:,na), bg )
   END DO
   !
   !
 END SUBROUTINE symtensor_zstar
   !-------------------------------------------------------------------------
 SUBROUTINE crys_to_cart( matr, bg )
   !-----------------------------------------------------------------------
   !! Crystal to cartesian axis conversion.
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(INOUT) :: matr(3,3)
   REAL(DP), INTENT(IN) :: bg(3,3)
   !! Axis conversion matrix
   !
   ! ... local variables
   !
   REAL(DP) :: work(3,3)
   INTEGER :: i,j,k,l
   !
   work(:,:) = 0.0_dp
   DO i = 1, 3
      DO j = 1, 3
         DO k = 1, 3
            DO l = 1, 3
               work(i,j) = work(i,j) + &
                           matr(k,l) * bg(i,k) * bg(j,l)
            END DO
         END DO
      END DO
   END DO
   matr(:,:) = work(:,:)
   !
 END SUBROUTINE crys_to_cart
   !-------------------------------------------------------------------------
 SUBROUTINE cart_to_crys( matr, at )
   !-----------------------------------------------------------------------
   !! Cartesian to crystal axis conversion.
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(INOUT) :: matr(3,3)
   REAL(DP), INTENT(IN) :: at(3,3)
   !! Axis conversion matrix
   !
   ! ... local variables
   !
   REAL(DP) :: work(3,3)
   INTEGER :: i,j,k,l
   !
   work(:,:) = 0.0_dp
   DO i = 1, 3
      DO j = 1, 3
         DO k = 1, 3
            DO l = 1, 3
               work(i,j) = work(i,j) + matr(k,l) * at(k,i) * at(l,j)
            END DO
         END DO
      END DO
   END DO
   !
   matr(:,:) = work(:,:)
   !
 END SUBROUTINE cart_to_crys


subroutine generate_zsimple_base(nat, zxx, nx)
   implicit none
   integer,intent(in) :: nat
   integer,intent(out) :: nx
   real(DP),intent(out) :: zxx(3,3,nat, 6*nat)
   integer :: i,j,n

   zxx = 0._dp
   nx = 0
   DO n = 1,nat
      DO i =1,3
         nx = nx+1
         zxx(i,i,n,nx) = 1._dp
      ENDDO
      DO i =1,2
      DO j =i+1,3
         nx = nx+1
         zxx(i,j,n,nx) = dsqrt(0.5_dp)
         zxx(j,i,n,nx) = dsqrt(0.5_dp)
      ENDDO
      ENDDO
   ENDDO
   IF(nx>6*nat) CALL errore("gen_sbase", "too many zeu?", 1)
end subroutine

subroutine generate_zstar_base(nat, zxx, z0, nx)
   implicit none
   integer,intent(in) :: nat
   integer,intent(out) :: nx
   real(DP),intent(in)  :: z0(3,3,nat, 6*nat)
   real(DP),intent(out) :: zxx(3,3,nat, 6*nat)
   integer :: i,j,n

   ! Only create basis for elements which were non-zero to begin with
   zxx = 0._dp
   nx = 0
   DO n = 1,nat
      DO i =1,3
         IF(ABS(z0(i,i,n,nx))>0._dp)THEN
            nx = nx+1
            zxx(i,i,n,nx) = 1._dp
         ENDIF
      ENDDO
      DO i =1,2
      DO j =i+1,3
         IF(ABS(z0(i,j,n,nx))>0._dp)THEN
            nx = nx+1
            zxx(i,j,n,nx) = dsqrt(0.5_dp)
            zxx(j,i,n,nx) = dsqrt(0.5_dp)
         ENDIF
      ENDDO
      ENDDO
   ENDDO
   IF(nx>6*nat) CALL errore("gen_sbase", "too many zeu?", 1)
end subroutine


! Generate a base for the space of dynamical zstar that respect the
! crystal symmetry
!---------------------------------------------------------------------
subroutine find_zstar_symm_base(rank, zbasis, nat, at, bg, nsym, irt, s, z0, method )
!---------------------------------------------------------------------

  USE kinds,     ONLY : DP
  implicit none
!
!   first the dummy variables
  integer,INTENT(out) :: rank
  real(DP),ALLOCATABLE,INTENT(out) :: zbasis(:,:,:,:)
  integer,INTENT(in)  :: nat
  REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
  integer,INTENT(in)  :: nsym
  integer,INTENT(in)  :: irt(48,nat), s(3,3,48)
  real(dp),INTENT(in),optional :: z0(3,3,nat)
  CHARACTER(len=*),INTENT(in) :: method 

! input: the q point
  integer :: i, j, k, nx, jx, na, nb, nb1, nb2, nb3, nb4

  real(DP),allocatable :: zxx(:,:,:,:)
  real(DP)    :: norzxx, sq_norzxx_m1
  real(DP),parameter :: eps_base = 1.d-8

  ALLOCATE(zxx(3,3,nat,6*nat))

   IF (method=="zstar") THEN
      IF(.not. present(z0)) CALL errore("generate d2 base", 'u0 is required with "mu"', 1)
      call generate_zstar_base(nat, zxx, z0, nx)
   ELSEIF (method=="simple") THEN
      call generate_zsimple_base(nat, zxx, nx)
   ELSE
      CALL errore("generate zstar base", 'Unknown method (can only be: mu, simple)', 1)
   ENDIF
   nb1 = nx ! save for printing

   ! symmetrize each zstar
   DO i = 1,nx
      CALL symtensor_zstar( nat, zxx(:,:,:,i), at, bg, nsym, s, irt )
   ENDDO

   ! Some zstar can be zero at this point, we throw them away
   jx = 0
   DO i = 1,nx
     norzxx = dotprodzstar(nat, zxx(:,:,:, i),zxx(:,:,:, i))
     IF(norzxx > eps_base)THEN
      jx = jx+1
      sq_norzxx_m1 = 1._dp / DSQRT(norzxx)
      zxx(:,:,:,jx) = zxx(:,:,:, i) * sq_norzxx_m1
      CALL enforce_adj(nat, zxx(:,:,:,jx))
      !print*, i, jx, norzxx
     ENDIF
     !WRITE(*,'(2i4,100(3(3f4.1,x),2x))') i, jx, zxx(:,:,:,jx)
   ENDDO
   !ioWRITE(stdout,'(2x,a,i8)') "Number of non-zero symmetric zstar: ", jx
   nx = jx
   nb2 = nx ! save for printing

   ! Graham-Schmidt
   jx = 0
   DO i = 1, nx
     DO j = 1, jx !i-1 ! I only orthogonalize w.r.t the non-zero zstar that I have found
         norzxx = dotprodzstar(nat, zxx(:,:,:,j), zxx(:,:,:,j))
         zxx(:,:,:,i) = zxx(:,:,:,i) - zxx(:,:,:,j) * dotprodzstar(nat, zxx(:,:,:,i), zxx(:,:,:,j))
     ENDDO
     norzxx = dotprodzstar(nat, zxx(:,:,:,i), zxx(:,:,:,i))
     IF(norzxx > eps_base)THEN
       jx = jx+1
       sq_norzxx_m1 = 1._dp / DSQRT(norzxx)
       zxx(:,:,:,jx) = zxx(:,:,:,i) * sq_norzxx_m1
       CALL enforce_adj(nat, zxx(:,:,:,jx))
       IF(jx<i) zxx(:,:,:,i) = 0._dp
       !print*, jx, i, norzxx
     ENDIF
     !WRITE(*,'(2i4,100(3(3f4.1,x),2x))') i, jx, zxx(:,:,:,jx)
   ENDDO
   !ioWRITE(stdout,'(2x,a,i8)') "Number of orthonormal zstar:", jx
   nx = jx
   nb3 = nx ! save for printing

  ! Purge zstar that have zero projection on the provided initial zstar
   IF(present(z0)) THEN
      jx = 0
      DO i = 1, nx
      IF( ABS(dotprodzstar(nat, z0, zxx(:,:,:,i))) > 1.d-6 ) THEN
         jx = jx+1
         IF(jx<i) zxx(:,:,:,jx) = zxx(:,:,:,i)
      ENDIF
      ENDDO
      nx = jx
   ENDIF
   !ioWRITE(stdout,'(2x,a,i8)') "Number of purged orthonormal zstar:", jx
   nb4 = nx ! save for printing

  WRITE(stdout,'(2x,a,4i8)') &
      "  Z^star basis (initial/symmetrized/orthogonal/purged) : ",    &
      nb1, nb2, nb3, nb4
  rank = nx
  ALLOCATE(zbasis(3,3,nat,rank))
  zbasis(:,:,:,1:rank) = zxx(:,:,:,1:rank)

  DEALLOCATE(zxx)

   return
end subroutine find_zstar_symm_base

subroutine enforce_adj(n,a)
  use kinds, only : dp
  implicit none
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(3,3,n)
  integer :: i,j,k
  do k = 1,n
   do i = 1, 2
   do j = i+1, 3
      a(i,j,k) = 0.5_dp * (a(i,j,k)+ (a(j,i,k)))
      a(j,i,k) =  a(i,j,k)
   enddo
   enddo
  enddo
end subroutine

function dotprodzstar(n,a,b) result(r)
  use kinds, only : dp
  implicit none
  integer,intent(in) :: n
  real(dp),intent(in) :: a(3,3,n), b(3,3,n)
  real(dp) :: r
  integer :: i,j,k
  r = 0._dp
!$OMP PARALLELDO DEFAULT(shared) PRIVATE(i,j,k) REDUCTION(+:r)
  do i = 1,n
    do j = 1,3
      r = r+a(j,j,i)*b(j,j,i)
      do k = j+1,3
         r = r + 0.5_dp*(a(j,k,i)*b(j,k,i) + a(k,j,i)*b(k,j,i))
      enddo
    enddo
  enddo
!$OMP END PARALLELDO
  !
end function


! Generate a base for the space of dynamical zstar that respect the
! crystal symmetry
!---------------------------------------------------------------------
subroutine recompose_zstar(nat, rank, zbasis, zstar_coef, zstar)
   !---------------------------------------------------------------------
   USE kinds,     ONLY : DP
   implicit none
   !
   !   first the dummy variables
   integer,INTENT(in)  :: nat
   integer,INTENT(in) :: rank
   real(DP),INTENT(in) :: zbasis(3,3,nat,rank)
   REAL(DP),INTENT(in) :: zstar_coef(rank)
   real(dp),INTENT(out) :: zstar(3,3,nat)

   ! input: the q point
   integer :: i

   zstar = 0._dp
!$OMP PARALLELDO DEFAULT(shared) PRIVATE(i) REDUCTION(+:zstar)
   DO i = 1, rank
      zstar = zstar + zstar_coef(i)*zbasis(:,:,:,i)
   ENDDO
!$OMP END PARALLELDO

end subroutine

!---------------------------------------------------------------------
subroutine zstar_to_supercell(nat, nat_sc, zstar, zstar_sc)
   !---------------------------------------------------------------------
   USE kinds,     ONLY : DP
   implicit none
   !
   !   first the dummy variables
   integer,INTENT(in)  :: nat, nat_sc
   real(dp),INTENT(in) :: zstar(3,3,nat)
   real(dp),INTENT(out) :: zstar_sc(3,3,nat_sc)

   ! input: the q point
   integer :: i

   DO i = 1,nat_sc,nat
      zstar_sc(:,:,i:i+nat-1) = zstar
   ENDDO

end subroutine

END MODULE

