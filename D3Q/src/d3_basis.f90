!
! Copyright (C) 2001-2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE d3_basis
  USE kinds, only: DP
  !
  TYPE d3_pattern_type
    COMPLEX(DP),ALLOCATABLE :: u(:,:)
  END TYPE d3_pattern_type
  ! One of the above for each q vector
!  TYPE(d3_pattern_type) :: patq(-3:3) ! symmetry of each q individually (the same as phonon at that q)
  TYPE(d3_pattern_type),ALLOCATABLE :: patq(:) ! symmetry of each q individually (the same as phonon at that q)
  !
CONTAINS
!-----------------------------------------------------------------------
pure SUBROUTINE allocate_d3_pattern(nat, pattern)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nat
  TYPE(d3_pattern_type),INTENT(INOUT) :: pattern
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  ALLOCATE( pattern%u(3*nat, 3*nat) )
  pattern%u = 0._dp
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE allocate_d3_pattern
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE deallocate_d3_pattern(pattern)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(d3_pattern_type),INTENT(INOUT) :: pattern
  DEALLOCATE( pattern%u )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE deallocate_d3_pattern
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE d3_cart2pat(d3in, nat, u1, u2, u3, d3out)
  !-----------------------------------------------------------------------
  !
  !   Rotates third derivative of the dynamical basis from cartesian axis
  !   to the basis of the modes
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout)        :: d3in(3*nat, 3*nat, 3*nat)  ! d3 in cartesian basis
  COMPLEX(DP),OPTIONAl,INTENT(out) :: d3out(3*nat, 3*nat, 3*nat) ! d3 in pattern basis
  COMPLEX(DP),INTENT(in) :: u1(3*nat, 3*nat), u2(3*nat, 3*nat), u3(3*nat, 3*nat) ! patterns
  INTEGER,INTENT(in) :: nat
  !
  COMPLEX(DP) :: work
  COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
  INTEGER :: a, b, c, i, j, k
  !
  ALLOCATE(d3tmp(3*nat, 3*nat, 3*nat))
  !
  DO k = 1, 3*nat
  DO j = 1, 3*nat
  DO i = 1, 3*nat
      !
      work = (0._dp, 0._dp)
      !
      DO c = 1, 3 * nat
      DO b = 1, 3 * nat
      DO a = 1, 3 * nat
          work = work + d3in(a, b, c) * &
                        u1(i,a) * u2(j,b) * u3(k,c)
      ENDDO
      ENDDO
      ENDDO
      !
      d3tmp(i, j, k) = work
      !
  ENDDO
  ENDDO
  ENDDO
  !
  IF(present(d3out)) THEN
    d3out = d3tmp
  ELSE
    d3in  = d3tmp
  ENDIF
  DEALLOCATE(d3tmp)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_cart2pat
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE d3_pat2cart(d3in, nat, u1, u2, u3, d3out)
  !-----------------------------------------------------------------------
  !
  !   Rotates third derivative of the dynamical basis from the basis
  !   of modes to cartesian axis
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout)       :: d3in(3*nat, 3*nat, 3*nat)  ! d3 in pattern basis
  COMPLEX(DP),OPTIONAL,INTENT(out):: d3out(3*nat, 3*nat, 3*nat) ! d3 in cartesian basis
  COMPLEX(DP),INTENT(in) :: u1(3*nat, 3*nat), u2(3*nat, 3*nat), u3(3*nat, 3*nat) ! patterns
  INTEGER,INTENT(in) :: nat
  !
  COMPLEX(DP) :: work
  COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
  INTEGER :: i, j, k, a, b, c
  !
  ALLOCATE(d3tmp(3*nat, 3*nat, 3*nat))
  !
  DO i = 1, 3*nat
  DO j = 1, 3*nat
  DO k = 1, 3*nat
      !
      work = (0._dp, 0._dp)
      !
      DO a = 1, 3*nat
      DO b = 1, 3*nat
      DO c = 1, 3*nat
        work = work + d3in(a, b, c) * &
                CONJG( u1(i,a)*u2(j,b)*u3(k,c) )
      ENDDO
      ENDDO
      ENDDO
      !
      d3tmp(i,j,k) = work
      !
  ENDDO
  ENDDO
  ENDDO
  !
  IF(present(d3out)) THEN
    d3out = d3tmp
  ELSE
    d3in  = d3tmp
  ENDIF
  DEALLOCATE(d3tmp)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_pat2cart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE d3_3idx_2_6idx(nat, d3, phi)
  !-----------------------------------------------------------------------
  !
  !   Changes the index from (3*nat)^3 to 3^3 times nat^3
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(in)  :: d3(3*nat, 3*nat, 3*nat)   ! d3 in cartesian basis
  COMPLEX(DP),INTENT(out) :: phi(3,3,3, nat,nat,nat)    ! d3 in pattern basis
  INTEGER,INTENT(in) :: nat
  !
  INTEGER :: i, j, k
  INTEGER :: at_i, at_j, at_k
  INTEGER :: pol_i, pol_j, pol_k
  !
  DO k = 1, 3*nat
    at_k  = 1 + (k-1)/3
    pol_k = k - 3*(at_k-1)
    !
    DO j = 1, 3*nat
      at_j  = 1 + (j-1)/3
      pol_j = j - 3*(at_j-1)
      !
      DO i = 1, 3*nat
        at_i  = 1 + (i-1)/3
        pol_i = i - 3*(at_i-1)
        !
        phi(pol_i, pol_j, pol_k, at_i, at_j, at_k) &
          = d3(i,j,k)
        !
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_3idx_2_6idx
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE d3_6idx_2_3idx(nat, phi, d3)
  !-----------------------------------------------------------------------
  !
  !   Changes the index from 3^3 times nat^3 to (3*nat)^3 
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(out)  :: d3(3*nat, 3*nat, 3*nat)   ! d3 in cartesian basis
  COMPLEX(DP),INTENT(in) :: phi(3,3,3, nat,nat,nat) ! d3 in pattern basis
  INTEGER,INTENT(in) :: nat
  !
  INTEGER :: i, j, k
  INTEGER :: at_i, at_j, at_k
  INTEGER :: pol_i, pol_j, pol_k
  !
  DO k = 1, 3*nat
    at_k  = 1 + (k-1)/3
    pol_k = k - 3*(at_k-1)
    !
    DO j = 1, 3*nat
      at_j  = 1 + (j-1)/3
      pol_j = j - 3*(at_j-1)
      !
      DO i = 1, 3*nat
        at_i  = 1 + (i-1)/3
        pol_i = i - 3*(at_i-1)
        !
        d3(i,j,k) &
          = phi(pol_i, pol_j, pol_k, at_i, at_j, at_k)
        !
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_6idx_2_3idx
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure SUBROUTINE d3_cart2crys(nat, at, bg, phi)
  !-----------------------------------------------------------------------
  !
  !   Rotates third derivative of the dynamical basis from the basis
  !   of modes to cartesian axis
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout) :: phi(3,3,3, nat,nat,nat)  ! phi in pattern basis
  REAL(DP),INTENT(in)       :: at(3,3), bg(3,3)         ! axis, and reciprocal axis
  INTEGER,INTENT(in)        :: nat
  !
  INTEGER :: na, nb, nc
  !
  DO na = 1, nat
     DO nb = 1, nat
        DO nc = 1, nat
           CALL trntnsc_3 (phi(:,:,:, na, nb, nc), at, bg, - 1)
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_cart2crys
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
pure SUBROUTINE d3_crys2cart(nat, at, bg, phi)
  !-----------------------------------------------------------------------
  !
  !   Rotates third derivative of the dynamical basis from the basis
  !   of modes to cartesian axis
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout) :: phi(3,3,3, nat,nat,nat)  ! phi in pattern basis
  REAL(DP),INTENT(in)       :: at(3,3), bg(3,3)         ! axis, and reciprocal axis
  INTEGER,INTENT(in)        :: nat
  !
  INTEGER :: na, nb, nc
  !
  DO na = 1, nat
     DO nb = 1, nat
        DO nc = 1, nat
           CALL trntnsc_3 (phi(:,:,:, na, nb, nc), at, bg, + 1)
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_crys2cart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
pure subroutine trntnsc_3 (phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a COMPLEX third order tensor
  !(like the derivative of the dynamical matrix)
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
  USE kinds, only : DP
  implicit none

  integer,intent(in) :: iflg
  ! input: gives the versus of the trans.

  complex (DP), intent(inout) :: phi (3, 3, 3)
  ! inp/out: the matrix to transform

  real (DP), intent(in) :: at (3, 3), bg (3, 3)
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice

  integer :: i, j, k, l, m, n
  !
  !  counters on polarizations
  !
  complex (DP) :: wrk (3, 3, 3)
  ! a work array

  if (iflg.gt.0) then
     !
     ! forward transformation (crystal to cartesian axis)
     !

!      call zcopy (27, phi, 1, wrk, 1)
     wrk = phi
     do m = 1, 3
        do i = 1, 3
           do j = 1, 3
              phi (m, i, j) = (0.d0, 0.d0)
              do n = 1, 3
                 do k = 1, 3
                    do l = 1, 3
                       phi (m, i, j) = phi (m, i, j) + wrk (n, k, l) * bg (i, k) &
                            * bg (j, l) * bg (m, n)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  else
     !
     ! backward transformation (cartesian to crystal axis)
     !
     do m = 1, 3
        do i = 1, 3
           do j = 1, 3
              wrk (m, i, j) = (0.d0, 0.d0)
              do n = 1, 3
                 do k = 1, 3
                    do l = 1, 3
                       wrk (m, i, j) = wrk (m, i, j) + phi (n, k, l) * at (k, i) &
                            * at (l, j) * at (n, m)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
!      call zcopy (27, wrk, 1, phi, 1)
     phi = wrk
  endif
  return


end subroutine trntnsc_3

!-----------------------------------------------------------------------
END MODULE d3_basis
!-----------------------------------------------------------------------
