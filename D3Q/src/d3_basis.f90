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
  TYPE(d3_pattern_type) :: patq(-3:3) ! symmetry of each q individually (the same as phonon at that q)
  !
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE allocate_d3_pattern(nat, pattern)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nat
  TYPE(d3_pattern_type),INTENT(INOUT) :: pattern
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  ALLOCATE( pattern%u(3*nat, 3*nat) )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE allocate_d3_pattern
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE deallocate_d3_pattern(pattern)
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
SUBROUTINE d3_cart2pat(d3in, nat, u1, u2, u3, d3out)
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
  DO i = 1, 3*nat
  DO j = 1, 3*nat
  DO k = 1, 3*nat
      !
      work = (0._dp, 0._dp)
      !
      DO a = 1, 3 * nat
      DO b = 1, 3 * nat
      DO c = 1, 3 * nat
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
SUBROUTINE d3_pat2cart(d3in, nat, u1, u2, u3, d3out)
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
SUBROUTINE d3_3idx_2_6idx(nat, d3, phi)
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
SUBROUTINE d3_6idx_2_3idx(nat, phi, d3)
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
SUBROUTINE d3_cart2crys(nat, at, bg, phi)
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
SUBROUTINE d3_crys2cart(nat, at, bg, phi)
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
END MODULE d3_basis
!-----------------------------------------------------------------------
