!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE rotate_d3
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE rotate_and_add_d3(phi, phi2, nat, isym, s, invs, irt, rtau, sxq)
!-----------------------------------------------------------------------
!  Rotates a third order matrix (phi) in crystal coordinates according
!  to the specified symmetry operation and add the rotated matrix
!  to phi2.  phi is left unmodified.
!
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi
  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER,INTENT(in) :: nat, isym, s (3, 3, 48), invs (48), irt (48, nat)
  ! number of atoms in the unit cell
  ! index of the symm.op.
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom

  COMPLEX(DP),INTENT(in)   :: phi (3, 3, 3, nat, nat, nat)
  COMPLEX(DP),INTENT(inout):: phi2 (3, 3, 3, nat, nat, nat)
  ! the input d3dyn.mat.
  ! in crystal coordinates
  ! the rotated d3dyn.mat
  ! in crystal coordinates
  REAL(DP),INTENT(in) :: rtau (3, 48, nat), sxq(3,3)
  ! for each atom and rotation gives
  ! the R vector involved
  ! the rotated q involved in this sym.op
  !
  !  local variables
  !
  INTEGER :: na, nb, nc, sna, snb, snc, ism1, i, j, k, l, m, n, ipol
  ! counters on atoms
  ! indices of rotated atoms
  ! index of the inverse symm.op.
  ! generic counters
  REAL(DP) :: arg
  ! argument of the phase
  COMPLEX(DP) :: phase, work

!   write(*,*)  ' pippo'
!   write(*,*)  'isym:', isym
!   write(*,*)  'invs(isym):', invs(isym)
!   do na = 1, nat
!   write(*,*) ' na:', na, 'sna:', irt(isym,na)
!   write(*,*) ' rtau', rtau(:,isym,na)
!   enddo
!   print*,"sxq1:", sxq(:,1)
!   print*,"sxq2:", sxq(:,2)
!   print*,"sxq3:", sxq(:,3)
!   print*, "s:", s(:,:,invs(isym))

  ism1 = invs(isym)
  !
  DO nc = 1, nat
     snc = irt(isym,nc)
     DO na = 1, nat
        sna = irt(isym,na)
        DO nb = 1, nat
           snb = irt(isym,nb)
           arg = 0._dp
           DO ipol = 1,3
              arg =  arg - ( sxq(ipol,1) * rtau(ipol, isym, nc) + &
                             sxq(ipol,2) * rtau(ipol, isym, na) + &
                             sxq(ipol,3) * rtau(ipol, isym, nb) ) * tpi
           ENDDO
           !
           phase = CMPLX(cos(arg),-sin(arg), kind=DP)
           DO m = 1, 3
              DO i = 1, 3
                 DO j = 1, 3
                    work = (0._dp, 0._dp)
                    DO k = 1, 3
                       DO l = 1, 3
                          DO n = 1, 3
                             work = work &
                                  + s(m,n,ism1) * s(i,k,ism1) * s(j,l,ism1) &
                                  * phi(n,k,l,nc,na,nb) * phase
                          ENDDO
                       ENDDO
                    ENDDO
                    phi2(m,i,j,snc,sna,snb) = phi2(m,i,j,snc,sna,snb) + work
                 ENDDO
              ENDDO

           ENDDO
        ENDDO
     ENDDO

  ENDDO
  RETURN
END SUBROUTINE rotate_and_add_d3

END MODULE rotate_d3
