!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE write_d3dyn_ascii
  CONTAINS
!-----------------------------------------------------------------------
subroutine write_d3dyn_XXX (xq1,xq2,xq3, phi, nat, iudyn)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell
  complex (DP) :: phi (3, 3, 3, nat, nat, nat)
  !  derivative of the dynamical matrix
  real (DP) :: xq1(3), xq2(3), xq3(3)
  ! the q vector
  !
  ! local variables
  !
  integer :: na, nb, nc, icar, jcar, kcar, i
  ! counters on atoms
  ! cartesian coordinate counters
  ! generic counter
  write (iudyn, "(/,5x,'Third derivative in cartesian axes',/)")
  write (iudyn, "(5x,'q1 = ( ',3f14.9,' ) ')") (xq1(icar), icar = 1, 3)
  write (iudyn, "(5x,'q2 = ( ',3f14.9,' ) ')") (xq2(icar), icar = 1, 3)
  write (iudyn, "(5x,'q3 = ( ',3f14.9,' ) ')") (xq3(icar), icar = 1, 3)
  do i = 1, 3 * nat
!      if (wrmode (i) ) then
        !
        write (iudyn, '(/,12x,"modo:",i5,/)') i
        nc = (i - 1) / 3 + 1
        kcar = i - 3 * (nc - 1)
        do na = 1, nat
           do nb = 1, nat
              write (iudyn, '(2i3)') na, nb
              do icar = 1, 3
                 write (iudyn, '(3(2f13.6,3x))') &
                 (phi (kcar, icar, jcar, nc, na, nb) , jcar = 1, 3)
              enddo
           enddo
        enddo
        !
!      endif
  enddo

  return

end subroutine write_d3dyn_XXX

!-----------------------------------------------------------------------
subroutine zero_d3dyn_XXX (phi, nat, thresh)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !              write(*,*) "dioporco"

  integer :: nat
  ! number of atom in the unit cell
  complex (DP) :: phi (3, 3, 3, nat, nat, nat)
  REAL(DP) :: thresh
  !
  ! local variables
  !
  integer :: na, nb, nc, icar, jcar, kcar
  ! counters on atoms
  ! cartesian coordinate counters
  ! generic counter
  do nc = 1, nat
  do nb = 1, nat
  do na = 1, nat
      do kcar = 1, 3
      do jcar = 1, 3
      do icar = 1, 3
            IF( ABS(phi (icar, jcar, kcar, na, nb, nc)) < thresh) THEN 
              phi (icar, jcar, kcar, na, nb, nc) = 0._dp
            ELSEIF( ABS(DBLE(phi (icar, jcar, kcar, na, nb, nc))) < thresh) THEN 
              phi (icar, jcar, kcar, na, nb, nc) = &
                (0._dp,1._dp) * DIMAG(phi (icar, jcar, kcar, na, nb, nc))
            ELSEIF( ABS(DIMAG(phi (icar, jcar, kcar, na, nb, nc))) < thresh) THEN 
              phi (icar, jcar, kcar, na, nb, nc) = &
                DBLE(phi (icar, jcar, kcar, na, nb, nc))
            ENDIF
      enddo
      enddo
      enddo
  enddo
  enddo
  enddo

  return

end subroutine zero_d3dyn_XXX

END MODULE 