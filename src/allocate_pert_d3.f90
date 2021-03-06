!
! Copyright (C) 2001-2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_pert_d3()
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities depending on the 
  ! maximum number of perturbations
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  !
  USE modes,      ONLY : npertx, t, tmq
  !
  IMPLICIT NONE
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  ALLOCATE (t(npertx, npertx, 48, 3*nat))
  ALLOCATE (tmq(npertx, npertx, 3*nat))
  !
  RETURN
END SUBROUTINE allocate_pert_d3
