!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE allocate_d3_module
CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE deallocate_d3()
  !-----------------------------------------------------------------------
  USE becmod,       ONLY : deallocate_bec_type, becp
  USE d3_efermi_shift, ONLY : ef_sh
  !
  IMPLICIT NONE
  !
  CALL deallocate_bec_type(becp)
  IF(allocated(ef_sh)) DEALLOCATE(ef_sh)
  !
  !-----------------------------------------------------------------------
END SUBROUTINE deallocate_d3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE allocate_d3()
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the third
  ! derivative of the total energy
  !
  ! Not much stuff left, as almost everything is now allocated on-place when used
  USE ions_base,    ONLY : nat
  USE uspp,         ONLY : nkb
  USE pwcom,        ONLY : nbnd, degauss, lgauss
  USE d3_efermi_shift, ONLY : ef_sh
  USE becmod, ONLY : allocate_bec_type, becp

  IMPLICIT NONE

  CALL allocate_bec_type( nkb, nbnd, becp )
  IF (lgauss) ALLOCATE (ef_sh( 3*nat))

  RETURN
END SUBROUTINE allocate_d3

END MODULE allocate_d3_module
