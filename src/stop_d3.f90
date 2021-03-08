!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE stop_d3_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE stop_d3 ()
!-----------------------------------------------------------------------
!
!    This routine closes all files before stopping
!    flag is no longer used
!
  USE eqv,              ONLY : dmuxc
  USE mp_global,        ONLY : mp_global_end
  USE d3_grid,          ONLY : d3_triplets

  IMPLICIT NONE
  !
  CALL clean_pw(.true.)
  !
  ! from d3_setup_q_independent
  DEALLOCATE(dmuxc)
  ! from d3_grid:
  DEALLOCATE(d3_triplets)

  CALL print_clock_d3()

  CALL mp_global_end()

  STOP
  RETURN
END SUBROUTINE stop_d3
END MODULE stop_d3_module
