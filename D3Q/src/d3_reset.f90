!
! Copyright (C) 2011 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE d3_reset_module
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_reset(print_clock, cleanup)
  !-----------------------------------------------------------------------
  ! From d3_Setup:
  USE kplus3q,            ONLY : reset_kplus3q, kplusq
  USE uspp,               ONLY : nlcc_any
  USE lr_symm_base,       ONLY : rtau
  USE d3_basis,           ONLY : patq, deallocate_d3_pattern
  USE d3_symmetry,        ONLY : symq, sym_gamma, deallocate_d3_symmetry
  USE d3_iofiles,         ONLY : closefild3
  USE allocate_d3_module, ONLY : deallocate_d3
  USE d3com,              ONLY : d3v, d3c
  USE d3_shuffle,         ONLY : d3_reset_permutations
  !
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: print_clock, cleanup
  INTEGER :: iq
  !
  ! form d3toten:
  CALL d3_reset_permutations()
  ! from d3_setup: 
  CALL reset_kplus3q(cleanup) !<-- also resets q_special_cases
  DEALLOCATE(rtau)
  NULLIFY(sym_gamma)
  DO iq = 1,3
    CALL deallocate_d3_pattern(patq(iq))
    CALL deallocate_d3_pattern(patq(-iq))
    CALL deallocate_d3_symmetry(symq(iq))
    !
  ENDDO
  !
  DEALLOCATE(patq)
  !
  ! From d3_init:
!   DO iq = -3,3
!     IF (kplusq(iq)%lstored) THEN
!       DEALLOCATE(d3v(iq)%loc)
!       IF ( nlcc_any ) DEALLOCATE(d3c(iq)%drc)
!     ELSE
!       NULLIFY(d3v(iq)%loc)
!       NULLIFY(d3c(iq)%drc)
!     ENDIF
!   ENDDO
  ! FIXME: the following lines may cause a memory leak
  DO iq = -3,3
    IF (kplusq(iq)%lstored) THEN
      IF(associated(d3v(iq)%loc)) DEALLOCATE(d3v(iq)%loc)
      IF ( nlcc_any .and. associated(d3c(iq)%drc)) DEALLOCATE(d3c(iq)%drc)
    ENDIF
    NULLIFY(d3v(iq)%loc)
    NULLIFY(d3c(iq)%drc)
  ENDDO
  !
  DEALLOCATE(d3v)
  DEALLOCATE(d3c)
  !
  CALL closefild3(cleanup)
  !
  CALL deallocate_d3()
  !
  IF(print_clock) CALL print_clock_d3_short()
  !
  !-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_reset_module
!-----------------------------------------------------------------------

