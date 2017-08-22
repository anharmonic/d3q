!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE d3_init_module
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_init
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : eps8, degspin, pi
  USE io_global,     ONLY : stdout
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_world,      ONLY : world_comm
  USE mp,            ONLY : mp_max, mp_min
  USE ions_base,     ONLY : ntyp => nsp
  USE gvect,         ONLY : ngm, g 
  USE cell_base,     ONLY : tpiba2, omega
  USE uspp_param,    ONLY : upf
  USE atom,          ONLY : msh, rgrid
  USE uspp,          ONLY : nlcc_any
  USE d3com,         ONLY : d3v, d3c
  USE mp,            ONLY : mp_barrier
  USE kplus3q,       ONLY : write_igkq_d3, nksq, kplusq, tot_nks, nbnd_max
  USE pwcom,         ONLY : lsda, ef
  USE control_lr,    ONLY : alpha_pv, nbnd_occ
  USE klist,         ONLY : lgauss, nelec
  USE wvfct,         ONLY : et, nbnd
  USE d3_debug,      ONLY : dbg_full_bands

  IMPLICIT NONE
  ! local vars:
  INTEGER :: ik, nt
  INTEGER :: iq, ibnd, ipol
  !
  CHARACTER(len=7),PARAMETER :: sub = 'd3_init'
  !
  ! ==============================================================================================
  ! Compute the number of occupated bands for each k point
  !
  CALL start_clock('d3_init')
  !
  IF(ALLOCATED(nbnd_occ))  DEALLOCATE(nbnd_occ)
  !
  CALL setup_nbnd_occ()
  !
  IF(dbg_full_bands)THEN
    WRITE(stdout,*) "WARNING!",&
                    "Solving the Sternheimer equation with all empty bands."
    nbnd_occ(:) = nbnd
  ENDIF
  !
  IF(lgauss)THEN
    nbnd_max = nbnd
  ELSE
    nbnd_max = NINT(nelec/degspin)
  ENDIF
  !
  CALL setup_alpha_pv()
  !
  WRITE(stdout,'(5x,"alpha_pv:",1f12.4)') alpha_pv
  !
  ! ===========================================================================================
  !
  ! NOTE: allocation has to be done in advance, or in the following loop pointers
  !       could point to yet-unallocated arrays
  ALLOCATE(d3v(-3:3))
  ALLOCATE(d3c(-3:3))
  !
  DO iq = -3,3
    NULLIFY(d3v(iq)%loc)
    NULLIFY(d3c(iq)%drc)
    IF (kplusq(iq)%lstored) THEN
      ALLOCATE(d3v(iq)%loc(ngm, ntyp))
      d3v(iq)%loc = 0._dp
      IF ( nlcc_any ) THEN
        ALLOCATE(d3c(iq)%drc(ngm, ntyp))
        d3c(iq)%drc = 0._dp
      ENDIF
    ENDIF
  ENDDO
  !
  DO iq = -3,3
    IF (kplusq(iq)%lstored) THEN
      !
      IF ( nlcc_any ) CALL set_drhoc(kplusq(iq)%xq, d3c(iq)%drc)
      !
      DO nt = 1, ntyp
          !
          IF (upf(nt)%tcoulombp) then
            CALL setlocq_coul( kplusq(iq)%xq, upf(nt)%zp, tpiba2, ngm, g, omega,d3v(iq)%loc(:,nt) )
          ELSE
            CALL setlocq( kplusq(iq)%xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
                      upf(nt)%vloc, upf(nt)%zp, tpiba2, ngm, g, omega, &
                      d3v(iq)%loc(:,nt) )
          END IF
          !
      END DO
      !
    ELSE
      !
      IF(.not. ASSOCIATED(d3v(kplusq(iq)%copy_of)%loc)) &
        CALL errore('d3_init', 'wrong order of init', 1)
      d3v(iq)%loc => d3v(kplusq(iq)%copy_of)%loc
      d3c(iq)%drc => d3c(kplusq(iq)%copy_of)%drc
      !
    ENDIF
  ENDDO
  !
  ! END: code imported from phq_init.
  !
  ! Writes to file the plane-wave ordering of each k (including k+q_x) point.
  !
  IF(nksq<=0) CALL errore('d3_init', 'nksq is zero', 1)
  !
  ! FIXME: The next subroutine compute the k+q+G lists for all q and k points,
  ! this has already been done somewhere during the nscf calculation
  ! hence doing it again here is in principle dangerous (different sort would
  ! cause wrong results)
  CALL write_igkq_d3()
  !
  CALL mp_barrier(world_comm)
  !
  CALL stop_clock('d3_init')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE d3_init
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_init_module
!-----------------------------------------------------------------------
