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
  USE mp_global,     ONLY : inter_pool_comm
  USE mp,            ONLY : mp_max, mp_min
  USE ions_base,     ONLY : ntyp => nsp
  USE pwcom,         ONLY : ngm, g, tpiba2, omega
  USE uspp_param,    ONLY : upf
  USE atom,          ONLY : msh, rgrid
  USE phcom,         ONLY : nlcc_any
  USE d3com,         ONLY : d3v, d3c
  USE mp,            ONLY : mp_barrier
  USE kplus3q,       ONLY : write_igkq_d3, nksq, kplusq, tot_nks

  USE pwcom,         ONLY : lsda, degauss, nelec, ngauss, &
                         ef, et, nbnd, xk
  USE phcom,         ONLY : nbnd_occ, alpha_pv

  IMPLICIT NONE
  ! local vars:
  INTEGER :: ik, nt
  INTEGER :: iq, ibnd, ipol
  REAL(DP) :: xmax, fac, TARGET, &  ! for nbnd_occ
              emin, emax            ! for alpha_pv
  ! parameters:
  REAL(DP),PARAMETER ::      small = 6.9626525973374d-5
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
  CHARACTER(len=7),PARAMETER :: sub = 'd3_init'
  !
  ! ==============================================================================================
  ! Compute the number of occupated bands for each k point
  !
  CALL start_clock('d3_init')
  !
  nbnd_occ = 0
  IF (degauss /= 0.d0) THEN
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     !
     ! - limit appropriated for gaussian broadening (used for all ngauss)
     !
     xmax = SQRT ( - LOG (SQRT (pi) * small) )
     !
     ! - limit appropriated for Fermi-Dirac
     !
     IF (ngauss == - 99) THEN
        fac = 1.d0 / SQRT (small)
        xmax = 2.d0 * LOG (0.5d0 * (fac + SQRT (fac * fac - 4.0d0) ) )
     ENDIF
     TARGET = ef + xmax * degauss
     DO ik = 1, tot_nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) < TARGET) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) == nbnd) &
             WRITE( stdout, '(5x,/,"Possibly too few bands at point ", &
             & i4,3f10.5)') ik,  (xk (ipol, ik) , ipol = 1, 3)
        IF (nbnd_occ (ik) >  nbnd) &
             WRITE( stdout, '(5x,/,"WARNING! Definitely too few bands at point ", &
             & i4,3f10.5)') ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
      WRITE(stdout,'(5x,a,i4)') "WARNING!! nbnd_occ == nbnd", nbnd
      nbnd_occ = nbnd
  ELSE
     IF (lsda) CALL infomsg (sub, 'occupation numbers probably wrong')
     nbnd_occ(1:tot_nks) = NINT(nelec/degspin)
  ENDIF
  !
  ! ==============================================================================================
  ! Computes alpha_pv
  !
  emin = et(1, 1)
  DO ik = 1, nksq
     DO ibnd = 1, nbnd_occ(ik)
        emin = MIN (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
  !
  emax = et(1, 1)
  DO ik = 1, nksq
     DO ibnd = 1, nbnd_occ(ik)
        emax = MAX (emax, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the maximum across pools
  CALL mp_max( emax, inter_pool_comm )
  !
  alpha_pv = 2._dp * (emax - emin)
  ! avoid zero value for alpha_pv
  alpha_pv = MAX (alpha_pv, 1.e-2_dp)
  WRITE(stdout,'(5x,"alpha_pv:",1f12.4)') alpha_pv
  !
  ! ==============================================================================================
  !
  ! NOTE: allocation has to be done in advance, or in the following loop pointers
  !       could point to yet-unallocated arrays
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
  DO ik = 1, nksq
    ! Note that this have already been done somewhere in the nscf calculation
    ! but here we want it on separate files for k, k+/-q1, q2 q3.
    ! Furthermore, I have no clue about *where* the nscf code (from run_nscf_d3)
    ! does it
    CALL write_igkq_d3(ik)
  ENDDO
  !
  CALL mp_barrier()
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