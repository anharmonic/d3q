!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! If enough RAM is available, pre-read all the Right dpsi in order to
! only read them once
#define __LOTTA_MEM
MODULE dpsi1dpsi2dv3_module
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE dpsi1dpsi2dv3(iq_rgt,iq_dH,iq_lft, d3dyn, order)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE io_global,  ONLY : stdout
  USE mp_pools,   ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE pwcom,      ONLY : npwx, nbnd, degauss, ngauss, et, ef
  USE units_lr,   ONLY : lrdwf 
  USE d3com,      ONLY : eps_delta
  USE kplus3q,    ONLY : kplusq, nksq, q_names, q_names2, q_sum_rule, nbnd_max
  USE d3_iofiles, ONLY : iu_dwfc, iu_psi_dH_psi, iu_dpsi_dH_psi, lrdpdvp, lrpdqvp
  USE control_lr, ONLY : nbnd_occ

  IMPLICIT NONE
  INTEGER,INTENT(in) :: iq_lft,iq_rgt,iq_dH
  COMPLEX(DP),VOLATILE,INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  LOGICAL,INTENT(in),OPTIONAL :: order
  !
#ifdef __LOTTA_MEM
  COMPLEX(DP),ALLOCATABLE,TARGET :: dpsi_rgt_buff(:,:,:)
  COMPLEX(DP),POINTER            :: dpsi_rgt(:,:)
#else
  COMPLEX(DP),ALLOCATABLE :: dpsi_rgt(:,:)
#endif
  COMPLEX(DP),ALLOCATABLE :: dpsi_lft(:,:), psi_dH_psi(:,:)
  COMPLEX(DP),ALLOCATABLE :: dpsi_dpsi(:,:)
  COMPLEX(DP) :: wrk, wrk2
  COMPLEX(DP),EXTERNAL :: ZDOTC

  INTEGER :: ik, ibnd, jbnd, ios
  REAL(DP) :: deltae, degaussm1 = 0._dp
  REAL(DP),ALLOCATABLE :: wg_lft(:), wg_rgt(:), dwg(:)
  COMPLEX(DP),ALLOCATABLE :: ps_lr(:,:), ps_rl(:,:)
  COMPLEX (DP), ALLOCATABLE :: d3dyn_tmp (:,:,:) ! workspace
  !
  LOGICAL :: lmetal
  CHARACTER(len=13),PARAMETER :: sub = 'dpsi1dpsi2dv3'
  !
  !INTEGER,ALLOCATABLE :: igk_dummy(:)
  INTEGER :: ik_lft,ik_rgt, npw_gamma
  INTEGER :: nrec_lft,nrec_rgt,nrec_lr,nrec_rl,nrec_dH
  !
  ! smearing and smearing derivative
  REAL(DP), EXTERNAL :: wgauss, w0gauss
  !
  ! perturbation indexes (see later)
  INTEGER,VOLATILE,TARGET  :: nu(3)
  INTEGER,VOLATILE,POINTER :: nu_l, nu_r, nu_h
  LOGICAL :: order_internal
  !
  CALL start_clock('dpsi1dpsi2dv3')
  !
  lmetal = (degauss /= 0._dp)
  !
#ifdef __LOTTA_MEM
  ALLOCATE(dpsi_rgt_buff(npwx,nbnd,3*nat))
  NULLIFY (dpsi_rgt)
#else
  ALLOCATE( dpsi_rgt(npwx, nbnd) )
#endif
  ALLOCATE( dpsi_lft(npwx, nbnd) )
  ALLOCATE( psi_dH_psi(nbnd, nbnd) )
  ALLOCATE( dpsi_dpsi(nbnd, nbnd) )
  !ALLOCATE( igk_dummy( npwx) )
  !
  IF (lmetal) THEN
     degaussm1 = 1._dp / degauss
     ALLOCATE(wg_lft(nbnd), wg_rgt(nbnd), dwg(nbnd))
     ALLOCATE(ps_lr( nbnd, nbnd) )
     ALLOCATE(ps_rl( nbnd, nbnd) )
  ENDIF
  !
  ALLOCATE(d3dyn_tmp(3*nat, 3*nat, 3*nat))
  d3dyn_tmp (:,:,:) = (0._dp, 0._dp)
  !
  !
  ! We have to integrate on the Brillouin zone the following quantity:
  ! <d_(-q_lft) psi_(k+q_lft)| d_(q_rgt) psi_(k-q_rgt)> <psi_(k-q_rgt)|d_(q_dH) H|psi_(k+q_lft)>
  !                                               dH -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! The necessary dH term is pre-computed (for each k point) in solve_linter_d3q and
  ! saved in unit iu_psi_dH_psi(iu_lft,iq_dH).
  ! The wavefunction derivative |d_(q_rgt) psi_(k-q_rgt)> is stored in iu_dwfc(-iq_rgt,iq_rgt),
  ! while <d_(-q_lft) psi_(k+q_lft)| is the c.c. of what's in iu_dwfc(iq_lft,-iq_lft).
  !
  WRITE(stdout, '(7x,"<d_",a3," psi_k",a3,"| d_",a2," psi_k",a3,"> <psi_k",a3,"|d_",a2," H|psi_k",a3,">")') &
      q_names(-iq_lft),ADJUSTL(q_names2(iq_lft)), ADJUSTL(q_names(iq_rgt)),&
      q_names2(-iq_rgt), q_names2(-iq_rgt),ADJUSTL(q_names(iq_dH)),ADJUSTL(q_names2(iq_lft))

  IF(iq_lft /= -q_sum_rule(iq_rgt, iq_dH)) &
    CALL errore(sub, "Invalid choice of iq's (sum coditions not respected)", 1)
  !
  ! nu(1) will contain the index of the mode associated with the perturbation
  ! at +/- q1, this may be the lft or rgt wavefunctions or the potential; in any case
  ! the perturbationpattern associated to q1 MUST correspond to the first index of D3.
  ! Accordingly, the pattern associated with q2 goes to the second index, q3 with the third.
  order_internal = .true.
  IF(present(order)) &
    order_internal = order
  IF(order_internal) THEN
    nu_l => nu(ABS(iq_lft)); nu_r => nu(ABS(iq_rgt)); nu_H => nu(ABS(iq_dH))
  ELSE
    ! in this case, we do not use the correct ordering (useful for debugging)
    nu_l => nu(1); nu_r => nu(2); nu_H => nu(3)
  ENDIF
  ! check for consistency:
  nu = (/ -1,0,1 /)
  IF(nu_l == nu_r .or. nu_H == nu_r .or. nu_H == nu_l) &
    CALL errore(sub, "Invalid choice of iq's (repeated)", 2)
  ! Rewind the units of igk, they can actually be the same, but 2 rewinds do no harm:
!  REWIND(kplusq(0)%iunigkq)
  !
!  write(stdout,'(a,i4)') "r psidvpsi:", iu_psi_dH_psi(iq_lft,iq_dH)
  !
  KPOINTS_LOOP : &
  DO ik = 1, nksq
    ! Get the k-points corresponding to k+q_lft and k-q_rgt
    ik_lft = kplusq( iq_lft)%ikqs(ik)
    ik_rgt = kplusq(-iq_rgt)%ikqs(ik)
    !
    ! read the number of plane-waves, because the dwfc that we use here are in the
    ! form |d_q psi_k-q> we take the planewaves from k+q-q which is just k.
    ! we do not actuallyt need the igk (it has to be read anyway)
    !READ (kplusq(0)%iunigkq, iostat = ios) npw_gamma, igk_dummy
    npw_gamma = kplusq(0)%ngkq(ik)
    !
    ! In the metal case we need to pre-compute some weights
    IF (lmetal) THEN
      wg_lft = 0._dp
      wg_rgt = 0._dp
      dwg    = 0._dp
      DO ibnd = 1, nbnd_max !_occ(ik_lft)
        wg_lft(ibnd) = wgauss ((ef - et(ibnd, ik_lft))*degaussm1, ngauss)
        dwg(ibnd)    = w0gauss((ef - et(ibnd, ik_lft))*degaussm1, ngauss) * degaussm1
      ENDDO
      DO ibnd = 1, nbnd_max !_occ(ik_rgt)
        wg_rgt(ibnd) = wgauss ((ef - et(ibnd, ik_rgt))*degaussm1, ngauss)
      ENDDO
    ENDIF
    !
#ifdef __LOTTA_MEM
    DO nu_r = 1, 3*nat
      nrec_rgt = (nu_r - 1) * nksq + ik
      CALL davcio(dpsi_rgt_buff(:,:,nu_r), lrdwf, iu_dwfc(-iq_rgt,iq_rgt), nrec_rgt, -1)
    ENDDO
#endif
    !
    LEFT_DWFC_LOOP : &
    DO nu_l = 1, 3*nat
      ! read |d_(-q_lft) psi_(k+q_lft)> (for all bands) from file:
      nrec_lft = (nu_l - 1) * nksq + ik
      CALL davcio(dpsi_lft, lrdwf, iu_dwfc(iq_lft,-iq_lft), nrec_lft, -1)
      !
      RIGHT_DWFC_LOOP : &
      DO nu_r = 1, 3*nat
        ! read |d_(q_rgt) psi_(k-q_rgt)> (for all bands) from file
        ! or copy it from dpsi_lft if q's are equal and we're doing the same
        ! perturbation:
#ifdef __LOTTA_MEM
        dpsi_rgt => dpsi_rgt_buff(:,:,nu_r)
#else
        nrec_rgt = (nu_r - 1) * nksq + ik
        IF( kplusq(iq_lft)%lsame(-iq_rgt) .and. nu_l == nu_r ) THEN
          dpsi_rgt = dpsi_lft
        ELSE
          CALL davcio(dpsi_rgt, lrdwf, iu_dwfc(-iq_rgt,iq_rgt), nrec_rgt, -1)
        ENDIF
#endif
        !
        ! In the metal case read some more <dpsi|dV|psi> terms
        IF (lmetal) THEN
          nrec_lr = nu_r + (nu_l-1)*3*nat + (ik-1)*(3*nat)**2
          CALL davcio(ps_lr, lrdpdvp, iu_dpsi_dH_psi(-iq_lft, iq_rgt), nrec_lr, -1)
          !
          nrec_rl = nu_l + (nu_r-1)*3*nat + (ik-1)*(3*nat)**2
          CALL davcio(ps_rl, lrdpdvp, iu_dpsi_dH_psi(iq_rgt, -iq_lft), nrec_rl, -1)
        ENDIF
        !
        ! precompute <dpsi|dpsi> for every band at this k-point, left and right perturbation
        DO ibnd = 1, nbnd_max !_occ(ik_lft) !nbnd
        DO jbnd = 1, nbnd_max !_occ(ik_rgt) !nbnd
            dpsi_dpsi(jbnd,ibnd) &
              = ZDOTC(npw_gamma, dpsi_lft(:,ibnd), 1, dpsi_rgt(:,jbnd), 1)
        ENDDO
        ENDDO
        !
        PSI_DH_PSI_LOOP : &
        DO nu_h = 1, 3*nat
          ! read the <psi|dH|psi> from file:
          nrec_dH = nu_h + (ik - 1) * 3*nat
          CALL davcio (psi_dH_psi, lrpdqvp, iu_psi_dH_psi(iq_lft,iq_dH), nrec_dH, -1)
          !
          wrk   = (0._dp, 0._dp)  ! <-- collected outside pools
          wrk2  = (0._dp, 0._dp)  ! <-- collected inside AND outside pools
          DO ibnd = 1, nbnd_max !_occ(ik_lft) !nbnd
          DO jbnd = 1, nbnd_max !_occ(ik_rgt) !nbnd
            IF (lmetal) THEN
              deltae = et(ibnd, ik_lft) - et(jbnd, ik_rgt)
              IF (ABS(deltae) > eps_delta) THEN
                  wrk = wrk + psi_dH_psi(jbnd, ibnd)  * &
                              (wg_lft(ibnd) * ps_lr(ibnd, jbnd) - &
                               wg_rgt(jbnd) * CONJG(ps_rl (jbnd, ibnd)) ) / deltae
              ELSE
                 wrk = wrk - psi_dH_psi(jbnd, ibnd) * dwg(ibnd) * ps_lr(ibnd, jbnd)
                 !
                 wrk2 = wrk2 - psi_dH_psi(jbnd, ibnd) &
                               * wg_lft(ibnd) * dpsi_dpsi(jbnd,ibnd)
              ENDIF
            ELSE
              ! insulators:
              wrk2 = wrk2 - psi_dH_psi(jbnd, ibnd) * dpsi_dpsi(jbnd,ibnd)
            ENDIF
          ENDDO
          ENDDO
          !
          CALL mp_sum(wrk2, intra_pool_comm )
          !
          d3dyn_tmp (nu(1), nu(2), nu(3)) = d3dyn_tmp (nu(1), nu(2), nu(3)) &
                                           + (wrk+wrk2) * kplusq(iq_lft)%wk(ik)
          !
        ENDDO &
        PSI_DH_PSI_LOOP
        !
      ENDDO &
      RIGHT_DWFC_LOOP 
    ENDDO &
    LEFT_DWFC_LOOP
  ENDDO &
  KPOINTS_LOOP
  !
#ifdef __MPI
  CALL mp_sum( d3dyn_tmp, inter_pool_comm )
#endif
  ! Add dpsidvdpsi contribution to the D3 matrix we got in input
  d3dyn = d3dyn + d3dyn_tmp
  !
  !DEALLOCATE(igk_dummy)
#ifdef __LOTTA_MEM
  DEALLOCATE(dpsi_lft, dpsi_rgt_buff)
  NULLIFY(dpsi_rgt)
#else
  DEALLOCATE(dpsi_lft, dpsi_rgt)
#endif
  DEALLOCATE(dpsi_dpsi)
  IF (lmetal) THEN
     DEALLOCATE( ps_lr, ps_rl )
     DEALLOCATE(wg_lft, wg_rgt, dwg)
  ENDIF
  DEALLOCATE(d3dyn_tmp)
  !
  CALL stop_clock('dpsi1dpsi2dv3')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE dpsi1dpsi2dv3
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE dpsi1dpsi2dv3_module
!-----------------------------------------------------------------------
