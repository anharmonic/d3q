!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE dpsi1dv2dpsi3_module
  !
#define __LOTTA_MEM
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE dpsi1dv2dpsi3 (iq_rgt,iq_dv,iq_lft,d3dyn) !, order)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE units_lr,   ONLY : lrdwf
  USE control_lr, ONLY : nbnd_occ
  !
  USE ions_base,  ONLY : nat
  USE fft_base,   ONLY : dfftp
  USE mp_pools,   ONLY : inter_pool_comm, intra_pool_comm
  USE io_global,  ONLY : stdout
  USE mp,         ONLY : mp_sum
  USE pwcom,      ONLY : lgauss, degauss, et, ef, ngauss, xk
  USE wvfct,      ONLY : npwx, nbnd
  USE uspp,       ONLY : nkb
  USE d3_basis,   ONLY : patq
  USE kplus3q,    ONLY : nksq, kplusq, q_sum_rule, q_names, nbnd_max
  USE d3_iofiles, ONLY : iu_dwfc
  USE uspp_init,  ONLY : init_us_2
  ! D3 subroutines called:
  USE dvdpsi_module
  USE dq_vscf_module

  IMPLICIT NONE
  ! the indexes of the q vectors that will be associated to <dpsi|, dH and |dpsi> respectively 
  INTEGER,INTENT(in) :: iq_lft, iq_dv, iq_rgt
  ! the 3rd order dynamical matrix (it is NOT initialized to zero)
  COMPLEX(DP),VOLATILE,INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  ! if present AND false, the elements of D3 will NOT be ordered (see later)
  !LOGICAL,OPTIONAL,INTENT(in) :: order
  !
  CHARACTER(len=13),PARAMETER :: sub = 'dpsi1dv2dpsi3'
  COMPLEX (DP), ALLOCATABLE :: dvloc(:,:)  ! derivative of local potential
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:)! |dH|dpsi>
  !
  COMPLEX (DP), ALLOCATABLE :: d3dyn_tmp (:,:,:) ! workspace
  ! non-local part of potential (beta projectors):
  COMPLEX(DP),POINTER :: vkb_rgt(:,:) => null(), vkb_lft(:,:) => null()
  ! description of the wavefunctions
#ifdef __LOTTA_MEM
  COMPLEX(DP),ALLOCATABLE,TARGET :: dpsi_lft_buff(:,:,:)
  COMPLEX(DP),ALLOCATABLE,TARGET :: dpsi_rgt_buff(:,:,:)
#endif
  COMPLEX (DP), POINTER :: dpsi_lft (:,:), dpsi_rgt (:,:)
  INTEGER,POINTER :: igk_lft(:) =>null(), igk_rgt(:) => null()
  INTEGER :: npw_lft, npw_rgt, nu_x
  ! band occupations with smearing (metals only)
  REAL(DP) :: degaussm1
  REAL(DP),ALLOCATABLE :: wga(:)
  ! perturbation indexes (see later)
  INTEGER,VOLATILE,TARGET  :: nu(3)
  INTEGER,VOLATILE,POINTER :: nu_r, nu_v, nu_l
  ! internal switches
  LOGICAL :: not_lsame!, order_internal
  ! Various counters, i/o auxiliary:
  INTEGER :: ik, ik_lft, ik_rgt, ik_gam, ibnd, nrec, ios
  COMPLEX (DP) :: wrk  ! workspace
  ! external functions:
  REAL(DP),EXTERNAL :: wgauss   ! gaussian weight
  COMPLEX(DP),EXTERNAL :: zdotc ! standard BLAS L1
  !
  CALL start_clock('dpsi1dv2dpsi3')
  !
  ! Initialization of wavefunction files, k+q grids, etc.
  not_lsame = .not.kplusq(-iq_lft)%lsame(iq_rgt)
  !
  ALLOCATE(igk_rgt(npwx))
  ALLOCATE(vkb_rgt(npwx, nkb))
  IF (not_lsame) THEN
    ALLOCATE(igk_lft(npwx))
    ALLOCATE(vkb_lft(npwx, nkb))
  ELSE
    igk_lft => igk_rgt
    vkb_lft => vkb_rgt
  ENDIF
  !
#ifdef __LOTTA_MEM
  ALLOCATE(dpsi_lft_buff(npwx,nbnd,3*nat))
  NULLIFY (dpsi_lft)
  ALLOCATE(dpsi_rgt_buff(npwx,nbnd,3*nat))
  NULLIFY (dpsi_rgt)
#else
  ALLOCATE(dpsi_lft(npwx, nbnd))
  ALLOCATE(dpsi_rgt(npwx, nbnd))
#endif
  ALLOCATE(dvpsi(npwx, nbnd))
  !
  ALLOCATE(dvloc(dfftp%nnr,3*nat))
  ALLOCATE(wga(nbnd))
  ALLOCATE(d3dyn_tmp(3*nat, 3*nat, 3*nat))
  d3dyn_tmp (:,:,:) = (0._dp, 0._dp)
  !
  !
  WRITE(stdout, '(7x,"<d^",a3," psi_k|d^",a3," V|d^",a3," psi_k>")') &
    q_names(-iq_lft), q_names(iq_dv), q_names(iq_rgt)
  IF(iq_lft /= -q_sum_rule(iq_rgt, iq_dv)) &
    CALL errore(sub, "Invalid choice of iq's (sum coditions not respected)", 1)
  !
  ! nu(1) will contain the index of the mode associated with the perturbation
  ! at +/- q1, this may be the lft or rgt wavefunctions or the potential; in any case
  ! the perturbationpattern associated to q1 MUST correspond to the first index of D3.
  ! Accordingly, the pattern associated with q2 goes to the second index, q3 with the third.
  !order_internal = .true.
  !IF(present(order)) &
  !  order_internal = order
  !IF(order_internal) THEN
    nu_l => nu(iq_lft)
    nu_v => nu(iq_dv)
    nu_r => nu(iq_rgt)
  !ELSE
  !  ! in this case, we do not use the correct ordering (useful for debugging)
  !  nu_l => nu(1)
  !  nu_v => nu(2)
  !  nu_r => nu(3)
  !ENDIF
  ! check for consistency:
  nu = (/ -1,0,1 /)
  IF(nu_l == nu_r .or. nu_v == nu_r .or. nu_v == nu_l) &
    CALL errore(sub, "Invalid choice of iq's (repeated)", 2)
  !
  ! pre-compute dvloc
  DO nu_v = 1,3*nat
    CALL dq_vscf(nu_v, dvloc(:,nu_v), kplusq(iq_dv)%xq, iq_dv, patq(iq_dv)%u)
  ENDDO
  !
  !REWIND (unit = kplusq(-iq_lft)%iunigkq)
  !IF(not_lsame) REWIND (unit = kplusq(iq_rgt)%iunigkq)
  !
  K_POINTS_LOOP : & !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  DO ik = 1, nksq
    !
    ik_rgt = kplusq( iq_rgt)%ikqs(ik)  ! global k-point index of k+q_rgt 
    ik_lft = kplusq(-iq_lft)%ikqs(ik) ! ... and of k+q_lft
    ik_gam = kplusq(0)%ikqs(ik)
    !
    npw_rgt = kplusq(iq_rgt)%ngkq(ik)
    igk_rgt = kplusq(iq_rgt)%igkq(:,ik)
    npw_lft = kplusq(-iq_lft)%ngkq(ik)
    ! igk_lft => igk_rgt : igk_lft is automatically a copy of igk_rgt
    IF(not_lsame) igk_lft = kplusq(-iq_lft)%igkq(:,ik)
    
    !READ (kplusq(iq_rgt)%iunigkq, iostat = ios) npw_rgt,igk_rgt
    !CALL errore (sub, 'reading iunigk (right)', ABS(ios) )
    !
    IF (not_lsame) THEN
      !READ (kplusq(-iq_lft)%iunigkq, iostat = ios) npw_lft, igk_lft
      !CALL errore (sub, 'reading iunigk (left)', ABS(ios) )
      !
      CALL init_us_2 (npw_rgt, igk_rgt, xk(1, ik_rgt), vkb_rgt)
      CALL init_us_2 (npw_lft, igk_lft, xk(1, ik_lft), vkb_lft)
    ELSE
      IF(ik_rgt /= ik_lft) CALL errore(sub, 'ik_lft /= ik_rgt', 1)
      !
      CALL init_us_2 (npw_rgt, igk_rgt, xk(1, ik_rgt), vkb_rgt)
      ! vkb_lft => vkb_rgt : vkb_lft is automatically a copy of vkb_rgt
    ENDIF
#ifdef __LOTTA_MEM
    DO nu_x = 1, 3 * nat
      ! pre-read "left" wavefunctions, this becomes a bottleneck for more than 3 atoms
      nrec = (nu_x - 1) * nksq + ik
      CALL davcio (dpsi_lft_buff(:,:,nu_x), lrdwf, iu_dwfc(0, -iq_lft), nrec, -1)
      ! pre-read "right" wavefunction, this is a bottleneck at 5 atoms
      CALL davcio (dpsi_rgt_buff(:,:,nu_x), lrdwf, iu_dwfc(0,  iq_rgt), nrec, -1)
    ENDDO
#endif
    !
    wga=0._dp
    IF (lgauss) THEN
      degaussm1 = 1._dp / degauss
      !
      DO ibnd = 1, nbnd_max !_occ(ik_gam)
        wga(ibnd) = kplusq(0)%wk(ik) * wgauss( (ef - et(ibnd, ik_gam))*degaussm1, ngauss )
      ENDDO
    ELSE
      !wga(1:nbnd_occ(ik_gam)) = kplusq(0)%wk(ik)
      wga(1:nbnd_max) = kplusq(0)%wk(ik)
    ENDIF
    !
    PERT_V_LOOP : & !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO nu_v = 1, 3 * nat
      ! compute the 1st order of the potential
!       CALL dq_vscf(nu_v, dvloc, kplusq(iq_dv)%xq, iq_dv, patq(iq_dv)%u)
      !
      PERT_R_LOOP : & !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO nu_r = 1, 3 * nat
        ! read the "right" 1st order wavefunctions (in dpsi_rgt)
#ifdef __LOTTA_MEM
        dpsi_rgt => dpsi_rgt_buff(:,:,nu_r)
#else
        nrec = (nu_r - 1) * nksq + ik
        CALL davcio (dpsi_rgt, lrdwf, iu_dwfc(0, iq_rgt), nrec, -1)
#endif
        !
        ! compute dH|dpsi> (non local part of dH included)
        CALL dvdpsi (nu_v, patq(iq_dv)%u, kplusq(iq_dv)%xq, dvloc(:,nu_v), &
                     npw_rgt,igk_rgt, npw_lft,igk_lft, vkb_rgt,vkb_lft, dpsi_rgt, dvpsi)
        !
        PERT_L_LOOP : & !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        DO nu_l = 1, 3 * nat
          !
          ! read the "left" 1st order wavefunctions (in dpsi_lft)
#ifdef __LOTTA_MEM
          dpsi_lft => dpsi_lft_buff(:,:,nu_l)
#else
          IF(not_lsame .or. nu_r /= nu_l) THEN
            nrec = (nu_l - 1) * nksq + ik
            CALL davcio (dpsi_lft, lrdwf, iu_dwfc(0, -iq_lft), nrec, -1)
          ELSE
            ! if we are doing the same mode at identical q's, we can just copy from memory
            dpsi_lft = dpsi_rgt
          ENDIF
#endif
          !
          ! Compute <dpsi|dH|dpsi>, including also the weight of the k-point and the smearing
          wrk = (0._dp, 0._dp)
          DO ibnd = 1, nbnd_max !_occ(ik_gam)
            wrk = wrk + wga(ibnd) * &
                  ZDOTC(npw_lft, dpsi_lft(:, ibnd), 1, dvpsi(:, ibnd), 1)
          ENDDO
#ifdef __MPI
          CALL mp_sum(wrk, intra_pool_comm)
#endif
          d3dyn_tmp (nu(1), nu(2), nu(3)) = d3dyn_tmp (nu(1), nu(2), nu(3)) + wrk
          !
        ENDDO & !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PERT_L_LOOP
        !
      ENDDO & !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      PERT_R_LOOP
      !
    ENDDO & !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    PERT_V_LOOP
    !
  ENDDO & !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  K_POINTS_LOOP
  !
#ifdef __MPI
  CALL mp_sum( d3dyn_tmp, inter_pool_comm )
#endif
  ! Add dpsidvdpsi contribution to the D3 matrix we got in input
  d3dyn = d3dyn + d3dyn_tmp
  !
  ! Clean-up
  DEALLOCATE(d3dyn_tmp)
#ifdef __LOTTA_MEM
  DEALLOCATE(dpsi_lft_buff, dpsi_rgt_buff)
  NULLIFY(dpsi_lft, dpsi_rgt)
#else
  DEALLOCATE(dpsi_lft, dpsi_rgt)
#endif
  DEALLOCATE(igk_rgt, vkb_rgt)
  IF(not_lsame) &
    DEALLOCATE (igk_lft, vkb_lft)
  DEALLOCATE(dvloc, dvpsi)
  DEALLOCATE(wga)
  !
  CALL stop_clock('dpsi1dv2dpsi3')
  !
  RETURN
END SUBROUTINE dpsi1dv2dpsi3
!
END MODULE dpsi1dv2dpsi3_module
