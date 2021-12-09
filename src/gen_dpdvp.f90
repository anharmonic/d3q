!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE dpsi1dv2psi_module
!-----------------------------------------------------------------------
!
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE gen_dpsi1dv2psi
  !-----------------------------------------------------------------------
  USE KINDS,      ONLY : DP
  USE d3_iofiles, ONLY : iu_dpsi_dH_psi, iudpdvp
  USE pwcom,      ONLY : lgauss
  USE kplus3q,    ONLY : kplusq, q_names, q_names2
  USE io_global,  ONLY : stdout
  USE d3_restart, ONLY : done_dwfc
  !USE davcio_debug
  !
  IMPLICIT NONE
  INTEGER :: iq_dp, iq_v, iq_dp_x, iq_v_x
  INTEGER :: i, unit
  LOGICAL :: done(-3:3,-3:3)
  INTEGER,PARAMETER :: ncases = 12
  INTEGER,PARAMETER :: cases(2,ncases) = &
        RESHAPE( (/ &
           1,-2,   1,-3, &
           2,-1,   2,-3, &
           3,-1,   3,-2, &
          -1, 2,  -1, 3, &
          -2, 1,  -2, 3, &
          -3, 1,  -3, 2 &
         /), (/2,ncases/) )

  IF (.not. lgauss) RETURN
  !
  CALL start_clock('gen_dpsi1dv2psi')
  !
  done = .false.
  DO i = 1, ncases
    !printed=.false.
    iq_dp = cases(1,i)
    iq_dp_x = kplusq(iq_dp)%copy_of
    !
    iq_v  = cases(2,i)
    iq_v_x = kplusq(iq_v)%copy_of
    !
    IF(.not.done(iq_dp_x,iq_v_x)) THEN
      IF(.not.done_dwfc(-iq_dp,iq_dp)) THEN
        WRITE(stdout,'(7x,"Skipping: < Pc d^",a," psi_k",a,"| d",a," V | psi_k",a," >")') &
          TRIM(q_names(-iq_dp)), TRIM(q_names2(iq_dp)), TRIM(q_names(iq_v)), TRIM(q_names2(-iq_v))
        WRITE(stdout,'(9x,"--> Wavefunction derivative not available (special case: q_x==Gamma)")')
      ELSE
        done(iq_dp_x,iq_v_x) = .true.
        done(iq_dp,iq_v) = .true.
        !
        unit = iu_dpsi_dH_psi(iq_dp, iq_v)
        CALL dpsi1dv2psi(unit, iq_dp, iq_v)
      ENDIF
    ELSE
      done(iq_dp,iq_v) = .true.
      WRITE(stdout,'(7x,"Skipping: < Pc d^",a," psi_k",a,"| d",a," V | psi_k",a," >")') &
        TRIM(q_names(-iq_dp)), TRIM(q_names2(iq_dp)), TRIM(q_names(iq_v)), TRIM(q_names2(-iq_v))
      WRITE(stdout,'(9x,"--> already done as < Pc d^",a," psi_k",a,"| d",a," V | psi_k",a," >")') &
        TRIM(q_names(-iq_dp_x)), TRIM(q_names2(iq_dp_x)), TRIM(q_names(iq_v_x)), TRIM(q_names2(-iq_v_x))
    ENDIF
    !
  ENDDO
  !
  !printed=.false.
  WRITE(stdout,'(7x,a)') 'Additional term to precompute for valence contribution:'
  ! IMPORTANT! The second parameter in the subroutine call defines which is "q" and
  ! which is "-q" of the  two non-zero q vectors this choice MUST be consistent
  ! with the iq_p in d3_valence_ij2!!!!
  IF(kplusq(1)%lgamma .and. kplusq(2)%lgamma .and. kplusq(3)%lgamma) THEN
    WRITE(stdout,'(9x,a)') 'none necessary for Gamma-only'
  ELSE IF(kplusq(1)%lgamma) THEN
    IF(done_dwfc(0,2)) CALL dpsi1dv2psi_gamma(iudpdvp(2), 2 )
    IF(done_dwfc(0,3)) CALL dpsi1dv2psi_gamma(iudpdvp(3), 3 )
  ELSE IF(kplusq(2)%lgamma) THEN
    IF(done_dwfc(0,1)) CALL dpsi1dv2psi_gamma(iudpdvp(1), 1 )
    IF(done_dwfc(0,3)) CALL dpsi1dv2psi_gamma(iudpdvp(3), 3 )
  ELSE IF(kplusq(3)%lgamma) THEN
    IF(done_dwfc(0,1)) CALL dpsi1dv2psi_gamma(iudpdvp(1), 1 )
    IF(done_dwfc(0,2)) CALL dpsi1dv2psi_gamma(iudpdvp(2), 2 )
  ENDIF
  !
  CALL stop_clock('gen_dpsi1dv2psi')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE gen_dpsi1dv2psi
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE dpsi1dv2psi(unit, iq_dp, iq_v)
  !-----------------------------------------------------------------------
  ! It comoputes the terms:
  !  < d^q_dp psi_(k-q_dp) | d^q_v V | psi_(k-q_v) >
  ! where:
  ! q_dp,q_v = 1,2,3,-1,-2,-3
  ! abs(q_v) /= abs(q_dp)
  ! The calculation is repeated for each perturbation and each band. 
  ! The resulting matrices are saved to unit in record
  !
  !
  USE kinds,       ONLY : DP
  USE klist,       ONLY : xk 
  USE units_lr,    ONLY : iuwfc, lrwfc, lrdwf
  USE ions_base,   ONLY : nat
  USE fft_base,    ONLY : dfftp
  USE wvfct,       ONLY : nbnd, npwx
  USE uspp,        ONLY : nkb
  USE mp_pools,    ONLY : intra_pool_comm
  USE mp,          ONLY : mp_sum
  USE kplus3q,     ONLY : kplusq, nksq, q_names, q_names2
  USE d3_basis,    ONLY : patq
  USE d3_iofiles,  ONLY : iu_dwfc, lrdpdvp
  USE io_global,   ONLY : stdout
  USE uspp_init,   ONLY : init_us_2
  ! D3 subroutines called:
  USE dvdpsi_module
  USE dq_vscf_module
  USE d3_restart,  ONLY : done_dwfc
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: unit
  INTEGER,INTENT(IN) :: iq_dp, iq_v
  !
  COMPLEX(DP),ALLOCATABLE :: dvloc(:,:), dpsidvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: vkb_gam(:,:), vkb_mqv(:,:) ! @ k, @ k-q_v
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), psi(:,:)
  INTEGER,ALLOCATABLE :: igk_gam(:), igk_mqv(:)
  INTEGER :: ik, ik_gam, ik_mqv
  INTEGER :: npw_gam, npw_mqv
  INTEGER :: nu_dp, nu_v
  INTEGER :: ibnd, jbnd
  INTEGER :: ios, nrec
  CHARACTER(len=11),PARAMETER:: sub = 'dpsi1dv2psi'
  !
  COMPLEX(DP),EXTERNAL :: ZDOTC
  !
  CALL start_clock('dpsi1dv2psi')
  !
  ALLOCATE(dvloc(dfftp%nnr,3*nat))
  ALLOCATE(dpsidvpsi(nbnd, nbnd))
  ALLOCATE(igk_gam(npwx), igk_mqv(npwx))
  ALLOCATE(vkb_gam(npwx, nkb))
  ALLOCATE(vkb_mqv(npwx, nkb))
  ALLOCATE(dvpsi(npwx, nbnd))
  ALLOCATE( dpsi(npwx, nbnd))
  ALLOCATE(  psi(npwx, nbnd))
  !
  IF(.not.done_dwfc(-iq_dp,iq_dp)) &
    CALL errore(sub, 'The necessary dwfc is not available!', 100+10*iq_dp+iq_v)
  !
  WRITE(stdout,'(7x,"< Pc d^",a," psi_k",a,"| d",a," V | psi_k",a," >")') &
    TRIM(q_names(-iq_dp)), TRIM(q_names2(iq_dp)), TRIM(q_names(iq_v)), TRIM(q_names2(-iq_v))
  !
  ! Pre-compute dvloc for every perturbation
  ! NOTE: computing it inside the k-point loop does not work with pools!
  DO nu_v = 1,3*nat
     CALL dq_vscf(nu_v, dvloc(:,nu_v), kplusq(iq_v)%xq, iq_v, patq(iq_v)%u)
  ENDDO
  !
!   REWIND (unit = kplusq( 0)%iunigkq)
!   REWIND (unit = kplusq(-iq_v)%iunigkq)
  !
  KPOINTS : &
  DO ik = 1, nksq
    !
    !READ(kplusq(0)%iunigkq, iostat = ios) npw_gam, igk_gam
!     IF(ios /=0) CALL errore(sub, 'Cannot read igk @ gamma', ABS(ios))
    ik_gam = kplusq( 0)%ikqs(ik)
    npw_gam = kplusq(0)%ngkq(ik)
    igk_gam = kplusq(0)%igkq(:,ik)
    !
    ik_mqv = kplusq(-iq_v)%ikqs(ik)
    npw_mqv = kplusq(-iq_v)%ngkq(ik)
    igk_mqv = kplusq(-iq_v)%igkq(:,ik)
!     IF(.not. kplusq(-iq_v)%lsame(0)) THEN
!       READ(kplusq(-iq_v)%iunigkq, iostat = ios) npw_mqv, igk_mqv
!       IF(ios /=0) CALL errore(sub, 'Cannot read igk @ -iq_v', ABS(ios))
!     ELSE
!       npw_mqv = npw_gam ; igk_mqv = igk_gam
!     ENDIF
    !
    !
    CALL init_us_2(npw_gam, igk_gam, xk(:, ik_gam), vkb_gam)
    CALL init_us_2(npw_mqv, igk_mqv, xk(:, ik_mqv), vkb_mqv)
    !
    CALL davcio(psi, lrwfc, iuwfc, ik_mqv, -1)
!     WRITE(20000*(mpime+1)+unit+100*ik,'(i5)') ik
!     WRITE(20000*(mpime+1)+unit+100*ik,'(20f12.6)') psi(1:10,1:3)
    !
    NU_V_PERT : &
    DO nu_v = 1, 3 * nat
      !
!       CALL dq_vscf(nu_v, dvloc, kplusq(iq_v)%xq, iq_v, patq(iq_v)%u)
      CALL dvdpsi(nu_v, patq(iq_v)%u, kplusq(iq_v)%xq, dvloc(:,nu_v), &
                  npw_mqv, igk_mqv, npw_gam, igk_gam, vkb_mqv, vkb_gam, psi, dvpsi)
      !
      DO nu_dp = 1, 3 * nat
        !
        nrec = (nu_dp-1)*nksq + ik
        CALL davcio(dpsi, lrdwf, iu_dwfc(-iq_dp,iq_dp), nrec, -1)
        !
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
              dpsidvpsi(ibnd, jbnd) = ZDOTC(npw_gam, dpsi(:,ibnd), 1, dvpsi(:,jbnd), 1)
          ENDDO
        ENDDO
#ifdef __MPI
        CALL mp_sum( dpsidvpsi, intra_pool_comm )
#endif
        nrec = nu_v + (nu_dp-1)*3*nat + (ik-1)*(3*nat)**2
        CALL davcio(dpsidvpsi, lrdpdvp, unit, nrec, +1)
!         write(10000+unit, '(3i5)') nu_v, nu_dp, ik
!         write(10000+unit, '(200f10.4)') dpsidvpsi
        !
      ENDDO
      !
     ENDDO &
    NU_V_PERT
    !
  ENDDO &
  KPOINTS
  !
  DEALLOCATE(dpsidvpsi)
  DEALLOCATE(igk_gam, igk_mqv)
  DEALLOCATE(vkb_gam, vkb_mqv)
  DEALLOCATE(dvloc, dvpsi, dpsi, psi)
  !
  CALL stop_clock('dpsi1dv2psi')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE dpsi1dv2psi
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE dpsi1dv2psi_gamma(unit, iq)
  !-----------------------------------------------------------------------
  ! It computes the terms:
  !  < d^-q psi_k | d^q V | psi_k >
  !
  USE kinds,       ONLY : DP
  USE klist,       ONLY : xk
  USE units_lr,    ONLY : iuwfc, lrwfc, lrdwf
  USE ions_base,   ONLY : nat
  USE fft_base,    ONLY : dfftp
  USE wvfct,       ONLY : nbnd, npwx
  USE uspp,        ONLY : nkb
  USE mp_pools,    ONLY : intra_pool_comm
  USE mp,          ONLY : mp_sum
  USE kplus3q,     ONLY : kplusq, nksq, q_names
  USE d3_basis,    ONLY : patq
  USE d3_iofiles,  ONLY : iu_dwfc, lrdpdvp
  USE io_global,   ONLY : stdout
  USE uspp_init,   ONLY : init_us_2
  ! D3 subroutines called:
  USE dvdpsi_module
  USE dq_vscf_module
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: unit
  INTEGER,INTENT(IN) :: iq
  !
  COMPLEX(DP),ALLOCATABLE :: dvloc(:,:), dpsidvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: vkb_g(:,:), vkb_q(:,:) ! @ k, @ k-q_v
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), psi(:,:)
  INTEGER,ALLOCATABLE :: igk_g(:), igk_q(:)
  INTEGER :: ik, ik_g, ik_q
  INTEGER :: npw_g, npw_q
  INTEGER :: nu_dp, nu_v
  INTEGER :: ibnd, jbnd
  INTEGER :: ios, nrec
  CHARACTER(len=17),PARAMETER:: sub = 'dpsi1dv2psi_gamma'
  !
  COMPLEX(DP),EXTERNAL :: ZDOTC
  !
  ALLOCATE(dvloc(dfftp%nnr,3*nat))
  ALLOCATE(dpsidvpsi(nbnd, nbnd))
  ALLOCATE(igk_g(npwx), igk_q(npwx))
  ALLOCATE(vkb_g(npwx, nkb), vkb_q(npwx, nkb))
  ALLOCATE(dvpsi(npwx, nbnd))
  ALLOCATE(dpsi(npwx, nbnd))
  ALLOCATE(psi(npwx, nbnd))
  !
  WRITE(stdout,'(7x,"< Pc dpsi_k/du(",a,")| dH/du(",a,") | psi_k >")') &
    TRIM(q_names(-iq)), TRIM(q_names(iq))
  !
  ! Pre-compute dvloc for every perturbation
  ! NOTE: computing it inside the k-point loop does not work with pools!
  DO nu_v = 1,3*nat
     CALL dq_vscf(nu_v, dvloc(:,nu_v), kplusq(iq)%xq, iq, patq(iq)%u)
  ENDDO
  !
!   REWIND (unit = kplusq( 0)%iunigkq)
!   REWIND (unit = kplusq(iq)%iunigkq)
  !
  KPOINTS : &
  DO ik = 1, nksq
    !
    ik_g  = kplusq(0)%ikqs(ik)
    npw_g = kplusq(0)%ngkq(ik)
    igk_g = kplusq(0)%igkq(:,ik)
!     READ(kplusq(0)%iunigkq, iostat = ios) npw_g, igk_g
!     IF(ios /=0) CALL errore(sub, 'Cannot read igk @ gamma', ABS(ios))
    !
    ik_q  = kplusq(iq)%ikqs(ik)
    npw_q = kplusq(iq)%ngkq(ik)
    igk_q = kplusq(iq)%igkq(:,ik)
!     IF(.not. kplusq(iq)%lgamma) THEN
!       READ(kplusq(iq)%iunigkq, iostat = ios) npw_q, igk_q
!       IF(ios /=0) CALL errore(sub, 'Cannot read igk @ iq', ABS(ios))
!     ELSE
!       npw_q = npw_g ; igk_q = igk_g
!     ENDIF
    !
    !
    CALL init_us_2 (npw_g, igk_g, xk(1, ik_g), vkb_g)
    CALL init_us_2 (npw_q, igk_q, xk(1, ik_q), vkb_q)
    !
    CALL davcio(psi, lrwfc, iuwfc, ik_g, -1)
    !
    NU_V_PERT : &
    DO nu_v = 1, 3 * nat
      !
!       CALL dq_vscf(nu_v, dvloc(:,nu_v), kplusq(iq)%xq, iq, patq(iq)%u)
      CALL dvdpsi(nu_v, patq(iq)%u, kplusq(iq)%xq, dvloc(:,nu_v), &
                  npw_g, igk_g, npw_q, igk_q, vkb_g, vkb_q, psi, dvpsi)
      !
      DO nu_dp = 1, 3 * nat
        !
        nrec = (nu_dp-1)*nksq + ik
        CALL davcio(dpsi, lrdwf, iu_dwfc(0,iq), nrec, -1)
        !
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
              dpsidvpsi(ibnd, jbnd) = &
                        ZDOTC(npw_q, dpsi(:,ibnd), 1, dvpsi(:,jbnd), 1)
          ENDDO
        ENDDO
#ifdef __MPI
        CALL mp_sum( dpsidvpsi, intra_pool_comm )
#endif
        nrec = nu_v + (nu_dp-1)*3*nat + (ik-1)*(3*nat)**2
        CALL davcio (dpsidvpsi, lrdpdvp, unit, nrec, + 1)
        !
      ENDDO
      !
     ENDDO &
    NU_V_PERT
    !
  ENDDO &
  KPOINTS
  !
  DEALLOCATE(dvloc)
  DEALLOCATE(dpsidvpsi)
  DEALLOCATE(igk_g, igk_q)
  DEALLOCATE(vkb_g)
  DEALLOCATE(vkb_q)
  DEALLOCATE(dvpsi)
  DEALLOCATE(dpsi)
  DEALLOCATE(psi)
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE dpsi1dv2psi_gamma
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE dpsi1dv2psi_module
!-----------------------------------------------------------------------
