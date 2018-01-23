!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE d3_nlcc_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_nlcc_0(d3dyn)
  !-----------------------------------------------------------------------
  !
  ! It calculates contribution due to non-linear-core-correction
  ! The variation of the density with respect to the perturbation must
  ! be corrected before calling this routine:
  ! while reading the variation of the density on unit iuaux and iud0rho
  ! it assumes it is the total density, i.e. sum of valence + core.
  !
  USE ions_base,  ONLY : nat, ityp, tau
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE gvect,      ONLY : ngm, g !, nl
  USE cell_base,  ONLY : tpiba2, tpiba, omega
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces, ONLY: fwfft
  USE uspp,       ONLY : nlcc_any
  USE d3com,      ONLY : d3c !, npert_i, npert_f
  !
  USE scf,        ONLY : rho, rho_core
  USE d3_basis,   ONLY : patq
  USE mp,         ONLY : mp_sum
  USE mp_world,   ONLY : world_comm
  USE mp_pools,   ONLY : intra_pool_comm

  implicit none
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)

  INTEGER :: na, nta, ig=0, i_cart, j_cart, k_cart, na_i, na_j, &
       na_k, nu_i, nu_j, nu_k, na_icart, nb_jcart, nc_kcart

  COMPLEX (DP) ::  work, expf
  COMPLEX (DP),ALLOCATABLE :: drc_exp(:,:), d3dyn0(:,:,:)

  REAL(DP),ALLOCATABLE :: vxc(:)
  REAL(DP) :: vtxc, etxc, arg, pref
  COMPLEX(DP),ALLOCATABLE :: vxc_g(:)
  !
  IF (.not. nlcc_any) RETURN
  !
  CALL start_clock('d3_nlcc_0')
  !
  ALLOCATE(drc_exp(ngm, nat))
  ALLOCATE(d3dyn0(3*nat, 3*nat, 3*nat))

  d3dyn0   = (0._dp, 0._dp)
  ! calculate v_x in real space than take it to g-space
  ! FIXME: drc_exp is dummy variable, v_xc should be changed instead
  !        This trick will BADLY FAIL for meta (non-local) functional!!!
  drc_exp = 0._dp
  ALLOCATE(vxc(dfftp%nnr))
  !
  CALL v_xc( rho, rho_core, drc_exp, etxc, vtxc, vxc )
  !
  ALLOCATE(vxc_g(dfftp%nnr))
  !
  FORALL(ig=1:dfftp%nnr) vxc_g(ig) = CMPLX(vxc(ig),0._dp,kind=DP)
  DEALLOCATE(vxc)
  !
  CALL fwfft('Rho', vxc_g, dfftp)
  !
  drc_exp  = (0._dp, 0._dp)
  DO na = 1, nat
     nta = ityp (na)
     DO ig = 1, ngm
        arg = - tpi * SUM(g(:,ig)*tau(:,na))
!         arg = - tpi * (g(1,ig)*tau(1,na)+g(2,ig)*tau(2,na)+g(3,ig)*tau(3,na))  !SUM(g(:,ig)*tau(:,na))
        expf = CMPLX(COS(arg), SIN(arg),kind=DP)
        drc_exp(ig, na) = expf * d3c(0)%drc(ig, nta)
     ENDDO
  ENDDO
  !
  pref = omega * tpiba2 * tpiba
  !
  DO na_i = 1, 3 * nat !npert_i, npert_f
    na = (na_i - 1) / 3 + 1
    i_cart = na_i - 3 * (na - 1)
    DO j_cart = 1, 3
      na_j = j_cart + 3 * (na - 1)
      DO k_cart = 1, 3
        na_k = k_cart + 3 * (na - 1)
        !
        work = (0._dp, 0._dp)
        DO ig = 1, ngm
          work = work + (0._dp, 1._dp) * g(i_cart,ig)*g(j_cart,ig)*g(k_cart,ig) &
                        * CONJG(vxc_g(dfftp%nl(ig))) * drc_exp (ig, na)
        ENDDO
        !
        d3dyn0 (na_i, na_j, na_k) = pref*work
        !
      ENDDO
    ENDDO
  ENDDO
  !
  CALL mp_sum(d3dyn0, intra_pool_comm)
  !
  DEALLOCATE(vxc_g, drc_exp)
  !
  ! rotate to the basis of modes
  DO nu_k = 1, 3 * nat
    DO nu_i = 1, 3 * nat
      DO nu_j = 1, 3 * nat
        work = (0._dp, 0._dp)
        DO nc_kcart = 1, 3 * nat
        DO na_icart = 1, 3 * nat
        DO nb_jcart = 1, 3 * nat
            work = work + patq(1)%u(nc_kcart, nu_k) &
                         *patq(2)%u(na_icart, nu_i) &
                         *patq(3)%u(nb_jcart, nu_j) &
                         *d3dyn0(nc_kcart, na_icart, nb_jcart)
        ENDDO
        ENDDO
        ENDDO
        !
        d3dyn(nu_k, nu_i, nu_j) = d3dyn(nu_k, nu_i, nu_j) + work
        !
      ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE(d3dyn0)
  !
  CALL stop_clock('d3_nlcc_0')
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_nlcc_0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_nlcc_123(iq_drho, iq_cci, iq_ccj, d3dyn)
  !-----------------------------------------------------------------------
  !
  ! It calculates contribution due to non-linear-core-correction
  ! The variation of the density with respect to the perturbation must
  ! be corrected before calling this routine:
  ! while reading the variation of the density on unit iuaux and iud0rho
  ! it assumes it is the total density, i.e. sum of valence + core.
  !
  USE ions_base,  ONLY : nat, ityp, tau
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE gvect,      ONLY : ngm, g !, nl
  USE cell_base,  ONLY : tpiba2, omega
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces, ONLY: fwfft
  USE uspp,       ONLY : nlcc_any
  USE eqv,        ONLY : dmuxc
  USE d3com,      ONLY : d3c, npert_i, npert_f
  !
  USE mp_pools,   ONLY : inter_pool_comm, intra_pool_comm, my_pool_id, npool
  USE mp,         ONLY : mp_sum, mp_bcast
  USE io_global,  ONLY : stdout, ionode_id
  USE d3_iofiles, ONLY : read_drho !iu_drho_cc_q, davcio_drho2
  USE d3_basis,   ONLY : patq
  USE kplus3q,    ONLY : kplusq, q_names

  IMPLICIT NONE
  COMPLEX(DP),VOLATILE,INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  INTEGER,INTENT(in) :: iq_cci, iq_ccj, iq_drho

  INTEGER :: na, nta, ig, ir, i_cart, j_cart, &
             na_icart, na_jcart

  REAL(DP) :: arg  ! argument of the phase factor
  REAL(DP) :: pref ! integration prefactor

  COMPLEX (DP) ::  expf, work
  COMPLEX (DP), ALLOCATABLE :: drc_exp (:,:), drho_dmu(:), d3dyn123 (:,:)

  INTEGER,VOLATILE,POINTER :: nu_cci=>null(), nu_ccj=>null(), nu_drho=>null()
  INTEGER,VOLATILE,TARGET  :: nu(3) = (/ 0, 0, 0 /)
  CHARACTER(len=11),PARAMETER :: sub='d3_nlcc_123'

  IF (.not.nlcc_any) RETURN
  !
  CALL start_clock('d3_nlcc_123')
  !
  pref = omega * tpiba2
  !
  WRITE(stdout, '(7x,"nlcc for drho d2V @ ",a2,1x,a2,1x,a2)') &
      ADJUSTL(q_names(iq_drho)),ADJUSTL(q_names(iq_cci)), ADJUSTL(q_names(iq_ccj))
  !
  nu_drho => nu(ABS(iq_drho));
  nu_cci  => nu(ABS(iq_cci));
  nu_ccj  => nu(ABS(iq_ccj))
  ! check for consistency:
  nu = (/ -1,0,1 /)
  IF(nu_drho == nu_cci .or. nu_drho == nu_ccj .or. nu_cci == nu_ccj) &
    CALL errore(sub, "Invalid choice of iq's (repeated)", 2)
  !
  ALLOCATE(d3dyn123(3*nat, 3*nat))
  !
  d3dyn123 = (0._dp, 0._dp)
  !
  FIRST_POOL_ONLY : &
  IF ( my_pool_id == 0 ) THEN
    !
    ALLOCATE(drho_dmu(dfftp%nnr))
    ALLOCATE(drc_exp(ngm, nat))
    drc_exp  = (0._dp, 0._dp)
    pref = 0.5_dp * omega * tpiba2
    !
    DO na = 1, nat
      nta = ityp (na)
      DO ig = 1, ngm
          arg  = - tpi * SUM( (g(:, ig)+kplusq(iq_drho)%xq(:))*tau(:,na) )
          expf = CMPLX(COS(arg), SIN(arg), kind=DP)
          drc_exp (ig, na) = expf*d3c(iq_drho)%drc(ig, nta)
      ENDDO
    ENDDO
    !
    RHO_PERTURBATION : &
    DO nu_drho = npert_i, npert_f !1, 3*nat
      !
      CALL read_drho(drho_dmu, iq_drho, nu_drho, with_core=.true., pool_only=.true.)
      !
      FORALL(ir=1:dfftp%nnr) drho_dmu(ir) = drho_dmu(ir) * dmuxc(ir, 1, 1)
      !
      CALL fwfft('Rho', drho_dmu, dfftp)
      !
      DO na = 1, nat
        DO i_cart = 1, 3
          na_icart = i_cart + 3 * (na - 1)
          !
          DO j_cart = 1, 3
            na_jcart = j_cart + 3 * (na - 1)
            !
            work = (0._dp, 0._dp)
            DO ig = 1, ngm
                work = work - CONJG(drho_dmu(dfftp%nl(ig))) * drc_exp(ig, na) &
                            *(g(i_cart, ig) + kplusq(iq_drho)%xq(i_cart)) &
                            *(g(j_cart, ig) + kplusq(iq_drho)%xq(j_cart))
            ENDDO
            !
            d3dyn123(na_icart, na_jcart) = pref * CONJG(work)
            !
          ENDDO
        ENDDO
      ENDDO
      !
#ifdef __MPI
      CALL mp_sum ( d3dyn123, intra_pool_comm )
      !CALL mp_sum ( d3dyn123, inter_pool_comm )
#endif
      ! The dynamical matrix was computed in cartesian axis and now we put
      ! it on the basis of the modes
      !
      DO nu_cci = 1, 3*nat
      DO nu_ccj = 1, 3*nat
        !
        work = (0._dp, 0._dp)
        !
        DO na = 1, nat
          DO i_cart = 1,3
            na_icart = 3*(na-1) + i_cart
            !
            DO j_cart = 1,3
              na_jcart = 3*(na-1) + j_cart
              !
              work = work &
                  +  patq(iq_cci)%u(na_icart,nu_cci) &
                   * patq(iq_ccj)%u(na_jcart,nu_ccj) &
                   * d3dyn123(na_icart,na_jcart)
            ENDDO
          ENDDO
        ENDDO
        !
        d3dyn(nu(1), nu(2), nu(3)) = d3dyn(nu(1), nu(2), nu(3)) &
                                    + work
      ENDDO
      ENDDO
        !
    ENDDO &
    RHO_PERTURBATION
    !
    DEALLOCATE(drho_dmu)
    DEALLOCATE(drc_exp)
    !
  ENDIF FIRST_POOL_ONLY
  !
  IF ( npool > 1 ) CALL mp_bcast( d3dyn, ionode_id, inter_pool_comm )
  !
  DEALLOCATE(d3dyn123)
  !
  CALL stop_clock('d3_nlcc_123')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_nlcc_123
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_nlcc_module
!-----------------------------------------------------------------------
