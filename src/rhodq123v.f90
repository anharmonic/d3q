!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE rhodq123v_module
CONTAINS
!----------------------------------------------------------------------
SUBROUTINE rhodq123v(d3dyn)
  !-----------------------------------------------------------------------
  !
  !  This routine calculates the electronic term: <psi|V"'|psi>
  !  of the third order dynamical matrix.
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : tpi
  USE ions_base,      ONLY : nat, ityp, ntyp => nsp, tau
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE uspp,           ONLY : dvan, nkb
  USE scf,            ONLY : rho
  USE gvect,          ONLY : g, ngm, igtongl !, nl
  USE wvfct,          ONLY : npwx, nbnd, wg
  USE vlocal,         ONLY : vloc
  USE klist,          ONLY : xk
  USE cell_base,      ONLY : omega, tpiba, tpiba2
  USE uspp_param,     ONLY : nh, nhm
  USE mp_pools,       ONLY : inter_pool_comm, intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE units_lr,       ONLY : iuwfc, lrwfc
  USE qpoint,         ONLY : nksq
  USE control_lr,     ONLY : nbnd_occ
  USE d3com,          ONLY : npert_i, npert_f
  USE kplus3q,        ONLY : kplusq
  USE d3_basis,       ONLY : patq
  USE d3_debug,       ONLY : dbgwrite_d3dyn
  USE uspp_init,      ONLY : init_us_2
  !
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)

  INTEGER :: icart, jcart, kcart, na_i, na_j, na_k, na, ng, ir, nt, &
             ik, ikk, ig, ibnd, ios, npw
  INTEGER :: ih, jh, ikb, ijkb0 ! indexes for projectors
  ! counters

  REAL(DP) :: gtau, fac, pref, tpiba3
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  COMPLEX(DP) :: zdotc, aux, aux2
  COMPLEX(DP),ALLOCATABLE :: d3dynwrk(:,:,:), d3dynwrk2(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: d3dynpat(:,:,:), d3dynpat2(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: rhog(:), work1 (:,:), work2 (:,:), work3 (:), alpha(:,:)
  COMPLEX(DP),ALLOCATABLE :: vkb(:,:), psi(:,:)
  INTEGER,ALLOCATABLE     :: igk(:)
  !
  CALL start_clock('rhodq123v')
  !
  ALLOCATE(rhog(dfftp%nnr))
  ALLOCATE(d3dynwrk (3*nat, 3*nat, 3*nat))
  ALLOCATE(d3dynwrk2(3*nat, 3*nat, 3*nat))
  !
  d3dynwrk = (0._dp, 0._dp)
  !
  FORALL(ir=1:dfftp%nnr) rhog(ir) = CMPLX(rho%of_r(ir, 1), 0._dp, kind=DP)
  !
  CALL fwfft('Rho', rhog, dfftp)
  !
  !     Contribution deriving from the local part of the potential
  !
  tpiba3 = tpiba2 * tpiba
  pref   = tpiba3 * omega
  !
  DO na_i = npert_i, npert_f
     na = (na_i - 1) / 3 + 1
     icart = na_i - 3 * (na - 1)
     DO jcart = 1, 3
        na_j = 3 * (na - 1) + jcart
        DO kcart = 1, 3
           na_k = 3 * (na - 1) + kcart
           DO ng = 1, ngm
              gtau = tpi * SUM(g(:, ng)*tau(:, na))
              !
              fac = pref * vloc (igtongl (ng), ityp (na)) &
                   * (  DBLE(rhog(dfftp%nl(ng))) * SIN(gtau)  &
                      +AIMAG(rhog(dfftp%nl(ng))) * COS(gtau) )
              !
              d3dynwrk (na_i, na_j, na_k) = d3dynwrk (na_i, na_j, na_k) &
                    + fac * g(icart, ng) * g(jcart, ng) * g(kcart, ng)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(rhog)
  !
  CALL mp_sum( d3dynwrk, intra_pool_comm )
  !
  !     Non local Kleinman-Bylander potential contribution
  !
  ALLOCATE(work1(npwx, 3))
  ALLOCATE(work2(npwx, 3))
  ALLOCATE(work3(npwx))
  ALLOCATE(vkb(npwx, nkb))
  ALLOCATE(igk(npwx))
  ALLOCATE(psi(npwx,nbnd))
  ALLOCATE(alpha(8,nhm))
  !
  d3dynwrk2 = (0._dp, 0._dp)
  !REWIND(unit = kplusq(0)%iunigkq)
  !
  DO ik = 1, nksq
     !READ (kplusq(0)%iunigkq, iostat = ios) npw, igk
     !CALL errore ('d3vrho', 'reading igk', abs (ios))
     ikk = kplusq(0)%ikqs(ik)
     npw = kplusq(0)%ngkq(ik)
     igk = kplusq(0)%igkq(:,ik)
     !
     CALL davcio (psi, lrwfc, iuwfc, ikk, -1)
     CALL init_us_2 (npw, igk, xk(1, ikk), vkb)

     DO kcart = 1, 3
        DO icart = 1, 3
           DO jcart = 1, 3
              DO ibnd = 1, nbnd_occ(ikk)
                 !
                 FORALL(ig=1:npw) work3(ig) = psi(ig, ibnd) * tpiba3 &
                                   * (xk(icart,ikk)+g(icart, igk(ig))) &
                                   * (xk(jcart,ikk)+g(jcart, igk(ig))) &
                                   * (xk(kcart,ikk)+g(kcart, igk(ig)))
                 !
                 FORALL(ig=1:npw) work2(ig, 1) = psi(ig, ibnd) * tpiba2 &
                                   * (xk(icart,ikk)+g(icart, igk(ig))) &
                                   * (xk(jcart,ikk)+g(jcart, igk(ig)))
                 !
                 FORALL(ig=1:npw) work2(ig, 2) = psi(ig, ibnd) * tpiba2 &
                                   * (xk(jcart,ikk)+g(jcart, igk(ig))) &
                                   * (xk(kcart,ikk)+g(kcart, igk(ig)))
                 !
                 FORALL(ig=1:npw) work2(ig, 3) = psi(ig, ibnd) * tpiba2 &
                                   * (xk(kcart,ikk)+g(kcart, igk(ig))) &
                                   * (xk(icart,ikk)+g(icart, igk(ig))) 
                 !
                 FORALL(ig=1:npw) work1(ig, 1) = psi(ig, ibnd) * tpiba &
                                   * (xk(kcart,ikk)+g(kcart, igk(ig)))
                 FORALL(ig=1:npw) work1(ig, 2) = psi(ig, ibnd) * tpiba &
                                   * (xk(icart,ikk)+g(icart, igk(ig)))
                 FORALL(ig=1:npw) work1(ig, 3) = psi(ig, ibnd) * tpiba &
                                   * (xk(jcart,ikk)+g(jcart, igk(ig)))
                 !
                 ijkb0 = 0
                 !
                 TYPES_LOOP : &
                 DO nt = 1, ntyp
                    ATOMS_LOOP : &
                    DO na = 1, nat
                       !
                       CORRECT_TYPE : &
                       IF (ityp (na) == nt) THEN
                          na_k = 3 * (na - 1) + kcart
                          na_i = 3 * (na - 1) + icart
                          na_j = 3 * (na - 1) + jcart
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             alpha(1, ih) = ZDOTC(npw, work3,       1, vkb(1,ikb),  1)
                             alpha(3, ih) = ZDOTC(npw, work1(1, 1), 1, vkb(1,ikb),  1)
                             alpha(5, ih) = ZDOTC(npw, work1(1, 2), 1, vkb(1,ikb),  1)
                             alpha(7, ih) = ZDOTC(npw, work1(1, 3), 1, vkb(1,ikb),  1)
                             !
                             alpha(2, ih) = ZDOTC(npw, vkb(1,ikb),  1, psi(1,ibnd), 1)
                             alpha(4, ih) = ZDOTC(npw, vkb(1,ikb),  1, work2(1, 1), 1)
                             alpha(6, ih) = ZDOTC(npw, vkb(1,ikb),  1, work2(1, 2), 1)
                             alpha(8, ih) = ZDOTC(npw, vkb(1,ikb),  1, work2(1, 3), 1)
                             
                          ENDDO
                          !
                          CALL mp_sum ( alpha, intra_pool_comm )
                          !
                          DO ih = 1, nh(nt)
                             DO jh = 1, nh(nt)
                               !
                               d3dynwrk2 (na_k, na_i, na_j) = d3dynwrk2 (na_k, na_i, na_j) - &
                                          2 * dvan(ih,jh,nt) * wg(ibnd, ikk) * &
                                    AIMAG(alpha(1,ih)*alpha(2,jh) + alpha(3,ih)*alpha(4,jh) +&
                                          alpha(5,ih)*alpha(6,jh) + alpha(7,ih)*alpha(8,jh))
                             ENDDO
                          ENDDO
                          !
                          ijkb0 = ijkb0 + nh(nt)
                          !
                       ENDIF CORRECT_TYPE
                    ENDDO ATOMS_LOOP
                 ENDDO TYPES_LOOP
                 !
              ENDDO !ibnd
           ENDDO !jcart
        ENDDO !ikart
     ENDDO !kcart 
  ENDDO !ik
  !
  DEALLOCATE(alpha)
  DEALLOCATE(work1, work2, work3)
  DEALLOCATE(psi)
  DEALLOCATE(vkb)
  DEALLOCATE(igk)
  !
  CALL mp_sum( d3dynwrk2, inter_pool_comm )
  !
  !   The dynamical matrix was computed in cartesian axis and now we put
  !   it on the basis of the modes
  !
  ALLOCATE( d3dynpat(3*nat, 3*nat, 3*nat))
  ALLOCATE(d3dynpat2(3*nat, 3*nat, 3*nat))
  d3dynpat  = (0._dp, 0._dp)
  d3dynpat2 = (0._dp, 0._dp)
  DO na_k = npert_i, npert_f
  DO na_i = 1, 3 * nat
  DO na_j = 1, 3 * nat
    aux  = (0._dp, 0._dp)
    aux2 = (0._dp, 0._dp)
    DO kcart = 1, 3 * nat
        DO icart = 1, 3 * nat
          DO jcart = 1, 3 * nat
              aux  = aux  + d3dynwrk(kcart, icart, jcart) &
                              * patq(1)%u(kcart, na_k) &
                              * patq(2)%u(icart, na_i) &
                              * patq(3)%u(jcart, na_j)
              aux2 = aux2 + d3dynwrk2(kcart, icart, jcart) &
                              * patq(1)%u(kcart, na_k) &
                              * patq(2)%u(icart, na_i) &
                              * patq(3)%u(jcart, na_j)
          ENDDO
        ENDDO
    ENDDO
    d3dynpat(na_k, na_i, na_j)  = aux
    d3dynpat2(na_k, na_i, na_j) = aux2
  ENDDO
  ENDDO
  ENDDO
  !
  CALL dbgwrite_d3dyn(d3dynpat,  'rd3v.1', 1)
  CALL dbgwrite_d3dyn(d3dynpat2, 'rd3v.2', 1)
  !
  d3dyn (:,:,:) = d3dyn (:,:,:) +  d3dynpat(:,:,:)+d3dynpat2(:,:,:)
  !
  DEALLOCATE(d3dynwrk2, d3dynwrk)
  DEALLOCATE(d3dynpat2, d3dynpat)
  !
  CALL stop_clock('rhodq123v')
  !
  RETURN
  !
  !----------------------------------------------------------------------
END SUBROUTINE rhodq123v
!-----------------------------------------------------------------------
!
END MODULE rhodq123v_module
