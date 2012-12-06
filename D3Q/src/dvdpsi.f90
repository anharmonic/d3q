!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dvdpsi_module
contains
!-----------------------------------------------------------------------
SUBROUTINE dvdpsi (nu_i, u_x, xq_dv, dvloc, npw_rgt, igk_rgt, npw_lft, &
                   igk_lft, vkb_rgt, vkb_lft, dpsi, dvpsi)
  !-----------------------------------------------------------------------
  !
  ! Receives in input the variation of the local part of the KS-potential
  ! and calculates dV(xq_dv)_KS*dpsi in G_space, for all bands
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE cell_base,  ONLY : tpiba
  USE fft_base,   ONLY : dfftp, dffts
  USE fft_interfaces, ONLY: invfft, fwfft
  USE gvect,      ONLY : g
  USE gvecs,      ONLY : nls
  USE wvfct,      ONLY : nbnd, npwx !, igk !, npw
  USE uspp,       ONLY : nkb, dvan
  USE uspp_param, ONLY : nh
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  INTEGER,INTENT(in)     :: nu_i                   !the mode under consideration
  COMPLEX(DP),INTENT(in) :: u_x(3*nat, 3*nat)      !the displacement patterns of the modes
  REAL(DP),INTENT(in)    :: xq_dv(3)               ! coordinates of the q point describing the perturbation
  INTEGER,INTENT(in)     :: npw_rgt, npw_lft             ! number of planewaves in vkb_rgt and vkb_lft respectively
  INTEGER,INTENT(in)     :: igk_rgt(npwx), igk_lft(npwx) ! wavefunctions mapping for vkb_rgt and vkb_lft respectively
  COMPLEX(DP),INTENT(in) :: vkb_rgt(npwx,nkb), vkb_lft(npwx,nkb) !non-local part of the potential at k and k+xq_dv
  COMPLEX(DP),INTENT(in) :: dvloc (dfftp%nnr)           ! local part of the KS potential
  COMPLEX(DP),INTENT(in) :: dpsi (npwx, nbnd)      ! the d-wavefunction at k
  COMPLEX(DP),INTENT(out):: dvpsi (npwx, nbnd)     ! variation of the KS potential applied to dpsi
  !
  ! Local variables
  !
  INTEGER :: na, mu, ig, ir, ibnd, nt, ih, jh, ijkb0, ikb, jkb
  ! the transformation modes patterns
  COMPLEX(DP),ALLOCATABLE :: aux (:), ps (:,:), wrk (:)
  ! work space
  COMPLEX(DP),EXTERNAL:: zdotc
  !
  ALLOCATE(aux(dfftp%nnr))
  ALLOCATE(ps(2, nbnd))
  ALLOCATE(wrk(npwx))
  !
  ! Apply the local potential in real space, then come back to g space
  !
  DO ibnd = 1, nbnd
     aux = (0._dp, 0._dp)
     ! order correctly the plane waves in dpsi
     FORALL(ig=1:npw_rgt) aux(nls(igk_rgt(ig))) = dpsi(ig, ibnd)
     ! take dpsi to real space
     CALL invfft('Wave', aux, dffts)
     ! compute vloc*dpsi in real-space (this changes the periodicity)
     FORALL(ir=1:dffts%nnr) aux(ir) = aux (ir) * dvloc (ir)
     ! take vloc*dpsi back to g-space
     CALL fwfft('Wave', aux, dffts)
     ! reset the order of plane waves according to the new periodicity
     FORALL(ig=1:npw_lft) dvpsi(ig, ibnd) = aux(nls(igk_lft(ig)))
  ENDDO
  !
  !    Now the contribution of the non local part in the KB form
  !
  ijkb0=0
  DO nt = 1, ntyp
     DO na = 1, nat
        CORRECT_TYPE : &
        IF (ityp(na)==nt) THEN
           mu = 3 * (na - 1)
           !
           IH_LOOP : DO ih = 1, nh(nt)
              ikb = ijkb0 + ih
              !
              JH_LOOP : DO jh = 1, nh(nt)
               jkb = ijkb0 + jh
               !
               IF ( ABS (u_x (mu + 1, nu_i) ) + ABS (u_x (mu + 2, nu_i) ) + &
                    ABS (u_x (mu + 3, nu_i) ) > 1.0d-12 ) THEN
                  !
                  ! first term: sum_l v_l beta_l(k+q+G) \sum_G' beta^*_l(k+G') (iG'*u) psi
                  !
                  DO ig = 1, npw_rgt
                     wrk (ig) = vkb_rgt(ig,jkb) * &
                           CONJG( (0._dp,1._dp)*tpiba * &
                                 (g(1, igk_rgt(ig) ) * u_x(mu+1, nu_i) + &
                                  g(2, igk_rgt(ig) ) * u_x(mu+2, nu_i) + &
                                  g(3, igk_rgt(ig) ) * u_x(mu+3, nu_i) ) )
                  ENDDO
                  DO ibnd = 1, nbnd
                     ps(1,ibnd) = dvan(ih,jh,nt) * ZDOTC(npw_rgt, wrk,           1, dpsi(1,ibnd), 1)
                     ps(2,ibnd) = dvan(ih,jh,nt) * ZDOTC(npw_rgt, vkb_rgt(1,jkb),1, dpsi(1,ibnd), 1)
                  ENDDO
                  !
                  ! when build is serial this call does nothing, we leave it there
                  !
                  CALL mp_sum ( ps, intra_pool_comm )
                  !
                  ! second term: sum_l E_l(-i(q+G)*u) beta_l(k+q+G)\sum_G'beta^*_l(k+G')ps
                  !
                  DO ig = 1, npw_lft
                     wrk (ig) = vkb_lft(ig,ikb) * (0._dp,-1._dp)*tpiba * &
                           ( (g(1, igk_lft(ig) ) + xq_dv(1) ) * u_x (mu+1, nu_i) +&
                           (g(2, igk_lft(ig) ) + xq_dv(2) ) * u_x (mu+2, nu_i) +&
                           (g(3, igk_lft(ig) ) + xq_dv(3) ) * u_x (mu+3, nu_i) )
                  ENDDO
                  !
                  DO ibnd = 1, nbnd
                     CALL ZAXPY(npw_lft,ps(1,ibnd),vkb_lft(1,ikb),1,dvpsi(1,ibnd),1)
                     CALL ZAXPY(npw_lft,ps(2,ibnd),wrk,           1,dvpsi(1,ibnd),1)
                  ENDDO
               ENDIF
               !
              ENDDO JH_LOOP
           ENDDO  IH_LOOP
           !
           ijkb0 = ijkb0 + nh(nt)
           !
        END IF CORRECT_TYPE
     END DO
  END DO
  DEALLOCATE (wrk)
  DEALLOCATE (ps)
  DEALLOCATE (aux)
  RETURN
END SUBROUTINE dvdpsi

END MODULE dvdpsi_module
