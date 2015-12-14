!
! Copyright(C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE d3_exc_gc_module
  CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_exc_gc(iq_i, iq_j, iq_k, d3dyn)
  !-----------------------------------------------------------------------
  !
  !    wat20100930:
  !    Calculates the gradient corrections to the derivative of the
  !    dynamical matrix
  !
  USE ions_base,    ONLY : nat
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : omega, alat
  USE scf,          ONLY : rho, rho_core
  USE lsda_mod,     ONLY : nspin
  USE funct,        ONLY : dft_is_gradient
  USE fft_base,     ONLY : dfftp
  USE gvect,        ONLY : g, ngm, nl
  USE qpoint,       ONLY : xq
  USE gc_ph,        ONLY : grho, dvxc_sr, dvxc_ss
  USE gc_d3,        ONLY : dvxc_rrr, dvxc_srr, &
                           dvxc_ssr, dvxc_sss
  USE io_global,    ONLY : ionode_id
  USE mp_global,    ONLY : inter_pool_comm, my_pool_id, &
                           npool, intra_pool_comm
  USE mp,           ONLY : mp_bcast, mp_sum 
  USE d3_iofiles,   ONLY : read_drho
  USE kplus3q,      ONLY : kplusq
!  USE q0modes,   ONLY : q0mode33
  USE units_ph,     ONLY : lrdrho, iudrho
!  USE units_d3,  ONLY : iud0rho
!  USE thirdorder,ONLY : d3dyn
!  USE d3aux,     ONLY : d3dyn_aux9
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: iq_i, iq_j, iq_k
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  !
  INTEGER :: ir, is, ipert, jpert, kpert, crd, nspin0
  COMPLEX(DP), ALLOCATABLE :: drho_dipert(:,:), drho_djpert(:,:), &
                              drho_dkpert(:,:), gdrdi(:,:,:), gdrdj(:,:,:), &
                              gr_gdrdi(:), gr_gdrdj(:), aux1(:), &
                              aux2(:,:,:), d3dyn1(:,:,:)
  REAL(DP), ALLOCATABLE :: d2muxc(:)
  COMPLEX(DP) :: aux
  REAL(DP) :: rhotot !user_qpoint(3), 
  REAL(DP),EXTERNAL :: d2mxc
  !
!  IF( .NOT. dft_is_gradient() ) RETURN
  CALL start_clock('d3_exc_gc')
  !
  nspin0=nspin
  !
  IF( .NOT. nspin0 == 1) &
     CALL errore('d3_exc_gc', ' Gradient corrections implemented for nspin = 1 only ! ', 1 )
  !
  ALLOCATE(d2muxc( dfftp%nnr))    
  !
  ALLOCATE(drho_dipert( dfftp%nnr, nspin0 ))
  ALLOCATE(drho_djpert( dfftp%nnr, nspin0 ))
  ALLOCATE(drho_dkpert( dfftp%nnr, nspin0 ))
  ALLOCATE(gdrdi      ( 3, dfftp%nnr, nspin0))
  ALLOCATE(gdrdj      ( 3, dfftp%nnr, nspin0))
  ALLOCATE(gr_gdrdi   ( dfftp%nnr ))
  ALLOCATE(gr_gdrdj   ( dfftp%nnr ))
  ALLOCATE(aux1       ( dfftp%nnr ))
  ALLOCATE(aux2       ( 3, dfftp%nnr, nspin0 ))
  ALLOCATE(d3dyn1     ( 3*nat, 3*nat, 3*nat )) 
  !
  d2muxc = 0._dp
  DO ir = 1, dfftp%nnr
  rhotot = rho%of_r(ir, 1) + rho_core (ir)
  IF (rhotot >   1.e-30_dp) d2muxc (ir) =  d2mxc( rhotot)
  IF (rhotot < - 1.e-30_dp) d2muxc (ir) = -d2mxc(-rhotot)
  ENDDO
  !
  !
  d3dyn1(:,:,:) = 0._dp
  !
  IPERT_LOOP : &
  DO ipert = 1, 3*nat
    !
    drho_dipert(:,:) = 0._dp
    gdrdi(:,:,:) = 0._dp
    gr_gdrdi(:) = 0._dp
    !
    CALL read_drho(drho_dipert, iq_i, ipert, with_core=.true., pool_only = .true.)
    !
    DO is = 1, nspin0
        !
        CALL qgradient (kplusq(iq_i)%xq, dfftp%nnr, drho_dipert(:, is), ngm, g, nl, &
                        alat, gdrdi(:,:,is) )
        !
        FORALL(ir = 1:dfftp%nnr) 
            gr_gdrdi(ir) = DOT_PRODUCT(grho(:,ir,is), gdrdi(:,ir,is))
        ENDFORALL
        !
    ENDDO
    !
    JPERT_LOOP : &
    DO jpert = 1, 3*nat
        !
        drho_djpert(:,:) = 0._dp
        gdrdj(:,:,:)     = 0._dp
        gr_gdrdj(:)      = 0._dp
        aux1(:)          = 0._dp
        aux2(:,:,:)      = 0._dp
        !
        CALL read_drho(drho_djpert, iq_j, jpert, with_core=.true., pool_only = .true.)
        !
        DO is = 1, nspin0
            !
            CALL qgradient (kplusq(iq_j)%xq, dfftp%nnr, drho_djpert(:, is), ngm, g, nl, &
                            alat, gdrdj(:,:,is) )
            !
            FORALL(ir = 1:dfftp%nnr) 
                gr_gdrdj(ir) = SUM(grho(:,ir,is) * gdrdj(:,ir,is))
            ENDFORALL
            !
        ENDDO
        !
        IF(nspin0 == 1) THEN
            !
            DO ir = 1, dfftp%nnr
              !
              DO crd = 1, 3
                !
                aux2(crd,ir,1) = &
                        - dvxc_sr(ir,1,1) * &
                          ( drho_dipert(ir,1) * gdrdj(crd,ir,1) &
                           +drho_djpert(ir,1) * gdrdi(crd,ir,1) &
                          ) &
                        - dvxc_ss(ir,1,1) * &
                          ( gr_gdrdj(ir) * gdrdi(crd,ir,1) &
                           +gr_gdrdi(ir) * gdrdj(crd,ir,1) &
                           +SUM(gdrdi(:,ir,1) * gdrdj(:,ir,1)) * grho(crd,ir,1)  &
                           ) &
                        - dvxc_srr(ir,1,1) * & 
                            drho_dipert(ir,1) * drho_djpert(ir,1) * grho(crd,ir,1) &
                        - dvxc_ssr(ir,1,1) * &
                          ( gr_gdrdj(ir) * drho_dipert(ir,1) + &
                            gr_gdrdi(ir) * drho_djpert(ir,1) ) * grho(crd,ir,1) & 
                        - dvxc_sss(ir,1,1) * &
                              gr_gdrdi(ir) * gr_gdrdj(ir) * grho(crd,ir,1)
                !
              ENDDO
              !
            ENDDO
            !
        ENDIF
        !
        DO is = 1, nspin0
            !    
            CALL qgrad_dot (kplusq(iq_j)%xq, dfftp%nnr, aux2(:,:,is), ngm, g, nl, alat, aux1)
            !
        ENDDO
        !
        DO ir = 1, dfftp%nnr
            !
            IF(nspin0 == 1) THEN
                !
                aux1(ir) = aux1(ir) + &
                    dvxc_rrr(ir,1,1) * drho_dipert(ir,1) * drho_djpert(ir,1) &
                  + dvxc_srr(ir,1,1) * &
                          ( gr_gdrdj(ir) * drho_dipert(ir,1) &
                           +gr_gdrdi(ir) * drho_djpert(ir,1) ) &
                  + dvxc_sr(ir,1,1) * &
                      SUM( gdrdi(:,ir,1)*  gdrdj(:,ir,1) ) &
                  + dvxc_ssr(ir,1,1) * gr_gdrdi(ir) * gr_gdrdj(ir)
                !
            ENDIF
            !
        ENDDO
        !
        KPERT_LOOP : &
        DO kpert = 1, 3*nat
            !
            CALL read_drho(drho_dkpert, iq_k, kpert, with_core=.true., pool_only = .true.)
            !
            !IF(nspin0 == 1) THEN
                !
                aux = 0._dp
                DO ir = 1, dfftp%nnr
                !
                aux = aux &!+ d2muxc(ir)  * drho_dipert(ir,1) * drho_djpert(ir,1) * drho_dkpert(ir,1)  &
                          + aux1(ir) * drho_dkpert(ir,1)
                !
                ENDDO
                !
            !ENDIF              
            !
            CALL mp_sum( aux, intra_pool_comm )
            !
            d3dyn1(ipert,jpert,kpert) = omega * aux / REAL(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, kind=DP)
            !
        ENDDO KPERT_LOOP
        !
    ENDDO JPERT_LOOP
    !
     !ENDIF
     !
  ENDDO IPERT_LOOP
  !
  IF( npool > 1 ) CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
  !
  d3dyn = d3dyn+d3dyn1
  !
  !d3dyn = d3dyn  + d3dyn1
  !d3dyn_aux9 = d3dyn1
  !
  DEALLOCATE(d2muxc)
  DEALLOCATE(drho_dipert)
  DEALLOCATE(drho_djpert)
  DEALLOCATE(drho_dkpert)
  DEALLOCATE(gdrdi)
  DEALLOCATE(gdrdj)
  DEALLOCATE(gr_gdrdi)
  DEALLOCATE(gr_gdrdj)
  DEALLOCATE(aux1)
  DEALLOCATE(aux2)
  DEALLOCATE(d3dyn1) 
  !
  CALL stop_clock('d3_exc_gc')
  RETURN
  !
END SUBROUTINE d3_exc_gc
END MODULE d3_exc_gc_module
