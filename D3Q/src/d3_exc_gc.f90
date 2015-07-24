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
SUBROUTINE d3_exc_gc(d3dyn)
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
!  USE q0modes,   ONLY : q0mode
  USE units_ph,     ONLY : lrdrho, iudrho
!  USE units_d3,  ONLY : iud0rho
!  USE thirdorder,ONLY : d3dyn
!  USE d3aux,     ONLY : d3dyn_aux9
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  !
  INTEGER :: ir, is, ipert, jpert, kpert, crd, nspin0
  COMPLEX(DP), ALLOCATABLE :: drho_dipert(:,:), drho_djpert(:,:), &
                              drho_dkpert(:,:), gdrdi(:,:,:), gdrdj(:,:,:), &
                              gr_gdrdi(:), gr_gdrdj(:), aux1(:), &
                              aux2(:,:,:), d3dyn1(:,:,:)
  REAL(DP), ALLOCATABLE :: d2muxc(:)
  COMPLEX(DP) :: aux
  REAL(DP) :: d2mxc, rhotot !user_qpoint(3), 
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
  !d2muxc(:) = 0.d0
  !   
  DO ir = 1, dfftp%nnr
     !
     rhotot = rho%of_r(ir, 1) + rho_core(ir)
     !
     IF(rhotot > 1.d-30) THEN
        d2muxc(ir) = d2mxc(rhotot)
     ELSE IF(rhotot < - 1.d-30) THEN
        d2muxc(ir) = - d2mxc( - rhotot)
     ELSE
        d2muxc(ir) = 0._dp
     ENDIF
     !
  ENDDO
  !
  !
  !
  d3dyn1(:,:,:) =(0.d0, 0.d0)
  !
  IPERT_LOOP : &
  DO ipert = 1, 3*nat
     !
     !IF(q0mode(ipert) ) THEN
    ! 
    drho_dipert(:,:) =(0.d0, 0.d0)
    gdrdi(:,:,:) =(0.d0, 0.d0)
    gr_gdrdi(:) =(0.d0, 0.d0)
    !
    !LP CALL davcio_drho(drho_dipert, lrdrho, iud0rho, ipert, - 1)
    CALL read_drho(drho_dipert, 1, ipert, with_core=.true.)
    !
    !  CASE : xq = 0
    !
    !user_qpoint(:) = 0.d0
    !
    DO is = 1, nspin0
        !
!         CALL qgradient(kplusq(1)%xq, nrx1, nrx2, nrx3, nr1, nr2, nr3, dfftp%nnr, &
!             drho_dipert(:, is), ngm, g, nl, alat, gdrdi(:, :, is) )
        CALL qgradient (kplusq(1)%xq, dfftp%nnr, drho_dipert(:, is), ngm, g, nl, &
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
        drho_djpert(:,:) =(0.d0, 0.d0)
        gdrdj(:,:,:) =(0.d0, 0.d0)
        gr_gdrdj(:) =(0.d0, 0.d0)
        aux1(:) =(0.d0, 0.d0)
        aux2(:,:,:) =(0.d0, 0.d0)
        !
        !LP CALL davcio_drho(drho_djpert, lrdrho, iudrho, jpert, - 1)
        !LP drho_djpert(:,:) = CONJG( drho_djpert(:,:) )
        CALL read_drho(drho_djpert, 2, jpert, with_core=.true.)
        !
        !
        !  CASE : xq = -q
        !
        !user_qpoint(:) = -xq(:)
        !
        DO is = 1, nspin0
            !
!             CALL qgradient(kplusq(2)%xq, nrx1, nrx2, nrx3, nr1, nr2, nr3, dfftp%nnr, &
!                 drho_djpert(1, is), ngm, g, nl, alat, gdrdj(1, 1, is) )
            CALL qgradient (kplusq(2)%xq, dfftp%nnr, drho_djpert(:, is), ngm, g, nl, &
                            alat, gdrdj(:,:,is) )
            !
            FORALL(ir = 1:dfftp%nnr) 
                gr_gdrdj(ir) = DOT_PRODUCT(grho(:,ir,is), gdrdj(:,ir,is))
            ENDFORALL
!             DO ir = 1, dfftp%nnr
!                 !
!                 gr_gdrdj(ir) = grho(1,ir,is) * gdrdj(1,ir,is) + &
!                                grho(2,ir,is) * gdrdj(2,ir,is) + &
!                                grho(3,ir,is) * gdrdj(3,ir,is)
!                 !
!             ENDDO
            !
        ENDDO
        !
        IF(nspin0 == 1) THEN
            !
            DO ir = 1, dfftp%nnr
              !
              DO crd = 1, 3
                !
                aux2(crd,ir,1) = - dvxc_sr(ir,1,1) * &
                       ( drho_dipert(ir,1) * gdrdj(crd,ir,1) + &
                         drho_djpert(ir,1) * gdrdi(crd,ir,1) ) - &
                        dvxc_ss(ir,1,1) * &
                       ( gr_gdrdj(ir) * gdrdi(crd,ir,1) + &
                         gr_gdrdi(ir) * gdrdj(crd,ir,1) + &
                       ( gdrdi(1,ir,1) * gdrdj(1,ir,1) + &
                         gdrdi(2,ir,1) * gdrdj(2,ir,1) + &
                         gdrdi(3,ir,1) * gdrdj(3,ir,1) ) * grho(crd,ir,1) ) - &
                        dvxc_srr(ir,1,1) * & 
                            drho_dipert(ir,1) * drho_djpert(ir,1) * grho(crd,ir,1) - &
                        dvxc_ssr(ir,1,1) * &
                       ( gr_gdrdj(ir) * drho_dipert(ir,1) + &
                            gr_gdrdi(ir) * drho_djpert(ir,1) ) * grho(crd,ir,1) - & 
                        dvxc_sss(ir,1,1) * &
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
            !user_qpoint(:) = -xq(:)
            !    
            CALL qgrad_dot (kplusq(2)%xq, dfftp%nnr, aux2(:,:,is), ngm, g, nl, alat, aux1)
            !LP call qgrad_dot(kplusq(2)%xq, nrx1, nrx2, nrx3, nr1, nr2, nr3, dfftp%nnr, &
            !LP    aux2(1, 1, is), ngm, g, nl, alat, aux1)
            !
        ENDDO
        !
        DO ir = 1, dfftp%nnr
            !
            IF(nspin0 == 1) THEN
                !
                aux1(ir) = aux1(ir) + &
                dvxc_rrr(ir,1,1) * drho_dipert(ir,1) * drho_djpert(ir,1) + &
                dvxc_srr(ir,1,1) * &
                       ( gr_gdrdj(ir) * drho_dipert(ir,1) + &
                         gr_gdrdi(ir) * drho_djpert(ir,1) ) + &
                dvxc_sr(ir,1,1) * &
                       ( gdrdi(1,ir,1) * gdrdj(1,ir,1) + &
                         gdrdi(2,ir,1) * gdrdj(2,ir,1) + &
                         gdrdi(3,ir,1) * gdrdj(3,ir,1) ) + &
                dvxc_ssr(ir,1,1) * gr_gdrdi(ir) * gr_gdrdj(ir)
                !
            ENDIF
            !
        ENDDO
        !
        KPERT_LOOP : &
        DO kpert = 1, 3*nat
            !
            !  CASE xq = q 
            !
            !user_qpoint(:) = xq(:)
            !
            !CALL davcio_drho(drho_dkpert, lrdrho, iudrho, kpert, - 1)
            CALL read_drho(drho_dkpert, 3, kpert, with_core=.true.)
            !
            IF(nspin0 == 1) THEN
                !
                aux = CMPLX(0.d0, 0.d0)
                !
                DO ir = 1, dfftp%nnr
                !
                aux = aux + d2muxc(ir) * drho_dipert(ir,1) * drho_djpert(ir,1) * drho_dkpert(ir,1) + &
                        aux1(ir) * drho_dkpert(ir,1)
                !
                ENDDO
                !
            ENDIF              
            !
            CALL mp_sum( aux, intra_pool_comm )
            !
            d3dyn1(ipert,jpert,kpert) = omega * aux /(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
            !
        ENDDO KPERT_LOOP
        !
    ENDDO JPERT_LOOP
    !
     !ENDIF
     !
  ENDDO IPERT_LOOP
  !
  IF( npool /= 1 ) CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
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
