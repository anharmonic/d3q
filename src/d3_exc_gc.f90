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
  !USE funct,        ONLY : dft_is_gradient
  !USE dft_par_mod,  ONLY: isgradient
  USE xc_lib,       ONLY : xclib_dft_is
  USE fft_base,     ONLY : dfftp
  USE gvect,        ONLY : g, ngm !, nl
  USE qpoint,       ONLY : xq
  USE gc_lr,        ONLY : grho, dvxc_s, dvxc_sr, dvxc_ss
  USE gc_d3,        ONLY : dvxc_rrr, dvxc_srr, &
                           dvxc_ssr, dvxc_sss
  USE io_global,    ONLY : ionode_id, stdout
  USE mp_pools,     ONLY : inter_pool_comm, my_pool_id, &
                           npool, intra_pool_comm
  USE mp,           ONLY : mp_bcast, mp_sum 
  USE d3_iofiles,   ONLY : read_drho
  USE kplus3q,      ONLY : kplusq, q_names
  USE units_ph,     ONLY : lrdrho, iudrho
  USE mp_world, ONLY : mpime
!  USE dgradcorr_module, ONLY : dgradcorr, qgradient, qgrad_dot
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: iq_i, iq_j, iq_k
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  !
  INTEGER :: ir, is, ipert, jpert, kpert, crd, nspin0
  COMPLEX(DP), ALLOCATABLE :: drho_dipert(:), drho_djpert(:), &
                              drho_dkpert(:), gdrdi(:,:), gdrdj(:,:), &
                              gr_gdrdi(:), gr_gdrdj(:), aux1(:), &
                              aux2(:,:), d3dyn1(:,:,:)
  REAL(DP), ALLOCATABLE :: d2muxc(:)
  COMPLEX(DP) :: aux
  REAL(DP) :: rhotot, domega !user_qpoint(3), 
  REAL(DP), PARAMETER :: epsr = 1.0d-5, epsg = 1.0d-10
  REAL(DP),EXTERNAL :: d2mxc
  !
  IF( .NOT. xclib_dft_is('gradient') ) RETURN
  CALL start_clock('d3_exc_gc')
  !
  nspin0=nspin
  domega = omega/ (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  !
  IF( .NOT. nspin0 == 1) &
     CALL errore('d3_exc_gc', ' Gradient corrections implemented for nspin = 1 only ! ', 1 )
  !
  !
  ALLOCATE(d3dyn1( 3*nat, 3*nat, 3*nat )) 
  d3dyn1(:,:,:) = 0._dp
  !
  ! This term is computed only on the first pool, better parallelisation is possible
  ! but this method keeps IO at a reasonable level
  ONLY_FIRST_POOL : &
  IF ( my_pool_id == 0 ) THEN
    !
    ALLOCATE(drho_dipert( dfftp%nnr ))
    ALLOCATE(drho_djpert( dfftp%nnr ))
    ALLOCATE(drho_dkpert( dfftp%nnr ))
    ALLOCATE(gdrdi      ( 3, dfftp%nnr))
    ALLOCATE(gdrdj      ( 3, dfftp%nnr))
    ALLOCATE(gr_gdrdi   ( dfftp%nnr ))
    ALLOCATE(gr_gdrdj   ( dfftp%nnr ))
    ALLOCATE(aux1       ( dfftp%nnr ))
    ALLOCATE(aux2       ( 3, dfftp%nnr ))
    !
    ALLOCATE (d2muxc(dfftp%nnr))
    d2muxc = 0._dp
    DO ir = 1, dfftp%nnr
      rhotot = rho%of_r(ir, 1) + rho_core (ir)
      IF (rhotot >   1.e-30_dp) d2muxc (ir) =  d2mxc( rhotot)
      IF (rhotot < - 1.e-30_dp) d2muxc (ir) = -d2mxc(-rhotot)
    ENDDO
    !
    WRITE(stdout,'(7x,3a)') q_names(iq_i), q_names(iq_j), q_names(iq_k)
    !
    !
    IPERT_LOOP : &
    DO ipert = 1, 3*nat
      !
      CALL read_drho(drho_dipert, iq_i, ipert, with_core=.true., pool_only = .true.)
!      CALL qgradient(kplusq(iq_i)%xq, dfftp%nnr, drho_dipert, ngm, g, dfftp%nl, &
!                    alat, gdrdi )
      CALL fft_qgradient (dfftp, drho_dipert, kplusq(iq_i)%xq, g, gdrdi)
      FORALL(ir = 1:dfftp%nnr) gr_gdrdi(ir) = DOT_PRODUCT(grho(:,ir,1), gdrdi(:,ir))
      !
      JPERT_LOOP : &
      DO jpert = 1, 3*nat
          !
          CALL read_drho(drho_djpert, iq_j, jpert, with_core=.true., pool_only = .true.)
!          CALL qgradient (kplusq(iq_j)%xq, dfftp%nnr, drho_djpert, ngm, g, dfftp%nl, &
!                          alat, gdrdj )
          CALL fft_qgradient (dfftp, drho_djpert, kplusq(iq_j)%xq, g, gdrdj)
          FORALL(ir = 1:dfftp%nnr) gr_gdrdj(ir) = SUM(grho(:,ir,1) * gdrdj(:,ir))
          !
          DO ir = 1, dfftp%nnr
            !aux2(:,ir,1) = 0._dp
            IF( ABS(rho%of_r(ir,1)) > epsr .and. SUM(grho(:,ir,1)**2)>epsg )THEN
              DO crd = 1, 3
                !
                aux2(crd,ir) = &
                        - dvxc_sr(ir,1,1) * &
                          ( drho_dipert(ir) * gdrdj(crd,ir) &
                          +drho_djpert(ir) * gdrdi(crd,ir) &
                          ) &
                        - dvxc_ss(ir,1,1) * &
                          ( gr_gdrdj(ir) * gdrdi(crd,ir) &
                          +gr_gdrdi(ir) * gdrdj(crd,ir) &
                          +SUM(gdrdi(:,ir) * gdrdj(:,ir)) * grho(crd,ir,1)  &
                          ) &
                        - dvxc_srr(ir,1,1) * & 
                            drho_dipert(ir) * drho_djpert(ir) * grho(crd,ir,1) &
                        - dvxc_ssr(ir,1,1) * &
                          ( gr_gdrdj(ir) * drho_dipert(ir) + &
                            gr_gdrdi(ir) * drho_djpert(ir) ) * grho(crd,ir,1) & 
                        - dvxc_sss(ir,1,1) * &
                              gr_gdrdi(ir) * gr_gdrdj(ir) * grho(crd,ir,1)
              ENDDO
            ELSE
              aux2(:,ir) = (0._dp, 0._dp)
            ENDIF
          ENDDO
          !
!          CALL qgrad_dot(kplusq(-iq_k)%xq, dfftp%nnr, aux2, ngm, g, dfftp%nl, alat, aux1)
          CALL fft_qgraddot(dfftp, aux2, kplusq(-iq_k)%xq, g, aux1)
          !
          DO ir = 1, dfftp%nnr
            !
            aux1(ir) = aux1(ir) + &
                dvxc_rrr(ir,1,1) * drho_dipert(ir) * drho_djpert(ir) &
              + dvxc_srr(ir,1,1) * &
                      ( gr_gdrdj(ir) * drho_dipert(ir) &
                      +gr_gdrdi(ir) * drho_djpert(ir) ) &
              + dvxc_sr(ir,1,1) * SUM( gdrdi(:,ir)*  gdrdj(:,ir) ) &
              + dvxc_ssr(ir,1,1) * gr_gdrdi(ir) * gr_gdrdj(ir)
            !
          ENDDO
          !
          KPERT_LOOP : &
          DO kpert = 1, 3*nat
              !
              CALL read_drho(drho_dkpert, iq_k, kpert, with_core=.true., pool_only = .true.)
              !
  !             aux = SUM(aux1(:) * drho_dkpert(:))
              !
              ! LDA contribution temporarily moved to another subroutine for debugging
              aux = 0._dp
              DO ir = 1, dfftp%nnr
                aux = aux + aux1(ir) * drho_dkpert(ir) &
                          + d2muxc(ir)*drho_dipert(ir)*drho_djpert(ir)*drho_dkpert(ir)
                          
              ENDDO
              !
              !CALL mp_sum( aux, intra_pool_comm )
              !
              d3dyn1(ipert,jpert,kpert) = domega * aux
              !
          ENDDO KPERT_LOOP
          !
      ENDDO JPERT_LOOP
      !
      !ENDIF
      !
    ENDDO IPERT_LOOP
    !
    DEALLOCATE(drho_dipert, drho_djpert, drho_dkpert)
    DEALLOCATE(gdrdi, gdrdj)
    DEALLOCATE(gr_gdrdi,gr_gdrdj)
    DEALLOCATE(aux1, aux2)
    DEALLOCATE(d2muxc)
    !
    ! SUM inside pool, better to do it once here than (3*nat)**3 times above
    CALL mp_sum( d3dyn1,  intra_pool_comm )
    !
  ELSE
    d3dyn1 = 0._dp
  ENDIF ONLY_FIRST_POOL
  !
  ! Let everybody get the correct matrix 
  ! (this is not strictly needed as only cpu 1 writes)
  IF( npool > 1 ) CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
  !
  d3dyn = d3dyn+d3dyn1
  DEALLOCATE(d3dyn1)
  !
  CALL stop_clock('d3_exc_gc')
  RETURN
  !
END SUBROUTINE d3_exc_gc
END MODULE d3_exc_gc_module
