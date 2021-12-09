!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! #define __OLD_NONCOLIN_GGA
MODULE gc_d3
  !
  USE kinds, ONLY: DP
  !
  USE gc_lr,            ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  !
  REAL(DP), ALLOCATABLE :: &
       dvxc_rrr(:,:,:),          &! dfftp%nnr, nspin, nspin), 
       dvxc_srr(:,:,:),          &! dfftp%nnr, nspin, nspin),
       dvxc_ssr(:,:,:),          &! dfftp%nnr, nspin, nspin), 
       dvxc_sss(:,:,:)            ! dfftp%nnr, nspin, nspin),
  !
  CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE setup_d3gc
    !-----------------------------------------------------------------------
    ! wat20100930:
    ! Allocate and setup all variable needed in the gradient correction case
    !
    !
    !
    USE constants,            ONLY : e2
    USE gvect,                ONLY : ngm, g !, nl
    USE lsda_mod,             ONLY : nspin
    !USE spin_orb,             ONLY : domag
    USE scf,                  ONLY : rho, rho_core, rhog_core
    USE noncollin_module,     ONLY : noncolin,domag
    USE kinds,                ONLY : DP
    !USE funct,                ONLY : dft_is_gradient
    !USE dft_par_mod,  ONLY: isgradient
    USE xc_lib,               ONLY : xclib_dft_is
    USE qe_drivers_gga,       ONLY : gcxc
!     USE gc_ph,            ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  !   USE gc_d3,            ONLY : dvxc_rrr, dvxc_srr, &
  !                                dvxc_ssr, dvxc_sss
    USE uspp,                 ONLY : nlcc_any
    USE fft_base,             ONLY : dfftp
    USE fft_interfaces,       ONLY : fwfft
    USE wavefunctions,        ONLY : psic
    USE io_global,            ONLY : stdout
    USE noncollin_module,     ONLY : nspin_gga
    USE qe_drivers_d_gga,     ONLY : dgcxc_unpol, d3gcxc
    !USE xc_gga,               ONLY : xc_gcx
    !
    IMPLICIT NONE
    !
    INTEGER :: ir, is, nspin0
    REAL(DP) :: grho2 (2), fac, grhox(1,3,1), rhox(1,1), &
        sx(1), sc(1), v1x(1), v2x(1), v1c(1), v2c(1), &
        vrrx(1), vsrx(1), vssx(1), vrrc(1), vsrc(1), vssc(1), &
        dvxcrrr, dvxcsrr, dvxcssr, dvxcsss, &
        vrrrx, vsrrx, vssrx, vsssx, &
        vrrrc, vsrrc, vssrc, vsssc
    REAL(DP), ALLOCATABLE :: rho_tot_r(:,:)
    COMPLEX(DP), ALLOCATABLE :: rho_tot_g(:,:)
    REAL (DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
    !
    IF ( .NOT. xclib_dft_is('gradient') ) RETURN

    WRITE(stdout, '(5x,a)') "Setting up GGA 2nd derivative"
    grho2 = 0._dp
    
    nspin0=nspin
    !
    IF ( .NOT. nspin0 == 1) &
        CALL errore ('setup_d3gc', ' Gradient corrections implemented for nspin = 1 only ! ', 1 )
    !
    ALLOCATE (grho    (  3    , dfftp%nnr   , nspin0))
    ALLOCATE (rho_tot_r  (  dfftp%nnr , nspin0))
    ALLOCATE (rho_tot_g  (  ngm , nspin0))
    ALLOCATE (dvxc_rr (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_sr (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_ss (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_s  (  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_rrr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_srr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_ssr(  dfftp%nnr , nspin0 , nspin0))
    ALLOCATE (dvxc_sss(  dfftp%nnr , nspin0 , nspin0))
    !
    grho    (:,:,:) = 0._dp
    dvxc_s  (:,:,:) = 0._dp
    dvxc_rr (:,:,:) = 0._dp
    dvxc_sr (:,:,:) = 0._dp
    dvxc_ss (:,:,:) = 0._dp
    dvxc_rrr(:,:,:) = 0._dp
    dvxc_srr(:,:,:) = 0._dp
    dvxc_ssr(:,:,:) = 0._dp
    dvxc_sss(:,:,:) = 0._dp
    !
    fac = 1._dp / DBLE (nspin0)
    IF (noncolin.AND.domag) THEN
      call errore('setup_d3gc',' domag not implemented',1)
    ELSE
      !
      IF (.NOT. nlcc_any) THEN
        DO is = 1, nspin0
          rho_tot_r(:,is) = rho%of_r(:,is)
        ENDDO
      ELSE
        DO is = 1, nspin0
          rho_tot_r(:,is) = fac * rho_core(:)  + rho%of_r(:,is)
        ENDDO
      ENDIF
      !
      !
      DO is = 1, nspin0
        psic(:) = rho_tot_r(:,is)
        CALL fwfft ('Rho', psic, dfftp)
        rho_tot_g(:,is) = psic(dfftp%nl(:))
        !CALL gradrho( dfftp%nnr, rho_tot_g(1,is), ngm, g, dfftp%nl, grho(:,:,is) )
        CALL fft_gradient_g2r( dfftp, rho_tot_g(:,is), g, grho(:,:,is) )
      ENDDO
      !
    END IF
    
    !WHERE(ABS(grho)>1.d+32) grho = 0._dp
    
    DO ir = 1, dfftp%nnr
      !print*, rho%of_r(ir,1), grho (1,ir,1), grho(2,ir,1), grho(3,ir,1)
      grho2(1) = grho (1,ir,1)**2 + grho(2,ir,1)**2 + grho(3,ir,1)**2
      grhox(1,1,1) = grho (1,ir,1)
      grhox(1,2,1) = grho (2,ir,1)
      grhox(1,3,1) = grho (3,ir,1)
      rhox = rho_tot_r(ir,1)
      IF (nspin0 == 1) THEN
          IF (ABS(rho_tot_r(ir, 1) ) > epsr .AND. grho2 (1) > epsg) THEN
            CALL gcxc( 1, rho_tot_r(ir,1), grho2, sx, sc, v1x, v2x, v1c, v2c)
            CALL dgcxc_unpol (1, rho_tot_r(ir,1), grho2, &
                              vrrx, vsrx, vssx, vrrc, vsrc, vssc)
            dvxc_rr (ir, 1, 1) = e2 * (vrrx(1) + vrrc(1))
            dvxc_sr (ir, 1, 1) = e2 * (vsrx(1) + vsrc(1))
            dvxc_ss (ir, 1, 1) = e2 * (vssx(1) + vssc(1))
            dvxc_s  (ir, 1, 1) = e2 * (v2x(1) + v2c(1))
            CALL d3gcxc(rho_tot_r(ir,1), grho2(1), vrrrx, vsrrx, vssrx, vsssx, &
                                        vrrrc, vsrrc, vssrc, vsssc )
            !
            dvxc_rrr(ir, 1, 1) = e2 * (vrrrx + vrrrc)
            dvxc_srr(ir, 1, 1) = e2 * (vsrrx + vsrrc)
            dvxc_ssr(ir, 1, 1) = e2 * (vssrx + vssrc)
            dvxc_sss(ir, 1, 1) = e2 * (vsssx + vsssc)
            !
          ENDIF
      ELSE
          CALL errore('setup_d3gc',' nspin>1 not implemented',1)
      ENDIF
      !
    ENDDO
    
!     WRITE(10006, '(2f12.6)') dvxc_rrr
!     WRITE(10007, '(2f12.6)') dvxc_srr
!     WRITE(10008, '(2f12.6)') dvxc_ssr
!     WRITE(10009, '(2f12.6)') dvxc_sss
    
    IF (noncolin.AND.domag) &
      CALL errore('setup_d3gc',' domag not implemented',1)
    !
    DEALLOCATE(rho_tot_r, rho_tot_g)
    !
    RETURN
    !
  END SUBROUTINE setup_d3gc

END MODULE gc_d3
