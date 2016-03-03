!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE dq_vscf_module

  INTERFACE dq_vscf
!     MODULE PROCEDURE dq_vscf_nuovo
    MODULE PROCEDURE dq_vscf_vecchio
  END INTERFACE dq_vscf


  CONTAINS
  
!   SUBROUTINE dq_vscf(nu_i, dvloc, xq_x, iq_x, u_x)
!     USE kinds,            ONLY : DP
!     USE ions_base,        ONLY : nat    
!     USE fft_base,         ONLY : dfftp
!     implicit none
!     INTEGER,INTENT(in)     :: nu_i ! index of the mode (1...3*nat)
!     INTEGER,INTENT(in)     :: iq_x ! index of the q vector
!     REAL(DP),INTENT(in)    :: xq_x (3)! input: coordinates of the q point
!     COMPLEX(DP),INTENT(in) :: u_x(3*nat, 3*nat)! displacement patterns for iq
!     COMPLEX(DP),INTENT(out) :: dvloc(dfftp%nnr)
!     !
!     COMPLEX(DP) :: dvloc_nuovo(dfftp%nnr)
!     COMPLEX(DP) :: dvloc_vecchio(dfftp%nnr)
!     INTEGER :: i
! 
!     CALL dq_vscf_nuovo(nu_i,   dvloc_nuovo, xq_x, iq_x, u_x)
!     CALL dq_vscf_vecchio(nu_i, dvloc_vecchio, xq_x, iq_x, u_x)
!     !
!     
!     DO i=1,dfftp%nnr
!       WRITE(10001, '(i6,2(2f12.6,4x),f12.6)') i,dvloc_nuovo(i), dvloc_vecchio(i), &
!       ABS(dvloc_nuovo(i)-dvloc_vecchio(i))*1000
!     ENDDO
!     !
!     dvloc = dvloc_nuovo
!     !
!   END SUBROUTINE dq_vscf
  
  
  !-----------------------------------------------------------------------
  SUBROUTINE dq_vscf_nuovo(nu_i, dvloc, xq_x, iq_x, u_x)
    !-----------------------------------------------------------------------
    !
    !   It reads the variation of the charge density from a file and
    !   calculates the variation of the local part of the variation of the
    !   K-S potential.
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : fpi, tpi, e2
    USE ions_base,        ONLY : nat, ityp, tau
    USE cell_base,        ONLY : tpiba, tpiba2, alat
    USE gvect,            ONLY : ngm, g, nl
    USE gvecs,            ONLY :  doublegrid
    USE fft_base,         ONLY : dfftp
    USE fft_interfaces,   ONLY : fwfft, invfft
    USE uspp_param,       ONLY : upf
    USE uspp,             ONLY : nlcc_any
    USE eqv,              ONLY : dmuxc 
    USE d3com,            ONLY : d3c, d3v
    USE d3_iofiles,       ONLY : read_drho!, addcore_d3
    USE scf,              ONLY : rho, rho_core
    USE gc_lr,            ONLY : grho
    USE gc_d3,            ONLY : dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
    USE d3com,            ONLY : d3c
    USE noncollin_module, ONLY : nspin_lsda, nspin_gga
    USE pwcom,            ONLY : nspin
    USE funct,            ONLY : dft_is_gradient
    !
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nu_i ! index of the mode (1...3*nat)
    INTEGER,INTENT(in)     :: iq_x ! index of the q vector
    REAL(DP),INTENT(in)    :: xq_x (3)! input: coordinates of the q point
    COMPLEX(DP),INTENT(in) :: u_x(3*nat, 3*nat)! displacement patterns for iq
    COMPLEX(DP),INTENT(out) :: dvloc(dfftp%nnr)
    !
    ! Local variables
    INTEGER :: ig, ir, mu, na, nt, is, is1
    REAL(DP) :: qg2, gtau, fac
    COMPLEX(DP) ::  guexp    
    COMPLEX(DP), ALLOCATABLE :: drho_tot (:), rho_tot (:), dvloc_g(:)
    !
    ! After 8 years, I still do not know why we need 3 spin variables
    IF (nspin/=1)      CALL errore('dq_vscf', 'nspin /= 1 not implemented', 1)
    IF (nspin_lsda/=1) CALL errore('dq_vscf','nspin_lsda /= 1 not implemented',1)
    IF (nspin_gga/=1)  CALL errore('dq_vscf','nspin_gga /= 1 not implemented',1)
  
    ALLOCATE ( drho_tot (dfftp%nnr) )
    ALLOCATE ( dvloc_g (dfftp%nnr) )
    ALLOCATE ( rho_tot (dfftp%nnr) )
 
    dvloc (:) = (0.d0, 0.d0)
    CALL read_drho(drho_tot, iq_x, nu_i, with_core=.true.)
    !CALL davcio_drho (drho_tot, lrdrho, iudrho_x, nu_i, - 1)
  
    fac = 1.d0 / DBLE (nspin_lsda)
    if (nlcc_any) then
      DO is = 1, nspin_lsda
          rho_tot(:) = rho%of_r(:, is) + fac * rho_core (:)
      ENDDO
    ENDIF
    !
    ! Add Exchange-Correlation contribution in real space
    dvloc(:) = drho_tot(:) * dmuxc(:,1,1)

    ! When spin will be implemented, do this:
!     DO is = 1, nspin
!       DO is1 = 1, nspin
!           DO ir = 1, nrxx
!             dvloc(ir,is) = dvloc(ir,is) + dmuxc(ir,is,is1) * drho_tot(ir,is1)
!           ENDDO
!       ENDDO
!     ENDDO
    !
    ! add gradient correction to xc, NB: if nlcc is true we need to add here
    ! its contribution. grho contains already the core charge
    !
    IF ( dft_is_gradient() ) THEN
       !PRINT*, "********* dgradcorr is in"
       CALL dgradcorr(rho_tot, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq_x, &
                      drho_tot, dfftp%nnr, nspin, nspin_gga, nl, ngm, g, alat, &
                      dvloc)
    ENDIF
    !
    ! Remove the core correction, if necessary
    IF(nlcc_any)THEN
      !CALL read_drho(drho_tot, iq_x, nu_i, with_core=.false.)
      CALL addcore_ofq(xq_x, nu_i, u_x, d3c(iq_x)%drc, drho_tot, -1)
      !CALL addcore_d3(xq_x, u_x, nu_i, d3c(iq_x)%drc, drho_tot, -1)
    ENDIF
    !
    ! copy the total (up+down) delta rho in drho_tot(*,1) and go to G-space
    !if (nspin == 2) drho_tot(:,1) = drho_tot(:,1) + drho_tot(:,2) 
    !
    !CALL cft3 (drho_tot(:), nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
    CALL fwfft('Dense', drho_tot, dfftp) 
    !
!     DO is = 1, nspin_lsda
      !
      ! Transforms the potential to G space
      !
      !CALL cft3 (dvloc, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
      !CALL fwfft('Dense', dvloc, dfftp) 
      dvloc_g = 0._dp
      !
      ! Hartree contribution in G space
      DO ig = 1, ngm
          qg2 = (g(1,ig)+xq_x(1))**2 + (g(2,ig)+xq_x(2))**2 + (g(3,ig)+xq_x(3))**2
          IF (qg2 > 1.d-8) THEN
            dvloc_g(nl(ig)) = dvloc_g(nl(ig)) + &
                                e2 * fpi * drho_tot(nl(ig)) / (tpiba2 * qg2)
          ENDIF
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! This is the derivative of the local potential, in G space
      DO na = 1, nat
          mu = 3 * (na - 1)
          IF ( (abs(u_x(mu+1,nu_i))) + abs(u_x(mu+2,nu_i)) + abs(u_x(mu+3,nu_i)) > 1.0d-12) THEN
            nt = ityp (na)
            DO ig = 1, ngm
                gtau = tpi * DOT_PRODUCT( g(:,ig)+xq_x(:), tau(:,na) )
                guexp = tpiba * ( (g(1,ig) + xq_x(1)) * u_x(mu+1,nu_i) + &
                                  (g(2,ig) + xq_x(2)) * u_x(mu+2,nu_i) + &
                                  (g(3,ig) + xq_x(3)) * u_x(mu+3,nu_i) ) &
                     * CMPLX(0.d0,-1.d0) * CMPLX(cos(gtau),-sin(gtau))
                        !* CMPLX(SIN(gtau),-COS(gtau))
                dvloc_g(nl(ig)) = dvloc_g(nl(ig)) + d3v(iq_x)%loc(ig,nt) * guexp
            ENDDO
          ENDIF
      ENDDO
      !
      !  Transforms back to real space
      !
      !CALL cft3 (dvloc, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
      CALL  invfft('Dense', dvloc_g, dfftp)
      !
      ! Add all the terms
      dvloc = dvloc + dvloc_g
      !
      if (doublegrid) CALL cinterpolate (dvloc, dvloc, - 1)
!     ENDDO


    DEALLOCATE (drho_tot)
    DEALLOCATE (rho_tot)
  
    RETURN
  END SUBROUTINE dq_vscf_nuovo
  
!-----------------------------------------------------------------------
SUBROUTINE dq_vscf_vecchio(nu_i, dvloc, xq_x, iq_x, u_x)
  !-----------------------------------------------------------------------
  !
  !   It reads the variation of the charge density from a file and
  !   calculates the variation of the local part of the variation of the
  !   K-S potential.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : fpi, tpi, e2
  USE ions_base,   ONLY : nat, ityp, tau
  USE cell_base,   ONLY : alat
  USE gvect,       ONLY : ngm, g, nl
  USE gvecs,       ONLY :  doublegrid
  USE fft_base,    ONLY : dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE cell_base,   ONLY : tpiba, tpiba2
  USE uspp_param,  ONLY : upf
  USE uspp,        ONLY : nlcc_any
  USE eqv,         ONLY : dmuxc 
  USE d3com,       ONLY : d3c, d3v
  USE d3_iofiles,  ONLY : read_drho
  USE gc_lr,            ONLY : dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s, grho
    USE scf,          ONLY : rho
    USE funct,        ONLY : dft_is_gradient
    USE noncollin_module, ONLY : nspin_gga, nspin_mag
  !
  IMPLICIT NONE
  INTEGER,INTENT(in)     :: nu_i ! mode under consideration
  ! index of the perturbation
  INTEGER,INTENT(in)     :: iq_x ! index of the q vector
  ! index of q point
  REAL(DP),INTENT(in)    :: xq_x (3)
  ! input: coordinates of the q point
  COMPLEX(DP),INTENT(in) :: u_x(3*nat, 3*nat)
  ! displacement patterns for iq
  COMPLEX(DP),INTENT(out) :: dvloc(dfftp%nnr)
  ! output: local part of the variation
  !         of the K_S potential
  !
  ! Local variables
  INTEGER :: ig, mu, na, nt ! counters
  REAL (DP) :: qg2, gtau
  ! the modulus of (q+G)^2
  ! auxiliary variable: g*tau
  COMPLEX(DP) ::  guexp   ! auxiliary variable: g*u*exp(gtau)
  COMPLEX(DP),PARAMETER :: mii = (0._dp, -1._dp)
  COMPLEX(DP),ALLOCATABLE :: aux1(:), aux2(:)
  !
  CALL start_clock('dq_vscf')
  !
  ALLOCATE(aux1(dfftp%nnr))
  ALLOCATE(aux2(dfftp%nnr))
  !
  CALL read_drho(aux2, iq_x, nu_i, with_core=.false.)
  !
  dvloc(:) = aux2(:) * dmuxc(:,1,1)
  CALL fwfft('Dense', aux2, dfftp) 

  aux1 = (0._dp, 0._dp)
  DO ig = 1, ngm
     qg2 = (g(1,ig)+xq_x(1))**2 + (g(2,ig)+xq_x(2))**2 + (g(3,ig)+xq_x(3))**2
     IF (qg2 > 1.d-8) THEN
        aux1(nl(ig)) = e2 * fpi * aux2(nl(ig)) / (tpiba2 * qg2)
     ENDIF
  ENDDO

  IF (nlcc_any) aux2= (0._dp, 0._dp)
  DO na = 1, nat
     mu = 3 * (na - 1)
     IF (ABS(u_x(mu+1,nu_i)) + ABS(u_x(mu+2,nu_i)) + &
         ABS(u_x(mu+3,nu_i)) > 1.d-12) THEN
        nt = ityp (na)
        DO ig = 1, ngm

           gtau = tpi * ( (g(1,ig) + xq_x(1)) * tau(1,na) + &
                          (g(2,ig) + xq_x(2)) * tau(2,na) + &
                          (g(3,ig) + xq_x(3)) * tau(3,na) )

           guexp = tpiba * ( (g(1,ig) + xq_x(1)) * u_x(mu+1,nu_i) + &
                             (g(2,ig) + xq_x(2)) * u_x(mu+2,nu_i) + &
                             (g(3,ig) + xq_x(3)) * u_x(mu+3,nu_i) ) &
                         * mii * EXP( mii * gtau)
           aux1 (nl(ig)) = aux1 (nl(ig)) + d3v(iq_x)%loc(ig,nt) * guexp
           IF (upf(nt)%nlcc) THEN
              aux2 (nl(ig)) = aux2 (nl(ig)) + d3c(iq_x)%drc(ig,nt) * guexp
           END IF
        ENDDO
     ENDIF

  ENDDO
  !
  CALL  invfft('Dense', aux1, dfftp)

  dvloc(:) = dvloc(:) + aux1 (:)
  IF (nlcc_any) THEN
     CALL invfft('Dense', aux2, dfftp)
     dvloc (:) = dvloc(:) + aux2 (:) * dmuxc(:,1,1)
  ENDIF
  !
  IF ( dft_is_gradient() ) THEN
    CALL read_drho(aux2, iq_x, nu_i, with_core=.true.)
    
!     CALL dgradcorr(rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq_x, &
!                    aux2, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvloc)
  ENDIF

  IF (doublegrid) call cinterpolate (dvloc, dvloc, - 1)
  !
  !
  DEALLOCATE (aux1, aux2)
  !
  CALL stop_clock('dq_vscf')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE dq_vscf_vecchio
!-----------------------------------------------------------------------  
  
  

END MODULE dq_vscf_module
