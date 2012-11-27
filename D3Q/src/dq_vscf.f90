!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dq_vscf_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE dq_vscf(nu_i, dvloc, xq_x, iq_x, u_x)
  !-----------------------------------------------------------------------
  !
  !   It reads the variation of the charge density from a file and
  !   calculates the variation of the local part of the variation of the
  !   K-S potential.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : fpi, tpi, e2
  USE ions_base,   ONLY : nat, ityp, tau
  USE gvect,       ONLY : ngm, g, nl
  USE gvecs,       ONLY :  doublegrid
  USE fft_base,    ONLY : dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE cell_base,   ONLY : tpiba, tpiba2
  USE uspp_param,  ONLY : upf
  USE nlcc_ph,     ONLY : nlcc_any
  USE eqv,         ONLY : dmuxc 
  USE d3com,       ONLY : d3c, d3v
  USE d3_iofiles,  ONLY : read_drho
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
  IF (doublegrid) call cinterpolate (dvloc, dvloc, - 1)
  !
  !
  DEALLOCATE (aux1, aux2)
  !
  CALL stop_clock('dq_vscf')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE dq_vscf
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
END MODULE dq_vscf_module
!-----------------------------------------------------------------------
