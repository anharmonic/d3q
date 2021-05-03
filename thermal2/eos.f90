!
!
! Copyright (C) 2003-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Adapted from origina PW/src/ev.f90 written by Eyvaz Isaev
! Dept of Physics, Chemistry and Biology (IFM), Linkoping University, Sweden
!
! Lorenzo Paulatto, 2019
!
! <<^V^\\=========================================//-//-//========//O\\//


MODULE eos
  USE kinds,     ONLY : DP
  USE constants, ONLY : bohr_radius_angs, ry_kbar
  REAL(DP), PARAMETER :: gpa_kbar = 10.0_dp
  INTEGER, PARAMETER  :: maxpar=4

  CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE eqstate(istat,npar,par,npt,v0,etot,efit,emin,chisq)
!-----------------------------------------------------------------------
! istat:
!  1: birch 1st order
!  2: birch 3rd order
!  3: keane       
!  4: murnaghan
!
      IMPLICIT NONE
      INTEGER, INTENT(in)  :: istat, npar, npt
      REAL(DP), INTENT(in) :: par(maxpar)
      REAL(DP), INTENT(out):: chisq
      REAL(DP), INTENT(in) :: v0(npt),etot(npt)
      REAL(DP),INTENT(out) :: efit(npt), emin
      INTEGER :: i
      REAL(DP) :: k0, dk0, d2k0, c0, c1, x, vol0, ddk
!
      vol0 = par(1)
      k0   = par(2)/ry_kbar ! converts k0 to Ry atomic units...
      dk0  = par(3)
      ! Only used by istat 2 and 3:
      d2k0 = par(4)*ry_kbar ! and d2k0/dp2 to (Ry a.u.)^(-1)
!
      IF(istat==1.or.istat==2) THEN
         IF(istat==1) THEN
            c0 = 0.0d0
         ELSE
            c0 = ( 9.d0*k0*d2k0 + 9.d0*dk0**2-63.d0*dk0+143.d0 )/48.d0
         ENDIF
         c1 = 3.d0*(dk0-4.d0)/8.d0
         DO i=1,npt
            x = vol0/v0(i)
            efit(i) = 9.d0*k0*vol0*( (-0.5d0+c1-c0)*x**(2.d0/3.d0)/2.d0 &
                         +( 0.50-2.d0*c1+3.d0*c0)*x**(4.d0/3.d0)/4.d0 &
                         +(       c1-3.d0*c0)*x**(6.d0/3.d0)/6.d0 &
                         +(            c0)*x**(8.d0/3.d0)/8.d0 &
                         -(-1.d0/8.d0+c1/6.d0-c0/8.d0) )
         ENDDO
      ELSE IF( istat==3 .or. istat==4) THEN
         IF(istat==3) THEN
            ddk = dk0 + k0*d2k0/dk0
         ELSE
            ddk = dk0
         ENDIF
         DO i=1,npt
            efit(i) = - k0*dk0/ddk*vol0/(ddk-1.d0) &
            + v0(i)*k0*dk0/ddk**2*( (vol0/v0(i))**ddk/(ddk-1.d0)+1.d0) &
            - k0*(dk0-ddk)/ddk*( v0(i)*log(vol0/v0(i)) + v0(i)-vol0 )
         ENDDO
      ELSE
        CALL errore("eqstate","wrong equation of state",1)
      ENDIF
!
!      emin = equilibrium energy obtained by minimizing chi**2
!
      emin = 0.0d0
      DO i = 1,npt
         emin = emin + etot(i)-efit(i)
      ENDDO
      emin = emin/npt
!
      chisq = 0.0d0
      DO i = 1,npt
          efit(i) = efit(i)+emin
          chisq   = chisq + (etot(i)-efit(i))**2
      ENDDO
      chisq = chisq/npt
!
      RETURN
    END SUBROUTINE eqstate
!
!

!-----------------------------------------------------------------------
      SUBROUTINE find_minimum2(istat,npar,par,npt,v0,etot,efit,emin,chisq)
!-----------------------------------------------------------------------
!
      USE random_numbers, ONLY : randy
      USE lmdif_module, ONLY : lmdif0
      IMPLICIT NONE
      INTEGER ,INTENT(in)  :: istat, npar, npt
      REAL(DP),INTENT(in)  :: v0(npt),etot(npt)
      REAL(DP),INTENT(out) :: efit(npt), emin
      REAL(DP),INTENT(out) :: par(maxpar)
      REAL(DP),INTENT(out) :: chisq
      !
      REAL(DP) :: ediff(npt)

      INTEGER :: i

      par(1) = v0(npt/2)
      par(2) = 500.0d0
      par(3) = 5.0d0
      par(4) = -0.01d0
      !
      CALL lmdif0(EOSDIFF, npt, npar, par, ediff, 1.d-12, i)
      !
      CALL eqstate(istat,npar,par,npt,v0,etot,efit,emin,chisq)
      

      CONTAINS

      ! This subroutine is passed to LMDIF to be minimized
      ! LMDIF takes as input the difference between f_fit and f_real
      !       and computes the chi^2 internally.
      SUBROUTINE EOSDIFF(m_, n_, par_, f_, i_)
            IMPLICIT NONE
              INTEGER,INTENT(in)  :: m_, n_
              INTEGER,INTENT(inout)   :: i_
              REAL(DP),INTENT(in)    :: par_(n_)
              REAL(DP),INTENT(out)   :: f_(m_)
              REAL(DP) :: chisq_
              !
              CALL eqstate(istat,n_,par_,npt,v0,etot,efit,emin,chisq_)
              f_ = etot - efit
            END SUBROUTINE
   
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE find_minimum(istat,npar,par,npt,v0,etot,efit,emin,chisq)
!-----------------------------------------------------------------------
!
!     Very Stupid Minimization
!
      USE random_numbers, ONLY : randy
      IMPLICIT NONE
      INTEGER ,INTENT(in)  :: istat, npar, npt
      REAL(DP),INTENT(in)  :: v0(npt),etot(npt)
      REAL(DP),INTENT(out) :: efit(npt), emin
      REAL(DP),INTENT(out) :: par(maxpar)
      REAL(DP),INTENT(out) :: chisq
      
      INTEGER  :: n,j,i
      INTEGER, PARAMETER :: nmin=10, nseek=10000
      REAL(DP) :: parnew(maxpar), &
                  parmin(maxpar), deltapar(maxpar), parmax(maxpar)
      REAL(DP) :: chinew
!
!      various initializations
!
!       print*, "istat", istat
!       print*, "npar", npar
!       print*, "npt", npt
!       print*, "v0", v0
!       print*, "etot", etot

      par(1) = v0(npt/2)
      par(2) = 500.0d0
      par(3) = 5.0d0
      par(4) = -0.01d0
!
      parmin(1) = 0.0d0
      parmin(2) = 0.0d0
      parmin(3) = 1.0d0
      parmin(4) = -1.0d0
!
      parmax(1) = 100000.d0
      parmax(2) = 100000.d0
      parmax(3) = 15.0d0
      parmax(4) = 0.0d0
!
      deltapar(1) = 1.0d0
      deltapar(2) = 100.d0
      deltapar(3) = 1.0d0
      deltapar(4) = 0.01d0

      chisq = 1.0d30
      chinew= 1.0d30
      CALL eqstate(istat,npar,par,npt,v0,etot,efit,emin,chisq)
      DO j = 1,nmin
         DO i = 1,nseek
            DO n = 1,npar
10             parnew(n) = par(n) + (0.5d0 - randy())*deltapar(n)
               IF(parnew(n)>parmax(n) .or. parnew(n)<parmin(n)) &
                  GOTO 10
            ENDDO
!
            CALL eqstate(istat,npar,parnew,npt,v0,etot,efit,emin,chinew)
!
            IF(chinew<chisq) THEN
               DO n = 1,npar
                  par(n) = parnew(n)
               ENDDO
               chisq = chinew
            ENDIF
         ENDDO
         DO n = 1,npar
            deltapar(n) = deltapar(n)/10.d0
         ENDDO
      ENDDO
!
      CALL eqstate(istat,npar,par,npt,v0,etot,efit,emin,chisq)
!
!       print*, "par", par
!       print*, "efit", efit
!       print*, "emin", emin
!       print*, "chisq", chisq
!       stop 0
      RETURN
    END SUBROUTINE find_minimum

END MODULE eos
! <<^V^\\=========================================//-//-//========//O\\//


