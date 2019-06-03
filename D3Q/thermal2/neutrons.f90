! Written by Lorenzo Paulatto (2018) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! This module contains the required subroutines to compute the effective cross section 
! of inelastic neutron scattering experiments
MODULE neutrons

  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info

 CONTAINS

  REAL(DP) FUNCTION neutron_form_factor(aname, amass)
  USE constants, ONLY : BOHR_RADIUS_SI
  IMPLICIT NONE
    ! Data from http://www.ncnr.nist.gov/resources/n-lengths/
    ! This software is open source, do not complain about the 
    ! table being incomplete: fill it yourself.
    !
    ! internal units are in femto-meters (aka fm, fermi)
    ! output is in bohr
    !
    CHARACTER(len=*),INTENT(in)  :: aname
    REAL(DP),OPTIONAL,INTENT(in) :: amass
    !
    REAL(DP),PARAMETER :: FERMI_TO_BOHR = 1.d-15/BOHR_RADIUS_SI
    !
    TYPE ffac
      CHARACTER(len=2) :: aname
      REAL(DP)         :: amass
      REAL(DP)         :: coh_b  ! coherent scattering length (fm)
      REAL(DP)         :: coh_xs ! coherent scattering cross section (barn)
    END TYPE ffac
    !
    INTEGER :: i_iso, i_type
    INTEGER,PARAMETER :: n_iso=6, n_types=100
    TYPE(ffac) :: factors(0:10,n_types)
    !
    ! Initialize
    factors = ffac( "", 0._dp, 0._dp, 0._dp)
    ! Hydrogen:
!     factors(0,1) = ffac( "H", 1.008_dp, -3.7390_dp)
!    ` factors(1,1) = ffac( "H", 1._dp,-3.7406_dp)
!     factors(2,1) = ffac( "D", 2._dp, 6.671_dp)
!     factors(3,1) = ffac( "T", 3._dp, 4.792_dp)
!     ! Palladium
!     factors(0,46) = ffac( "Pd", 106.42_dp, 5.91_dp)
!     factors(1,46) = ffac( "Pd", 102._dp, 7.7_dp)
!     factors(2,46) = ffac( "Pd", 104._dp, 7.7_dp)
!     factors(3,46) = ffac( "Pd", 105._dp, 5.5_dp)
!     factors(4,46) = ffac( "Pd", 106._dp, 6.4_dp)
!     factors(5,46) = ffac( "Pd", 108._dp, 4.1_dp)
!     factors(6,46) = ffac( "Pd", 110._dp, 7.7_dp)
    !
    ! Selenium
    factors(0,34) = ffac("Se", 78.971_dp, 7.970_dp, 7.98_dp)
    factors(1,34) = ffac("Se", 74_dp,     0.8_dp,   0.1_dp)
    factors(2,34) = ffac("Se", 76_dp,     12.2_dp,  18.7_dp)
    factors(3,34) = ffac("Se", 77_dp,     8.25_dp,  8.6_dp)
    factors(4,34) = ffac("Se", 78_dp,     8.24_dp,  8.5_dp)
    factors(5,34) = ffac("Se", 80_dp,     7.48_dp,  7.03_dp)
    factors(6,34) = ffac("Se", 82_dp,     6.34_dp,  5.05_dp)
    !
    ! Bismuth
    factors(0,83) = ffac("Bi", 208.9804_dp, 8.532_dp, 9.148_dp)
    factors(1,83) = ffac("Bi", 209._dp,     8.532_dp, 9.148_dp)
    !
    ! Tellurium
    factors(0,52) = ffac( "Te", 127.6_dp, 5.80_dp, 4.23_dp)
    factors(1,52) = ffac( "Te", 120._dp, 5.30_dp, 3.5_dp)
    factors(2,52) = ffac( "Te", 122._dp, 3.80_dp, 1.8_dp)
    factors(3,52) = ffac( "Te", 123._dp, 0.00_dp, 0.52_dp) ! -CMPLX(0.05, 0.11)
    factors(4,52) = ffac( "Te", 124._dp, 7.96_dp, 8._dp)
    factors(5,52) = ffac( "Te", 125._dp, 5.02_dp, 3.17_dp)
    factors(6,52) = ffac( "Te", 126._dp, 5.56_dp, 3.88_dp)
    factors(7,52) = ffac( "Te", 128._dp, 5.89_dp, 4.36_dp)
    factors(8,52) = ffac( "Te", 130._dp, 6.02_dp, 4.55_dp)
    !
    ! Lead
    factors(0,82) = ffac( "Pb", 207.2_dp, 9.405_dp, 11.115_dp)
    factors(1,82) = ffac( "Pb", 204._dp, 9.90_dp, 12.3_dp )
    factors(2,82) = ffac( "Pb", 206._dp, 9.22_dp, 10.68_dp)
    factors(3,82) = ffac( "Pb", 207._dp, 9.28_dp, 10.82_dp)
    factors(4,82) = ffac( "Pb", 208._dp, 9.50_dp, 11.34_dp)
    
    IF(.not. present(amass)) THEN
      WRITE(*,*) "Looking for ", TRIM(aname)
      !
      ! Deuterium an Tritium are special:
      IF( aname == "D" ) THEN
        neutron_form_factor = FERMI_TO_BOHR*factors(2,1)%coh_b
        RETURN
      ELSEIF( aname == "T" ) THEN
        neutron_form_factor = FERMI_TO_BOHR*factors(3,1)%coh_b
        RETURN
      ENDIF
      ! Scan with only atomic symbol:
      DO i_type = 1,n_types
        IF(aname == factors(0,i_type)%aname) THEN
          neutron_form_factor = FERMI_TO_BOHR*factors(0,i_type)%coh_b
          RETURN
        ENDIF
      ENDDO
      !
    ELSE
      !
      WRITE(*,*) "Looking for ", TRIM(aname), "with mass", amass
      ! Scan with both symbol and mass:
      DO i_type = 1,n_types
      DO i_iso  = 1,n_iso
        IF( aname == factors(i_iso,i_type)%aname .and. &
           ABS(amass-factors(i_iso,i_type)%amass)<1.e-6_dp ) THEN
          neutron_form_factor = FERMI_TO_BOHR*factors(i_iso,i_type)%coh_b
          RETURN
        ENDIF
      ENDDO
      ENDDO
      !
    ENDIF
    !
    CALL errore("neutron_form_factor", "could not find '"//aname//"'",1)
    !
  END FUNCTION neutron_form_factor
  !
  ! Compute convoluting function
  SUBROUTINE neutron_function_convolution(w, xq, xg, elastic_peak_norm, &
               elastic_peak_width,T, S, fc2, &
               res_frac, res0_cmm1, ne, ee, spf)
    USE functions,      ONLY : f_bose, f_ngauss, f_lorentz
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE fc2_interpolate,     ONLY : freq_phq_safe
    IMPLICIT NONE
    REAL(DP) :: elastic_peak_norm
    REAL(DP) :: elastic_peak_width
    REAL(DP),INTENT(in) :: w, xq(3), xg(3), T
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(in) :: res_frac  ! as a fraction of the scattering energy
    REAL(DP),INTENT(in) :: res0_cmm1 ! resolution at zero energy, in cmm1
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ee(ne)
    REAL(DP),INTENT(inout) :: spf(ne,S%nat3)
    !
    INTEGER     :: nu, ie, je, je_range, je_range_dw, je_range_up
    !
    !REAL(DP) :: freq(S%nat3), xqq(3)
    !COMPLEX(DP) :: zz(S%nat3, S%nat3)
    !
    REAL(DP) :: G(S%nat3), F, Sq(ne,S%nat3)
    REAL(DP) :: sigma, sigma_p, sigma_m, norm, new_norm
    REAL(DP) :: dee(ne), peak0(ne)
    REAL(DP),PARAMETER :: FWHM_TO_SIGMA = 0.42466_dp !=( 2 *  sqrt(2._dp *log(2._dp)) )**-1
    !
    ! prepare phonon frequencies and patterns
!     WRITE(*,*)
!     WRITE(*,'(a,3f12.6)') "Convolution q-point: ",  xq
!     WRITE(*,'(a,3f12.6)') "            g-point: ",  xg

    !xqq = xq+xg
    !write(*,"(3f12.6)") xqq
    !
    !CALL freq_phq_safe(xqq, S, fc2, freq, zz)
    !
    ! Put spf in its normalized form
    ! (currently it is in 1/cmm1^2, we get it to 1/cmm1 and convert those to Ry)
    ! there is also a factor pi/2 that's lost somewhere in lw.x, fix it in the future
    DO nu=1,S%nat3
      spf(:,nu) = spf(:,nu)*ee*RY_TO_CMM1*2/pi
    ENDDO

    ! Compute the form factor in a more general form than really necessary
    ! for Simpson integration
    dee(1:ne-1) = ee(2:ne)-ee(1:ne-1)
    dee(ne) = dee(ne-1)
    DO nu=1,S%nat3
     CALL simpson(ne, spf(:,nu), dee, norm)
      !print*, "norm", nu, norm
      IF(ABS(norm-1._dp)>.1_dp) WRITE(*,*) "Strange norm", nu, norm
    ENDDO
    !
    !
    G = neutron_cross_section(xq, xg, S, fc2)
    !
    ! Doing the convolution with a gaussian in the least efficient way possible
    IF(res_frac>0._dp .or. res0_cmm1>0._dp )THEN
      Sq = 0._dp
      DO ie=1,ne
      !
      sigma = MAX(ABS(res_frac*ee(ie)), ABS(res0_cmm1)/RY_TO_CMM1)*FWHM_TO_SIGMA
      !WRITE(stdout,*) "Convolution with experimental width:", sigma, " cm^-1"
      je_range = 5*sigma/dee(ie)
      je_range_dw = MIN(je_range, ie-1)
      je_range_up = MIN(je_range, ne-ie)
      
        DO je = ie-je_range_dw, ie+je_range_up
          !
          DO nu = 1,S%nat3
            F = f_ngauss(ee(ie)-ee(je),sigma)
            !
            ! I use spf*e which seems to be the correctly normalizable form
  !         Sq(je,nu)=Sq(je,nu)+de*ee(ie)*spf(ie,nu)*F
            Sq(je,nu)=Sq(je,nu)+dee(ie)*spf(ie,nu)*F
            !
          ENDDO
        ENDDO
      ENDDO
    ELSE
      Sq = spf
    ENDIF


    !
    !
    ! Add the elastic peak at zero, note that this could be done before or after
    ! convoluting with the exp resolution, both approaches are legit, but this is more 
    ! flexible as it allows narrower elastic peaks
    !
    IF(elastic_peak_norm/=0._dp .and. elastic_peak_width /= 0)THEN
      peak0 = 0._dp
      DO ie = 1,ne
        !IF(ee(ie)/=0._dp) THEN
          !peak0(ie) = f_bose(ee(ie),T)!/ee(ie)
          peak0(ie) = f_lorentz(ee(ie), elastic_peak_width)
        !ENDIF
        !write(999, '(2e25.15)') ee(ie), peak0(ie)
      ENDDO
      !
      CALL simpson(ne, peak0, dee, norm)
      print*, "norm0", norm
      peak0 = elastic_peak_norm*peak0/norm
      
      DO nu = 1, S%nat3
        !spf(:,nu) = spf(:,nu) + peak0(:)/S%nat3
        sq(:,nu) = sq(:,nu) + peak0(:)
      ENDDO
!       DO nu=1,S%nat3
!         CALL simpson(ne, spf(:,nu), dee, norm)
!         print*, "norm2", nu, norm
!       ENDDO
    ENDIF    
    !     DO nu=1,S%nat3
!       CALL simpson(ne, RY_TO_CMM1*ee*spf(:,nu)*2/pi, dee, norm)
!       print*, "norm3", nu, norm
!     ENDDO
    !
    !
    ! Apply the structure factor
    DO nu = 1,S%nat3
      Sq(:,nu) = G(nu) * Sq(:,nu)
    ENDDO
    !

    
    ! renomalize to more or less the same as before, to make comparison easier
!     norm = SUM(Sq)
!     Sq = Sq/norm/pi
    !norm = MAXVAL(abs(Sq))      
    !Sq = Sq/norm
    spf = sq
    !
!     WRITE(999, '(2(a,3f12.6))') "# xq = ", xq, "      xg = ", xg
!     DO ie = 1, ne
!       WRITE(999,'(2f12.6,100e20.10e4)') w, ee(ie)*RY_TO_CMM1, SUM(Sq(ie,:)), Sq(ie,:)
!     ENDDO
!     WRITE(999,*)
    
    !
  END SUBROUTINE neutron_function_convolution 

  FUNCTION neutron_cross_section(xq, xg, S, fc2) RESULT(G)
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE fc2_interpolate,     ONLY : freq_phq_safe
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: xq(3), xg(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    !
    COMPLEX(DP) :: cscal, cscal2, csum
    INTEGER     :: nu, ia, idx, it
    !
    REAL(DP) :: freq(S%nat3), xqq(3), ff(S%ntyp)
    COMPLEX(DP) :: zz(S%nat3, S%nat3)
    !
    REAL(DP) :: G(S%nat3)
    !
    ! prepare phonon frequencies and patterns
    xqq = xq+xg
    !write(*,"(3f12.6)") xqq
    !
    DO it = 1, S%ntyp
      ff(it) = neutron_form_factor(S%atm(it))/S%celldm(1)
    ENDDO
    !
    CALL freq_phq_safe(xqq, S, fc2, freq, zz)
    !
    DO nu=1,S%nat3
      G(nu)=0._dp
      csum =(0._dp,0._dp)
      DO ia=1,S%nat
          idx=(ia-1)*3
          ! 2 \pi i q \dot \tau
          cscal  = CMPLX(0._dp, tpi*SUM(xqq*S%tau(:,ia)), kind=DP) 
          !print*, xqq*S%tau(:,ia)
          ! 2 \pi q \dot zz^{ia}_nu
          cscal2 = tpi* ( xqq(1)*zz(idx+1,nu) &
                         +xqq(2)*zz(idx+2,nu) &
                         +xqq(3)*zz(idx+3,nu) )
!           cscal2 = 1._dp
          ! \sum F_{ia} exp(\pi i q \dot \tau) 2 \pi q \dot zz^{ia}_nu / sqrt(M_{ia})
          csum = csum+ff(S%ityp(ia))*EXP(cscal)* cscal2 !/DSQRT(S%amass(S%ityp(ia)))
      ENDDO
      G(nu) = ABS(csum)**2
    ENDDO
    !
    !print*, G
  END FUNCTION
  
 END MODULE
 
 
 
