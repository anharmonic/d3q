!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Triple licenced under the CeCILL licence v 2.1:
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see:
!  <http://www.gnu.org/copyleft/gpl.txt>
!  and under the MIT licence see:
!  <https://opensource.org/licenses/mit-license.php>
!
! The isotope data comes from:
! http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some
! It is probably not copyrightable.
!
! The data has been formatted with two vim macros and a bit of patience.
!
! Atomic Weights and Isotopic Compositions for elements up to Z=99 (too much already!)
!
! NOTE: For Artificial Elements there is no natural isotope concentration! 
!       You can only select a specific isotope.
!
! NOTE: Deuterium and Tritium can be requested as "D" and "T". H will return 
!       the natural isotope concentration of Hydrogen.
!
! NOTE: Some radioactive isotopes have zero concentration
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE nist_isotopes_db
    USE kinds, ONLY : DP
#include "mpi_thermal.h"  
    
    PUBLIC :: element, search_nist, compute_gs
    ! >>>>>>>>
    PRIVATE !<<<<<<<<<
! <<^V^\\=========================================//-//-//========//O\\//
    TYPE isotope
      INTEGER :: anum  ! atomic number ...
      REAL(DP):: amass ! ... mass ...
      REAL(DP):: aconc ! ... and concentration of different isotopes

    END TYPE isotope
    INTERFACE isotope
      MODULE PROCEDURE new_isotope
    END INTERFACE isotope
    !
    TYPE element
      CHARACTER(LEN=2)     :: aname ! name of element
      INTEGER              :: z     !
      INTEGER              :: nisot = 0 ! number of isotopes
      TYPE(isotope),ALLOCATABLE :: isot(:)
      CONTAINS
        procedure :: assign_element
        generic   :: assignment(=) => assign_element
#ifdef __INTEL
        final     :: destroy_element
#endif        
    END TYPE
    INTERFACE element
      MODULE PROCEDURE new_element
    END INTERFACE element
  !
  CONTAINS
  !
  TYPE(element) FUNCTION new_element(aname,z,nisot) RESULT(self)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: aname
    INTEGER,INTENT(in)          :: z, nisot
    self%nisot = nisot
    self%aname = aname
    self%z     = z
    ALLOCATE(self%isot(self%nisot) )
    RETURN
  END FUNCTION new_element
  !
  SUBROUTINE destroy_element(self)
    IMPLICIT NONE
    TYPE(element),INTENT(inout) :: self
    IF(allocated(self%isot)) DEALLOCATE(self%isot)
    RETURN
  END SUBROUTINE destroy_element
  !
  SUBROUTINE assign_element(to, from)
    IMPLICIT NONE
    CLASS(element),INTENT(out) :: to
    CLASS(element),INTENT(in) :: from
    INTEGER :: i
    !
    to%aname = from%aname
    to%z     = from%z
    to%nisot = from%nisot
    ALLOCATE(to%isot(to%nisot))
    DO i = 1,to%nisot
      to%isot(i)%anum  = from%isot(i)%anum
      to%isot(i)%amass = from%isot(i)%amass
      to%isot(i)%aconc = from%isot(i)%aconc
    ENDDO
    RETURN
  END SUBROUTINE
  !
  TYPE(isotope) FUNCTION new_isotope(i,mass,conc) RESULT(self)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: i
    REAL(DP),INTENT(in):: mass, conc ! get them in single precision, for simplicity
                                 ! there are a just couple concentrations given with more than 7 digits
    self%anum = i
    self%amass = DBLE(mass)
    self%aconc = DBLE(conc)
    RETURN
  END FUNCTION new_isotope
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the standard atomic mass (gm) and its standard deviation (gs)
  SUBROUTINE compute_gs(gm, gs, aname, anum, nisot, amass, aconc)
    IMPLICIT NONE
    REAL(DP),INTENT(out) :: gm, gs
    !
    CHARACTER(len=*),INTENT(in) :: aname
    INTEGER,INTENT(in)    :: anum  ! atomic number: = 0 take the natural isotope occurence
                                   !                > 0 take the mass of that specific isotope
    INTEGER,INTENT(in) :: nisot    ! number of isotopes provided in input (set = 0 if automatic, use natural values)
    
    REAL(DP),INTENT(inout),ALLOCATABLE,OPTIONAL:: amass(:) ! isotope atomic mass provided in input
    REAL(DP),INTENT(inout),ALLOCATABLE,OPTIONAL:: aconc(:) ! isotope concentration provided in input
    !
    ! method 1: 
    !
    ! Set anum=0 and nisot=0 in input, receive the natural occurrence of isotopes in output
    ! amass and aconc, if present, will be allocated an filled ccoredingly.
    !
    ! method 2:
    ! Set anum>0, the mass for that specific isotope will be returned, gs will be zero, 
    ! amass and aconc, if present, will be allocated to size 1 and contain gm and 1.0 respecively
    !
    ! method 3:
    ! Set nisot > 0, in this case it is assumed that amass and aconc are already allocated 
    ! and contain the masses and concentration of the isotopes. gs and gm are computed and returned
    !
    ! This subroutine is a bit ugly, but it works, and it does not have to be nice as it is
    ! not really critical or anything. It is also long to compile.
    !
    TYPE(element) :: elem
    !
    INTEGER :: i
    REAL(DP),ALLOCATABLE :: amass_(:) ! isotope atomic mass
    REAL(DP),ALLOCATABLE :: aconc_(:) ! isotope concentration
    !
    NISOT_IF : &
    IF(nisot>0)THEN
      IF(.not.PRESENT(amass) .or. .not.PRESENT(aconc))&
        CALL errore("compute_gs", "You need to provide amass and aconc for nisot>0",1)
      IF(.not.ALLOCATED(amass) .or. .not.allocated(aconc))&
        CALL errore("compute_gs", "You need to provide amass and aconc for nisot>0",2)
      IF( ABS(SUM(aconc)-1._dp)>1.e-4_dp) THEN
        !CALL errore("compute_gs", "the sum of isotope concentrations is not 100%", 3)
        ioWRITE(*,'(5x,3a)') "WARNING: isotopes of ", TRIM(aname), " did not sum to 100%: renormalized"
        aconc = aconc/SUM(aconc)
      ENDIF
      !
      ioWRITE(*,'(5x,2a)') "Element from input: ", aname
      ioWRITE(*,'(8x,a,i3,a)') "number of isotopes:", nisot, ". Mass, concentration:"
      DO i = 1, nisot
        ioWRITE(*,'(10x,2f12.6)') amass(i), aconc(i)
      ENDDO
        
      gm = SUM( aconc*amass )
      gs = SUM( aconc*(amass-gm)**2 ) / gm**2
      
    ELSE NISOT_IF
      
      elem = search_nist(aname)
      ioWRITE(*,'(5x,2a)') "Found in NIST: ", elem%aname
    
      IF(anum>0)THEN
        !
        gs = 0._dp
        gm = -1._dp
        !
        DO i = 1, elem%nisot
          IF(elem%isot(i)%anum == anum) gm = elem%isot(i)%amass
        ENDDO
        IF(gm<0._dp) CALL errore("compute_gs", "Requested isotope not found", 5)
        !
      ELSE
        !
        ALLOCATE(amass_(elem%nisot))
        ALLOCATE(aconc_(elem%nisot))
        DO i = 1, elem%nisot
          amass_(i) = elem%isot(i)%amass
          aconc_(i) = elem%isot(i)%aconc
        ENDDO
        IF( ABS(SUM(aconc_)-1._dp)>1.e-4_dp) &
               CALL errore("compute_gs", "You probably chose an artificial element: "&
                           //"natural concentration of isotopes is undefined!", 10)
        !
!         IF(present(amass)) amass = amass_ ! should be allocated automatically, fortran 2003
!         IF(present(aconc)) aconc = aconc_ ! idem
        !
        gm = SUM( aconc_*amass_ )
        gs = SUM( aconc_*(amass_-gm)**2 ) / gm**2
        !
        ioWRITE(*,'(8x,a,i3,a)') "number of natural isotopes:", elem%nisot, ". Mass, concentration:"
        DO i = 1, elem%nisot
          ioWRITE(*,'(10x,2f12.6)') amass_(i), aconc_(i)
        ENDDO
      ENDIF
    ENDIF &
    NISOT_IF
    !
    ioWRITE(*,'(8x,a,2f12.6)')  "average mass, relative standard dev.:", gm, gs

  END SUBROUTINE compute_gs
  ! \/o\________\\\_________________________________________/^>
  
  TYPE(element) FUNCTION search_nist(name) RESULT(self)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: name
    TYPE(element) :: H, D, T, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, &
    Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, &
    Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, &
    Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, &
    Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es
    !
    !_________________________________________________________________________
      H = new_element("H", 1, 3)
        H%isot(1) = new_isotope(1, 1.00782503207_dp , 0.999885_dp )
        H%isot(2) = new_isotope(2, 2.0141017778_dp ,  0.000115_dp )
        H%isot(3) = new_isotope(3, 3.016049277_dp ,   0.000000_dp )
      D = new_element("D", 1, 1)
        D%isot(1) = new_isotope(2, 2.0141017778_dp , 1.0_dp )
      T = new_element("T", 1, 1)
        T%isot(1) = new_isotope(3, 3.0160492777_dp , 1.0_dp )
    !_________________________________________________________________________
      He = new_element("He" , 2 , 2)
        He%isot(1) = new_isotope(3, 3.0160293191_dp ,  0.00000134_dp )
        He%isot(2) = new_isotope(4, 4.00260325415_dp , 0.99999866_dp )
    !_________________________________________________________________________
      Li= new_element("Li", 3, 2)
        Li%isot(1) = new_isotope(6, 6.015122795_dp , 0.0759_dp )
        Li%isot(2) = new_isotope(7, 7.01600455_dp ,  0.9241_dp )
    !_________________________________________________________________________
      Be= new_element("Be", 4, 1)
        Be%isot(1) = new_isotope(9, 9.0121822_dp , 1.0000_dp )
    !_________________________________________________________________________
      B= new_element("B", 5, 2)
        B%isot(1) = new_isotope(10, 10.0129370_dp , 0.199_dp )
        B%isot(2) = new_isotope(11, 11.0093054_dp , 0.801_dp )
    !_________________________________________________________________________
      C= new_element("C", 6, 3)
        C%isot(1) = new_isotope(12, 12.0000000_dp ,    0.9893_dp )
        C%isot(2) = new_isotope(13, 13.0033548378_dp , 0.0107_dp )
        C%isot(3) = new_isotope(14, 14.003241989_dp ,  0.0_dp )
    !_________________________________________________________________________
      N= new_element("N", 7, 2)
        N%isot(1) = new_isotope(14, 14.0030740048_dp , 0.99636_dp )
        N%isot(2) = new_isotope(15, 15.0001088982_dp , 0.00364_dp )
    !_________________________________________________________________________
      O= new_element("O", 8, 3)
        O%isot(1) = new_isotope(16, 15.99491461956_dp , 0.99757_dp )
        O%isot(2) = new_isotope(17, 16.99913170_dp ,    0.00038_dp )
        O%isot(3) = new_isotope(18, 17.9991610_dp ,     0.00205_dp )
    !_________________________________________________________________________
      F= new_element("F", 9, 1)
        F%isot(1) = new_isotope(19, 18.99840322_dp , 1.0000_dp )
    !_________________________________________________________________________
    Ne= new_element("Ne", 10, 3)
        Ne%isot(1) = new_isotope(20, 19.9924401754_dp , 0.9048_dp )
        Ne%isot(2) = new_isotope(21, 20.99384668_dp ,   0.0027_dp )
        Ne%isot(3) = new_isotope(22, 21.991385114_dp ,  0.0925_dp )
    !_________________________________________________________________________
    Na= new_element("Na", 11, 1)
        Na%isot(1) = new_isotope(23, 22.9897692809_dp , 1.0000_dp )
    !_________________________________________________________________________
    Mg= new_element("Mg", 12, 3)
        Mg%isot(1) = new_isotope(24, 23.985041700_dp , 0.7899_dp )
        Mg%isot(2) = new_isotope(25, 24.98583692_dp ,  0.1000_dp )
        Mg%isot(3) = new_isotope(26, 25.982592929_dp , 0.1101_dp )
    !_________________________________________________________________________
    Al= new_element("Al", 13, 1)
        Al%isot(1) = new_isotope(27, 26.98153863_dp , 1.0000_dp )
    !_________________________________________________________________________
    Si= new_element("Si", 14, 3)
        Si%isot(1) = new_isotope(28, 27.9769265325_dp , 0.92223_dp )
        Si%isot(2) = new_isotope(29, 28.976494700_dp ,  0.04685_dp )
        Si%isot(3) = new_isotope(30, 29.97377017_dp ,   0.03092_dp )
    !_________________________________________________________________________
    P= new_element("P", 15, 1)
        P%isot(1) = new_isotope(31, 30.97376163_dp , 1.0000_dp )
    !_________________________________________________________________________
    S= new_element("S", 16, 4)
        S%isot(1) = new_isotope(32, 31.97207100_dp , 0.9499_dp )
        S%isot(2) = new_isotope(33, 32.97145876_dp , 0.0075_dp )
        S%isot(3) = new_isotope(34, 33.96786690_dp , 0.0425_dp )
        S%isot(4) = new_isotope(36, 35.96708076_dp , 0.0001_dp )
    !_________________________________________________________________________
    Cl= new_element("Cl", 17, 2)
        Cl%isot(1) = new_isotope(35, 34.96885268_dp , 0.7576_dp )
        Cl%isot(2) = new_isotope(37, 36.96590259_dp , 0.2424_dp )
    !_________________________________________________________________________
    Ar= new_element("Ar", 18, 3)
        Ar%isot(1) = new_isotope(36, 35.967545106_dp ,  0.003365_dp )
        Ar%isot(2) = new_isotope(38, 37.9627324_dp ,    0.000632_dp )
        Ar%isot(3) = new_isotope(40, 39.9623831225_dp , 0.996003_dp )
    !_________________________________________________________________________
    K= new_element("K", 19, 3)
        K%isot(1) = new_isotope(39, 38.96370668_dp , 0.932581_dp )
        K%isot(2) = new_isotope(40, 39.96399848_dp , 0.000117_dp )
        K%isot(3) = new_isotope(41, 40.96182576_dp , 0.067302_dp )
    !_________________________________________________________________________
    Ca= new_element("Ca", 20, 6)
        Ca%isot(1) = new_isotope(40, 39.96259098_dp , 0.96941_dp )
        Ca%isot(2) = new_isotope(42, 41.95861801_dp , 0.00647_dp )
        Ca%isot(3) = new_isotope(43, 42.9587666_dp ,  0.00135_dp )
        Ca%isot(4) = new_isotope(44, 43.9554818_dp ,  0.02086_dp )
        Ca%isot(5) = new_isotope(46, 45.9536926_dp ,  0.00004_dp )
        Ca%isot(6) = new_isotope(48, 47.952534_dp ,   0.00187_dp )
    !_________________________________________________________________________
    Sc= new_element("Sc", 21, 1)
        Sc%isot(1) = new_isotope(45, 44.9559119_dp , 1.0000_dp )
    !_________________________________________________________________________
    Ti= new_element("Ti", 22, 5)
        Ti%isot(1) = new_isotope(46, 45.9526316_dp , 0.0825_dp )
        Ti%isot(2) = new_isotope(47, 46.9517631_dp , 0.0744_dp )
        Ti%isot(3) = new_isotope(48, 47.9479463_dp , 0.7372_dp )
        Ti%isot(4) = new_isotope(49, 48.9478700_dp , 0.0541_dp )
        Ti%isot(5) = new_isotope(50, 49.9447912_dp , 0.0518_dp )
    !_________________________________________________________________________
    V= new_element("V", 23, 2)
        V%isot(1) = new_isotope(50, 49.9471585_dp , 0.00250_dp )
        V%isot(2) = new_isotope(51, 50.9439595_dp , 0.99750_dp )
    !_________________________________________________________________________
    Cr= new_element("Cr", 24, 4)
        Cr%isot(1) = new_isotope(50, 49.9460442_dp , 0.04345_dp )
        Cr%isot(2) = new_isotope(52, 51.9405075_dp , 0.83789_dp )
        Cr%isot(3) = new_isotope(53, 52.9406494_dp , 0.09501_dp )
        Cr%isot(4) = new_isotope(54, 53.9388804_dp , 0.02365_dp )
    !_________________________________________________________________________
    Mn= new_element("Mn", 25, 1)
        Mn%isot(1) = new_isotope(55, 54.9380451_dp , 1.0000_dp )
    !_________________________________________________________________________
    Fe= new_element("Fe", 26, 4)
        Fe%isot(1) = new_isotope(54, 53.9396105_dp , 0.05845_dp )
        Fe%isot(2) = new_isotope(56, 55.9349375_dp , 0.91754_dp )
        Fe%isot(3) = new_isotope(57, 56.9353940_dp , 0.02119_dp )
        Fe%isot(4) = new_isotope(58, 57.9332756_dp , 0.00282_dp )
    !_________________________________________________________________________
    Co= new_element("Co", 27, 1)
        Co%isot(1) = new_isotope(59, 58.9331950_dp , 1.0000_dp )
    !_________________________________________________________________________
    Ni= new_element("Ni", 28, 5)
        Ni%isot(1) = new_isotope(58, 57.9353429_dp , 0.680769_dp )
        Ni%isot(2) = new_isotope(60, 59.9307864_dp , 0.262231_dp )
        Ni%isot(3) = new_isotope(61, 60.9310560_dp , 0.011399_dp )
        Ni%isot(4) = new_isotope(62, 61.9283451_dp , 0.036345_dp )
        Ni%isot(5) = new_isotope(64, 63.9279660_dp , 0.009256_dp )
    !_________________________________________________________________________
    Cu= new_element("Cu", 29, 2)
        Cu%isot(1) = new_isotope(63, 62.9295975_dp , 0.6915_dp )
        Cu%isot(2) = new_isotope(65, 64.9277895_dp , 0.3085_dp )
    !_________________________________________________________________________
    Zn= new_element("Zn", 30, 5)
        Zn%isot(1) = new_isotope(64, 63.9291422_dp , 0.48268_dp )
        Zn%isot(2) = new_isotope(66, 65.9260334_dp , 0.27975_dp )
        Zn%isot(3) = new_isotope(67, 66.9271273_dp , 0.04102_dp )
        Zn%isot(4) = new_isotope(68, 67.9248442_dp , 0.19024_dp )
        Zn%isot(5) = new_isotope(70, 69.9253193_dp , 0.00631_dp )
    !_________________________________________________________________________
    Ga= new_element("Ga", 31, 2)
        Ga%isot(1) = new_isotope(69, 68.9255736_dp , 0.60108_dp )
        Ga%isot(2) = new_isotope(71, 70.9247013_dp , 0.39892_dp )
    !_________________________________________________________________________
    Ge= new_element("Ge", 32, 5)
        Ge%isot(1) = new_isotope(70, 69.9242474_dp , 0.2038_dp )
        Ge%isot(2) = new_isotope(72, 71.9220758_dp , 0.2731_dp )
        Ge%isot(3) = new_isotope(73, 72.9234589_dp , 0.0776_dp )
        Ge%isot(4) = new_isotope(74, 73.9211778_dp , 0.3672_dp )
        Ge%isot(5) = new_isotope(76, 75.9214026_dp , 0.0783_dp )
    !_________________________________________________________________________
    As= new_element("As", 33, 1)
        As%isot(1) = new_isotope(75, 74.9215965_dp , 1.0000_dp )
    !_________________________________________________________________________
    Se= new_element("Se", 34, 6)
        Se%isot(1) = new_isotope(74, 73.9224764_dp , 0.0089_dp )
        Se%isot(2) = new_isotope(76, 75.9192136_dp , 0.0937_dp )
        Se%isot(3) = new_isotope(77, 76.9199140_dp , 0.0763_dp )
        Se%isot(4) = new_isotope(78, 77.9173091_dp , 0.2377_dp )
        Se%isot(5) = new_isotope(80, 79.9165213_dp , 0.4961_dp )
        Se%isot(6) = new_isotope(82, 81.9166994_dp , 0.0873_dp )
    !_________________________________________________________________________
    Br= new_element("Br", 35, 2)
        Br%isot(1) = new_isotope(79, 78.9183371_dp , 0.5069_dp )
        Br%isot(2) = new_isotope(81, 80.9162906_dp , 0.4931_dp )
    !_________________________________________________________________________
    Kr= new_element("Kr", 36, 6)
        Kr%isot(1) = new_isotope(78, 77.9203648_dp ,  0.00355_dp )
        Kr%isot(2) = new_isotope(80, 79.9163790_dp ,  0.02286_dp )
        Kr%isot(3) = new_isotope(82, 81.9134836_dp ,  0.11593_dp )
        Kr%isot(4) = new_isotope(83, 82.914136_dp ,   0.11500_dp )
        Kr%isot(5) = new_isotope(84, 83.911507_dp ,   0.56987_dp )
        Kr%isot(6) = new_isotope(86, 85.91061073_dp , 0.17279_dp )
    !_________________________________________________________________________
    Rb= new_element("Rb", 37, 2)
        Rb%isot(1) = new_isotope(85, 84.911789738_dp , 0.7217_dp )
        Rb%isot(2) = new_isotope(87, 86.909180527_dp , 0.2783_dp )
    !_________________________________________________________________________
    Sr= new_element("Sr", 38, 4)
        Sr%isot(1) = new_isotope(84, 83.913425_dp ,  0.0056_dp )
        Sr%isot(2) = new_isotope(86, 85.9092602_dp , 0.0986_dp )
        Sr%isot(3) = new_isotope(87, 86.9088771_dp , 0.0700_dp )
        Sr%isot(4) = new_isotope(88, 87.9056121_dp , 0.8258_dp )
    !_________________________________________________________________________
    Y= new_element("Y", 39, 1)
        Y%isot(1) = new_isotope(89, 88.9058483_dp , 1.0000_dp )
    !_________________________________________________________________________
    Zr= new_element("Zr", 40, 5)
        Zr%isot(1) = new_isotope(90, 89.9047044_dp , 0.5145_dp )
        Zr%isot(2) = new_isotope(91, 90.9056458_dp , 0.1122_dp )
        Zr%isot(3) = new_isotope(92, 91.9050408_dp , 0.1715_dp )
        Zr%isot(4) = new_isotope(94, 93.9063152_dp , 0.1738_dp )
        Zr%isot(5) = new_isotope(96, 95.9082734_dp , 0.0280_dp )
    !_________________________________________________________________________
    Nb= new_element("Nb", 41, 1)
        Nb%isot(1) = new_isotope(93, 92.9063781_dp , 1.0000_dp )
    !_________________________________________________________________________
    Mo= new_element("Mo", 42, 7)
        Mo%isot(1) = new_isotope(92, 91.906811_dp ,  0.1477_dp )
        Mo%isot(2) = new_isotope(94, 93.9050883_dp , 0.0923_dp )
        Mo%isot(3) = new_isotope(95, 94.9058421_dp , 0.1590_dp )
        Mo%isot(4) = new_isotope(96, 95.9046795_dp , 0.1668_dp )
        Mo%isot(5) = new_isotope(97, 96.9060215_dp , 0.0956_dp )
        Mo%isot(6) = new_isotope(98, 97.9054082_dp , 0.2419_dp )
        Mo%isot(7) = new_isotope(100, 99.907477_dp , 0.0967_dp )
    !_________________________________________________________________________
    Tc= new_element("Tc", 43, 3)
        Tc%isot(1) = new_isotope(97, 96.906365_dp ,  0.0_dp )
        Tc%isot(2) = new_isotope(98, 97.907216_dp ,  0.0_dp )
        Tc%isot(3) = new_isotope(99, 98.9062547_dp , 0.0_dp )
    !_________________________________________________________________________
    Ru= new_element("Ru", 44, 7)
        Ru%isot(1) = new_isotope(96, 95.907598_dp ,    0.0554_dp )
        Ru%isot(2) = new_isotope(98, 97.905287_dp ,    0.0187_dp )
        Ru%isot(3) = new_isotope(99, 98.9059393_dp ,   0.1276_dp )
        Ru%isot(4) = new_isotope(100, 99.9042195_dp ,  0.1260_dp )
        Ru%isot(5) = new_isotope(101, 100.9055821_dp , 0.1706_dp )
        Ru%isot(6) = new_isotope(102, 101.9043493_dp , 0.3155_dp )
        Ru%isot(7) = new_isotope(104, 103.905433_dp ,  0.1862_dp )
    !_________________________________________________________________________
    Rh= new_element("Rh", 45, 1)
        Rh%isot(1) = new_isotope(103, 102.905504_dp , 1.0000_dp )
    !_________________________________________________________________________
    Pd= new_element("Pd", 46, 6)
        Pd%isot(1) = new_isotope(102, 101.905609_dp , 0.0102_dp )
        Pd%isot(2) = new_isotope(104, 103.904036_dp , 0.1114_dp )
        Pd%isot(3) = new_isotope(105, 104.905085_dp , 0.2233_dp )
        Pd%isot(4) = new_isotope(106, 105.903486_dp , 0.2733_dp )
        Pd%isot(5) = new_isotope(108, 107.903892_dp , 0.2646_dp )
        Pd%isot(6) = new_isotope(110, 109.905153_dp , 0.1172_dp )
    !_________________________________________________________________________
    Ag= new_element("Ag", 47, 2)
        Ag%isot(1) = new_isotope(107, 106.905097_dp , 0.51839_dp )
        Ag%isot(2) = new_isotope(109, 108.904752_dp , 0.48161_dp )
    !_________________________________________________________________________
    Cd= new_element("Cd", 48, 8)
        Cd%isot(1) = new_isotope(106, 105.906459_dp ,  0.0125_dp )
        Cd%isot(2) = new_isotope(108, 107.904184_dp ,  0.0089_dp )
        Cd%isot(3) = new_isotope(110, 109.9030021_dp , 0.1249_dp )
        Cd%isot(4) = new_isotope(111, 110.9041781_dp , 0.1280_dp )
        Cd%isot(5) = new_isotope(112, 111.9027578_dp , 0.2413_dp )
        Cd%isot(6) = new_isotope(113, 112.9044017_dp , 0.1222_dp )
        Cd%isot(7) = new_isotope(114, 113.9033585_dp , 0.2873_dp )
        Cd%isot(8) = new_isotope(116, 115.904756_dp ,  0.0749_dp )
    !_________________________________________________________________________
    In= new_element("In", 49, 2)
        In%isot(1) = new_isotope(113, 112.904058_dp , 0.0429_dp )
        In%isot(2) = new_isotope(115, 114.903878_dp , 0.9571_dp )
    !_________________________________________________________________________
    Sn= new_element("Sn", 50, 10)
        Sn%isot(1)  = new_isotope(112, 111.904818_dp ,  0.0097_dp )
        Sn%isot(2)  = new_isotope(114, 113.902779_dp ,  0.0066_dp )
        Sn%isot(3)  = new_isotope(115, 114.903342_dp ,  0.0034_dp )
        Sn%isot(4)  = new_isotope(116, 115.901741_dp ,  0.1454_dp )
        Sn%isot(5)  = new_isotope(117, 116.902952_dp ,  0.0768_dp )
        Sn%isot(6)  = new_isotope(118, 117.901603_dp ,  0.2422_dp )
        Sn%isot(7)  = new_isotope(119, 118.903308_dp ,  0.0859_dp )
        Sn%isot(8)  = new_isotope(120, 119.9021947_dp , 0.3258_dp )
        Sn%isot(9)  = new_isotope(122, 121.9034390_dp , 0.0463_dp )
        Sn%isot(10) = new_isotope(124, 123.9052739_dp , 0.0579_dp )
    !_________________________________________________________________________
    Sb= new_element("Sb", 51, 2)
        Sb%isot(1) = new_isotope(121, 120.9038157_dp , 0.5721_dp )
        Sb%isot(2) = new_isotope(123, 122.9042140_dp , 0.4279_dp )
    !_________________________________________________________________________
    Te= new_element("Te", 52, 8)
        Te%isot(1) = new_isotope(120, 119.904020_dp ,  0.0009_dp )
        Te%isot(2) = new_isotope(122, 121.9030439_dp , 0.0255_dp )
        Te%isot(3) = new_isotope(123, 122.9042700_dp , 0.0089_dp )
        Te%isot(4) = new_isotope(124, 123.9028179_dp , 0.0474_dp )
        Te%isot(5) = new_isotope(125, 124.9044307_dp , 0.0707_dp )
        Te%isot(6) = new_isotope(126, 125.9033117_dp , 0.1884_dp )
        Te%isot(7) = new_isotope(128, 127.9044631_dp , 0.3174_dp )
        Te%isot(8) = new_isotope(130, 129.9062244_dp , 0.3408_dp )
    !_________________________________________________________________________
    I= new_element("I", 53, 1)
        I%isot(1) = new_isotope(127, 126.904473_dp , 1.0000_dp )
    !_________________________________________________________________________
    Xe= new_element("Xe", 54, 9)
        Xe%isot(1) = new_isotope(124, 123.9058930_dp , 0.000952_dp )
        Xe%isot(2) = new_isotope(126, 125.904274_dp ,  0.000890_dp )
        Xe%isot(3) = new_isotope(128, 127.9035313_dp , 0.019102_dp )
        Xe%isot(4) = new_isotope(129, 128.9047794_dp , 0.264006_dp )
        Xe%isot(5) = new_isotope(130, 129.9035080_dp , 0.040710_dp )
        Xe%isot(6) = new_isotope(131, 130.9050824_dp , 0.212324_dp )
        Xe%isot(7) = new_isotope(132, 131.9041535_dp , 0.269086_dp )
        Xe%isot(8) = new_isotope(134, 133.9053945_dp , 0.104357_dp )
        Xe%isot(9) = new_isotope(136, 135.907219_dp ,  0.088573_dp )
    !_________________________________________________________________________
    Cs= new_element("Cs", 55, 1)
        Cs%isot(1) = new_isotope(133, 132.905451933_dp , 1.0000_dp )
    !_________________________________________________________________________
    Ba= new_element("Ba", 56, 7)
        Ba%isot(1) = new_isotope(130, 129.9063208_dp , 0.00106_dp )
        Ba%isot(2) = new_isotope(132, 131.9050613_dp , 0.00101_dp )
        Ba%isot(3) = new_isotope(134, 133.9045084_dp , 0.02417_dp )
        Ba%isot(4) = new_isotope(135, 134.9056886_dp , 0.06592_dp )
        Ba%isot(5) = new_isotope(136, 135.9045759_dp , 0.07854_dp )
        Ba%isot(6) = new_isotope(137, 136.9058274_dp , 0.11232_dp )
        Ba%isot(7) = new_isotope(138, 137.9052472_dp , 0.71698_dp )
    !_________________________________________________________________________
    La= new_element("La", 57, 2)
        La%isot(1) = new_isotope(138, 137.907112_dp , 0.00090_dp )
        La%isot(2) = new_isotope(139, 138.9063533_dp , 0.99910_dp )
    !_________________________________________________________________________
    Ce= new_element("Ce", 58, 4)
        Ce%isot(1) = new_isotope(136, 135.907172_dp ,  0.00185_dp )
        Ce%isot(2) = new_isotope(138, 137.905991_dp ,  0.00251_dp )
        Ce%isot(3) = new_isotope(140, 139.9054387_dp , 0.88450_dp )
        Ce%isot(4) = new_isotope(142, 141.909244_dp ,  0.11114_dp )
    !_________________________________________________________________________
    Pr= new_element("Pr", 59, 1)
        Pr%isot(1) = new_isotope(141, 140.9076528_dp , 1.0000_dp )
    !_________________________________________________________________________
    Nd= new_element("Nd", 60, 7)
        Nd%isot(1) = new_isotope(142, 141.9077233_dp , 0.272_dp )
        Nd%isot(2) = new_isotope(143, 142.9098143_dp , 0.122_dp )
        Nd%isot(3) = new_isotope(144, 143.9100873_dp , 0.238_dp )
        Nd%isot(4) = new_isotope(145, 144.9125736_dp , 0.083_dp )
        Nd%isot(5) = new_isotope(146, 145.9131169_dp , 0.172_dp )
        Nd%isot(6) = new_isotope(148, 147.916893_dp , 0.057_dp )
        Nd%isot(7) = new_isotope(150, 149.920891_dp , 0.056_dp )
    !_________________________________________________________________________
    Pm= new_element("Pm", 61, 2)
        Pm%isot(1) = new_isotope(145, 144.912749_dp ,  0.0_dp )
        Pm%isot(2) = new_isotope(147, 146.9151385_dp , 0.0_dp )
    !_________________________________________________________________________
    Sm= new_element("Sm", 62, 7)
        Sm%isot(1) = new_isotope(144, 143.911999_dp ,  0.0307_dp )
        Sm%isot(2) = new_isotope(147, 146.9148979_dp , 0.1499_dp )
        Sm%isot(3) = new_isotope(148, 147.9148227_dp , 0.1124_dp )
        Sm%isot(4) = new_isotope(149, 148.9171847_dp , 0.1382_dp )
        Sm%isot(5) = new_isotope(150, 149.9172755_dp , 0.0738_dp )
        Sm%isot(6) = new_isotope(152, 151.9197324_dp , 0.2675_dp )
        Sm%isot(7) = new_isotope(154, 153.9222093_dp , 0.2275_dp )
    !_________________________________________________________________________
    Eu= new_element("Eu", 63, 2)
        Eu%isot(1) = new_isotope(151, 150.9198502_dp , 0.4781_dp )
        Eu%isot(2) = new_isotope(153, 152.9212303_dp , 0.5219_dp )
    !_________________________________________________________________________
    Gd= new_element("Gd", 64, 7)
        Gd%isot(1) = new_isotope(152, 151.9197910_dp , 0.0020_dp )
        Gd%isot(2) = new_isotope(154, 153.9208656_dp , 0.0218_dp )
        Gd%isot(3) = new_isotope(155, 154.9226220_dp , 0.1480_dp )
        Gd%isot(4) = new_isotope(156, 155.9221227_dp , 0.2047_dp )
        Gd%isot(5) = new_isotope(157, 156.9239601_dp , 0.1565_dp )
        Gd%isot(6) = new_isotope(158, 157.9241039_dp , 0.2484_dp )
        Gd%isot(7) = new_isotope(160, 159.9270541_dp , 0.2186_dp )
    !_________________________________________________________________________
    Tb= new_element("Tb", 65, 1)
        Tb%isot(1) = new_isotope(159, 158.9253468_dp , 1.0000_dp )
    !_________________________________________________________________________
    Dy= new_element("Dy", 66, 7)
        Dy%isot(1) = new_isotope(156, 155.924283_dp ,  0.00056_dp )
        Dy%isot(2) = new_isotope(158, 157.924409_dp ,  0.00095_dp )
        Dy%isot(3) = new_isotope(160, 159.9251975_dp , 0.02329_dp )
        Dy%isot(4) = new_isotope(161, 160.9269334_dp , 0.18889_dp )
        Dy%isot(5) = new_isotope(162, 161.9267984_dp , 0.25475_dp )
        Dy%isot(6) = new_isotope(163, 162.9287312_dp , 0.24896_dp )
        Dy%isot(7) = new_isotope(164, 163.9291748_dp , 0.28260_dp )
    !_________________________________________________________________________
    Ho= new_element("Ho", 67, 1)
        Ho%isot(1) = new_isotope(165, 164.9303221_dp , 1.0000_dp )
    !_________________________________________________________________________
    Er= new_element("Er", 68, 6)
        Er%isot(1) = new_isotope(162, 161.928778_dp ,  0.00139_dp )
        Er%isot(2) = new_isotope(164, 163.929200_dp ,  0.01601_dp )
        Er%isot(3) = new_isotope(166, 165.9302931_dp , 0.33503_dp )
        Er%isot(4) = new_isotope(167, 166.9320482_dp , 0.22869_dp )
        Er%isot(5) = new_isotope(168, 167.9323702_dp , 0.26978_dp )
        Er%isot(6) = new_isotope(170, 169.9354643_dp , 0.14910_dp )
    !_________________________________________________________________________
    Tm= new_element("Tm", 69, 1)
        Tm%isot(1) = new_isotope(169, 168.9342133_dp , 1.0000_dp )
    !_________________________________________________________________________
    Yb= new_element("Yb", 70, 7)
        Yb%isot(1) = new_isotope(168, 167.933897_dp ,  0.0013_dp )
        Yb%isot(2) = new_isotope(170, 169.9347618_dp , 0.0304_dp )
        Yb%isot(3) = new_isotope(171, 170.9363258_dp , 0.1428_dp )
        Yb%isot(4) = new_isotope(172, 171.9363815_dp , 0.2183_dp )
        Yb%isot(5) = new_isotope(173, 172.9382108_dp , 0.1613_dp )
        Yb%isot(6) = new_isotope(174, 173.9388621_dp , 0.3183_dp )
        Yb%isot(7) = new_isotope(176, 175.9425717_dp , 0.1276_dp )
    !_________________________________________________________________________
    Lu= new_element("Lu", 71, 2)
        Lu%isot(1) = new_isotope(175, 174.9407718_dp , 0.9741_dp )
        Lu%isot(2) = new_isotope(176, 175.9426863_dp , 0.0259_dp )
    !_________________________________________________________________________
    Hf= new_element("Hf", 72, 6)
        Hf%isot(1) = new_isotope(174, 173.940046_dp ,  0.0016_dp )
        Hf%isot(2) = new_isotope(176, 175.9414086_dp , 0.0526_dp )
        Hf%isot(3) = new_isotope(177, 176.9432207_dp , 0.1860_dp )
        Hf%isot(4) = new_isotope(178, 177.9436988_dp , 0.2728_dp )
        Hf%isot(5) = new_isotope(179, 178.9458161_dp , 0.1362_dp )
        Hf%isot(6) = new_isotope(180, 179.9465500_dp , 0.3508_dp )
    !_________________________________________________________________________
    Ta= new_element("Ta", 73, 2)
        Ta%isot(1) = new_isotope(180, 179.9474648_dp , 0.00012_dp )
        Ta%isot(2) = new_isotope(181, 180.9479958_dp , 0.99988_dp )
    !_________________________________________________________________________
    W= new_element("W", 74, 5)
        W%isot(1) = new_isotope(180, 179.946704_dp ,  0.0012_dp )
        W%isot(2) = new_isotope(182, 181.9482042_dp , 0.2650_dp )
        W%isot(3) = new_isotope(183, 182.9502230_dp , 0.1431_dp )
        W%isot(4) = new_isotope(184, 183.9509312_dp , 0.3064_dp )
        W%isot(5) = new_isotope(186, 185.9543641_dp , 0.2843_dp )
    !_________________________________________________________________________
    Re= new_element("Re", 75, 2)
        Re%isot(1) = new_isotope(185, 184.9529550_dp , 0.3740_dp )
        Re%isot(2) = new_isotope(187, 186.9557531_dp , 0.6260_dp )
    !_________________________________________________________________________
    Os= new_element("Os", 76, 7)
        Os%isot(1) = new_isotope(184, 183.9524891_dp , 0.0002_dp )
        Os%isot(2) = new_isotope(186, 185.9538382_dp , 0.0159_dp )
        Os%isot(3) = new_isotope(187, 186.9557505_dp , 0.0196_dp )
        Os%isot(4) = new_isotope(188, 187.9558382_dp , 0.1324_dp )
        Os%isot(5) = new_isotope(189, 188.9581475_dp , 0.1615_dp )
        Os%isot(6) = new_isotope(190, 189.9584470_dp , 0.2626_dp )
        Os%isot(7) = new_isotope(192, 191.9614807_dp , 0.4078_dp )
    !_________________________________________________________________________
    Ir= new_element("Ir", 77, 2)
        Ir%isot(1) = new_isotope(191, 190.9605940_dp , 0.373_dp )
        Ir%isot(2) = new_isotope(193, 192.9629264_dp , 0.627_dp )
    !_________________________________________________________________________
    Pt= new_element("Pt", 78, 6)
        Pt%isot(1) = new_isotope(190, 189.959932_dp ,  0.00014_dp )
        Pt%isot(2) = new_isotope(192, 191.9610380_dp , 0.00782_dp )
        Pt%isot(3) = new_isotope(194, 193.9626803_dp , 0.32967_dp )
        Pt%isot(4) = new_isotope(195, 194.9647911_dp , 0.33832_dp )
        Pt%isot(5) = new_isotope(196, 195.9649515_dp , 0.25242_dp )
        Pt%isot(6) = new_isotope(198, 197.967893_dp ,  0.07163_dp )
    !_________________________________________________________________________
    Au= new_element("Au", 79, 1)
        Au%isot(1) = new_isotope(197, 196.9665687_dp , 1.0000_dp )
    !_________________________________________________________________________
    Hg= new_element("Hg", 80, 7)
        Hg%isot(1) = new_isotope(196, 195.965833_dp ,  0.0015_dp )
        Hg%isot(2) = new_isotope(198, 197.9667690_dp , 0.0997_dp )
        Hg%isot(3) = new_isotope(199, 198.9682799_dp , 0.1687_dp )
        Hg%isot(4) = new_isotope(200, 199.9683260_dp , 0.2310_dp )
        Hg%isot(5) = new_isotope(201, 200.9703023_dp , 0.1318_dp )
        Hg%isot(6) = new_isotope(202, 201.9706430_dp , 0.2986_dp )
        Hg%isot(7) = new_isotope(204, 203.9734939_dp , 0.0687_dp )
    !_________________________________________________________________________
    Tl= new_element("Tl", 81, 2)
        Tl%isot(1) = new_isotope(203, 202.9723442_dp , 0.2952_dp )
        Tl%isot(2) = new_isotope(205, 204.9744275_dp , 0.7048_dp )
    !_________________________________________________________________________
    Pb= new_element("Pb", 82, 4)
        Pb%isot(1) = new_isotope(204, 203.9730436_dp , 0.014_dp )
        Pb%isot(2) = new_isotope(206, 205.9744653_dp , 0.241_dp )
        Pb%isot(3) = new_isotope(207, 206.9758969_dp , 0.221_dp )
        Pb%isot(4) = new_isotope(208, 207.9766521_dp , 0.524_dp )
    !_________________________________________________________________________
    Bi= new_element("Bi", 83, 1)
        Bi%isot(1) = new_isotope(209, 208.9803987_dp , 1.0000_dp )
    !_________________________________________________________________________
    Po= new_element("Po", 84, 2)
        Po%isot(1) = new_isotope(209, 208.9824304_dp , 0.0_dp )
        Po%isot(2) = new_isotope(210, 209.9828737_dp , 0.0_dp )
    !_________________________________________________________________________
    At= new_element("At", 85, 2)
        At%isot(1) = new_isotope(210, 209.987148_dp ,  0.0_dp )
        At%isot(2) = new_isotope(211, 210.9874963_dp , 0.0_dp )
    !_________________________________________________________________________
    Rn= new_element("Rn", 86, 3)
        Rn%isot(1) = new_isotope(211, 210.990601_dp ,  0.0_dp )
        Rn%isot(2) = new_isotope(220, 220.0113940_dp , 0.0_dp )
        Rn%isot(3) = new_isotope(222, 222.0175777_dp , 0.0_dp )
    !_________________________________________________________________________
    Fr= new_element("Fr", 87, 1)
        Fr%isot(1) = new_isotope(223, 223.0197359_dp , 0.0_dp )
    !_________________________________________________________________________
    Ra= new_element("Ra", 88, 4)
        Ra%isot(1) = new_isotope(223, 223.0185022_dp , 0.0_dp )
        Ra%isot(2) = new_isotope(224, 224.0202118_dp , 0.0_dp )
        Ra%isot(3) = new_isotope(226, 226.0254098_dp , 0.0_dp )
        Ra%isot(4) = new_isotope(228, 228.0310703_dp , 0.0_dp )
    !_________________________________________________________________________
    Ac= new_element("Ac", 89, 1)
        Ac%isot(1) = new_isotope(227, 227.0277521_dp , 0.0_dp )
    !_________________________________________________________________________
    Th= new_element("Th", 90, 2)
        Th%isot(1) = new_isotope(230, 230.0331338_dp , 0.0_dp )
        Th%isot(2) = new_isotope(232, 232.0380553_dp , 1.0000_dp )
    !_________________________________________________________________________
    Pa= new_element("Pa", 91, 1)
        Pa%isot(1) = new_isotope(231, 231.0358840_dp , 1.0000_dp )
    !_________________________________________________________________________
    U= new_element("U", 92, 5)
        U%isot(1) = new_isotope(233, 233.0396352_dp , 0.0_dp )
        U%isot(2) = new_isotope(234, 234.0409521_dp , 0.000054_dp )
        U%isot(3) = new_isotope(235, 235.0439299_dp , 0.007204_dp )
        U%isot(4) = new_isotope(236, 236.0455680_dp , 0.0_dp )
        U%isot(5) = new_isotope(238, 238.0507882_dp , 0.992742_dp )
    !_________________________________________________________________________
    Np= new_element("Np", 93, 2)
        Np%isot(1) = new_isotope(236, 236.046570_dp ,  0.0_dp )
        Np%isot(2) = new_isotope(237, 237.0481734_dp , 0.0_dp )
    !_________________________________________________________________________
    Pu= new_element("Pu", 94, 6)
        Pu%isot(1) = new_isotope(238, 238.0495599_dp , 0.0_dp )
        Pu%isot(2) = new_isotope(239, 239.0521634_dp , 0.0_dp )
        Pu%isot(3) = new_isotope(240, 240.0538135_dp , 0.0_dp )
        Pu%isot(4) = new_isotope(241, 241.0568515_dp , 0.0_dp )
        Pu%isot(5) = new_isotope(242, 242.0587426_dp , 0.0_dp )
        Pu%isot(6) = new_isotope(244, 244.064204_dp ,  0.0_dp )
    !_________________________________________________________________________
    Am= new_element("Am", 95, 2)
        Am%isot(1) = new_isotope(241, 241.0568291_dp , 0.0_dp )
        Am%isot(2) = new_isotope(243, 243.0613811_dp , 0.0_dp )
    !_________________________________________________________________________
    Cm= new_element("Cm", 96, 6)
        Cm%isot(1) = new_isotope(243, 243.0613891_dp , 0.0_dp )
        Cm%isot(2) = new_isotope(244, 244.0627526_dp , 0.0_dp )
        Cm%isot(3) = new_isotope(245, 245.0654912_dp , 0.0_dp )
        Cm%isot(4) = new_isotope(246, 246.0672237_dp , 0.0_dp )
        Cm%isot(5) = new_isotope(247, 247.070354_dp ,  0.0_dp )
        Cm%isot(6) = new_isotope(248, 248.072349_dp ,  0.0_dp )
    !_________________________________________________________________________
    Bk= new_element("Bk", 97, 2)
        Bk%isot(1) = new_isotope(247, 247.070307_dp ,  0.0_dp )
        Bk%isot(2) = new_isotope(249, 249.0749867_dp , 0.0_dp )
    !_________________________________________________________________________
    Cf= new_element("Cf", 98, 4)
        Cf%isot(1) = new_isotope(249, 249.0748535_dp , 0.0_dp )
        Cf%isot(2) = new_isotope(250, 250.0764061_dp , 0.0_dp )
        Cf%isot(3) = new_isotope(251, 251.079587_dp ,  0.0_dp )
        Cf%isot(4) = new_isotope(252, 252.081626_dp ,  0.0_dp )
    !_________________________________________________________________________
    Es= new_element("Es", 99, 1)
        Es%isot(1) = new_isotope(252, 252.082980_dp , 0.0_dp )
    !_________________________________________________________________________
    SELECT CASE ( TRIM(name) )
      CASE("H")
        self=H
      CASE("D")
        self=D
      CASE("T")
        self=T
      CASE("He")
        self=He
      CASE("Li")
        self=Li
      CASE("Be")
        self=Be
      CASE("B")
        self=B
      CASE("C")
        self=C
      CASE("N")
        self=N
      CASE("O")
        self=O
      CASE("F")
        self=F
      CASE("Ne")
        self=Ne
      CASE("Na")
        self=Na
      CASE("Mg")
        self=Mg
      CASE("Al")
        self=Al
      CASE("Si")
        self=Si
      CASE("P")
        self=P
      CASE("S")
        self=S
      CASE("Cl")
        self=Cl
      CASE("Ar")
        self=Ar
      CASE("K")
        self=K
      CASE("Ca")
        self=Ca
      CASE("Sc")
        self=Sc
      CASE("Ti")
        self=Ti
      CASE("V")
        self=V
      CASE("Cr")
        self=Cr
      CASE("Mn")
        self=Mn
      CASE("Fe")
        self=Fe
      CASE("Co")
        self=Co
      CASE("Ni")
        self=Ni
      CASE("Cu")
        self=Cu
      CASE("Zn")
        self=Zn
      CASE("Ga")
        self=Ga
      CASE("Ge")
        self=Ge
      CASE("As")
        self=As
      CASE("Se")
        self=Se
      CASE("Br")
        self=Br
      CASE("Kr")
        self=Kr
      CASE("Rb")
        self=Rb
      CASE("Sr")
        self=Sr
      CASE("Y")
        self=Y
      CASE("Zr")
        self=Zr
      CASE("Nb")
        self=Nb
      CASE("Mo")
        self=Mo
      CASE("Tc")
        self=Tc
      CASE("Ru")
        self=Ru
      CASE("Rh")
        self=Rh
      CASE("Pd")
        self=Pd
      CASE("Ag")
        self=Ag
      CASE("Cd")
        self=Cd
      CASE("In")
        self=In
      CASE("Sn")
        self=Sn
      CASE("Sb")
        self=Sb
      CASE("Te")
        self=Te
      CASE("I")
        self=I
      CASE("Xe")
        self=Xe
      CASE("Cs")
        self=Cs
      CASE("Ba")
        self=Ba
      CASE("La")
        self=La
      CASE("Ce")
        self=Ce
      CASE("Pr")
        self=Pr
      CASE("Nd")
        self=Nd
      CASE("Pm")
        self=Pm
      CASE("Sm")
        self=Sm
      CASE("Eu")
        self=Eu
      CASE("Gd")
        self=Gd
      CASE("Tb")
        self=Tb
      CASE("Dy")
        self=Dy
      CASE("Ho")
        self=Ho
      CASE("Er")
        self=Er
      CASE("Tm")
        self=Tm
      CASE("Yb")
        self=Yb
      CASE("Lu")
        self=Lu
      CASE("Hf")
        self=Hf
      CASE("Ta")
        self=Ta
      CASE("W")
        self=W
      CASE("Re")
        self=Re
      CASE("Os")
        self=Os
      CASE("Ir")
        self=Ir
      CASE("Pt")
        self=Pt
      CASE("Au")
        self=Au
      CASE("Hg")
        self=Hg
      CASE("Tl")
        self=Tl
      CASE("Pb")
        self=Pb
      CASE("Bi")
        self=Bi
      CASE("Po")
        self=Po
      CASE("At")
        self=At
      CASE("Rn")
        self=Rn
      CASE("Fr")
        self=Fr
      CASE("Ra")
        self=Ra
      CASE("Ac")
        self=Ac
      CASE("Th")
        self=Th
      CASE("Pa")
        self=Pa
      CASE("U")
        self=U
      CASE("Np")
        self=Np
      CASE("Pu")
        self=Pu
      CASE("Am")
        self=Am
      CASE("Cm")
        self=Cm
      CASE("Bk")
        self=Bk
      CASE("Cf")
        self=Cf
      CASE("Es")
        self=Es
      CASE DEFAULT
        CALL errore("search_nist", "element '"//TRIM(name)//"' not found.", 1)
    END SELECT
  END FUNCTION search_nist

END MODULE nist_isotopes_db


