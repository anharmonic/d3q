!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
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
MODULE nist_isotopes_database
    USE kinds, ONLY : DP
    
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
        final     :: destroy_element
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
    REAL,INTENT(in):: mass, conc ! get them in single precision, for simplicity
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
    TYPE(element) :: elem
    !
    INTEGER :: i, nisot_
    REAL(DP),ALLOCATABLE :: amass_(:) ! isotope atomic mass
    REAL(DP),ALLOCATABLE :: aconc_(:) ! isotope concentration
    !
    NISOT_IF : &
    IF(nisot>0)THEN
      IF(.not.PRESENT(amass) .or. .not.PRESENT(aconc))&
        CALL errore("compute_gs", "You need to provide amass and aconc for nisot>0",1)
      IF(.not.ALLOCATED(amass) .or. .not.allocated(aconc))&
        CALL errore("compute_gs", "You need to provide amass and aconc for nisot>0",2)
      IF( ABS(SUM(aconc)-1._dp)>1.e-4_dp) &
        CALL errore("compute_gs", "the sum of isotope concentrations is not 100%", 3)
        
      gm = SUM( aconc*amass )
      gs = SUM( aconc*(amass-gm)**2 ) / gm**2
      
    ELSE NISOT_IF
      
      elem = search_nist(aname)
      WRITE(*,'(5x,2a)') "Found in NIST: ", elem%aname
    
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
      ENDIF
    ENDIF &
    NISOT_IF
    !
    WRITE(*,'(8x,2a,i3)') "number of isotopes selected: ", elem%nisot
    DO i = 1, elem%nisot
      WRITE(*,'(10x,2f12.6)') amass_(i), aconc_(i)
    ENDDO

    WRITE(*,'(8x,a,2f12.6)')  "Average mass, relative standard dev.:", gm, gs

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
        H%isot(1) = new_isotope(1, 1.00782503207, 0.999885)
        H%isot(2) = new_isotope(2, 2.0141017778,  0.000115)
        H%isot(3) = new_isotope(3, 3.016049277,   0.000000)
      D = new_element("D", 1, 1)
        D%isot(1) = new_isotope(2, 2.0141017778, 1.0)
      T = new_element("T", 1, 1)
        T%isot(1) = new_isotope(3, 3.0160492777, 1.0)
    !_________________________________________________________________________
      He = new_element("He" , 2 , 2)
        He%isot(1) = new_isotope(3, 3.0160293191,  0.00000134)
        He%isot(2) = new_isotope(4, 4.00260325415, 0.99999866)
    !_________________________________________________________________________
      Li= new_element("Li", 3, 2)
        Li%isot(1) = new_isotope(6, 6.015122795, 0.0759)
        Li%isot(2) = new_isotope(7, 7.01600455,  0.9241)
    !_________________________________________________________________________
      Be= new_element("Be", 4, 1)
        Be%isot(1) = new_isotope(9, 9.0121822, 1.0000)
    !_________________________________________________________________________
      B= new_element("B", 5, 2)
        B%isot(1) = new_isotope(10, 10.0129370, 0.199)
        B%isot(2) = new_isotope(11, 11.0093054, 0.801)
    !_________________________________________________________________________
      C= new_element("C", 6, 3)
        C%isot(1) = new_isotope(12, 12.0000000,    0.9893)
        C%isot(2) = new_isotope(13, 13.0033548378, 0.0107)
        C%isot(3) = new_isotope(14, 14.003241989,  0.0)
    !_________________________________________________________________________
      N= new_element("N", 7, 2)
        N%isot(1) = new_isotope(14, 14.0030740048, 0.99636)
        N%isot(2) = new_isotope(15, 15.0001088982, 0.00364)
    !_________________________________________________________________________
      O= new_element("O", 8, 3)
        O%isot(1) = new_isotope(16, 15.99491461956, 0.99757)
        O%isot(2) = new_isotope(17, 16.99913170,    0.00038)
        O%isot(3) = new_isotope(18, 17.9991610,     0.00205)
    !_________________________________________________________________________
      F= new_element("F", 9, 1)
        F%isot(1) = new_isotope(19, 18.99840322, 1.0000)
    !_________________________________________________________________________
    Ne= new_element("Ne", 10, 3)
        Ne%isot(1) = new_isotope(20, 19.9924401754, 0.9048)
        Ne%isot(2) = new_isotope(21, 20.99384668,   0.0027)
        Ne%isot(3) = new_isotope(22, 21.991385114,  0.0925)
    !_________________________________________________________________________
    Na= new_element("Na", 11, 1)
        Na%isot(1) = new_isotope(23, 22.9897692809, 1.0000)
    !_________________________________________________________________________
    Mg= new_element("Mg", 12, 3)
        Mg%isot(1) = new_isotope(24, 23.985041700, 0.7899)
        Mg%isot(2) = new_isotope(25, 24.98583692,  0.1000)
        Mg%isot(3) = new_isotope(26, 25.982592929, 0.1101)
    !_________________________________________________________________________
    Al= new_element("Al", 13, 1)
        Al%isot(1) = new_isotope(27, 26.98153863, 1.0000)
    !_________________________________________________________________________
    Si= new_element("Si", 14, 3)
        Si%isot(1) = new_isotope(28, 27.9769265325, 0.92223)
        Si%isot(2) = new_isotope(29, 28.976494700,  0.04685)
        Si%isot(3) = new_isotope(30, 29.97377017,   0.03092)
    !_________________________________________________________________________
    P= new_element("P", 15, 1)
        P%isot(1) = new_isotope(31, 30.97376163, 1.0000)
    !_________________________________________________________________________
    S= new_element("S", 16, 4)
        S%isot(1) = new_isotope(32, 31.97207100, 0.9499)
        S%isot(2) = new_isotope(33, 32.97145876, 0.0075)
        S%isot(3) = new_isotope(34, 33.96786690, 0.0425)
        S%isot(4) = new_isotope(36, 35.96708076, 0.0001)
    !_________________________________________________________________________
    Cl= new_element("Cl", 17, 2)
        Cl%isot(1) = new_isotope(35, 34.96885268, 0.7576)
        Cl%isot(2) = new_isotope(37, 36.96590259, 0.2424)
    !_________________________________________________________________________
    Ar= new_element("Ar", 18, 3)
        Ar%isot(1) = new_isotope(36, 35.967545106,  0.003365)
        Ar%isot(2) = new_isotope(38, 37.9627324,    0.000632)
        Ar%isot(3) = new_isotope(40, 39.9623831225, 0.996003)
    !_________________________________________________________________________
    K= new_element("K", 19, 3)
        K%isot(1) = new_isotope(39, 38.96370668, 0.932581)
        K%isot(2) = new_isotope(40, 39.96399848, 0.000117)
        K%isot(3) = new_isotope(41, 40.96182576, 0.067302)
    !_________________________________________________________________________
    Ca= new_element("Ca", 20, 6)
        Ca%isot(1) = new_isotope(40, 39.96259098, 0.96941)
        Ca%isot(2) = new_isotope(42, 41.95861801, 0.00647)
        Ca%isot(3) = new_isotope(43, 42.9587666,  0.00135)
        Ca%isot(4) = new_isotope(44, 43.9554818,  0.02086)
        Ca%isot(5) = new_isotope(46, 45.9536926,  0.00004)
        Ca%isot(6) = new_isotope(48, 47.952534,   0.00187)
    !_________________________________________________________________________
    Sc= new_element("Sc", 21, 1)
        Sc%isot(1) = new_isotope(45, 44.9559119, 1.0000)
    !_________________________________________________________________________
    Ti= new_element("Ti", 22, 5)
        Ti%isot(1) = new_isotope(46, 45.9526316, 0.0825)
        Ti%isot(2) = new_isotope(47, 46.9517631, 0.0744)
        Ti%isot(3) = new_isotope(48, 47.9479463, 0.7372)
        Ti%isot(4) = new_isotope(49, 48.9478700, 0.0541)
        Ti%isot(5) = new_isotope(50, 49.9447912, 0.0518)
    !_________________________________________________________________________
    V= new_element("V", 23, 2)
        V%isot(1) = new_isotope(50, 49.9471585, 0.00250)
        V%isot(2) = new_isotope(51, 50.9439595, 0.99750)
    !_________________________________________________________________________
    Cr= new_element("Cr", 24, 4)
        Cr%isot(1) = new_isotope(50, 49.9460442, 0.04345)
        Cr%isot(2) = new_isotope(52, 51.9405075, 0.83789)
        Cr%isot(3) = new_isotope(53, 52.9406494, 0.09501)
        Cr%isot(4) = new_isotope(54, 53.9388804, 0.02365)
    !_________________________________________________________________________
    Mn= new_element("Mn", 25, 1)
        Mn%isot(1) = new_isotope(55, 54.9380451, 1.0000)
    !_________________________________________________________________________
    Fe= new_element("Fe", 26, 4)
        Fe%isot(1) = new_isotope(54, 53.9396105, 0.05845)
        Fe%isot(2) = new_isotope(56, 55.9349375, 0.91754)
        Fe%isot(3) = new_isotope(57, 56.9353940, 0.02119)
        Fe%isot(4) = new_isotope(58, 57.9332756, 0.00282)
    !_________________________________________________________________________
    Co= new_element("Co", 27, 1)
        Co%isot(1) = new_isotope(59, 58.9331950, 1.0000)
    !_________________________________________________________________________
    Ni= new_element("Ni", 28, 5)
        Ni%isot(1) = new_isotope(58, 57.9353429, 0.680769)
        Ni%isot(2) = new_isotope(60, 59.9307864, 0.262231)
        Ni%isot(3) = new_isotope(61, 60.9310560, 0.011399)
        Ni%isot(4) = new_isotope(62, 61.9283451, 0.036345)
        Ni%isot(5) = new_isotope(64, 63.9279660, 0.009256)
    !_________________________________________________________________________
    Cu= new_element("Cu", 29, 2)
        Cu%isot(1) = new_isotope(63, 62.9295975, 0.6915)
        Cu%isot(2) = new_isotope(65, 64.9277895, 0.3085)
    !_________________________________________________________________________
    Zn= new_element("Zn", 30, 5)
        Zn%isot(1) = new_isotope(64, 63.9291422, 0.48268)
        Zn%isot(2) = new_isotope(66, 65.9260334, 0.27975)
        Zn%isot(3) = new_isotope(67, 66.9271273, 0.04102)
        Zn%isot(4) = new_isotope(68, 67.9248442, 0.19024)
        Zn%isot(5) = new_isotope(70, 69.9253193, 0.00631)
    !_________________________________________________________________________
    Ga= new_element("Ga", 31, 2)
        Ga%isot(1) = new_isotope(69, 68.9255736, 0.60108)
        Ga%isot(2) = new_isotope(71, 70.9247013, 0.39892)
    !_________________________________________________________________________
    Ge= new_element("Ge", 32, 5)
        Ge%isot(1) = new_isotope(70, 69.9242474, 0.2038)
        Ge%isot(2) = new_isotope(72, 71.9220758, 0.2731)
        Ge%isot(3) = new_isotope(73, 72.9234589, 0.0776)
        Ge%isot(4) = new_isotope(74, 73.9211778, 0.3672)
        Ge%isot(5) = new_isotope(76, 75.9214026, 0.0783)
    !_________________________________________________________________________
    As= new_element("As", 33, 1)
        As%isot(1) = new_isotope(75, 74.9215965, 1.0000)
    !_________________________________________________________________________
    Se= new_element("Se", 34, 6)
        Se%isot(1) = new_isotope(74, 73.9224764, 0.0089)
        Se%isot(2) = new_isotope(76, 75.9192136, 0.0937)
        Se%isot(3) = new_isotope(77, 76.9199140, 0.0763)
        Se%isot(4) = new_isotope(78, 77.9173091, 0.2377)
        Se%isot(5) = new_isotope(80, 79.9165213, 0.4961)
        Se%isot(6) = new_isotope(82, 81.9166994, 0.0873)
    !_________________________________________________________________________
    Br= new_element("Br", 35, 2)
        Br%isot(1) = new_isotope(79, 78.9183371, 0.5069)
        Br%isot(2) = new_isotope(81, 80.9162906, 0.4931)
    !_________________________________________________________________________
    Kr= new_element("Kr", 36, 6)
        Kr%isot(1) = new_isotope(78, 77.9203648,  0.00355)
        Kr%isot(2) = new_isotope(80, 79.9163790,  0.02286)
        Kr%isot(3) = new_isotope(82, 81.9134836,  0.11593)
        Kr%isot(4) = new_isotope(83, 82.914136,   0.11500)
        Kr%isot(5) = new_isotope(84, 83.911507,   0.56987)
        Kr%isot(6) = new_isotope(86, 85.91061073, 0.17279)
    !_________________________________________________________________________
    Rb= new_element("Rb", 37, 2)
        Rb%isot(1) = new_isotope(85, 84.911789738, 0.7217)
        Rb%isot(2) = new_isotope(87, 86.909180527, 0.2783)
    !_________________________________________________________________________
    Sr= new_element("Sr", 38, 4)
        Sr%isot(1) = new_isotope(84, 83.913425,  0.0056)
        Sr%isot(2) = new_isotope(86, 85.9092602, 0.0986)
        Sr%isot(3) = new_isotope(87, 86.9088771, 0.0700)
        Sr%isot(4) = new_isotope(88, 87.9056121, 0.8258)
    !_________________________________________________________________________
    Y= new_element("Y", 39, 1)
        Y%isot(1) = new_isotope(89, 88.9058483, 1.0000)
    !_________________________________________________________________________
    Zr= new_element("Zr", 40, 5)
        Zr%isot(1) = new_isotope(90, 89.9047044, 0.5145)
        Zr%isot(2) = new_isotope(91, 90.9056458, 0.1122)
        Zr%isot(3) = new_isotope(92, 91.9050408, 0.1715)
        Zr%isot(4) = new_isotope(94, 93.9063152, 0.1738)
        Zr%isot(5) = new_isotope(96, 95.9082734, 0.0280)
    !_________________________________________________________________________
    Nb= new_element("Nb", 41, 1)
        Nb%isot(1) = new_isotope(93, 92.9063781, 1.0000)
    !_________________________________________________________________________
    Mo= new_element("Mo", 42, 7)
        Mo%isot(1) = new_isotope(92, 91.906811,  0.1477)
        Mo%isot(2) = new_isotope(94, 93.9050883, 0.0923)
        Mo%isot(3) = new_isotope(95, 94.9058421, 0.1590)
        Mo%isot(4) = new_isotope(96, 95.9046795, 0.1668)
        Mo%isot(5) = new_isotope(97, 96.9060215, 0.0956)
        Mo%isot(6) = new_isotope(98, 97.9054082, 0.2419)
        Mo%isot(7) = new_isotope(100, 99.907477, 0.0967)
    !_________________________________________________________________________
    Tc= new_element("Tc", 43, 3)
        Tc%isot(1) = new_isotope(97, 96.906365,  0.0)
        Tc%isot(2) = new_isotope(98, 97.907216,  0.0)
        Tc%isot(3) = new_isotope(99, 98.9062547, 0.0)
    !_________________________________________________________________________
    Ru= new_element("Ru", 44, 7)
        Ru%isot(1) = new_isotope(96, 95.907598,    0.0554)
        Ru%isot(2) = new_isotope(98, 97.905287,    0.0187)
        Ru%isot(3) = new_isotope(99, 98.9059393,   0.1276)
        Ru%isot(4) = new_isotope(100, 99.9042195,  0.1260)
        Ru%isot(5) = new_isotope(101, 100.9055821, 0.1706)
        Ru%isot(6) = new_isotope(102, 101.9043493, 0.3155)
        Ru%isot(7) = new_isotope(104, 103.905433,  0.1862)
    !_________________________________________________________________________
    Rh= new_element("Rh", 45, 1)
        Rh%isot(1) = new_isotope(103, 102.905504, 1.0000)
    !_________________________________________________________________________
    Pd= new_element("Pd", 46, 6)
        Pd%isot(1) = new_isotope(102, 101.905609, 0.0102)
        Pd%isot(2) = new_isotope(104, 103.904036, 0.1114)
        Pd%isot(3) = new_isotope(105, 104.905085, 0.2233)
        Pd%isot(4) = new_isotope(106, 105.903486, 0.2733)
        Pd%isot(5) = new_isotope(108, 107.903892, 0.2646)
        Pd%isot(6) = new_isotope(110, 109.905153, 0.1172)
    !_________________________________________________________________________
    Ag= new_element("Ag", 47, 2)
        Ag%isot(1) = new_isotope(107, 106.905097, 0.51839)
        Ag%isot(2) = new_isotope(109, 108.904752, 0.48161)
    !_________________________________________________________________________
    Cd= new_element("Cd", 48, 8)
        Cd%isot(1) = new_isotope(106, 105.906459,  0.0125)
        Cd%isot(2) = new_isotope(108, 107.904184,  0.0089)
        Cd%isot(3) = new_isotope(110, 109.9030021, 0.1249)
        Cd%isot(4) = new_isotope(111, 110.9041781, 0.1280)
        Cd%isot(5) = new_isotope(112, 111.9027578, 0.2413)
        Cd%isot(6) = new_isotope(113, 112.9044017, 0.1222)
        Cd%isot(7) = new_isotope(114, 113.9033585, 0.2873)
        Cd%isot(8) = new_isotope(116, 115.904756,  0.0749)
    !_________________________________________________________________________
    In= new_element("In", 49, 2)
        In%isot(1) = new_isotope(113, 112.904058, 0.0429)
        In%isot(2) = new_isotope(115, 114.903878, 0.9571)
    !_________________________________________________________________________
    Sn= new_element("Sn", 50, 10)
        Sn%isot(1)  = new_isotope(112, 111.904818,  0.0097)
        Sn%isot(2)  = new_isotope(114, 113.902779,  0.0066)
        Sn%isot(3)  = new_isotope(115, 114.903342,  0.0034)
        Sn%isot(4)  = new_isotope(116, 115.901741,  0.1454)
        Sn%isot(5)  = new_isotope(117, 116.902952,  0.0768)
        Sn%isot(6)  = new_isotope(118, 117.901603,  0.2422)
        Sn%isot(7)  = new_isotope(119, 118.903308,  0.0859)
        Sn%isot(8)  = new_isotope(120, 119.9021947, 0.3258)
        Sn%isot(9)  = new_isotope(122, 121.9034390, 0.0463)
        Sn%isot(10) = new_isotope(124, 123.9052739, 0.0579)
    !_________________________________________________________________________
    Sb= new_element("Sb", 51, 2)
        Sb%isot(1) = new_isotope(121, 120.9038157, 0.5721)
        Sb%isot(2) = new_isotope(123, 122.9042140, 0.4279)
    !_________________________________________________________________________
    Te= new_element("Te", 52, 8)
        Te%isot(1) = new_isotope(120, 119.904020,  0.0009)
        Te%isot(2) = new_isotope(122, 121.9030439, 0.0255)
        Te%isot(3) = new_isotope(123, 122.9042700, 0.0089)
        Te%isot(4) = new_isotope(124, 123.9028179, 0.0474)
        Te%isot(5) = new_isotope(125, 124.9044307, 0.0707)
        Te%isot(6) = new_isotope(126, 125.9033117, 0.1884)
        Te%isot(7) = new_isotope(128, 127.9044631, 0.3174)
        Te%isot(8) = new_isotope(130, 129.9062244, 0.3408)
    !_________________________________________________________________________
    I= new_element("I", 53, 1)
        I%isot(1) = new_isotope(127, 126.904473, 1.0000)
    !_________________________________________________________________________
    Xe= new_element("Xe", 54, 9)
        Xe%isot(1) = new_isotope(124, 123.9058930, 0.000952)
        Xe%isot(2) = new_isotope(126, 125.904274,  0.000890)
        Xe%isot(3) = new_isotope(128, 127.9035313, 0.019102)
        Xe%isot(4) = new_isotope(129, 128.9047794, 0.264006)
        Xe%isot(5) = new_isotope(130, 129.9035080, 0.040710)
        Xe%isot(6) = new_isotope(131, 130.9050824, 0.212324)
        Xe%isot(7) = new_isotope(132, 131.9041535, 0.269086)
        Xe%isot(8) = new_isotope(134, 133.9053945, 0.104357)
        Xe%isot(9) = new_isotope(136, 135.907219,  0.088573)
    !_________________________________________________________________________
    Cs= new_element("Cs", 55, 1)
        Cs%isot(1) = new_isotope(133, 132.905451933, 1.0000)
    !_________________________________________________________________________
    Ba= new_element("Ba", 56, 7)
        Ba%isot(1) = new_isotope(130, 129.9063208, 0.00106)
        Ba%isot(2) = new_isotope(132, 131.9050613, 0.00101)
        Ba%isot(3) = new_isotope(134, 133.9045084, 0.02417)
        Ba%isot(4) = new_isotope(135, 134.9056886, 0.06592)
        Ba%isot(5) = new_isotope(136, 135.9045759, 0.07854)
        Ba%isot(6) = new_isotope(137, 136.9058274, 0.11232)
        Ba%isot(7) = new_isotope(138, 137.9052472, 0.71698)
    !_________________________________________________________________________
    La= new_element("La", 57, 2)
        La%isot(1) = new_isotope(138, 137.907112, 0.00090)
        La%isot(2) = new_isotope(139, 138.9063533, 0.99910)
    !_________________________________________________________________________
    Ce= new_element("Ce", 58, 4)
        Ce%isot(1) = new_isotope(136, 135.907172,  0.00185)
        Ce%isot(2) = new_isotope(138, 137.905991,  0.00251)
        Ce%isot(3) = new_isotope(140, 139.9054387, 0.88450)
        Ce%isot(4) = new_isotope(142, 141.909244,  0.11114)
    !_________________________________________________________________________
    Pr= new_element("Pr", 59, 1)
        Pr%isot(1) = new_isotope(141, 140.9076528, 1.0000)
    !_________________________________________________________________________
    Nd= new_element("Nd", 60, 7)
        Nd%isot(1) = new_isotope(142, 141.9077233, 0.272)
        Nd%isot(2) = new_isotope(143, 142.9098143, 0.122)
        Nd%isot(3) = new_isotope(144, 143.9100873, 0.238)
        Nd%isot(4) = new_isotope(145, 144.9125736, 0.083)
        Nd%isot(5) = new_isotope(146, 145.9131169, 0.172)
        Nd%isot(6) = new_isotope(148, 147.916893, 0.057)
        Nd%isot(7) = new_isotope(150, 149.920891, 0.056)
    !_________________________________________________________________________
    Pm= new_element("Pm", 61, 2)
        Pm%isot(1) = new_isotope(145, 144.912749,  0.0)
        Pm%isot(2) = new_isotope(147, 146.9151385, 0.0)
    !_________________________________________________________________________
    Sm= new_element("Sm", 62, 7)
        Sm%isot(1) = new_isotope(144, 143.911999,  0.0307)
        Sm%isot(2) = new_isotope(147, 146.9148979, 0.1499)
        Sm%isot(3) = new_isotope(148, 147.9148227, 0.1124)
        Sm%isot(4) = new_isotope(149, 148.9171847, 0.1382)
        Sm%isot(5) = new_isotope(150, 149.9172755, 0.0738)
        Sm%isot(6) = new_isotope(152, 151.9197324, 0.2675)
        Sm%isot(7) = new_isotope(154, 153.9222093, 0.2275)
    !_________________________________________________________________________
    Eu= new_element("Eu", 63, 2)
        Eu%isot(1) = new_isotope(151, 150.9198502, 0.4781)
        Eu%isot(2) = new_isotope(153, 152.9212303, 0.5219)
    !_________________________________________________________________________
    Gd= new_element("Gd", 64, 7)
        Gd%isot(1) = new_isotope(152, 151.9197910, 0.0020)
        Gd%isot(2) = new_isotope(154, 153.9208656, 0.0218)
        Gd%isot(3) = new_isotope(155, 154.9226220, 0.1480)
        Gd%isot(4) = new_isotope(156, 155.9221227, 0.2047)
        Gd%isot(5) = new_isotope(157, 156.9239601, 0.1565)
        Gd%isot(6) = new_isotope(158, 157.9241039, 0.2484)
        Gd%isot(7) = new_isotope(160, 159.9270541, 0.2186)
    !_________________________________________________________________________
    Tb= new_element("Tb", 65, 1)
        Tb%isot(1) = new_isotope(159, 158.9253468, 1.0000)
    !_________________________________________________________________________
    Dy= new_element("Dy", 66, 7)
        Dy%isot(1) = new_isotope(156, 155.924283,  0.00056)
        Dy%isot(2) = new_isotope(158, 157.924409,  0.00095)
        Dy%isot(3) = new_isotope(160, 159.9251975, 0.02329)
        Dy%isot(4) = new_isotope(161, 160.9269334, 0.18889)
        Dy%isot(5) = new_isotope(162, 161.9267984, 0.25475)
        Dy%isot(6) = new_isotope(163, 162.9287312, 0.24896)
        Dy%isot(7) = new_isotope(164, 163.9291748, 0.28260)
    !_________________________________________________________________________
    Ho= new_element("Ho", 67, 1)
        Ho%isot(1) = new_isotope(165, 164.9303221, 1.0000)
    !_________________________________________________________________________
    Er= new_element("Er", 68, 6)
        Er%isot(1) = new_isotope(162, 161.928778,  0.00139)
        Er%isot(2) = new_isotope(164, 163.929200,  0.01601)
        Er%isot(3) = new_isotope(166, 165.9302931, 0.33503)
        Er%isot(4) = new_isotope(167, 166.9320482, 0.22869)
        Er%isot(5) = new_isotope(168, 167.9323702, 0.26978)
        Er%isot(6) = new_isotope(170, 169.9354643, 0.14910)
    !_________________________________________________________________________
    Tm= new_element("Tm", 69, 1)
        Tm%isot(1) = new_isotope(169, 168.9342133, 1.0000)
    !_________________________________________________________________________
    Yb= new_element("Yb", 70, 7)
        Yb%isot(1) = new_isotope(168, 167.933897,  0.0013)
        Yb%isot(2) = new_isotope(170, 169.9347618, 0.0304)
        Yb%isot(3) = new_isotope(171, 170.9363258, 0.1428)
        Yb%isot(4) = new_isotope(172, 171.9363815, 0.2183)
        Yb%isot(5) = new_isotope(173, 172.9382108, 0.1613)
        Yb%isot(6) = new_isotope(174, 173.9388621, 0.3183)
        Yb%isot(7) = new_isotope(176, 175.9425717, 0.1276)
    !_________________________________________________________________________
    Lu= new_element("Lu", 71, 2)
        Lu%isot(1) = new_isotope(175, 174.9407718, 0.9741)
        Lu%isot(2) = new_isotope(176, 175.9426863, 0.0259)
    !_________________________________________________________________________
    Hf= new_element("Hf", 72, 6)
        Hf%isot(1) = new_isotope(174, 173.940046,  0.0016)
        Hf%isot(2) = new_isotope(176, 175.9414086, 0.0526)
        Hf%isot(3) = new_isotope(177, 176.9432207, 0.1860)
        Hf%isot(4) = new_isotope(178, 177.9436988, 0.2728)
        Hf%isot(5) = new_isotope(179, 178.9458161, 0.1362)
        Hf%isot(6) = new_isotope(180, 179.9465500, 0.3508)
    !_________________________________________________________________________
    Ta= new_element("Ta", 73, 2)
        Ta%isot(1) = new_isotope(180, 179.9474648, 0.00012)
        Ta%isot(2) = new_isotope(181, 180.9479958, 0.99988)
    !_________________________________________________________________________
    W= new_element("W", 74, 5)
        W%isot(1) = new_isotope(180, 179.946704,  0.0012)
        W%isot(2) = new_isotope(182, 181.9482042, 0.2650)
        W%isot(3) = new_isotope(183, 182.9502230, 0.1431)
        W%isot(4) = new_isotope(184, 183.9509312, 0.3064)
        W%isot(5) = new_isotope(186, 185.9543641, 0.2843)
    !_________________________________________________________________________
    Re= new_element("Re", 75, 2)
        Re%isot(1) = new_isotope(185, 184.9529550, 0.3740)
        Re%isot(2) = new_isotope(187, 186.9557531, 0.6260)
    !_________________________________________________________________________
    Os= new_element("Os", 76, 7)
        Os%isot(1) = new_isotope(184, 183.9524891, 0.0002)
        Os%isot(2) = new_isotope(186, 185.9538382, 0.0159)
        Os%isot(3) = new_isotope(187, 186.9557505, 0.0196)
        Os%isot(4) = new_isotope(188, 187.9558382, 0.1324)
        Os%isot(5) = new_isotope(189, 188.9581475, 0.1615)
        Os%isot(6) = new_isotope(190, 189.9584470, 0.2626)
        Os%isot(7) = new_isotope(192, 191.9614807, 0.4078)
    !_________________________________________________________________________
    Ir= new_element("Ir", 77, 2)
        Ir%isot(1) = new_isotope(191, 190.9605940, 0.373)
        Ir%isot(2) = new_isotope(193, 192.9629264, 0.627)
    !_________________________________________________________________________
    Pt= new_element("Pt", 78, 6)
        Pt%isot(1) = new_isotope(190, 189.959932,  0.00014)
        Pt%isot(2) = new_isotope(192, 191.9610380, 0.00782)
        Pt%isot(3) = new_isotope(194, 193.9626803, 0.32967)
        Pt%isot(4) = new_isotope(195, 194.9647911, 0.33832)
        Pt%isot(5) = new_isotope(196, 195.9649515, 0.25242)
        Pt%isot(6) = new_isotope(198, 197.967893,  0.07163)
    !_________________________________________________________________________
    Au= new_element("Au", 79, 1)
        Au%isot(1) = new_isotope(197, 196.9665687, 1.0000)
    !_________________________________________________________________________
    Hg= new_element("Hg", 80, 7)
        Hg%isot(1) = new_isotope(196, 195.965833,  0.0015)
        Hg%isot(2) = new_isotope(198, 197.9667690, 0.0997)
        Hg%isot(3) = new_isotope(199, 198.9682799, 0.1687)
        Hg%isot(4) = new_isotope(200, 199.9683260, 0.2310)
        Hg%isot(5) = new_isotope(201, 200.9703023, 0.1318)
        Hg%isot(6) = new_isotope(202, 201.9706430, 0.2986)
        Hg%isot(7) = new_isotope(204, 203.9734939, 0.0687)
    !_________________________________________________________________________
    Tl= new_element("Tl", 81, 2)
        Tl%isot(1) = new_isotope(203, 202.9723442, 0.2952)
        Tl%isot(2) = new_isotope(205, 204.9744275, 0.7048)
    !_________________________________________________________________________
    Pb= new_element("Pb", 82, 4)
        Pb%isot(1) = new_isotope(204, 203.9730436, 0.014)
        Pb%isot(2) = new_isotope(206, 205.9744653, 0.241)
        Pb%isot(3) = new_isotope(207, 206.9758969, 0.221)
        Pb%isot(4) = new_isotope(208, 207.9766521, 0.524)
    !_________________________________________________________________________
    Bi= new_element("Bi", 83, 1)
        Bi%isot(1) = new_isotope(209, 208.9803987, 1.0000)
    !_________________________________________________________________________
    Po= new_element("Po", 84, 2)
        Po%isot(1) = new_isotope(209, 208.9824304, 0.0)
        Po%isot(2) = new_isotope(210, 209.9828737, 0.0)
    !_________________________________________________________________________
    At= new_element("At", 85, 2)
        At%isot(1) = new_isotope(210, 209.987148,  0.0)
        At%isot(2) = new_isotope(211, 210.9874963, 0.0)
    !_________________________________________________________________________
    Rn= new_element("Rn", 86, 3)
        Rn%isot(1) = new_isotope(211, 210.990601,  0.0)
        Rn%isot(2) = new_isotope(220, 220.0113940, 0.0)
        Rn%isot(3) = new_isotope(222, 222.0175777, 0.0)
    !_________________________________________________________________________
    Fr= new_element("Fr", 87, 1)
        Fr%isot(1) = new_isotope(223, 223.0197359, 0.0)
    !_________________________________________________________________________
    Ra= new_element("Ra", 88, 4)
        Ra%isot(1) = new_isotope(223, 223.0185022, 0.0)
        Ra%isot(2) = new_isotope(224, 224.0202118, 0.0)
        Ra%isot(3) = new_isotope(226, 226.0254098, 0.0)
        Ra%isot(4) = new_isotope(228, 228.0310703, 0.0)
    !_________________________________________________________________________
    Ac= new_element("Ac", 89, 1)
        Ac%isot(1) = new_isotope(227, 227.0277521, 0.0)
    !_________________________________________________________________________
    Th= new_element("Th", 90, 2)
        Th%isot(1) = new_isotope(230, 230.0331338, 0.0)
        Th%isot(2) = new_isotope(232, 232.0380553, 1.0000)
    !_________________________________________________________________________
    Pa= new_element("Pa", 91, 1)
        Pa%isot(1) = new_isotope(231, 231.0358840, 1.0000)
    !_________________________________________________________________________
    U= new_element("U", 92, 5)
        U%isot(1) = new_isotope(233, 233.0396352, 0.0)
        U%isot(2) = new_isotope(234, 234.0409521, 0.000054)
        U%isot(3) = new_isotope(235, 235.0439299, 0.007204)
        U%isot(4) = new_isotope(236, 236.0455680, 0.0)
        U%isot(5) = new_isotope(238, 238.0507882, 0.992742)
    !_________________________________________________________________________
    Np= new_element("Np", 93, 2)
        Np%isot(1) = new_isotope(236, 236.046570,  0.0)
        Np%isot(2) = new_isotope(237, 237.0481734, 0.0)
    !_________________________________________________________________________
    Pu= new_element("Pu", 94, 6)
        Pu%isot(1) = new_isotope(238, 238.0495599, 0.0)
        Pu%isot(2) = new_isotope(239, 239.0521634, 0.0)
        Pu%isot(3) = new_isotope(240, 240.0538135, 0.0)
        Pu%isot(4) = new_isotope(241, 241.0568515, 0.0)
        Pu%isot(5) = new_isotope(242, 242.0587426, 0.0)
        Pu%isot(6) = new_isotope(244, 244.064204,  0.0)
    !_________________________________________________________________________
    Am= new_element("Am", 95, 2)
        Am%isot(1) = new_isotope(241, 241.0568291, 0.0)
        Am%isot(2) = new_isotope(243, 243.0613811, 0.0)
    !_________________________________________________________________________
    Cm= new_element("Cm", 96, 6)
        Cm%isot(1) = new_isotope(243, 243.0613891, 0.0)
        Cm%isot(2) = new_isotope(244, 244.0627526, 0.0)
        Cm%isot(3) = new_isotope(245, 245.0654912, 0.0)
        Cm%isot(4) = new_isotope(246, 246.0672237, 0.0)
        Cm%isot(5) = new_isotope(247, 247.070354,  0.0)
        Cm%isot(6) = new_isotope(248, 248.072349,  0.0)
    !_________________________________________________________________________
    Bk= new_element("Bk", 97, 2)
        Bk%isot(1) = new_isotope(247, 247.070307,  0.0)
        Bk%isot(2) = new_isotope(249, 249.0749867, 0.0)
    !_________________________________________________________________________
    Cf= new_element("Cf", 98, 4)
        Cf%isot(1) = new_isotope(249, 249.0748535, 0.0)
        Cf%isot(2) = new_isotope(250, 250.0764061, 0.0)
        Cf%isot(3) = new_isotope(251, 251.079587,  0.0)
        Cf%isot(4) = new_isotope(252, 252.081626,  0.0)
    !_________________________________________________________________________
    Es= new_element("Es", 99, 1)
        Es%isot(1) = new_isotope(252, 252.082980, 0.0)
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

END MODULE nist_isotopes_database


