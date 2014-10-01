!
! Phys. Rev. B 75, 174508 – Published 14 May 2007 
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013-2014 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
MODULE sqom_program
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info
  !
  TYPE sqom_input_type
    !
    CHARACTER(len=16) :: calculation ! lw=linewidth, spf=spectral function
    CHARACTER(len=256) :: outdir
    CHARACTER(len=256) :: e_file     ! put this in fron of file names
                                     ! must match the k-points grid 
    !
    LOGICAL :: elastic_peak
    !
    CHARACTER(len=256) :: file_mat2
    LOGICAL            :: asr2
    !
    ! for spectral function:
    REAL(DP) :: qq(3)
    REAL(DP) :: neutron_resolution
  END TYPE sqom_input_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(input, S, fc2)
    USE io_global,      ONLY : stdout
    USE q_grid,         ONLY : q_grid_type, setup_path, setup_simple_grid
    USE constants,      ONLY : RY_TO_CMM1
    USE more_constants, ONLY : INVALID
    USE wrappers,       ONLY : f_mkdir_safe
    !
    IMPLICIT NONE
    !
    TYPE(sqom_input_type),INTENT(out) :: input
    TYPE(forceconst2_grid),INTENT(out) :: fc2
    TYPE(ph_system_info),INTENT(out)   :: S    
    !
    ! Input variable, and defaul values:
    CHARACTER(len=16)   :: calculation = "neutron" ! "spf"
    CHARACTER(len=256)  :: file_mat2  = INVALID ! no default
    CHARACTER(len=256)  :: e_file     = INVALID ! default: don't use it
    REAL(DP)            :: qq(3)       = 0._dp
    !
    LOGICAL             :: asr2 = .false.
    LOGICAL             :: elastic_peak = .false.
    !
    INTEGER             :: ne = -1, nq = -1
    !
    REAL(DP)  :: neutron_resolution = 0.1_dp ! fraction of the energy
    !
    NAMELIST  / sqominput / &
      calculation, &
      file_mat2, asr2, &
      qq, e_file, &
      elastic_peak, neutron_resolution
      
    WRITE(*,*) "Waiting for input"
    !
    READ(*, sqominput)
    WRITE(*, sqominput)
    !
    IF(TRIM(file_mat2) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    
    input%file_mat2             = file_mat2
    input%asr2                  = asr2
    input%elastic_peak          = elastic_peak
    input%e_file                = e_file
    input%qq                    = qq
    input%neutron_resolution    = neutron_resolution
    !
    CALL READ_DATA(input, s, fc2)
    !
  END SUBROUTINE READ_INPUT
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_DATA(input, S, fc2)
    USE iso_c_binding,  ONLY : c_int
    USE input_fc,       ONLY : same_system, read_fc2, &
                               aux_system, div_mass_fc2
    USE asr2_module,    ONLY : impose_asr2
    USE io_global,      ONLY : stdout
    IMPLICIT NONE
    !
    TYPE(sqom_input_type),INTENT(in)        :: input
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    TYPE(ph_system_info),INTENT(inout)   :: S
    !
    INTEGER(kind=c_int) :: kb
    !
    CALL read_fc2(input%file_mat2, S,  fc2)
    !
    CALL aux_system(S)
    !
    CALL memstat(kb)
    WRITE(stdout,*) "Reading : done."
    WRITE(stdout,*) "Memory used : ", kb/1000, "Mb"
    !
    IF(input%asr2) CALL impose_asr2(S%nat, fc2)
    CALL div_mass_fc2(S, fc2)
    !
  END SUBROUTINE READ_DATA
  !
  REAL(DP) FUNCTION neutron_form_factor(aname, amass)
  USE constants, ONLY : BOHR_RADIUS_SI
  IMPLICIT NONE
    ! Data from http://www.ncnr.nist.gov/resources/n-lengths/
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
      REAL(DP)         :: coh_b
    END TYPE ffac
    !
    INTEGER :: i_iso, i_type
    INTEGER,PARAMETER :: n_iso=6, n_types=46
    TYPE(ffac) :: factors(0:10,100)
    !
    ! Initialize
    factors = ffac( "", 0._dp, 0._dp)
    ! Hydrogen:
    factors(0,1) = ffac( "H", 1.008_dp, -3.7390_dp)
    factors(1,1) = ffac( "H", 1._dp,-3.7406_dp)
    factors(2,1) = ffac( "D", 2._dp, 6.671_dp)
    factors(3,1) = ffac( "T", 3._dp, 4.792_dp)
    ! Palladium
    factors(0,46) = ffac( "Pd", 106.42_dp, 5.91_dp)
    factors(1,46) = ffac( "Pd", 102._dp, 7.7_dp)
    factors(2,46) = ffac( "Pd", 104._dp, 7.7_dp)
    factors(3,46) = ffac( "Pd", 105._dp, 5.5_dp)
    factors(4,46) = ffac( "Pd", 106._dp, 6.4_dp)
    factors(5,46) = ffac( "Pd", 108._dp, 4.1_dp)
    factors(6,46) = ffac( "Pd", 110._dp, 7.7_dp)
    !
    IF(.not. present(amass)) THEN
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
  SUBROUTINE neutron_function_convolution(xq, xg, T, S, fc2, ff, neutron_lw, ne, ee, spf)
    USE functions,      ONLY : f_bose, f_ngauss
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE linewidth,      ONLY : freq_phq_safe
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: xq(3), xg(3), T
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(in) :: ff(S%ntyp)
    REAL(DP),INTENT(in) :: neutron_lw ! as a fraction of the scattering energy
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ee(ne)
    REAL(DP),INTENT(inout) :: spf(ne,S%nat3)
    !
    COMPLEX(DP) :: cscal, cscal2, csum
    INTEGER     :: nu, ia, idx, ie, je
    !
    REAL(DP) :: freq(S%nat3)
    COMPLEX(DP) :: zz(S%nat3, S%nat3)
    !
    REAL(DP) :: G(S%nat3), F, Sq(ne,S%nat3)
    REAL(DP) :: de, sigma, sigma_p, sigma_m, norm, xqq(3)
    REAL(DP),PARAMETER :: FWHM_TO_SIGMA = 0.42466_dp !=( 2 *  sqrt(2._dp *log(2._dp)) )**-1
    !
    ! prepare phonon frequencies and patterns
    WRITE(*,*)
    WRITE(*,'(a,3f12.6)') "Convolution q-point: ",  xq
    WRITE(*,'(a,3f12.6)') "            g-point: ",  xg

    xqq = xq+xg
    !
    CALL freq_phq_safe(xqq, S, fc2, freq, zz)
!     zz=transpose(zz)
    !
    WRITE(*,'(a,99e15.3)') "ff: ", ff
    DO nu = 1,S%nat3
      WRITE(*,'(a,99(3(2f8.3,3x)6x))') "zz: ", zz(:,nu)
    ENDDO
    !
    DO nu=1,S%nat3
      G(nu)=0._dp
      csum =(0._dp,0._dp)
      DO ia=1,S%nat
          idx=(ia-1)*3
          ! 2 \pi i q \dot \tau
          cscal  = CMPLX(0._dp, tpi*SUM(xq*S%tau(:,ia)), kind=DP) 
          ! 2 \pi q \dot zz^{ia}_nu
          cscal2 = tpi* ( xqq(1)*zz(idx+1,nu) &
                         +xqq(2)*zz(idx+2,nu) &
                         +xqq(3)*zz(idx+3,nu) )
!           cscal2 = 1._dp
          ! \sum F_{ia} exp(\pi i q \dot \tau) 2 \pi q \dot zz^{ia}_nu / sqrt(M_{ia})
          csum = csum+ff(S%ityp(ia))*EXP(cscal)* cscal2 !/SQRT(S%amass(S%ityp(ia)))
      ENDDO
      G(nu) = ABS(csum)**2
    ENDDO
    ! Renormalize G to make things easier
    G= G/MAXVAL(G)
    WRITE(*,'(a,99e15.3)') "Renormalized G: ", G
    !
    !
    ! Add the elastic peak at zero
    DO nu = 1, S%nat3
    DO ie = 1,ne
      IF(ee(ie)/=0._dp) THEN
        ! I use spf*e which seems to be the correctly normalizable form
        ! I also divide by omega, the two effects should cancel out
        spf(ie,nu) = spf(ie,nu) * (f_bose(ee(ie),T)+1._dp) !/ ee(ie)
      ENDIF
    ENDDO
    ENDDO
!     spf = Sq
    !
    Sq = 0._dp
    ! Note: "de" is pointless as I'm going to renormalize later, 
    ! leaving it here for future developement
    de = 1._dp !/ REAL(ne, kind=DP)
    !
    DO nu = 1,S%nat3
    DO ie=1,ne
      !
      DO je=1,ne
        !
        sigma = MAX(ABS(neutron_lw*ee(je)), 1.e-3_dp/RYTOEV)*FWHM_TO_SIGMA
        F = f_ngauss(ee(ie)-ee(je),sigma)
        !
        ! I use spf*e which seems to be the correctly normalizable form
!         Sq(je,nu)=Sq(je,nu)+de*ee(ie)*spf(ie,nu)*F
        Sq(je,nu)=Sq(je,nu)+de*spf(ie,nu)*F
        !
      ENDDO
    ENDDO
    ENDDO
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
    norm = MAXVAL(abs(Sq))      
    Sq = Sq/norm
    !
    WRITE(999, '(2(a,3f12.6))') "# xq = ", xq, "      xg = ", xg
    DO ie = 1, ne
      WRITE(999,'(f12.6,100e20.10e4)') ee(ie)*RY_TO_CMM1, SUM(Sq(ie,:)), Sq(ie,:)
    ENDDO
    WRITE(999,*)
    
    !
  END SUBROUTINE neutron_function_convolution
  !
  SUBROUTINE scan_spf(filename, nq, ne, nat3)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: filename
    INTEGER,INTENT(out)  :: nq, ne
    INTEGER,INTENT(in)   :: nat3
    INTEGER,PARAMETER :: u = 1001
    !
    CHARACTER(len=4096) :: buffer
    INTEGER :: ios
    LOGICAL :: reading, reading_ee
    !
    CHARACTER(len=8) :: a1, a2
    INTEGER  :: iq, iq_, ie, j
    REAL(DP) :: r1, r2, r3, r4(nat3), xq(3)
    !
    OPEN(unit=u, file=filename, status="OLD", action="READ")
    WRITE(*,*) "Scanning..."
    !
    iq = 0
    ie = 0
    !
    reading=.true.
    reading_ee = .true.
    !
    READING_LOOP : &
    DO WHILE (reading)
      READ(u, '(a4096)', iostat=ios) buffer
!       WRITE(*,*) TRIM(buffer)
      IF(ios/=0) EXIT READING_LOOP
      !
      READ(buffer,*,iostat=ios) a1, a2
      IF(ios/=0) CYCLE READING_LOOP
      
      IF(a1=="#" .and. a2 == "xq" )THEN
        a1=""; a2=""
        iq = iq+1
        READ(buffer,*,iostat=ios) a1, a2, iq_, xq
        IF(ios/=0) CALL errore("read_spf", "error reading xq "//TRIM(buffer), 1)
        WRITE(*,'(a,i6,3f12.6)') "Found q-point: ",  iq, xq
        !
        DO WHILE(reading_ee)
          READ(u, '(a4096)', iostat=ios) buffer
          IF(ios/=0) EXIT READING_LOOP
          !
          READ(buffer,*,iostat=ios) r1, r2, r3, (r4(j), j =1,nat3)
          IF(ios/=0) THEN
            READ(buffer,*,iostat=ios) a1, a2
            IF(ios==0.and.a1=="#" .and. a2 == "xq" ) iq = iq+1
            reading_ee = .false.
            CYCLE READING_LOOP
          ELSE
            ie = ie +1
          ENDIF
        ENDDO
        !
      ENDIF
    ENDDO &
    READING_LOOP 
    !
    WRITE(*,*) "Found: ", iq, ie
    nq = iq
    ne = ie
    !
    CLOSE(u)
    IF(iq/=nq) CALL errore("read_spf", "missing q-points", 1)
    !
  END SUBROUTINE scan_spf
  !
  SUBROUTINE read_spf(filename, nq, ne, nat3, xq, ee, spf)
    USE constants,      ONLY : RY_TO_CMM1
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: filename
    INTEGER,INTENT(in)  :: nq, ne, nat3
    INTEGER,PARAMETER :: u = 1001
    REAL(DP),INTENT(out) :: xq(3,nq)
    REAL(DP),INTENT(out) :: spf(ne,nat3,nq)
    REAL(DP),INTENT(out) :: ee(ne)
    !
    CHARACTER(len=4096) :: buffer
    INTEGER :: ios
    LOGICAL :: reading
    !
    CHARACTER(len=8) :: a1, a2
    INTEGER  :: iq, iq_, ie, j
    REAL(DP) :: r1, r2
    !
    OPEN(unit=u, file=filename, status="OLD", action="READ")
    WRITE(*,*) "Going to read:", nq, ne, nat3
    !
    iq = 0
    !
    reading=.true.
    READING_LOOP : &
    DO WHILE (reading)
      READ(u, '(a4096)', iostat=ios) buffer
!       WRITE(*,*) TRIM(buffer)
      IF(ios/=0) EXIT READING_LOOP
      !
      READ(buffer,*,iostat=ios) a1, a2
      IF(ios/=0) CYCLE READING_LOOP
      
      IF(a1=="#" .and. a2 == "xq" )THEN
        a1=""; a2=""
        iq = iq+1
        READ(buffer,*,iostat=ios) a1, a2, iq_, xq(:,iq)
        IF(ios/=0) CALL errore("read_spf", "error reading xq "//TRIM(buffer), 1)
        WRITE(*,'(a,i6,3f12.6)') "Read q-point: ",  iq, xq(:,iq)
        !
        DO ie = 1,ne
          READ(u, '(a4096)', iostat=ios) buffer
          IF(ios/=0) CALL errore("read_spf", "error reading spf line", ie)
          READ(buffer,*,iostat=ios) r1, ee(ie), r2, (spf(ie,j,iq), j =1,nat3)
          IF(ios/=0) CALL errore("read_spf", "error reading from spf line "//TRIM(buffer), ie)
        ENDDO
        !
      ENDIF
    ENDDO &
    READING_LOOP 
    !
    CLOSE(u)
    IF(iq/=nq) CALL errore("read_spf", "missing q-points", 1)
    !
    ee = ee/RY_TO_CMM1
    !
  END SUBROUTINE read_spf
  !
END MODULE sqom_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM sqom

  USE kinds,            ONLY : DP
  USE sqom_program
  USE input_fc,         ONLY : print_citations_linewidth
  USE q_grid,           ONLY : q_grid_type !, setup_simple_grid
  USE parameters,       ONLY : ntypx
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(ph_system_info)   :: S
  TYPE(sqom_input_type)  :: sqominput
  !
  REAL(DP),ALLOCATABLE :: xq(:,:)
  REAL(DP),ALLOCATABLE :: spf(:,:,:)
  REAL(DP),ALLOCATABLE :: ee(:)
  INTEGER :: nq, ne, i
  !
  REAL(DP) :: ff(ntypx)
  !
  CALL print_citations_linewidth()
  
!   print*, neutron_form_factor("H")
!   print*, neutron_form_factor("D")
!   print*, neutron_form_factor("T",3._dp)
!   print*, neutron_form_factor("Pd")
!   print*, neutron_form_factor("Cu")

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT(sqominput, S, fc2)
  S%amass(2) = 2*S%amass(2)
  !
  CALL scan_spf(sqominput%e_file, nq, ne, S%nat3)
  ALLOCATE(xq(3,nq))
  ALLOCATE(spf(ne,S%nat3,nq))
  ALLOCATE(ee(ne))
  CALL read_spf(sqominput%e_file, nq, ne, S%nat3, xq, ee, spf)
  !
  ! ...
  S%atm(2)="D"
  DO i = 1, S%ntyp
    ff(i) = neutron_form_factor(S%atm(i))/S%celldm(1)
  ENDDO
  !
  DO i = 1, nq
    CALL neutron_function_convolution(xq(:,i), sqominput%qq, 80._dp, S, fc2, ff, &
                                      sqominput%neutron_resolution, ne, ee, spf(:,:,i))
  ENDDO
  !
END PROGRAM sqom
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













