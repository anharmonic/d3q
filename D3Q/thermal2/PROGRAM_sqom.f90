!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Convolutes a phonon spectral function with a certain experimental linewidth,
! which is assumed proportional to the energy (as in neutron scattering experiments)
! Reads files produced by lw.x 
!
! NOTE: heavily work in progress, not user friendly.
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
#include "mpi_thermal.h"
  !
  TYPE sqom_input_type
    !
    CHARACTER(len=16) :: calculation ! lw=linewidth, spf=spectral function
    CHARACTER(len=256) :: outdir
    CHARACTER(len=256) :: postfix
    CHARACTER(len=256),ALLOCATABLE :: e_file(:)
    INTEGER :: n_files = 1
    !
    LOGICAL :: elastic_peak
    REAL(DP) :: T 
    !
    CHARACTER(len=256) :: file_mat2
    LOGICAL            :: asr2
    !
    ! for spectral function:
    LOGICAL  :: convolution
    REAL(DP) :: qq(3)
    REAL(DP) :: neutron_res_frac
    REAL(DP) :: neutron_res0_cmm1 ! resolution at zero energy (cm-1)
  END TYPE sqom_input_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(input, S, fc2)
    USE q_grids,        ONLY : q_grid, setup_path, setup_simple_grid
    USE constants,      ONLY : RY_TO_CMM1
    USE more_constants, ONLY : INVALID
    USE wrappers,       ONLY : f_mkdir_safe
    USE code_input,     ONLY : parse_command_line
    !
    IMPLICIT NONE
    !
    TYPE(sqom_input_type),INTENT(out) :: input
    TYPE(forceconst2_grid),INTENT(out) :: fc2
    TYPE(ph_system_info),INTENT(out)   :: S    
    !
    ! Input variable, and defaul values:
    CHARACTER(len=16)   :: calculation = "neutron"
    CHARACTER(len=256)  :: outdir = "./" 
    CHARACTER(len=256)  :: postfix = INVALID, aux
    CHARACTER(len=256)  :: file_mat2  = INVALID ! no default
    !CHARACTER(len=256)  :: e_file(:)     = INVALID ! default: don't use it
    REAL(DP)            :: qq(3)       = 0._dp
    !
    LOGICAL             :: asr2 = .false.
    LOGICAL             :: elastic_peak = .false.
    LOGICAL             :: convolution = .true.
    !
    INTEGER             :: ne = -1, nq = -1, n_files = 1
    !
    REAL(DP)  :: neutron_res_frac = 0.025_dp ! fraction of the energy
    REAL(DP)  :: neutron_res0_cmm1 = 0.4 ! base resolution in cmm1
    REAL(DP)  :: T = 300._dp ! temperature for elastic peak
    INTEGER :: i, input_unit
    CHARACTER(len=256)  :: input_file
    CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER,EXTERNAL :: find_free_unit
    !
    NAMELIST  / sqominput / &
      calculation, outdir, &
      file_mat2, asr2, &
      n_files, &
      convolution, &
      neutron_res_frac, neutron_res0_cmm1, &
      qq, &
      elastic_peak, T
      
    WRITE(*,*) "Waiting for input"
    !
    input_file="input.SQOM"
    CALL parse_command_line(input_file)
    IF(TRIM(input_file)=="-")THEN
      ioWRITE(stdout,'(2x,3a)') "Warning! Reading standard input will probably not work with MPI"
      input_unit = 5
    ELSE
      ioWRITE(stdout,'(2x,3a)') "Reading input file '", TRIM(input_file), "'"
      input_unit = find_free_unit()
      OPEN(unit=input_unit, file=input_file, status="OLD", action="READ")
    ENDIF
    !
    READ(input_unit, sqominput)
    WRITE(stdout, sqominput)
    !
    !
    ALLOCATE(input%e_file(n_files))
    DO i = 1,n_files
      READ(input_unit,'(a256)') input%e_file(i)
    ENDDO
    CLOSE(input_unit)

    IF(TRIM(file_mat2) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    IF(TRIM(postfix) == INVALID) THEN
      postfix = "_"//TRIM(int_to_char(NINT(qq(1))))//&
                     TRIM(int_to_char(NINT(qq(2))))//&
                     TRIM(int_to_char(NINT(qq(3))))
    ENDIF

    input%calculation         = calculation
    input%postfix             = postfix
    input%outdir              = outdir
    input%file_mat2           = file_mat2
    input%asr2                = asr2
    input%n_files             = n_files
    input%elastic_peak        = elastic_peak
    input%T                   = T
    input%convolution         = convolution
    input%qq                  = qq
    input%neutron_res_frac    = neutron_res_frac
    input%neutron_res0_cmm1   = neutron_res0_cmm1
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
    ioWRITE(stdout,*) "Reading : done."
    ioWRITE(stdout,*) "Memory used : ", kb/1000, "Mb"
    !
    IF(input%asr2) CALL impose_asr2("simple",S%nat, fc2)
    CALL div_mass_fc2(S, fc2)
    !
  END SUBROUTINE READ_DATA
  !
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
!     factors(1,1) = ffac( "H", 1._dp,-3.7406_dp)
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
  SUBROUTINE neutron_function_convolution(w, xq, xg, do_elastic, T, S, fc2, ff, &
      res_frac, res0_cmm1, ne, ee, spf)
    USE functions,      ONLY : f_bose, f_ngauss
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE fc2_interpolate,     ONLY : freq_phq_safe
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: do_elastic
    REAL(DP),INTENT(in) :: w, xq(3), xg(3), T
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(in) :: ff(S%ntyp)
    REAL(DP),INTENT(in) :: res_frac  ! as a fraction of the scattering energy
    REAL(DP),INTENT(in) :: res0_cmm1 ! resolution at zero energy, in cmm1
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
    REAL(DP) :: sigma, sigma_p, sigma_m, norm, xqq(3)
    REAL(DP) :: dee(ne), peak0(ne)
    REAL(DP),PARAMETER :: FWHM_TO_SIGMA = 0.42466_dp !=( 2 *  sqrt(2._dp *log(2._dp)) )**-1
    !
    ! prepare phonon frequencies and patterns
!     WRITE(*,*)
!     WRITE(*,'(a,3f12.6)') "Convolution q-point: ",  xq
!     WRITE(*,'(a,3f12.6)') "            g-point: ",  xg

    xqq = xq+xg
    !
    CALL freq_phq_safe(xqq, S, fc2, freq, zz)
    !
    ! Put spf in its normalized form
    ! (currently it is in 1/cmm1^2, we get it to 1/cmm1 and convert those to Ry)
    ! there is also a factor pi/2 that's lost somewhere in lw.x, fix it in the future
    DO nu=1,S%nat3
      spf(:,nu) = spf(:,nu)*ee*RY_TO_CMM1*2/pi
    ENDDO

    ! Compute the de factor in a more general form than really necessary
    ! for Simpson integration
    dee(1:ne-1) = ee(2:ne)-ee(1:ne-1)
    dee(ne) = dee(ne-1)
    !dee(1) = dee(1)
    DO nu=1,S%nat3
      CALL simpson(ne, spf(:,nu), dee, norm)
!       print*, "norm", nu, norm
    ENDDO
    !
    !print*, ff
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
    !G = G/SUM(G)
    ! Renormalize G makes things easier but is not consistent 
    ! among different q points (the inconsistency is not huge)
    !G= G/MAXVAL(G)
    ! norm = DSQRT(SUM(G**2))
!     IF(norm/=0._dp) THEN
!       G = G/norm
!     ENDIF
!    WRITE(*,'(a,3f12.6,3x,99e15.6)') "xq+Q/G", xqq, DSQRT(SUM(G**2)), SUM(G)
    !
    !
    ! Add the elastic peak at zero
    IF(do_elastic)THEN
      peak0 = 0._dp
      DO ie = 1,ne
      IF(ee(ie)/=0._dp) THEN
        peak0(ie) = f_bose(ee(ie),T)/ee(ie)
      ENDIF
!       write(999, '(2e25.15)') ee(ie), peak0(ie)
      ENDDO
!       CALL simpson(ne, peak0, dee, norm)
!       print*, "norm0", norm
      
      DO nu = 1, S%nat3
        spf(:,nu) = spf(:,nu) + peak0(:)/S%nat3
      ENDDO
!       DO nu=1,S%nat3
!         CALL simpson(ne, spf(:,nu), dee, norm)
!         print*, "norm2", nu, norm
!       ENDDO
    ENDIF
!     spf = Sq
    !
    !
    Sq = 0._dp
    !
    ! Doing the convolution with a gaussian in the least efficient way possible
    DO ie=1,ne
      !
      DO je=1,ne
        !
        DO nu = 1,S%nat3
          sigma = MAX(ABS(res_frac*ee(je)), res0_cmm1/RY_TO_CMM1)*FWHM_TO_SIGMA
          F = f_ngauss(ee(ie)-ee(je),sigma)
          !
          ! I use spf*e which seems to be the correctly normalizable form
  !         Sq(je,nu)=Sq(je,nu)+de*ee(ie)*spf(ie,nu)*F
          Sq(je,nu)=Sq(je,nu)+dee(ie)*spf(ie,nu)*F
          !
        ENDDO
      ENDDO
    ENDDO
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
        READ(u,*)
        READ(u,*)
        a1=""; a2=""
        iq = iq+1
        READ(buffer,*,iostat=ios) a1, a2, iq_, xq
        IF(ios/=0) CALL errore("read_spf", "error reading xq "//TRIM(buffer), 1)
!         WRITE(*,'(a,i6,3f12.6)') "Found q-point: ",  iq, xq
        !
        DO WHILE(reading_ee)
          READ(u, '(a4096)', iostat=ios) buffer
          IF(ios/=0) EXIT READING_LOOP
          !
          READ(buffer,*,iostat=ios) r1, r2, r3, (r4(j), j =1,nat3)
          IF(ios/=0) THEN
            READ(buffer,*,iostat=ios) a1, a2
            IF(ios==0.and.a1=="#" .and. a2 == "xq" ) THEN
              READ(u,*)
              READ(u,*)
              iq = iq+1
            ENDIF
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
  SUBROUTINE read_spf(filename, nq, ne, nat3, xq, w, ee, spf)
    USE constants,      ONLY : RY_TO_CMM1
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: filename
    INTEGER,INTENT(in)  :: nq, ne, nat3
    INTEGER,PARAMETER :: u = 1001
    REAL(DP),INTENT(out) :: xq(3,nq)
    REAL(DP),INTENT(out) :: w(nq) !path length
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
    WRITE(*,'(x,3a)') "Processing: '", TRIM(filename),"'"
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
        READ(u,*)
        READ(u,*)
        a1=""; a2=""
        iq = iq+1
        READ(buffer,*,iostat=ios) a1, a2, iq_, xq(:,iq)
        !IF(iq_/=iq) CALL errore("read_spf", "unexpected iq", 1)
        IF(ios/=0) CALL errore("read_spf", "error reading xq "//TRIM(buffer), 1)
!         WRITE(*,'(a,i6,3f12.6)') "Read q-point: ",  iq, xq(:,iq)
        !
        DO ie = 1,ne
          READ(u, '(a4096)', iostat=ios) buffer
          IF(ios/=0) CALL errore("read_spf", "error reading spf line", ie)
          READ(buffer,*,iostat=ios) w(iq), ee(ie), r2, (spf(ie,j,iq), j =1,nat3)
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
  USE more_constants,   ONLY : print_citations_linewidth
  USE q_grids,          ONLY : q_grid !, setup_simple_grid
  USE parameters,       ONLY : ntypx
  USE constants,        ONLY : RY_TO_CMM1
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi, mpi_bsum, &
                               mpi_broadcast, &
                               ionode, num_procs, my_id
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(ph_system_info)   :: S
  TYPE(sqom_input_type)  :: sqi
  !
  REAL(DP),ALLOCATABLE :: xq(:,:), w(:)
  REAL(DP),ALLOCATABLE :: spf(:,:,:), spf_merge(:,:,:)
  REAL(DP),ALLOCATABLE :: ee(:)
  INTEGER :: nq, ne, iq, ie, ifile, it
  !
  REAL(DP) :: ff(ntypx)
  !
  
!   print*, neutron_form_factor("H")
!   print*, neutron_form_factor("D")
!   print*, neutron_form_factor("T",3._dp)
!   print*, neutron_form_factor("Pd")
!   print*, neutron_form_factor("Cu")
  CALL start_mpi()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT(sqi, S, fc2)
  S%amass(2) = 2*S%amass(2)
  !
  IF(ionode)THEN
    CALL scan_spf(sqi%e_file(1), nq, ne, S%nat3)
    DO it = 1, S%ntyp
      ff(it) = neutron_form_factor(S%atm(it))/S%celldm(1)
    ENDDO
    !ff=ff/MAXVAL(ff)
  ENDIF
  CALL mpi_broadcast(nq)
  CALL mpi_broadcast(ne)
  CALL mpi_broadcast(S%ntyp, ff)
  
  ALLOCATE(xq(3,nq))
  ALLOCATE(w(nq))
  ALLOCATE(spf(ne,S%nat3,nq))
  ALLOCATE(ee(ne))
  !
  !
  ALLOCATE(spf_merge(ne,S%nat3,nq))
  spf_merge = 0._dp
  !
  !DO ifile = 1,sqi%n_files
  DO ifile = 1+my_id, sqi%n_files, num_procs
    !
    CALL read_spf(sqi%e_file(ifile), nq, ne, S%nat3, xq, w, ee, spf)
    !
    IF( sqi%convolution ) THEN
!       IF(ifile==1+my_id)THEN
!       !S%atm(2)="D"
!         DO it = 1, S%ntyp
!           ff(it) = neutron_form_factor(S%atm(it))/S%celldm(1)
!         ENDDO        
!       ENDIF
      DO iq = 1, nq
        CALL neutron_function_convolution(w(iq), xq(:,iq), sqi%qq, &
        sqi%elastic_peak, sqi%T, S, fc2, ff, &
        sqi%neutron_res_frac, sqi%neutron_res0_cmm1, &
        ne, ee, spf(:,:,iq))
      ENDDO
    ENDIF
    !
    spf_merge = spf_merge + spf
    !
  ENDDO
  !
  CALL mpi_bsum(ne,S%nat3,nq,spf_merge)
  spf_merge = spf_merge / DBLE(sqi%n_files)
  !
  IF(ionode)THEN
    OPEN(998, file=TRIM(sqi%e_file(1))//trim(sqi%postfix))
    DO iq = 1, nq
      WRITE(998, '(2(a,3f12.6))') "# xq = ", xq(:,iq)
      WRITE(998, '("#")')
      WRITE(998, '("#")')
      DO ie = 1, ne
        WRITE(998,'(2f12.6,100e20.10e4)') w(iq), ee(ie)*RY_TO_CMM1,&
                SUM(spf_merge(ie,:,iq)), spf_merge(ie,:,iq)
      ENDDO
      WRITE(998,*)
    ENDDO
    CLOSE(998)
  ENDIF
  !
  CALL stop_mpi()
  !
  CALL print_citations_linewidth()
  !
END PROGRAM sqom
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













