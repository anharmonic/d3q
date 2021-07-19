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
    !LOGICAL :: elastic_peak
    REAL(DP) :: elastic_peak_norm  = 0._dp
    REAL(DP) :: elastic_peak_width = 0._dp
    REAL(DP) :: T 
    !
    CHARACTER(len=256) :: file_mat2
    LOGICAL            :: asr2
    !
    ! for spectral function:
    LOGICAL  :: convolution
    REAL(DP) :: qq(3)
    REAL(DP) :: res_frac
    REAL(DP) :: res0_cmm1 ! resolution at zero energy (cm-1)
  END TYPE sqom_input_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT_SQOM(input, S, fc2)
    USE q_grids,        ONLY : q_grid, setup_path, setup_simple_grid
    USE constants,      ONLY : RY_TO_CMM1
    USE more_constants, ONLY : INVALID, write_conf
    USE clib_wrappers,       ONLY : f_mkdir_safe
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
    CHARACTER(len=256)  :: postfix = INVALID    ! add to the name of the first file to create out file
    CHARACTER(len=256)  :: file_mat2  = INVALID ! no default
    !CHARACTER(len=256)  :: e_file(:)     = INVALID ! default: don't use it
    REAL(DP)            :: qq(3)       = 0._dp
    !
    LOGICAL             :: asr2 = .false.
    !LOGICAL             :: elastic_peak = .false.
    REAL(DP)            :: elastic_peak_norm = 0._dp
    REAL(DP)            :: elastic_peak_width= 0._dp
    LOGICAL             :: convolution = .true.
    !
    INTEGER             :: ne = -1, nq = -1, n_files = 1
    !
    REAL(DP)  :: res_frac = 0.0_dp ! fraction of the energy
    REAL(DP)  :: res0_cmm1 = 0._dp ! base resolution in cmm1
    REAL(DP)  :: T = 300._dp ! temperature for elastic peak
    INTEGER :: i, input_unit
    CHARACTER(len=256)  :: input_file
    CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER,EXTERNAL :: find_free_unit
    !
    NAMELIST  / sqominput / &
      calculation, outdir, postfix, &
      file_mat2, asr2, &
      n_files, &
      convolution, &
      res_frac, res0_cmm1, &
      qq, &
      elastic_peak_norm, elastic_peak_width, T
      
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

    IF(TRIM(file_mat2) == INVALID) CALL errore('READ_INPUT_SQOM', 'Missing file_mat2', 1)

    IF(calculation=="neutron") THEN
      IF(TRIM(postfix) == INVALID) THEN
!         postfix = TRIM(int_to_char(NINT(qq(1))))//&
!                   TRIM(int_to_char(NINT(qq(2))))//&
!                   TRIM(int_to_char(NINT(qq(3))))
        postfix=TRIM(write_conf(1, 3, qq))//"_"&
              //TRIM(write_conf(2, 3, qq))//"_"&
              //TRIM(write_conf(3, 3, qq))
      ENDIF
      
    
    ELSE IF (calculation=="ixs") THEN
      IF(TRIM(postfix) == INVALID) postfix="ixs"
    ELSE
      CALL errore("read_input", 'unknown calculation', 1)
    ENDIF

    input%calculation         = calculation
    input%postfix             = postfix
    input%outdir              = outdir
    input%file_mat2           = file_mat2
    input%asr2                = asr2
    !
    input%n_files             = n_files
    input%elastic_peak_norm   = elastic_peak_norm
    input%elastic_peak_width  = elastic_peak_width
    input%T                   = T
    input%convolution         = convolution
    input%res_frac    = res_frac
    input%res0_cmm1   = res0_cmm1
    !
    CALL READ_DATA(input, s, fc2)
    CALL cryst_to_cart(1,qq,S%bg, +1)
    input%qq                  = qq
    WRITE(*,'(a,3f12.6)') "qq in 2pi/alat: ", qq
    !
  END SUBROUTINE READ_INPUT_SQOM
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_DATA(input, S, fc2)
    USE input_fc,       ONLY : same_system, read_fc2, &
                               aux_system, div_mass_fc2
    USE asr2_module,    ONLY : impose_asr2
    USE clib_wrappers,        ONLY : memstat
    IMPLICIT NONE
    !
    TYPE(sqom_input_type),INTENT(in)        :: input
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    TYPE(ph_system_info),INTENT(inout)   :: S
    !
    INTEGER :: kb
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
  !
  ! Compute convoluting function
  SUBROUTINE ixs_function_convolution_fft(S, res_cmm1, ne, ee, spf)
    USE fft_scalar,     ONLY : cft_1z
    USE functions,      ONLY : f_ngauss, f_psvoigt
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE fc2_interpolate,     ONLY : freq_phq_safe
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: res_cmm1
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ee(ne)
    REAL(DP),INTENT(inout) :: spf(ne,S%nat3)
    !
    INTEGER     :: nu, ie, je, je_range
    !
    REAL(DP) :: F, Sq(ne,S%nat3)
    COMPLEX(DP) :: aux(2*ne), auxg(2*ne)
    COMPLEX(DP), SAVE, ALLOCATABLE :: psvoigtg(:)
    LOGICAL, SAVE :: first = .true.
    REAL(DP) :: sigma, fwhm, e_avg
    REAL(DP) :: dee
    !
    dee = (ee(ne)-ee(1))/(ne-1)
    !
    IF(first)THEN
      first = .false.
      ALLOCATE(psvoigtg(2*ne))
      fwhm = res_cmm1/RY_TO_CMM1
      e_avg = (ee(ne)+ee(1))*0.5_dp
      DO ie = 1,ne
        aux(ie) = f_psvoigt(ee(ie)-e_avg,fwhm,0.3_dp)
      ENDDO
      aux(ne+1) = 0._dp
      CALL cft_1z(aux, 1, 2*ne, 1, -1, psvoigtg)
      psvoigtg = 2*ne*dee*psvoigtg
    ENDIF
    !
    DO nu = 1,S%nat3
      aux(1:ne) = spf(1:ne,nu)
      aux(ne+1:2*ne) = 0._dp
      CALL cft_1z(aux,  1, 2*ne, 1, -1, auxg)
      auxg = auxg*psvoigtg
      CALL cft_1z(auxg, 1, 2*ne, 1, +1, aux)
      Sq(1:ne,nu)  = DBLE(aux(ne/2:ne/2+ne-1))
      !Sq(ne/2:ne,nu) = dee*DBLE(aux(1:ne/2))
      !Sq(:,nu) = DBLE(aux)
    ENDDO
    !
    ! Put back and eventually add the weight
    spf = sq
    !
  END SUBROUTINE ixs_function_convolution_fft
  
  ! Compute convoluting function
  SUBROUTINE ixs_function_convolution(S, res_cmm1, ne, ee, spf)
    USE fft_scalar,     ONLY : cft_1z
    USE functions,      ONLY : f_psvoigt
    USE constants,      ONLY : tpi, RY_TO_CMM1, pi, K_BOLTZMANN_RY, RYTOEV
    USE fc2_interpolate,ONLY : freq_phq_safe
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: res_cmm1
    INTEGER,INTENT(in)  :: ne
    REAL(DP),INTENT(in) :: ee(ne)
    REAL(DP),INTENT(inout) :: spf(ne,S%nat3)
    !
    INTEGER  :: nu, ie, je, je_range_dw, je_range_up
    REAL(DP) :: F, dee, fwhm
    !
    COMPLEX(DP), SAVE, ALLOCATABLE :: psvoigt(:)
    INTEGER, SAVE :: je_range = -1
    LOGICAL, SAVE :: first = .true.
    !
    REAL(DP) :: Sq(ne,S%nat3)
    !
    dee = (ee(ne)-ee(1))/(ne-1)
    !
    IF(first)THEN
      first = .false.
      fwhm = res_cmm1/RY_TO_CMM1
      je_range = 30*fwhm/dee
      ALLOCATE(psvoigt(-je_range:je_range))
      !
      DO je = -je_range, je_range
        psvoigt(je) = dee*f_psvoigt(dee*je,fwhm,0.3_dp)
      ENDDO
      !
    ENDIF
    !
    Sq = 0._dp
    !
    ! Doing the convolution with a gaussian in the least efficient way possible
    DO nu = 1,S%nat3
      DO ie=1,ne
        je_range_dw = MIN(je_range, ie-1)
        je_range_up = MIN(je_range, ne-ie)
        DO je=-je_range_dw, je_range_up
          !
          ! I use spf*e which seems to be the correctly normalizable form
          Sq(ie,nu)=Sq(ie,nu)+spf(ie+je,nu)*psvoigt(je)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    ! Put back and eventually add the weight
    spf = sq
    !
  END SUBROUTINE ixs_function_convolution
  
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
!         WRITE(*,'(a,i6,3f12.6)') "Found q-point: ",  iq, xq
        !
        DO WHILE(reading_ee)
          SKIP_COMMENTS : DO
            READ(u, '(a4096)', iostat=ios) buffer
            IF(buffer(1:1) /= "#") EXIT SKIP_COMMENTS
          ENDDO SKIP_COMMENTS
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
          !READ(u, '(a4096)', iostat=ios) buffer
          SKIP_COMMENTS : DO
            READ(u, '(a4096)', iostat=ios) buffer
            IF(buffer(1:1) /= "#") EXIT SKIP_COMMENTS
          ENDDO SKIP_COMMENTS
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
  USE neutrons,         ONLY : neutron_function_convolution
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
  REAL(DP) :: ff(ntypx), norm_g
  !
  
!   print*, neutron_form_factor("H")
!   print*, neutron_form_factor("D")
!   print*, neutron_form_factor("T",3._dp)
!   print*, neutron_form_factor("Pd")
!   print*, neutron_form_factor("Cu")
  CALL start_mpi()

  ! READ_INPUT_SQOM also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT_SQOM(sqi, S, fc2)
  S%amass(2) = 2*S%amass(2)
  !
  IF(ionode)THEN
    CALL scan_spf(sqi%e_file(1), nq, ne, S%nat3)
!     IF(sqi%calculation == "neutron")THEN
!       DO it = 1, S%ntyp
!         ff(it) = neutron_form_factor(S%atm(it))/S%celldm(1)
!       ENDDO
!     ENDIF
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
    CALCULATION_TYPE : &
    IF(sqi%calculation == "neutron")THEN
      IF( sqi%convolution ) THEN
  !       IF(ifile==1+my_id)THEN
  !       !S%atm(2)="D"
  !         DO it = 1, S%ntyp
  !           ff(it) = neutron_form_factor(S%atm(it))/S%celldm(1)
  !         ENDDO        
  !       ENDIF
        DO iq = 1, nq
          CALL neutron_function_convolution(w(iq), xq(:,iq), sqi%qq, &
          sqi%elastic_peak_norm, sqi%elastic_Peak_width/RY_TO_CMM1, sqi%T,&
          S, fc2, sqi%res_frac, sqi%res0_cmm1, &
          ne, ee, spf(:,:,iq))
        ENDDO
      ENDIF
      !
    ELSE IF(sqi%calculation == "ixs")THEN
      DO iq = 1, nq
      !(w, xq, S, res_cmm1, ne, ee, spf)
        !CALL ixs_function_convolution(S,sqi%res0_cmm1, ne, ee, spf(:,:,iq))
        CALL ixs_function_convolution_fft(S,sqi%res0_cmm1, ne, ee, spf(:,:,iq))
      ENDDO
      !
    ELSE
      CALL errore("sqom", "not a valid type of calculation: "//TRIM(sqi%calculation),1)
    ENDIF CALCULATION_TYPE
    !
    spf_merge = spf_merge + spf
    !
  ENDDO
  !
  CALL mpi_bsum(ne,S%nat3,nq,spf_merge)
  spf_merge = spf_merge / DBLE(sqi%n_files)
  !
  norm_g = DSQRT(SUM(sqi%qq**2))
  IF(ionode)THEN
    OPEN(998, file=TRIM(sqi%e_file(1))//"_"//trim(sqi%postfix))
    DO iq = 1, nq
      WRITE(998, '(2(a,3f12.6))') "# xq = ", xq(:,iq)
      WRITE(998, '("#")')
      WRITE(998, '("#")')
      DO ie = 1, ne
        WRITE(998,'(2f12.6,100e20.10e3)') w(iq)+norm_g, ee(ie)*RY_TO_CMM1,&
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













