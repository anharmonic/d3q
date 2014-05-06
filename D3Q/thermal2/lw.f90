!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013-2014 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
MODULE more_constants
  USE kinds, ONLY : DP
  REAL(DP),PARAMETER :: RY_TO_JOULE =  0.5* 4.35974394e-18
  REAL(DP),PARAMETER :: RY_TO_SECOND = 2* 2.418884326505e-17
  REAL(DP),PARAMETER :: RY_TO_METER = 5.2917721092e-11
  REAL(DP),PARAMETER :: RY_TO_WATTMM1KM1 = RY_TO_JOULE / (RY_TO_SECOND * RY_TO_METER)
END MODULE more_constants

MODULE linewidth_program
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       forceconst3_grid, &
                       ph_system_info
                       
  REAL(DP) :: default_sigma = 10._dp
  
  TYPE lwinput_type
    !
    CHARACTER(len=16) :: calculation ! lw=linewidth, spf=spectral function
    CHARACTER(len=16) :: mode        ! "full" or "simple" spectral function 
    CHARACTER(len=256) :: outdir
    !
    LOGICAL :: exp_t_factor
    !
    CHARACTER(len=256) :: file_mat3
    CHARACTER(len=256) :: file_mat2
    LOGICAL            :: asr2
    !
    INTEGER            :: nconf
    REAL(DP),ALLOCATABLE :: T(:), sigma(:)
    ! for spectral function:
    INTEGER :: ne
    REAL(DP) :: e0, de
    !
    INTEGER :: nk(3)
  END TYPE lwinput_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(input, qpath, S, fc2, fc3)
    USE io_global,      ONLY : stdout
    USE q_grid,         ONLY : q_grid_type, setup_path
    USE constants,      ONLY : RY_TO_CMM1
    USE wrappers,       ONLY : f_mkdir_safe
    !
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(out) :: input
    TYPE(q_grid_type),INTENT(out)  :: qpath
    TYPE(forceconst2_grid),INTENT(out) :: fc2
    TYPE(forceconst3_grid),INTENT(out) :: fc3
    TYPE(ph_system_info),INTENT(out)   :: S    
    !
    ! Input variable, and defaul values:
    CHARACTER(len=3),PARAMETER :: INVALID = '///'
    CHARACTER(len=16)  ::  calculation = "lw" ! "spf"
    CHARACTER(len=256) :: file_mat3 = INVALID
    CHARACTER(len=256) :: file_mat2 = INVALID
    CHARACTER(len=256) :: outdir = './'
    LOGICAL            :: asr2 = .false.
    INTEGER            :: nconf = 1
    INTEGER            :: nq = -1
    INTEGER            :: nk(3) = (/-1, -1, -1/)
    LOGICAL            :: exp_t_factor = .false.
    !
    INTEGER  :: ne = -1
    REAL(DP) :: de = 1._dp, e0 = 0._dp
    !
    REAL(DP) :: xq(3), xq0(3)
    INTEGER  :: ios, ios2, i, naux
    CHARACTER(len=1024) :: line, word
    CHARACTER(len=16)   :: word2, word3
    CHARACTER(LEN=256), EXTERNAL :: TRIMCHECK
    !
    NAMELIST  / lwinput / &
      calculation, outdir, &
      file_mat2, file_mat3, asr2, &
      nconf, nq, nk, &
      ne, de, e0, exp_t_factor

    WRITE(*,*) "Waiting for input"
    !
    READ(*, lwinput)
    WRITE(*, lwinput)
    !
    IF(TRIM(file_mat2) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    IF(TRIM(file_mat3) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat3', 1)
    IF(ANY(nk<0)) CALL errore('READ_INPUT', 'Missing nk', 1)    
    IF(nq<0)      CALL errore('READ_INPUT', 'Missing nq', 1)    
    IF(nconf<0)   CALL errore('READ_INPUT', 'Missing nconf', 1)    
    
    input%file_mat2 = file_mat2
    input%file_mat3 = file_mat3
    input%outdir    = TRIMCHECK(outdir)
    input%asr2      = asr2
    input%nconf     = nconf
    input%nk        = nk
    input%exp_t_factor = exp_t_factor
    !
    ios = f_mkdir_safe(input%outdir)
    IF(ios>0) CALL errore('READ_INPUT', 'cannot create directory: "'//TRIM(input%outdir)//'"',1)
    !
    CALL READ_DATA (input, s, fc2, fc3)
    !
    READ(calculation,*,iostat=ios) input%calculation, input%mode
    IF(ios/=0) THEN
      input%calculation = calculation
      input%mode = "full"
    ENDIF
    !
    IF(TRIM(input%calculation) == 'spf' .and. ne < 0) &
      CALL errore('READ_INPUT', 'Missing ne for spf calculation', 1)    
    input%ne = ne
    input%de = de
    input%e0 = e0
    !
    ALLOCATE(input%T(nconf), input%sigma(nconf))
    !
    READ(*,'(a1024)', iostat=ios) line
    READ_CARDS : &
    DO WHILE (ios==0)
      !
      READ(line,*,iostat=ios2) word
      IF(ios2/=0) THEN
        READ(*,'(a1024)', iostat=ios) line
        CYCLE READ_CARDS 
      ENDIF
      !
      SELECT CASE (TRIM(word))
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE ("QPOINTS")
        WRITE(*,*) "Reading QPOINTS"
        !
        qpath%basis = 'cartesian'
        READ(line,*,iostat=ios) word2, word3
        IF(ios==0) qpath%basis = TRIM(word3(1:9))
        !
        QPOINT_LOOP : &
        DO i = 1, nq
          READ(*,'(a1024)', iostat=ios) line
          IF(ios/=0) CALL errore("READ_INPUT","Expecting q point: input error.", 1)

          READ(line,*, iostat=ios2) xq(1), xq(2), xq(3), naux
          IF(ios2==0) THEN
            IF(i>1) THEN
              xq0 = qpath%xq(:,qpath%nq)
            ELSE
              xq0 = 0._dp
            ENDIF
            IF(SUM((xq-xq0)**2)<1.d-6) CALL errore("READ_INPUT", "q-point repeated, please think for the planet.", 1)
            CALL setup_path(xq0, xq, naux, qpath)
            CYCLE QPOINT_LOOP
          ENDIF

          READ(line,*, iostat=ios2) xq(1), xq(2), xq(3)
          IF(ios2==0) THEN
            CALL setup_path(xq, xq, 1, qpath)
            CYCLE QPOINT_LOOP
          ENDIF
          
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting q point, got: '"//TRIM(line)//"'.", 1)
          EXIT QPOINT_LOOP
          
        ENDDO &
        QPOINT_LOOP
        
        WRITE(*,"(2x,a,i4,a)") "Read ", qpath%nq, " q-points, "//TRIM(qpath%basis)//" basis"
        !
        IF(TRIM(qpath%basis) == "crystal")  THEN
          CALL cryst_to_cart(qpath%nq,qpath%xq,S%bg, +1)
          qpath%basis = "cartesian"
          WRITE(*,"(4x,a,i4,a)") "q-points converted to cartesian basis"
        ENDIF
        WRITE(*,"(2x,3f12.6)") qpath%xq
        WRITE(*,*)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE ("CONFIGS")
        WRITE(*,*) "Reading CONFIGS", nconf
        DO i = 1,nconf
          READ(*,'(a1024)', iostat=ios2) line
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration: input error.", 1)
          READ(line,*,iostat=ios2) input%sigma(i), input%T(i)
          IF(ios2/=0) THEN
            READ(line,*,iostat=ios2) input%T(i)
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration, got: '"//TRIM(line)//"'.", 1)
            IF(i>1) THEN
              input%sigma(i) = input%sigma(i-1)
            ELSE
              input%sigma = default_sigma
            ENDIF
          ENDIF
        ENDDO
        WRITE(*,'(2x,a,/,100(8f9.1,/))') "Temperatures", input%T
        WRITE(*,'(2x,a,/,100(8f9.1,/))') "Sigma       ", input%sigma
        WRITE(*,*)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE DEFAULT
!         WRITE(*,*) "Skip:", TRIM(line)
      END SELECT

    READ(*,'(a1024)', iostat=ios) line
    word = ''
    ENDDO &
    READ_CARDS
    !
    !
  END SUBROUTINE READ_INPUT
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_DATA(input, s, fc2, fc3)
    USE iso_c_binding,  ONLY : c_int
    USE input_fc,       ONLY : same_system, read_fc2, read_fc3, &
                               aux_system, div_mass_fc2, div_mass_fc3
    USE asr2_module,    ONLY : impose_asr2
    USE io_global,      ONLY : stdout
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)        :: input
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    TYPE(forceconst3_grid),INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(inout)   :: s
    TYPE(ph_system_info) :: s3
    !
    INTEGER(kind=c_int) :: kb
    !
    CALL read_fc2(input%file_mat2, S,  fc2)
    CALL read_fc3(input%file_mat3, S3, fc3)
    IF(.not.same_system(S, S3)) THEN
      PRINT*, "WARNING! FC2 and FC3 systems DO NOT MATCH !!!"
    ENDIF
    !
    CALL aux_system(S)
    !
    CALL memstat(kb)
    WRITE(stdout,*) "Reading : done."
    WRITE(stdout,*) "Memory used : ", kb/1000, "Mb"
    !
    IF(input%asr2) CALL impose_asr2(S%nat, fc2)
    CALL div_mass_fc2(S, fc2)
    CALL div_mass_fc3(S, fc3)
    !
  END SUBROUTINE READ_DATA
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE LW_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE interp_fc,   ONLY : fftinterp_mat2, mat2_diag
    USE linewidth,   ONLY : linewidth_q, selfnrg_q, spectre_q, freq_phq
    USE constants,   ONLY : RY_TO_CMM1
!     USE ph_velocity
    USE q_grid,      ONLY : q_grid_type, setup_simple_grid
    USE functions,   ONLY : write_temperature
    USE nanoclock
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qpath
    !
    COMPLEX(DP),ALLOCATABLE :: D(:,:)
    REAL(DP) :: w2(S%nat3)
    REAL(DP) :: pl,dpl
    INTEGER :: i, it
    TYPE(q_grid_type) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP)   :: lw(S%nat3,input%nconf)
    !
    ALLOCATE(D(S%nat3, S%nat3))
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                                "lw.T"//TRIM(write_temperature(it,input%nconf,input%T))//&
                                "s"//TRIM(write_temperature(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
    ENDDO
    CALL print_header()
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
    DO i = 1,qpath%nq
      WRITE(*,'(i6,3f15.8)') i, qpath%xq(:,i)
      !
      CALL freq_phq(qpath%xq(:,i), S, fc2, w2, D)
        !
        ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
        ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
        ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
        ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
        !  => 2d = 2 sqrt(log(2) c => d = sqrt(log(2)) d = 0.83255 c
        !input%sigma = 0.83255 *input%sigma/RY_TO_CMM1
      IF (TRIM(input%mode) == "full") THEN
        ls = selfnrg_q(qpath%xq(:,i), input%nconf, input%T, input%sigma/RY_TO_CMM1, &
                        S, grid, fc2, fc3)
        DO it = 1,input%nconf
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,6f12.6,6x,6e15.5,6x,6e15.5,6x,6e15.5)') &
                i,pl,qpath%xq(:,i), w2*RY_TO_CMM1, -DIMAG(ls(:,it))*RY_TO_CMM1, DBLE(ls(:,it))*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE IF (TRIM(input%mode) == "real") THEN
        lw = linewidth_q(qpath%xq(:,i), input%nconf, input%T, input%sigma/RY_TO_CMM1, &
                        S, grid, fc2, fc3)
        DO it = 1,input%nconf
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,6f12.6,6x,6e15.5,6x,6e15.5,6x,6e15.5)') &
                i,pl,qpath%xq(:,i), w2*RY_TO_CMM1, lw(:,it)*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE
        CALL errore('LW_QBZ_LINE', 'wrong mode (real or full)', 1)
      ENDIF
      !
      !
      IF(i>1) dpl = SQRT(SUM( (qpath%xq(:,i-1)-qpath%xq(:,i))**2 ))
      pl = pl + dpl
    ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    DEALLOCATE(D)
    !
  END SUBROUTINE LW_QBZ_LINE
  !   
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE SPECTR_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE interp_fc,   ONLY : fftinterp_mat2, mat2_diag
    USE linewidth,   ONLY : spectre_q, simple_spectre_q, add_exp_t_factor
    USE constants,   ONLY : RY_TO_CMM1
!     USE ph_velocity
    USE q_grid,      ONLY : q_grid_type, setup_simple_grid
    USE functions,   ONLY : write_temperature
    USE nanoclock
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qpath
    !
    REAL(DP) :: pl,dpl
    INTEGER :: iq, it, ie
    TYPE(q_grid_type) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    !
    REAL(DP),PARAMETER :: inv2_RY_TO_CMM1 = 1/(RY_TO_CMM1)**2
    !
    REAL(DP),ALLOCATABLE :: ener(:), spectralf(:,:,:)
    !
    ALLOCATE(spectralf(input%ne,S%nat3,input%nconf))
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              "spf.T"//TRIM(write_temperature(it,input%nconf,input%T))//&
                              "s"//TRIM(write_temperature(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "# spectral function mode: ", input%mode
      WRITE(1000+it, *) "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      WRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      CALL flush_unit(1000+it)
    ENDDO
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
    DO iq = 1,qpath%nq
      WRITE(*,'(i6,3f15.8)') iq, qpath%xq(:,iq)
      DO it = 1,input%nconf
        WRITE(1000+it, *)
        WRITE(1000+it, *) "#  xq",  qpath%xq(:,iq)
      ENDDO
      !
      IF (TRIM(input%mode) == "full") THEN
        spectralf = spectre_q(qpath%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, &
                                  S, grid, fc2, fc3, input%ne, ener/RY_TO_CMM1)
      ELSE IF (TRIM(input%mode) == "simple") THEN
        spectralf = simple_spectre_q(qpath%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, &
                                  S, grid, fc2, fc3, input%ne, ener/RY_TO_CMM1)
      ELSE
        CALL errore("SPECTR_QBZ_LINE", 'unknown mode "'//TRIM(input%mode)//'"', 1)
      ENDIF
      !
      IF(input%exp_t_factor) CALL add_exp_t_factor(input%nconf, input%T, input%ne, S%nat3, ener/RY_TO_CMM1, spectralf)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          WRITE(1000+it, '(2f14.8,100e14.6)') &
                pl, ener(ie), SUM(spectralf(ie,:,it))*inv2_RY_TO_CMM1, spectralf(ie,:,it)*inv2_RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ENDDO
      !
    ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    DEALLOCATE(spectralf, ener)
    !
  END SUBROUTINE SPECTR_QBZ_LINE
  !   
  END MODULE linewidth_program


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE linewidth_program
  USE environment,      ONLY : environment_start, environment_end
  USE constants,        ONLY : RY_TO_CMM1
  USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : print_citations_linewidth
  USE q_grid,           ONLY : q_grid_type
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(forceconst3_grid) :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(lwinput_type)     :: lwinput
  TYPE(q_grid_type)      :: qpath
  
!   REAL(DP) :: xq(3)
!   INTEGER  :: ios
!   REAL(DP) :: M(3),K(3),G(3)
  
  CALL mp_world_start(world_comm)
  CALL environment_start('LW')
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT(lwinput, qpath, S, fc2, fc3)
  !
  IF(TRIM(lwinput%calculation) == "lw") THEN
    !
    CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE &
  IF(TRIM(lwinput%calculation) == "spf") THEN
    !
    CALL SPECTR_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE
    CALL errore("lw", "what else to do?", 1)
  ENDIF
  !
  CALL environment_end('LW')
  CALL mp_world_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













