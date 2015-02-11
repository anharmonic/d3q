!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
MODULE linewidth_program
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info
  USE sparse_fc,ONLY : forceconst3
                       
  REAL(DP) :: default_sigma = 10._dp
  
  ! NOTE: energies in the LWINPUT_TYPE structure are assumed to be in CM^-1
  !       in the rest of the code energies are in Rydberg!
  TYPE lwinput_type
    !
    CHARACTER(len=16) :: calculation ! lw=linewidth, spf=spectral function
    CHARACTER(len=16) :: mode        ! "full" or "simple" spectral function 
    CHARACTER(len=256) :: outdir
    CHARACTER(len=256) :: prefix     ! put this in fron of file names
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
    ! for final state:
    REAL(DP) :: e_initial
    REAL(DP) :: q_initial(3)
    LOGICAL  :: q_resolved
    REAL(DP) :: sigmaq
    ! for isotope contribution to lw
    LOGICAL  :: isotopic_disorder
    !
    INTEGER :: nk(3)
  END TYPE lwinput_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(input, qpath, S, fc2, fc3)
    USE io_global,      ONLY : stdout
    USE q_grid,         ONLY : q_grid_type, setup_path, setup_simple_grid
    USE constants,      ONLY : RY_TO_CMM1
    USE more_constants, ONLY : INVALID, MASS_DALTON_TO_RY
    USE wrappers,       ONLY : f_mkdir_safe
    USE sparse_fc,      ONLY : forceconst3
    USE nist_isotopes_database, ONLY : compute_gs
    !
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(out) :: input
    TYPE(q_grid_type),INTENT(out)  :: qpath
    TYPE(forceconst2_grid),INTENT(out) :: fc2
    CLASS(forceconst3),POINTER,INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(out)   :: S
    !
    ! Input variable, and defaul values:
    CHARACTER(len=16)  :: calculation = "lw" ! "spf"
    CHARACTER(len=256) :: file_mat3  = INVALID ! no default
    CHARACTER(len=256) :: file_mat2  = INVALID ! no default
    CHARACTER(len=256) :: prefix     = INVALID ! default: calculation.mode
    !
    CHARACTER(len=256) :: outdir = './'
    LOGICAL            :: asr2 = .false.
    INTEGER            :: nconf = 1
    INTEGER            :: nq = -1
    INTEGER            :: nk(3) = (/-1, -1, -1/)
    LOGICAL            :: exp_t_factor = .false.
    !
    INTEGER  :: ne = -1
    REAL(DP) :: de = 1._dp, e0 = 0._dp
    REAL(DP) :: e_initial = -1._dp
    REAL(DP) :: q_initial(3) = 0._dp
    LOGICAL  :: q_resolved = .false.
    REAL(DP) :: sigmaq = 0.1_dp
    !
    LOGICAL  :: isotopic_disorder = .true.
    REAL(DP),ALLOCATABLE :: isotopes_mass(:), isotopes_conc(:), auxm(:), auxs(:)
    INTEGER :: n_isotopes, atomic_N
    !
    REAL(DP) :: xq(3), xq0(3)
    INTEGER  :: ios, ios2, i, j, naux, nq1, nq2, nq3
    CHARACTER(len=1024) :: line, word
    CHARACTER(len=16)   :: word2, word3
    CHARACTER(LEN=256), EXTERNAL :: TRIMCHECK
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    !
    LOGICAL :: qpoints_ok=.false., configs_ok=.false., isotopes_ok=.false.
    !
    NAMELIST  / lwinput / &
      calculation, outdir, prefix, &
      file_mat2, file_mat3, asr2, &
      nconf, nq, nk, &
      ne, de, e0, e_initial, q_initial, q_resolved, sigmaq, exp_t_factor, &
      isotopic_disorder

    WRITE(*,*) "Waiting for input"
    !
    READ(*, lwinput)
    WRITE(*, lwinput)
    !
    IF(TRIM(file_mat2) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    IF(TRIM(file_mat3) == INVALID ) CALL errore('READ_INPUT', 'Missing file_mat3', 1)
    IF(ANY(nk<0)) CALL errore('READ_INPUT', 'Missing nk', 1)    
    IF(nconf<0)   CALL errore('READ_INPUT', 'Missing nconf', 1)    
    
    input%file_mat2 = file_mat2
    input%file_mat3 = file_mat3
    input%outdir    = TRIMCHECK(outdir)
    input%asr2      = asr2
    input%nconf     = nconf
    input%nk        = nk
    input%exp_t_factor = exp_t_factor
    input%isotopic_disorder = isotopic_disorder
    !
    ios = f_mkdir_safe(input%outdir)
    IF(ios>0) CALL errore('READ_INPUT', 'cannot create directory: "'//TRIM(input%outdir)//'"',1)
    !
    CALL READ_DATA(input, s, fc2, fc3)
    !
    READ(calculation,*,iostat=ios) input%calculation, input%mode
    IF(ios/=0) THEN
      input%calculation = calculation
      input%mode = "full"
    ENDIF
    WRITE(*,*) "calculation: ", input%calculation, input%mode
    !
    IF(nq<0.and.TRIM(input%calculation)/="grid")      CALL errore('READ_INPUT', 'Missing nq', 1)    
    !
    IF(TRIM(prefix)==INVALID)THEN
      input%prefix = TRIM(input%calculation)//"."//TRIM(input%mode)
    ELSE
      input%prefix = prefix
    ENDIF
    !
    IF(TRIM(input%calculation) == 'spf' .and. ne < 0) &
      CALL errore('READ_INPUT', 'Missing ne for spf calculation', 1)    
    input%ne = ne
    input%de = de
    input%e0 = e0
    !
    IF(TRIM(input%calculation) == 'final' .and. e_initial < 0) &
      CALL errore('READ_INPUT', 'Missing e_initial for final state calculation', 1)    
    input%e_initial = e_initial
    input%q_initial = q_initial
    ! if we also want the q of the final state:
    input%q_resolved = q_resolved
    input%sigmaq  = sigmaq
    !
    ! Read the configurations:
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
        IF(qpoints_ok) CALL errore("READ_INPUT", "Won't reads q-points twice", 1)
        qpoints_ok = .true.
        !
        IF(TRIM(input%calculation) == "grid")THEN
          READ(*,*,iostat=ios) nq1, nq2, nq3
          IF(ios/=0) CALL errore("READ_INPUT", "reading nq1, nq2, nq3 for q-grid calculation", 1)
          line=''
          WRITE(*,"(2x,a,i4,a)") "Read ", qpath%nq, " q-points, "//TRIM(qpath%basis)//" basis"
          !CALL setup_simple_grid(S, nq1,nq2,nq3, qpath) this is done later
          CYCLE READ_CARDS
        ENDIF
        WRITE(*,*) "Reading QPOINTS"
        !
        qpath%basis = 'cartesian'
        READ(line,*,iostat=ios) word2, word3
        IF(ios==0) qpath%basis = TRIM(word3(1:9))
        !
        QPOINT_LOOP : & ! ..............................................................
        DO i = 1, nq
          READ(*,'(a1024)', iostat=ios) line
          IF(ios/=0) CALL errore("READ_INPUT","Expecting q point: input error.", 1)
          !
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
          !
          READ(line,*, iostat=ios2) xq(1), xq(2), xq(3)
          IF(ios2==0) THEN
            CALL setup_path(xq, xq, 1, qpath)
            CYCLE QPOINT_LOOP
          ENDIF
          !
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting q point, got: '"//TRIM(line)//"'.", 1)
          EXIT QPOINT_LOOP
          !
        ENDDO &
        QPOINT_LOOP ! .................................................................
        
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
        IF(configs_ok) CALL errore("READ_INPUT", "Won't reads configs twice", 1)
        configs_ok = .true.
        !
        WRITE(*,*) "Reading CONFIGS", nconf
        DO i = 1,nconf
          READ(*,'(a1024)', iostat=ios2) line
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration: input error.", 1)
          !
          ! try to read sigma and temperature
          READ(line,*,iostat=ios2) input%sigma(i), input%T(i)
          ! If it fails, read just temperature
          IF(ios2/=0) THEN
            READ(line,*,iostat=ios2) input%T(i)
            ! If this fails, complain
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting configuration, got: '"//TRIM(line)//"'.", 1)
            ! reuse previous value of sigma if we read jus ttemperature
            IF(i>1) THEN
              input%sigma(i) = input%sigma(i-1)
            ELSE
              CALL errore("READ_INPUT","I need at least one value of sigma", 1)
            ENDIF
          ENDIF
        ENDDO
        WRITE(*,'(2x,a,/,100(8f9.1,/))') "Temperatures:", input%T
        WRITE(*,'(2x,a,/,100(8f9.1,/))') "Smearings:   ", input%sigma
        WRITE(*,*)

      CASE ("ISOTOPES")
        IF(isotopes_ok) CALL errore("READ_INPUT", "Won't reads isotopes twice", 1)
        isotopes_ok = .true.
        !
        WRITE(*,*) "Reading ISOTOPES", S%ntyp
        IF(.not.isotopic_disorder) WRITE(*,*) "WARNING! you did not se isotopic_disorder to true!"
        ALLOCATE(auxs(S%ntyp), auxm(S%ntyp))
        !
        ISOTOPE_TYPE_LOOP : &
        DO i = 1,S%ntyp
          READ(*,'(a1024)', iostat=ios2) line
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope: input error.", 1)
          !

!     DO it = 1, S%ntyp
!       CALL compute_gs(auxm(it), auxs(it), S%atm(it), 0, 0)
!     ENDDO
!     ! Rearrange the cross sections with the mode index, for simplicity
!     nu = 0
!     DO ia = 1, S%nat
!       it = S%ityp(ia)
!       DO alpha = 1,3
!         nu = nu+1
!         gm(nu)  = auxm(it)
!         gs2(nu) = auxs(it)
!       ENDDO
!     ENDDO
          ! Try to read isotope name and method
          READ(line,*,iostat=ios2) word2, word3
          !
          ! If this fails, complain
          IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope, got: '"//TRIM(line)//"'.", 1)
          IF(word2 /= S%atm(i)) &
            WRITE(*,"(a,i3,2a)") "WARNING: isotope name from input does not match FC file", i, word2, S%atm(i)
          !
          !
          IF (word3=="natural") THEN
            CALL compute_gs(auxm(i), auxs(i), word2, 0, 0)
            
          ELSE IF (word3=="N") THEN
            READ(line,*,iostat=ios2) word2, word3, atomic_N
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope atomic number: input error.", 2)
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, 0)
            
          ELSE IF (word3=="M") THEN
            n_isotopes = 1
            ALLOCATE(isotopes_mass(n_isotopes), isotopes_conc(n_isotopes))
            isotopes_conc = 1._dp
            READ(line,*,iostat=ios2) word2, word3, isotopes_mass(1)
            IF(ios2/=0) CALL errore("READ_INPUT","Expecting isotope atomic mass: input error.", 3)
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, n_isotopes, isotopes_mass, isotopes_conc)
            DEALLOCATE(isotopes_mass, isotopes_conc)
            
          ELSE IF (word3=="isotopes") THEN
            READ(line,*,iostat=ios2) word2, word3, n_isotopes
            ALLOCATE(isotopes_mass(n_isotopes), isotopes_conc(n_isotopes))
            DO j = 1,n_isotopes
              READ(*,*, iostat=ios2) isotopes_mass(j), isotopes_conc(j)
              IF(ios2/=0) CALL errore("READ_INPUT","REading isotope line: input error.", 4)
            ENDDO
            CALL compute_gs(auxm(i), auxs(i), word2, atomic_N, n_isotopes, isotopes_mass, isotopes_conc)
            DEALLOCATE(isotopes_mass, isotopes_conc)

          ELSE
            CALL errore("READ_INPUT","I do not understand this isotope choice '"//TRIM(word2)//"'.", 1)
          ENDIF
        ENDDO &
        ISOTOPE_TYPE_LOOP 
        !
        ! Replace atomic masses with those from input:
        S%amass = auxm * MASS_DALTON_TO_RY
        S%amass_variance = auxs
        DEALLOCATE(auxs, auxm)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CASE DEFAULT
        IF(TRIM(line) /= '') WRITE(*,*) "Skip:", TRIM(line)
      END SELECT
      !
      READ(*,'(a1024)', iostat=ios) line
      word = ''
      !
    ENDDO &
    READ_CARDS
    !
    IF(TRIM(input%calculation)=="grid") THEN
      IF(.not.qpoints_ok) THEN
        nq1=nk(1); nq2=nk(2); nq3=nk(3)
        qpoints_ok = .true.
      ENDIF
      CALL setup_simple_grid(S, nq1, nq2, nq3, qpath)
      input%prefix = TRIM(input%prefix)//&
                "."//TRIM(int_to_char(nq1))// &
                "."//TRIM(int_to_char(nq2))// &
                "."//TRIM(int_to_char(nq3))
    ENDIF
    !
    ! Set natural isotope concentration for every species, if not read from input
    IF(input%isotopic_disorder.and..not.isotopes_ok)THEN
      DO j = 1,S%ntyp
        CALL compute_gs(S%amass(i), S%amass_variance(i), S%atm(j), 0, 0)
      ENDDO
      S%amass = S%amass * MASS_DALTON_TO_RY
    ENDIF
    
    IF(.not.qpoints_ok) CALL errore("READ_INPUT", "I did not find QPOINTS card", 1)
    IF(.not.configs_ok) CALL errore("READ_INPUT", "I did not find CONFIGS card", 1)
    !
  END SUBROUTINE READ_INPUT
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_DATA(input, s, fc2, fc3)
    USE iso_c_binding,  ONLY : c_int
    USE input_fc,       ONLY : same_system, read_fc2, aux_system, div_mass_fc2
    USE asr2_module,    ONLY : impose_asr2
    USE io_global,      ONLY : stdout
    USE sparse_fc,      ONLY : read_fc3, forceconst3
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)        :: input
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    CLASS(forceconst3),POINTER,INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(inout)   :: s
    TYPE(ph_system_info) :: s3
    !
    INTEGER(kind=c_int) :: kb
    !
    CALL read_fc2(input%file_mat2, S,  fc2)
    fc3 => read_fc3(input%file_mat3, S3)
!     CALL read_fc3(input%file_mat3, S3, fc3)
    !
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
!     CALL div_mass_fc3(S, fc3)
    CALL fc3%div_mass(S)
    !
  END SUBROUTINE READ_DATA
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE LW_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag
    USE linewidth,      ONLY : linewidth_q, selfnrg_q, spectre_q, freq_phq
    USE constants,      ONLY : RY_TO_CMM1
!     USE ph_velocity
    USE q_grid,         ONLY : q_grid_type, setup_simple_grid
    USE more_constants, ONLY : write_conf
    USE timers
    USE sparse_fc,      ONLY : forceconst3
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qpath
    !
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)
    REAL(DP) :: pl,dpl
    INTEGER :: iq, it
    TYPE(q_grid_type) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP)   :: lw(S%nat3,input%nconf)
    REAL(DP)   :: sigma_ry(input%nconf)
    CHARACTER(len=15) :: f1, f2
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
      IF (TRIM(input%mode) == "full") THEN
        WRITE(1000+it, *) "# calculation of linewidth (gamma_n) [and lineshift (delta_n)]"
      ELSE
        WRITE(1000+it, *) "# calculation of linewidth (gamma_n)"
      ENDIF
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
    ENDDO
    ! Prepare formats to write out data
    WRITE(f1,'(i6,a)') S%nat3, "f12.6,6x,"
    WRITE(f2,'(i6,a)') S%nat3, "e15.5,6x,"
    !
    ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
    ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
    ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
    ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
    !  => 2d = 2 sqrt(log(2) c
    !      d = sqrt(log(2)) c
    !      d = 0.83255 c
    ! IMHO: you need to use a sigma that is 0.6 (=0.5/0.83255) times smaller when using
    ! linewidth_q than when using selfnrg_q in order to get the same values
    ! THIS FACTOR IS NOW INCLUDED DIRECTLY IN SUM_LINEWIDTH_MODES!
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    !
    DO iq = 1,qpath%nq
      WRITE(*,'(i6,3f15.8)') iq, qpath%xq(:,iq)
      !
      CALL freq_phq(qpath%xq(:,iq), S, fc2, w2, D)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      IF (TRIM(input%mode) == "full") THEN
        ls = selfnrg_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                        S, grid, fc2, fc3)
        DO it = 1,input%nconf
!           WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,6f12.6,6x,6e15.5,6x,6e15.5)') &
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,'//f1//f2//f2//'x)') &
                iq,pl,qpath%xq(:,iq), w2*RY_TO_CMM1, -DIMAG(ls(:,it))*RY_TO_CMM1, DBLE(ls(:,it))*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE IF (TRIM(input%mode) == "real") THEN
        lw = linewidth_q(qpath%xq(:,iq), input%nconf, input%T,  sigma_ry, &
                        S, grid, fc2, fc3)
        DO it = 1,input%nconf
!           WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,6f12.6,6x,6e15.5)') &
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,'//f1//f2//'x)') &
                iq,pl,qpath%xq(:,iq), w2*RY_TO_CMM1, lw(:,it)*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE
        CALL errore('LW_QBZ_LINE', 'wrong mode (real or full)', 1)
      ENDIF
      !
    ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    CALL t_freq%print() 
    CALL t_sum%print() 
    CALL t_fc3int%print() 
    CALL t_fc3m2%print() 
    CALL t_fc3rot%print() 
    !
  END SUBROUTINE LW_QBZ_LINE
  !   
  !  
  SUBROUTINE SPECTR_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag
    USE linewidth,      ONLY : spectre_q, simple_spectre_q, add_exp_t_factor, freq_phq
    USE constants,      ONLY : RY_TO_CMM1
    USE q_grid,         ONLY : q_grid_type, setup_simple_grid
    USE more_constants, ONLY : write_conf
    USE nanoclock
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(inout)  :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qpath
    !
    REAL(DP) :: pl,dpl
    INTEGER :: iq, it, ie
    TYPE(q_grid_type) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP)   :: sigma_ry(input%nconf)
    !
    REAL(DP),ALLOCATABLE :: ener(:), spectralf(:,:,:)
    !
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)
    ALLOCATE(spectralf(input%ne,S%nat3,input%nconf))
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    ener = ener/RY_TO_CMM1
    !
    ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
    ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
    ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
    ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
    !  => 2d = 2 sqrt(log(2) c
    !      d = sqrt(log(2)) c
    !      d = 0.83255 c
    ! IMHO: you need to use a sigma that is 0.6 (=0.5/0.83255) times smaller when using
    ! linewidth_q than when using selfnrg_q in order to get the same values
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                              "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "# spectral function mode: ", input%mode
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      WRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      CALL flush_unit(1000+it)
    ENDDO
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
    DO iq = 1,qpath%nq
      CALL freq_phq(qpath%xq(:,iq), S, fc2, w2, D)
      WRITE(*,'(i6,3f15.8,5x,6f12.6)') iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      !
      DO it = 1,input%nconf
        WRITE(1000+it, *)
        WRITE(1000+it, '(a,i6,3f15.8,5x,6f12.6)') "#  xq",  iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      ENDDO
      !
      IF (TRIM(input%mode) == "full") THEN
        spectralf = spectre_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                                  S, grid, fc2, fc3, input%ne, ener)
      ELSE IF (TRIM(input%mode) == "simple") THEN
        spectralf = simple_spectre_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                                  S, grid, fc2, fc3, input%ne, ener)
      ELSE
        CALL errore("SPECTR_QBZ_LINE", 'unknown mode "'//TRIM(input%mode)//'"', 1)
      ENDIF
      !
      IF(input%exp_t_factor) CALL add_exp_t_factor(input%nconf, input%T, input%ne, S%nat3, ener, spectralf)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          WRITE(1000+it, '(2f14.8,100e14.6)') &
                pl, ener(ie)*RY_TO_CMM1, SUM(spectralf(ie,:,it))/RY_TO_CMM1**2, spectralf(ie,:,it)/RY_TO_CMM1**2
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
  !  
  SUBROUTINE FINAL_STATE_LINE(input, qpath, S, fc2, fc3)
    USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag
    USE  final_state,   ONLY : final_state_q
    USE constants,      ONLY : RY_TO_CMM1
    USE q_grid,         ONLY : q_grid_type, setup_simple_grid
    USE more_constants, ONLY : write_conf
    USE nanoclock
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(inout)  :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qpath
    !
    REAL(DP) :: pl,dpl, e_inital_ry 
    INTEGER :: iq, it, ie
    TYPE(q_grid_type) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP) :: sigma_ry(input%nconf)
    !
    REAL(DP),PARAMETER :: inv2_RY_TO_CMM1 = 1/(RY_TO_CMM1)**2
    !
    REAL(DP),ALLOCATABLE :: ener(:), fstate(:,:,:)
    CHARACTER(len=256) :: filename
    !
    ALLOCATE(fstate(input%ne,S%nat3,input%nconf))
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    ! Convert to Rydberg:
    ener = ener/RY_TO_CMM1
    e_inital_ry = input%e_initial/RY_TO_CMM1
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    DO it = 1,input%nconf
      filename = TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                              "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out"
      OPEN(unit=1000+it, file=filename)
      WRITE(*,*) "opening ", TRIM(filename)
      WRITE(1000+it, *) "# final state decompositions, mode: ", input%mode
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      WRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      CALL flush_unit(1000+it)
    ENDDO
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
!     DO iq = 1,qpath%nq
!       WRITE(*,'(i6,3f15.8)') iq, qpath%xq(:,iq)
!       !
!       DO it = 1,input%nconf
!         WRITE(1000+it, *)
!         WRITE(1000+it, '(a,i6,3f15.8)') "#  xq",  iq, qpath%xq(:,iq)
!       ENDDO
      !
      fstate = final_state_q(input%q_initial, qpath, input%nconf, input%T, sigma_ry, &
                             S, grid, fc2, fc3, e_inital_ry, input%ne, ener, &
                             input%q_resolved, input%sigmaq, input%outdir, input%prefix)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          WRITE(1000+it, '(2f14.8,100e18.6e4)') &
                pl, ener(ie), SUM(fstate(ie,:,it)), fstate(ie,:,it)
          CALL flush_unit(1000+it)
        ENDDO
      ENDDO
      !
!     ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    DEALLOCATE(fstate, ener)
    !
  END SUBROUTINE FINAL_STATE_LINE
  !   
  END MODULE linewidth_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!               !               !               !               !               !               !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE linewidth_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : print_citations_linewidth
  USE q_grid,           ONLY : q_grid_type !, setup_simple_grid
  USE sparse_fc,        ONLY : forceconst3
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(lwinput_type)     :: lwinput
  TYPE(q_grid_type)      :: qpath

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT(lwinput, qpath, S, fc2, fc3)
  !
  IF(    TRIM(lwinput%calculation) == "lw"   &
    .or. TRIM(lwinput%calculation) == "grid" &
    ) THEN
    !
    CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE &
  IF(TRIM(lwinput%calculation) == "spf") THEN
    !
    CALL SPECTR_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
!   ELSE &
!   IF(TRIM(lwinput%calculation) == "grid") THEN
!     !
!     ! Compute the linewidth over the grid, will be reused later to do self-consistent linewidth
!     CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
  ELSE &
  IF(TRIM(lwinput%calculation) == "final") THEN
    !
    CALL FINAL_STATE_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE
    CALL errore("lw", "what else to do?", 1)
  ENDIF
  !
!   CALL environment_end('LW')
!   CALL mp_world_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

