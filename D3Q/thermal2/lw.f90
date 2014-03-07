!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1 
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
    CHARACTER(len=256) :: file_mat3
    CHARACTER(len=256) :: file_mat2
    LOGICAL            :: asr2
    INTEGER            :: nconf
    REAL(DP),ALLOCATABLE :: T(:), sigma(:)
    !
    INTEGER :: nk(3)
  END TYPE lwinput_type
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT(input, qpath)
    USE io_global,      ONLY : stdout
    USE q_grid,         ONLY : q_grid_type, setup_path
    USE constants,      ONLY : RY_TO_CMM1
    !
    IMPLICIT NONE
    !
    TYPE(lwinput_type),INTENT(out) :: input
    TYPE(q_grid_type),INTENT(out)  :: qpath
    !
    ! Input variable, and defaul values:
    CHARACTER(len=256) :: file_mat3 = '///'
    CHARACTER(len=256) :: file_mat2 = '///'
    LOGICAL            :: asr2 = .false.
    INTEGER            :: nconf = 1
    INTEGER            :: nq = -1
    INTEGER            :: nk(3) = (/-1, -1, -1/)
    !
    REAL(DP) :: xq(3), xq0(3)
    INTEGER  :: ios, ios2, i, naux
    CHARACTER(len=1024) :: line, word
    !
    NAMELIST  / lwinput / &
      file_mat2, file_mat3, asr2, &
      nconf, nq, nk

    WRITE(*,*) "Waiting for input"
    !
    READ(*, lwinput)
    
    IF(TRIM(file_mat2) == '///' ) CALL errore('READ_INPUT', 'Missing file_mat2', 1)
    IF(TRIM(file_mat3) == '///' ) CALL errore('READ_INPUT', 'Missing file_mat3', 1)
    IF(ANY(nk<0)) CALL errore('READ_INPUT', 'Missing nk', 1)    
    IF(nq<0)      CALL errore('READ_INPUT', 'Missing nq', 1)    
    IF(nconf<0)   CALL errore('READ_INPUT', 'Missing nconf', 1)    
    
    input%file_mat2 = file_mat2
    input%file_mat3 = file_mat3
    input%asr2      = asr2
    input%nconf     = nconf
    input%nk        = nk

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
      CASE ("QPOINTS")
        WRITE(*,*) "Reading QPOINTS"
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
          IF(ios/=0) CALL errore("READ_INPUT","Expecting q point, got: '"//TRIM(line)//"'.", 1)
          EXIT QPOINT_LOOP
          
        ENDDO &
        QPOINT_LOOP
        
        WRITE(*,"(2x,a,i4,a)") "Read ", qpath%nq, " q-points."
        WRITE(*,"(2x,3f12.6)") qpath%xq
        WRITE(*,*)
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
    USE linewidth,   ONLY : linewidth_q, lineshift_q
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
!     REAL(DP) :: lw(S%nat3)
    REAL(DP) :: pl,dpl
    INTEGER :: i, it
    TYPE(q_grid_type) :: grid
!     REAL(DP)   :: ratio(S%nat3)
!     INTEGER,PARAMETER :: nconf = 21
!     REAL(DP)   :: T(nconf), sigma
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    !
    ALLOCATE(D(S%nat3, S%nat3))
    !
!     WRITE(*,'(2x,a,/,100(12f6.1))') "Temperatures", input%T
!     WRITE(*,'(2x,a,/,100(12f6.1))') "Sigma       ", input%sigma

    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
!     T = (/ 300._dp, 500._dp /)
!     T = (/ 1._dp, 5._dp, 10._dp, 50._dp, 100._dp, 150._dp, 200._dp, &
!            250._dp, 300._dp, 350._dp, 400._dp, 450._dp, 500._dp, 550._dp, &
!            600._dp, 700._dp, 800._dp, 900._dp, 1000._dp, 1200._dp, 1500._dp /)
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file="lw.T"//TRIM(write_temperature(it,input%nconf,input%T))//&
                                "s"//TRIM(write_temperature(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
    ENDDO
!     CALL print_header()
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
    DO i = 1,qpath%nq
      WRITE(*,'(i6,3f15.8)') i, qpath%xq(:,i)
      CALL fftinterp_mat2(qpath%xq(:,i), S, fc2, D)
      CALL mat2_diag(S, D, w2)
      !
      WHERE(w2>0._dp)
        w2=SQRT(w2)
      ELSEWHERE
        w2= -SQRT(-w2)
      ENDWHERE
      !
        ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
        ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
        ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
        ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
        !  => 2d = 2 sqrt(log(2) c => d = sqrt(log(2)) d = 0.83255 c
        !input%sigma = 0.83255 *input%sigma/RY_TO_CMM1
      ls = lineshift_q(qpath%xq(:,i), input%nconf, input%T, 0.83255 *input%sigma/RY_TO_CMM1, &
                       S, grid, fc2, fc3)
      !
      DO it = 1,input%nconf
        WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,6f12.6,6x,6e15.5,6x,6e15.5)') &
              i,pl,qpath%xq(:,i), w2*RY_TO_CMM1, -DIMAG(ls(:,it))*RY_TO_CMM1, DBLE(ls(:,it))*RY_TO_CMM1
        CALL flush_unit(1000+it)
      ENDDO
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

  CALL READ_INPUT(lwinput, qpath)
  CALL READ_DATA (lwinput, s, fc2, fc3)
  !CALL impose_asr2(S%nat, fc2)
 
!   CALL QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/),&
!                    200, S, fc2)

!   G = (/ 0._dp, 0._dp, 0._dp /)
!   M = (/ 0.5_dp,      sqrt(3._dp)/6, 0._dp /)
!   K = (/ 1._dp/3._dp, sqrt(3._dp)/3, 0._dp /)
!   
!   CALL setup_path(G, M, 50, qpath)
!   CALL setup_path(M, K, 42, qpath)
!   CALL setup_path(K, G, 36, qpath)
  !
  CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
  !
  CALL environment_end('LW')
  CALL mp_world_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













