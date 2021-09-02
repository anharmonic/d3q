!
! Written by Lorenzo Paulatto (2017) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE cmdline_param_module

  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  ! 
  ! Each of the following functions will parse a single command line argument
  ! after parsing, the argument will be removed from command line, calling it
  ! again will cause an error.
  ! 
  PUBLIC :: cmdline_param_char
  ! Example: 
  !   char_var = cmdline_param_char('c' [, default="def value"] [,found])
  ! Return a character array with content "value" from argument '-c value',
  ! will stop if argument is not found, unless default is specified.
  ! I the optional logical argument "found" is present, it will be set to 
  ! false when the default value is used, to true if the argument is found.
  !
  PUBLIC :: cmdline_param_dble
  ! Example: 
  !   dble_var = cmdline_param_char('f' [, default=0.d0] [,found])
  ! Same as above, but tries to read a double-precision FP number
  !
  PUBLIC :: cmdline_param_int
  ! Example: 
  !   int_var = cmdline_param_char('i' [, default=1] [,found])
  ! Same as above, but tries to read an integer number
  !
  PUBLIC :: cmdline_param_logical
  ! Example:
  !    bool_var = cmdline_param_logical('b')
  ! Returns true if "-b" is found in the command line, false otherwise.
  !
  PUBLIC :: cmdline_check_exausted
  ! Check that there are no command line options left in the command line,
  ! issues an error message and stops if any is found.
  !
  PUBLIC :: cmdline_residual
  ! Returns all that is left in the command line after stripping the command name
  ! (i.e. the executable name, argument 0). It also calls cmdline_check_exausted
  !
  PUBLIC :: cmdline_to_namelist
  !  cmdline_to_namelist(nname, nunit)
  ! Creates and saves to nunit a file containing the namelist nname, the namelist will be
  ! generated from the "--" command line arguments, like this:
  !   --option value 
  ! becomes
  !   option = value
  ! The code detects if value can be read as a number or a logical it will stay as is, otherwise
  ! it will be put in quotes.  Setting array values like
  !   --nk 1,1,1
  ! Usually work, but it is a bit magic
  INTEGER :: length = -1
  CHARACTER(len=:),ALLOCATABLE :: command
  PUBLIC :: fgetpid
  !
  CONTAINS

  INTEGER FUNCTION fgetpid()
    USE random_numbers, ONLY : randy
!#if defined(__INTEL_COMPILER)
!    USE IFPORT, ONLY : getpid
!    fgetpid = INT(getpid(), kind=8)
!#elif defined (__GFORTRAN__)
!    fgetpid = INT(getpid(), kind=8)
!#else 
    fgetpid = INT(randy()*1.d+6)
!#endif
  END FUNCTION

  REAL(DP) FUNCTION cmdline_param_dble(switch, default, found)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    REAL(DP),OPTIONAL,INTENT(in) :: default
    LOGICAL,OPTIONAL,INTENT(out) :: found
    !
    CHARACTER(len=:),ALLOCATABLE :: char_value
    INTEGER :: ios
    REAL(DP) :: value
    !
    char_value = cmdline_param_char(switch, default=' ', found=found)
    READ(char_value,*,iostat=ios) value
    IF(ios/=0) THEN
      IF(present(default)) THEN
        value = default
      ELSE
        CALL errore("cmdline_param", "required parameter not found '-"//switch//"'", 1)
       ENDIF
    ENDIF
    DEALLOCATE(char_value)
    cmdline_param_dble = value
    RETURN
    !
  END FUNCTION
  !
  INTEGER FUNCTION cmdline_param_int(switch, default, found)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    INTEGER,OPTIONAL,INTENT(in) :: default
    LOGICAL,OPTIONAL,INTENT(out) :: found
    !
    CHARACTER(len=:),ALLOCATABLE :: char_value
    INTEGER :: ios
    INTEGER :: value
    !
    char_value = cmdline_param_char(switch, default=' ', found=found)
    READ(char_value,*,iostat=ios) value
    IF(ios/=0) THEN
      IF(present(default)) THEN
        value = default
      ELSE
        CALL errore("cmdline_param", "required parameter not found '-"//switch//"'", 1)
       ENDIF
    ENDIF
    DEALLOCATE(char_value)
    cmdline_param_int = value
    RETURN
    !
  END FUNCTION
  !
  LOGICAL FUNCTION cmdline_param_logical(switch)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    !
    CHARACTER(len=:),ALLOCATABLE :: char_value
    LOGICAL :: value
    !
    char_value = cmdline_param_char(switch, default=' ', found=value)
    DEALLOCATE(char_value)
    cmdline_param_logical = value
    RETURN
    !
  END FUNCTION
  !
  ! Throw away the command name (arg 0) and return the rest
  FUNCTION cmdline_residual()
    CHARACTER(len=:),ALLOCATABLE :: cmdline_residual
    INTEGER :: i
    !
    CALL cmdline_check_exausted()
    !
    DO i =1, length
      IF(command(i:i) == " ") EXIT
    ENDDO
    ALLOCATE(character(length-i+1) :: cmdline_residual)
    cmdline_residual = command(i:length)
    RETURN
  END FUNCTION
  !
  SUBROUTINE cmdline_check_exausted()
    CHARACTER(len=:),ALLOCATABLE :: cmdline_residual
    INTEGER :: i
    DO i =1, length-2
      IF(command(i:i+1) == " -") THEN
        CALL errore("cmdline_check_exausted", &
                    "unknown argument '"//command(i:i+2)//"', use '-h' for help",1)
      ENDIF
    ENDDO
  END SUBROUTINE
  !
  FUNCTION cmdline_param_char(switch, default, found)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    CHARACTER(len=*),OPTIONAL,INTENT(in) :: default
    CHARACTER(len=:),ALLOCATABLE :: cmdline_param_char
    LOGICAL,OPTIONAL,INTENT(out) :: found
    !
    INTEGER :: i, j, start_switch, start_arg, end_arg, len_arg, len_switch
    LOGICAL :: start, empty
    !
    ! On first call, allocate the space and read in the command line
    IF(length<0) THEN
      CALL GET_COMMAND(length=length)
      ALLOCATE(character(length) :: command)
      CALL GET_COMMAND(command=command)
    ENDIF
    !
    IF(present(found)) found = .false.
    !
    start = .true.
    empty = .false.
    start_arg = -1
    end_arg = -1
    !
    SCAN_COMMAND : &
    DO i = 1, length-1
      IF(start .and. command(i:i+1) == '-'//switch)THEN
        IF(present(found)) found = .true.
        start_switch = i
        !
        IF(start_switch+2>=length)THEN
            start_arg = 0
            empty = .true.
            EXIT SCAN_COMMAND
        ENDIF
        !
        SCAN_ARG : &
        DO j = i+2,length
          IF(command(j:j) /= ' ' .and. start_arg<0) THEN
            IF(command(j:j)=='-')THEN
              start_arg = 0
              empty = .true.
              EXIT SCAN_COMMAND
            ELSE
              start_arg = j
            ENDIF
          ENDIF
          IF(command(j:j) == ' ' .and. start_arg>0) THEN
            end_arg = j-1
            EXIT SCAN_ARG
          ENDIF
        ENDDO &
        SCAN_ARG
        IF(end_arg<0) end_arg = length
        EXIT SCAN_COMMAND
        !
      END IF
      IF(command(i:i) == ' ') start = .true.
    ENDDO &
    SCAN_COMMAND
    !
    IF(start_arg<0) THEN
      IF(present(default)) THEN
        ALLOCATE(character(len(default)) :: cmdline_param_char)
        cmdline_param_char = default
        RETURN
      ELSE
        CALL errore("cmdline_param", "required parameter not found '-"//switch//"'; &
                     &use '-h' for help", 1)
      ENDIF
    ENDIF
    IF(empty)THEN
      IF(.not.present(found))THEN
        CALL errore("cmdline_param", "missing value for parameter '-"//switch//"'; &
                    &use '-h' for help", 1)
      ENDIF
      !print*, ">>"//trim(command)//".."
      command(start_switch:length-2) = command(start_switch+2:length)
      command(length-1:length) = ''
      length = length-2
      ALLOCATE(character(0) :: cmdline_param_char)
      cmdline_param_char = ''
      !print*, ".."//trim(command(1:length))//"<<"
      RETURN
    ENDIF
    !
    len_arg = end_arg-start_arg+1
    ALLOCATE(character(len_arg) :: cmdline_param_char)
    cmdline_param_char = command(start_arg:end_arg)
    !
    !print*, ">>"//trim(command)//".."
    len_switch = end_arg-start_switch+1
    IF(end_arg<length)THEN
      command(start_switch:length-len_switch) = command(end_arg+1:length)
      command(length-len_switch+1:) = ''
    ELSE
      command(start_arg:length) = ''
    ENDIF
    length = length - len_switch
    !print*, ".."//trim(command(1:length))//"<<"
    !
  END FUNCTION
  !
  SUBROUTINE cmdline_to_namelist(nname, nunit)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: nname
    INTEGER,INTENT(in) :: nunit
    INTEGER :: nargs, i, ios, alen, astatus
    CHARACTER(len=512) :: argv, cdummy
    REAL(dp) :: rdummy
    LOGICAL  :: ldummy
    LOGICAL, EXTERNAL :: matches
    nargs = COMMAND_ARGUMENT_COUNT()

    WRITE(nunit,*) "&"//TRIM(ADJUSTL(nname))

    DO i = 1, nargs
      CALL GET_COMMAND_ARGUMENT(i,argv,alen,astatus)
      IF(astatus<0) CALL errore("cmdline_to_namelist","truncated command",1)
      IF(argv(1:2)=="--") THEN
        IF(alen>=3) WRITE(unit=nunit,fmt="(2a)",advance='no') "  ",TRIM(argv(3:))
        ! write an "=" if it is not already in the argument
        IF(.not.matches("=",argv)) WRITE(unit=nunit,fmt="(a)",advance='no') " = "
      ELSE IF (argv(1:1) == "-") THEN
        ! A simple way to get rid of single "-" options: they are put in the file as comments
        WRITE(unit=nunit,fmt="(3a)",advance='no') "  ! ",TRIM(argv)
        IF(i<nargs)THEN
          CALL GET_COMMAND_ARGUMENT(i+1,argv,alen,astatus)
          IF(argv(1:1) == "-") WRITE(nunit,*)
        ENDIF
      ELSE
        ! Try to determine if command is a number, a logical or a string
        ! unfortunately compilers can be very generous when readin, i.e.
        ! "file_02" could be read as a number and break this mechanism
        ! So I check if first character is not a number

        IF(.not.matches(argv(1:1),"01234567890-+")) THEN
          ! it could still be a logical type (.f. ot .false. or just false)  
          READ(argv,*,iostat=ios) ldummy
          IF(ios==0) THEN
            WRITE(nunit,*) TRIM(argv), ","
            CYCLE
          ELSE
            WRITE(nunit,*) '"'//TRIM(argv)//'"', ","
            CYCLE
          ENDIF
        ELSE
          ! here it is definitely a number
          WRITE(nunit,*) argv, ","
          CYCLE
        ENDIF
      ENDIF
    ENDDO

    WRITE(nunit,*) 
    WRITE(nunit,*) "/"
    WRITE(nunit,*) 
    CALL FLUSH(nunit)

  END SUBROUTINE
  !
END MODULE

