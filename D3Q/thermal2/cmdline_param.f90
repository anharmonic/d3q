!
MODULE cmdline_param_module

  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  

!   INTERFACE cmdline_param
!     MODULE PROCEDURE cmdline_param_char
!     MODULE PROCEDURE cmdline_param_dble
!     MODULE PROCEDURE cmdline_param_int
!   END INTERFACE
  !
  INTEGER :: length = -1
  CHARACTER(len=:),ALLOCATABLE :: command
  !
  CONTAINS

  REAL(DP) FUNCTION cmdline_param_dble(switch, default)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    REAL(DP),OPTIONAL,INTENT(in) :: default
    !
    CHARACTER(len=:),ALLOCATABLE :: char_value
    INTEGER :: ios
    REAL(DP) :: value
    !
    char_value = cmdline_param_char(switch, default=' ')
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
  INTEGER FUNCTION cmdline_param_int(switch, default)
    IMPLICIT NONE
    CHARACTER(len=1),INTENT(in)  :: switch
    INTEGER,OPTIONAL,INTENT(in) :: default
    !
    CHARACTER(len=:),ALLOCATABLE :: char_value
    INTEGER :: ios
    INTEGER :: value
    !
    char_value = cmdline_param_char(switch, default=' ')
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
        CALL errore("cmdline_param", "required parameter not found '-"//switch//"'", 1)
      ENDIF
    ENDIF
    IF(empty)THEN
      print*, ">>"//trim(command)//".."
      command(start_switch:length-2) = command(start_switch+2:length)
      command(length-1:length) = ''
      length = length-2
      ALLOCATE(character(0) :: cmdline_param_char)
      cmdline_param_char = ''
      print*, ".."//trim(command(1:length))//"<<"
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
END MODULE

