!
MODULE more_constants
  USE kinds, ONLY : DP
  REAL(DP),PARAMETER :: RY_TO_JOULE =  0.5* 4.35974394e-18
  REAL(DP),PARAMETER :: RY_TO_SECOND = 2* 2.418884326505e-17
  REAL(DP),PARAMETER :: RY_TO_METER = 5.2917721092e-11
  REAL(DP),PARAMETER :: RY_TO_WATTMM1KM1 = RY_TO_JOULE / (RY_TO_SECOND * RY_TO_METER)
  CHARACTER(len=3),PARAMETER :: INVALID = '///'
  
  CONTAINS
  !
  CHARACTER(len=256) &
  FUNCTION sigma_file_name(prefix, nq1, nq2, nq3, cT, csigma)
    IMPLICIT NONE
    CHARACTER(len=256),INTENT(in) :: prefix
    INTEGER,INTENT(in) :: nq1, nq2, nq3
    CHARACTER(len=6),INTENT(in) ::  cT, csigma
    !
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    !
    sigma_file_name= TRIM(prefix)//&
                "."//TRIM(int_to_char(nq1))// &
                "."//TRIM(int_to_char(nq2))// &
                "."//TRIM(int_to_char(nq3))// &
                "."//TRIM(cT)// &
                "."//TRIM(csigma)//".out"
    !
  END FUNCTION sigma_file_name
  !
  ! Write a number from a list using as many digits after the dot as the longest in the list.
  CHARACTER(len=6) &
  FUNCTION write_temperature(it,nt,T) RESULT(str)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: it, nt
    REAL(DP),INTENT(in) :: T(nt)
    
    INTEGER :: max_digit_left=0, max_digit_right=0, jt
    REAL(DP) :: Tx
    CHARACTER(len=64):: fmt1='', fmt2=''
    
    DO jt = 1,nt      
      max_digit_left = MAX( max_digit_left , CEILING(LOG10(T(jt))) )
      IF(T(jt)>0) THEN
        Tx = 1._dp/(T(jt)-INT(T(jt)))
        max_digit_right = MAX( max_digit_right , CEILING(LOG10(Tx)) )
      ENDIF
    ENDDO
    
    str=""
    WRITE(fmt1,'(i6)') max_digit_right
    fmt2 = "(1f6."//TRIM(ADJUSTL(fmt1))//")"
  !     print*, fmt2, max_digit_left, max_digit_right
    
    WRITE(str,fmt2) T(it)
    str=ADJUSTL(str)
    
  END FUNCTION
  !
  CHARACTER(len=6) &
  FUNCTION write_sigma(it,nat3,nt,sigma) RESULT(str)
    USE constants, ONLY : eps12
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: it, nt, nat3
    REAL(DP),INTENT(in) :: sigma(nat3,nt)
    REAL(DP) :: csigma(nt)
    INTEGER :: jt, is
    LOGICAL :: same
    CHARACTER(len=6),EXTERNAL :: int_to_char

    same = .true.
    DO is = 1,nat3
      same=same.and.(ABS(sigma(1,it)-sigma(is,it))<eps12)
    ENDDO
    IF(.not.same) THEN
        str="X"//TRIM(int_to_char(it))//"."
        RETURN
    ENDIF
    !
    csigma = 0._dp
    DO jt = 1,nt
      same = .true.
      DO is = 1,nat3
        same=same.and.(ABS(sigma(1,jt)-sigma(is,jt))<eps12)
      ENDDO
      IF(same) csigma(jt) = sigma(1,jt)
    ENDDO
    csigma = 2*csigma
    str = write_temperature(it,nt,csigma)
    
  END FUNCTION
  
END MODULE more_constants
!