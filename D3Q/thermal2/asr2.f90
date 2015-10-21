!
! Written by Lorenzo Paulatto (2013-2014) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
MODULE asr2_module
  USE kinds,    ONLY : DP
  
  CONTAINS
  SUBROUTINE impose_asr2(method, nat,fc)
    USE input_fc,       ONLY : forceconst2_grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    CHARACTER(len=*),INTENT(in) :: method
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    !
    INTEGER :: iR, a,b, i,j, nu,mu
    REAL(DP):: delta, norm
    !
    IF(method=="no")RETURN
    !
    DO a = 1,3
    DO b = 1,3
      DO i = 1,nat
      nu = 3*(i-1)+a
        !
        delta = 0._dp
        norm = 0._dp
        !
        DO j = 1,nat
        mu = 3*(j-1)+b
        DO iR = 1,fc%n_R
          delta = delta+fc%FC(mu,nu,iR)
          norm  = norm +fc%FC(mu,nu,iR)**2 
        ENDDO
        ENDDO
        
        IF(method=="simple")THEN
          mu = 3*(i-1)+b
          fc%FC(mu,nu,fc%i_0) = fc%FC(mu,nu,fc%i_0) - delta
        ELSEIF(method=="diagonal")THEN
          mu = 3*(i-1)+b
          norm=0._dp
          DO iR = 1,fc%n_R
            norm  = norm +fc%FC(mu,nu,iR)**2
          ENDDO
          DO iR = 1,fc%n_R
            fc%FC(mu,nu,iR) = fc%FC(mu,nu,iR) - delta*fc%FC(mu,nu,iR)**2/norm
          ENDDO
        ELSEIF(method=="spread")THEN
          norm = 0._dp
          DO j = 1,nat
          mu = 3*(j-1)+b
          DO iR = 1,fc%n_R
            norm  = norm +fc%FC(mu,nu,iR)**2 
          ENDDO
          ENDDO
          !  
          DO j = 1,nat
          mu = 3*(j-1)+b
          DO iR = 1,fc%n_R
            fc%FC(mu,nu,iR) = fc%FC(mu,nu,iR) - delta*fc%FC(mu,nu,iR)**2/norm
          ENDDO
          ENDDO
        ELSE
          CALL errore("impose_asr2", "unknown method "//TRIM(method), 1)
        ENDIF
       
      ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE impose_asr2
  ! \/o\________\\\_________________________________________/^>

END MODULE asr2_module






