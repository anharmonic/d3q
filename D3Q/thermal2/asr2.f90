!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
MODULE asr2_module
  USE kinds,    ONLY : DP
  
  CONTAINS
  SUBROUTINE impose_asr2(nat,fc)
    USE input_fc,       ONLY : forceconst2_grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    !
    INTEGER :: iR, a,b, i,j, nu,mu
    REAL(DP):: delta
    !
    !
    DO a = 1,3
    DO b = 1,3
      DO i = 1,nat
      nu = 3*(i-1)+a
        !
        delta = 0._dp
        !
        DO j = 1,nat
        mu = 3*(j-1)+b
        DO iR = 1,fc%n_R
          delta = delta+fc%FC(mu,nu,iR)
        ENDDO
        ENDDO

        mu = 3*(i-1)+b
        fc%FC(mu,nu,fc%i_0) = fc%FC(mu,nu,fc%i_0) - delta
        !
      ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE impose_asr2
  ! \/o\________\\\_________________________________________/^>

END MODULE asr2_module






