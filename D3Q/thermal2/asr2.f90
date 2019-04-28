!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE asr2_module
  USE kinds,    ONLY : DP
  
  CONTAINS
  SUBROUTINE impose_asr2(method, nat,fc,zeu)
    USE input_fc,       ONLY : forceconst2_grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    CHARACTER(len=*),INTENT(in) :: method
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    REAL(DP),INTENT(inout),OPTIONAL :: zeu(3,3,nat)
    !
    INTEGER :: iR, a,b, i,j, nu,mu, na
    REAL(DP):: delta, norm
    !
    IF(method=="no" .or. method=="none")RETURN
    !
    DO i = 1,nat
      DO a = 1,3
      DO b = 1,3
      nu = 3*(i-1)+a
        !
        delta = 0._dp
        !norm = 0._dp
        !
        DO j = 1,nat
        mu = 3*(j-1)+b
        DO iR = 1,fc%n_R
          delta = delta+fc%FC(mu,nu,iR)
          !norm  = norm +fc%FC(mu,nu,iR)**2 
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
          CALL errore("impose_asr2", "spread asr is buggy/wrong"//TRIM(method), 1)
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
    IF(present(zeu)) THEN
    print*, "doing zeu asr"
     do i=1,3
        do j=1,3
           delta=0._dp
           do na=1,nat
              delta = delta + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - delta/nat
           end do
        end do
     end do
    ENDIF
    !
  END SUBROUTINE impose_asr2
  ! \/o\________\\\_________________________________________/^>

END MODULE asr2_module






