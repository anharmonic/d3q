!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE gruneisen_module

  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info
                       
  PRIVATE
  REAL(DP),PARAMETER :: h = 1.e-7_dp
  
!  INTERFACE velocity
!    MODULE PROCEDURE gruneisen
!  END INTERFACE  
  !
  PUBLIC :: gruneisen
         
  CONTAINS
 !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION gruneisen(S0,Sp,Sm ,fc0,fcp,fcm, xq)
    USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    USE functions,        ONLY : rotate_d2, is_gamma
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc0,fcp,fcm
    TYPE(ph_system_info),INTENT(in)   :: S0, Sp, Sm
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: gruneisen(S0%nat3)
    !
    COMPLEX(DP) :: D2(S0%nat3,S0%nat3), U(S0%nat3,S0%nat3), W(S0%nat3,S0%nat3)
    REAL(DP)    :: w2p(S0%nat3), w2m(S0%nat3), w2(S0%nat3), grun(S0%nat3)
    !
    INTEGER :: nu
    REAL(DP) :: domega
    !
    domega = (Sp%omega-Sm%omega)
!    IF(  ABS(ABS(Sp%omega+Sm%omega)/2-S0%omega)>1.d-6 ) &
!      CALL errore("gruneisen","+ and - calculations do not average to the unperturbed one",1)
    !
    CALL fftinterp_mat2(xq, S0, fc0, U)
    CALL mat2_diag(S0%nat3, U, w2)
    WHERE(w2>=0._dp)
      w2 = DSQRT(w2)
    ELSEWHERE
      w2 = -DSQRT(-w2)
    END WHERE
    !
    !
    CALL fftinterp_mat2(xq, Sp, fcp, D2)
    W = rotate_d2(S0%nat3, D2, U)
    FORALL(nu = 1:S0%nat3) w2p(nu) = REAL(W(nu,nu),kind=DP)
    WHERE(w2p>=0._dp)
      w2p = DSQRT(w2p)
    ELSEWHERE
      w2p = -DSQRT(-w2p)
    END WHERE

    CALL fftinterp_mat2(xq, Sm, fcm, D2)
    W = rotate_d2(S0%nat3, D2, U)
    FORALL(nu = 1:S0%nat3) w2m(nu) = REAL(W(nu,nu),kind=DP)
    WHERE(w2m>=0._dp)
      w2m = DSQRT(w2m)
    ELSEWHERE
      w2m = -DSQRT(-w2m)
    END WHERE
    !
    grun = 0._dp
    WHERE (w2/=0._dp)
      grun = - S0%omega/w2 * (w2p-w2m)/domega
    ELSEWHERE
      grun = 0._dp
    END WHERE
    !
    IF(is_gamma(xq, S0%at)) grun(1:3) = 0._dp 
    gruneisen = grun
    !
  END FUNCTION gruneisen
  !
END MODULE gruneisen_module














