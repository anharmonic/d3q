!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!

MODULE casimir_linewidth
  !
  USE kinds, ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid
  !
  CONTAINS
  !
  ! Computes Casimir linewidth at a certain q-point, i.e. it computes the velocity
  ! and then calls the other function to compute the linewidth
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION casimir_linewidth_q(xq0, l_casimir, casimir_dir, S, fc2) &
  RESULT(lw)
    !
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE ph_velocity, ONLY : velocity_proj
    !
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq0(3)
    REAL(DP),INTENT(in) :: l_casimir
    REAL(DP),INTENT(in) :: casimir_dir(3)
    REAL(DP) :: vel(3,S%nat3)
    !
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(S%nat3)
    !
    vel = velocity_proj(S, fc2, xq0)
    lw = casimir_linewidth_vel(vel, l_casimir, casimir_dir, S%nat3)
     !
  END FUNCTION casimir_linewidth_q

  ! COmputes the Casimir linewidth from ph group velocity
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION casimir_linewidth_vel(vel, l_casimir, casimir_dir, nat3) &
  RESULT(lw)
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: l_casimir
    REAL(DP),INTENT(in) :: casimir_dir(3)
    REAL(DP),INTENT(in) :: vel(3,nat3)
    !
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(nat3)
    REAL(DP) :: inv_lcas
    INTEGER  :: nu
    REAL(DP),PARAMETER :: eps = 1.e-20_dp
    !
    ! Case 1: Velocity projected along casimir_dir
    IF (SUM(ABS(casimir_dir)) > eps) THEN
      inv_lcas = 1/ l_casimir / DSQRT(SUM(casimir_dir**2))
      DO nu = 1,nat3
          lw(nu) = inv_lcas * ABS(DOT_PRODUCT(casimir_dir,vel(:,nu)))
      ENDDO
    !
    ! Case 2: Modulus of velocity (as in PRB 87, 214303 (2013))
    ELSE
      inv_lcas = 1/ l_casimir
      DO nu = 1,nat3
          lw(nu) = inv_lcas * SQRT(SUM(vel(:,nu)**2))
      ENDDO
      
    ENDIF
          !
  END FUNCTION casimir_linewidth_vel
  !
END MODULE casimir_linewidth
!
