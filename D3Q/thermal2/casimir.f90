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
  FUNCTION casimir_linewidth_q(xq0, l_sample, sample_dir, S, fc2) &
  RESULT(lw)
    !
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE ph_velocity, ONLY : velocity
    !
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq0(3)
    REAL(DP),INTENT(in) :: l_sample
    REAL(DP),INTENT(in) :: sample_dir(3)
    REAL(DP) :: vel(3,S%nat3)
    !
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(S%nat3)
    !
    vel = velocity(S, fc2, xq0)
    lw = casimir_linewidth_vel(vel, l_sample, sample_dir, S%nat3)
     !
  END FUNCTION casimir_linewidth_q

  ! Returns the HALF width half maximum (note the 0.5 factor) of phonon modes
  ! due to Casimir scattering with grain boundary or sample boundary
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION casimir_linewidth_vel(vel, l_sample, sample_dir, nat3) &
  RESULT(lw)
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: l_sample
    REAL(DP),INTENT(in) :: sample_dir(3)
    REAL(DP),INTENT(in) :: vel(3,nat3)
    !
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(nat3)
    REAL(DP) :: inv_lcas
    INTEGER  :: nu
    REAL(DP),PARAMETER :: eps = 1.e-10_dp
    !
    ! Case 1: Velocity projected along sample_dir
    IF (SUM(ABS(sample_dir)) > eps) THEN
      inv_lcas = 0.5_dp / l_sample / DSQRT(SUM(sample_dir**2))
      DO nu = 1,nat3
          lw(nu) = inv_lcas * ABS(DOT_PRODUCT(sample_dir,vel(:,nu)))
      ENDDO
    !
    ! Case 2: Modulus of velocity (as in PRB 87, 214303 (2013))
    ELSE
      inv_lcas = 0.5_dp / l_sample
      DO nu = 1,nat3
          lw(nu) = inv_lcas * DSQRT(SUM(vel(:,nu)**2))
      ENDDO
    ENDIF
    !
  END FUNCTION casimir_linewidth_vel
  !
  ! Returns .true. if the phonon mean free path is longer than l_sample
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION mfp_scatter_vel(vel, inverse_lifetime, l_sample, sample_dir) &
  RESULT(scatter)
    !
    IMPLICIT NONE
    !
    ! inverse_lifetime should be 1/(2*gamma)
    REAL(DP),INTENT(in) :: vel(3)
    REAL(DP),INTENT(in) :: inverse_lifetime
    REAL(DP),INTENT(in) :: l_sample
    REAL(DP),INTENT(in) :: sample_dir(3)
    ! FUNCTION RESULT:
    LOGICAL :: scatter
    REAL(DP) :: effective_vel
    !
    REAL(DP),PARAMETER :: eps = 1.e-10_dp
    !
!     IF(inverse_lifetime==0._dp)THEN
!       scatter = .true.
!       RETURN
!     ENDIF
    !
    ! Case 1: Velocity projected along sample_dir
    IF (SUM(ABS(sample_dir)) > eps) THEN
      effective_vel = ABS(SUM( vel*sample_dir ))
    !
    ! Case 2: Modulus of velocity
    ELSE
      effective_vel = DSQRT(SUM( vel**2 ))
    ENDIF
    !
    IF (effective_vel > l_sample*inverse_lifetime) THEN
      scatter = .true.
    ELSE
      scatter = .false.
    ENDIF
    !
  END FUNCTION mfp_scatter_vel
  !
END MODULE casimir_linewidth
