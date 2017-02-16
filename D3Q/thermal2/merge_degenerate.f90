!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE merge_degenerate
  USE kinds,      ONLY : DP
  USE constants,  ONLY : RY_TO_CMM1
  USE timers,     ONLY : t_merged
#include "mpi_thermal.h"
  IMPLICIT NONE
  ! 
  !NOTE: pay attention that the number of modes must be the rightmost dimension
  !      of the array you wish to simmetrize! (or you will have to do a loop) 
  !
  PRIVATE
  ! Consider two bands degenerate when they differ by less than 10^-3 cm^-1
  REAL(DP),PARAMETER :: eps_freq = (1.e-3_dp/RY_TO_CMM1)
  LOGICAL,PARAMETER :: disable_merge = .false.
  LOGICAL,PARAMETER :: disable_merge_vel = .false.
  !
  PUBLIC :: merge_degen, merge_degenerate_velocity
  !
  INTERFACE merge_degen
    MODULE PROCEDURE merge_degenerate_real
    MODULE PROCEDURE merge_degenerate_cmplx
    MODULE PROCEDURE merge_degenerate_real_vec
    MODULE PROCEDURE merge_degenerate_cmplx_vec
    MODULE PROCEDURE merge_degenerate_real_mat
  END INTERFACE
  !
  CONTAINS

  ! Average the lw of degenerate modes
  SUBROUTINE merge_degenerate_real(nat3, lw, w)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    REAL(DP),INTENT(inout) :: lw(nat3)
    REAL(DP),INTENT(in)    :: w(nat3)
    !
    INTEGER :: i,j,k
    REAL(DP) :: avg_lw
    
    IF(disable_merge) RETURN
    timer_CALL t_merged%start()
    j=0
    DO i = 1, nat3-1
      ! Skip modes that have already been treated
      IF(i>j)THEN
        ! Find how many modes are degenerate with this one
        DO j = i+1,nat3
          ! As soon as you find one that is not degenerate, stop seeking
          !print*,i,j,ABS(w(i)-w(j))*RY_TO_CMM1,eps_freq
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        ! Go back to the highest degenerate one
        ! (NOTE: j will be nat3+1 if the last mode is degenerate)
        j=j-1
        ! Average the lw over the degenerate modes
        IF(j>i)THEN
          !print*, i,j
          avg_lw = SUM(lw(i:j))
          lw(i:j) = avg_lw/(j-i+1)
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_real
  !
  ! Same as merge_degenerate_real but for a complex quantity (i.e. self-energy)
  SUBROUTINE merge_degenerate_cmplx(nat3, lw, w)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    COMPLEX(DP),INTENT(inout) :: lw(nat3)
    REAL(DP),INTENT(in)       :: w(nat3)
    !
    INTEGER :: i,j,k
    COMPLEX(DP) :: avg_lw
    IF(disable_merge) RETURN
    timer_CALL t_merged%start()
    j=0
    DO i = 1, nat3-1
      ! Skip modes that have already been treated
      IF(i>j)THEN
        ! Find how many modes are degenerate with this one
        DO j = i+1,nat3
          ! As soon as you find one that is not degenerate, stop seeking
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        ! Go back to the highest degenerate one
        ! (NOTE: j will be nat3+1 if the last mode is degenerate)
        j=j-1
        ! Average the lw over the degenerate modes
        IF(j>i)THEN
          avg_lw = SUM(lw(i:j))
          lw(i:j) = avg_lw/(j-i+1)
          !print*,"found degen", i,j, avg_lw/(j-i+1)*RY_TO_CMM1
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_cmplx
  !
  ! Same as merge_degenerate_real, but acts on a vector
  SUBROUTINE merge_degenerate_real_vec(ndim, nat3, vec, w)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndim,nat3
    REAL(DP),INTENT(inout) :: vec(ndim,nat3)
    REAL(DP),INTENT(in)     :: w(nat3)
    !
    INTEGER :: i,j,k
    REAL(DP) :: avg_vec(ndim)
    IF(disable_merge) RETURN
    timer_CALL t_merged%start()
    j=0
    DO i = 1, nat3-1
      IF(i>j)THEN
        DO j = i+1,nat3
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        j=j-1
        IF(j>i)THEN
          avg_vec = 0._dp
          DO k = i,j
            avg_vec = avg_vec + vec(:,k)
          ENDDO
          avg_vec = avg_vec/(j-i+1)
          DO k = i,j
            vec(:,k) = avg_vec 
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_real_vec
  !
  ! Same as merge_degenerate_real_vec, but for a vector of complex numbers
  SUBROUTINE merge_degenerate_cmplx_vec(ndim, nat3, vec, w)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndim,nat3
    COMPLEX(DP),INTENT(inout) :: vec(ndim,nat3)
    REAL(DP),INTENT(in)       :: w(nat3)
    INTEGER :: i,j,k
    COMPLEX(DP) :: avg_vec(ndim)
    IF(disable_merge) RETURN
    timer_CALL t_merged%start()
    j=0
    DO i = 1, nat3-1
      IF(i>j)THEN
        DO j = i+1,nat3
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        j=j-1
        IF(j>i)THEN
          avg_vec = 0._dp
          DO k = i,j
            avg_vec = avg_vec + vec(:,k)
          ENDDO
          avg_vec = avg_vec/(j-i+1)
          DO k = i,j
            vec(:,k) = avg_vec 
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_cmplx_vec
  !
  ! Same as merge_degenerate_real, but acts on a vector
  SUBROUTINE merge_degenerate_real_mat(ndim,mdim, nat3, mat, w)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndim,mdim,nat3
    REAL(DP),INTENT(inout) :: mat(ndim,mdim,nat3)
    REAL(DP),INTENT(in)     :: w(nat3)
    INTEGER :: i,j,k
    REAL(DP) :: avg_mat(ndim,mdim)
    IF(disable_merge) RETURN
    timer_CALL t_merged%start()
    j=0
    DO i = 1, nat3-1
      IF(i>j)THEN
        DO j = i+1,nat3
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        j=j-1
        IF(j>i)THEN
          avg_mat = 0._dp
          DO k = i,j
            avg_mat = avg_mat + mat(:,:,k)
          ENDDO
          avg_mat = avg_mat/(j-i+1)
!          ioWRITE(stdout,*) "merging", avg_mat(1,1), mat(1,1,i:j)
          DO k = i,j
            mat(:,:,k) = avg_mat 
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_real_mat
  !
  ! Average the velocity coming from degenerate modes
  ! SPECIAL VERSION for VELOCITY ONLY: takes \omega^2 in input instead of \omega
  !                                    assumes vector to have length 3
  SUBROUTINE merge_degenerate_velocity(nat3, vel, w2)
    USE kinds, ONLY : DP
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    REAL(DP) :: vel(3,nat3)
    REAL(DP) :: w2(nat3)
    !
    INTEGER :: i,j,k
    REAL(DP) :: avg_vel(3), w(nat3)
    
    ! NOTE: checking |a^2-b^2| < epsilon^2 is a completely different
    !       thing than checking |a-b| < epsilon
    IF(disable_merge_vel) RETURN
    timer_CALL t_merged%start()
    w= DSQRT(ABS(w2))
    j=0
    DO i = 1, nat3-1
      ! Skip the modes that have already been treated
      IF(i>j)THEN
        ! Find how many modes are degenerate with this one
        DO j = i+1,nat3
          ! As soon as you find one that is not degenerate, stop seeking
          IF(ABS(w(i)-w(j))>eps_freq) EXIT
        ENDDO
        ! Go back to the highest degenerate one
        ! (NOTE: j will be nat3+1 if the last mode is degenerate)
        j=j-1
        ! Average the velocity over the degenerate modes
        IF(j>i)THEN
          avg_vel = 0._dp
          DO k = i,j
            avg_vel = avg_vel + vel(:,k)
          ENDDO
          avg_vel = avg_vel/(j-i+1)
          DO k = i,j
            vel(:,k) = avg_vel 
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    timer_CALL t_merged%stop()
  END SUBROUTINE merge_degenerate_velocity
  !
  END MODULE
