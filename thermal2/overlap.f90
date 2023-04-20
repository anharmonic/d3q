
MODULE overlap
  !
  USE kinds, ONLY : DP
#include "mpi_thermal.h"

  TYPE :: order_type
    COMPLEX(DP),ALLOCATABLE,PRIVATE :: e_prev(:,:)
    INTEGER,ALLOCATABLE     :: idx(:)
    REAL(DP) :: xq_ref(3)
    CONTAINS
      PROCEDURE :: set => set_idx
      PROCEDURE :: set_path => set_path_idx
      PROCEDURE :: set_ref => set_idx_reference
      PROCEDURE :: set_xq_ref => set_xq_ref
      PROCEDURE :: reset => reset_idx
  END TYPE
  !

  CONTAINS
  !
  SUBROUTINE reset_idx(self)
    IMPLICIT NONE
    CLASS(order_type),INTENT(inout) :: self
    DEALLOCATE(self%idx, self%e_prev)
    RETURN
  END SUBROUTINE
  
  ! Specialized version of set_idx, that also checks the path length:
  ! if the path length at this point is not the same than at the previous point,
  ! it means that there is a discontinuity, do not attempt to compute the overlap
  ! just keep the current order
  ! Same story when the length is reset to zero
  SUBROUTINE set_path_idx(self, nat3, w, e, iq, nq, lpath, xq) !RESULT(order_out)
    USE constants, ONLY : RY_TO_CMM1
    !USE merge_degenerate, ONLY : merge_degen
    IMPLICIT NONE
    CLASS(order_type),INTENT(inout) :: self
    !
    INTEGER,INTENT(in) :: nat3
    REAL(DP),INTENT(in) :: w(nat3)
    COMPLEX(DP),INTENT(in) :: e(nat3,nat3)
    !
    INTEGER,INTENT(in) :: iq,nq
    REAL(DP),INTENT(in) ::  lpath(nq), xq(3,nq)
    LOGICAL :: turning
    INTEGER :: i
    !
    ! Do not change the order of the bands, under any circumstance when:
    ! 1. the path has just turned
    ! 2. the path is jumping q->q+G (length is the same as previous point) 
    ! 3. the path has been reset (length is zero)
    ! 4. we're at Gamma
    IF(iq>1 .and. ALLOCATED(self%e_prev))THEN
    
      IF(iq<nq)THEN
        turning = SUM( ((xq(:,iq-1) - xq(:,iq))&
                      -(xq(:,iq) - xq(:,iq+1)))**2 ) > 1.d-6    
      ELSE
         turning = .true.
      ENDIF
      
      !
      IF(lpath(iq-1)==lpath(iq) .or. lpath(iq)==0._dp &
         .or. turning .or. SUM(ABS(xq(:,iq-1)))<1.d-4 ) THEN
        DO i = 1, nat3
          self%e_prev(:,i) = e(:,self%idx(i))
        ENDDO
        !print*, i, xq(:,i), "quick"
        ! quick return
        RETURN
      ENDIF
    ENDIF
    !
    ! otherwise, do the full overlap analysis:
    CALL set_idx(self, nat3, w, e) 
    !
  END SUBROUTINE
  !
  ! Special version of set idx which checks the overlap not at the q-point in question,
  ! but at a small q in the same direction q0 = epsilon * q/|q|. 
  ! This *could* ensure that continuity is consistent at crossings. 
  ! EXPERIMENT ! 
  SUBROUTINE set_idx_reference(self, nat3, xq, S, fc2) !RESULT(order_out)
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE fc2_interpolate,    ONLY : freq_phq, freq_phq_hat
    USE q_grids,            ONLY : q_grid, setup_path, q_grid_destroy

    !USE merge_degenerate, ONLY : merge_degen
    IMPLICIT NONE
    CLASS(order_type),INTENT(inout) :: self
    !
    INTEGER,INTENT(in) :: nat3
    ! REAL(DP),INTENT(in) :: w(nat3)
    ! COMPLEX(DP),INTENT(in) :: e(nat3,nat3)
    !
    REAL(DP),INTENT(in) ::  xq(3)
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    !
    REAL(DP) :: w0(S%nat3), xq0(3), normq
    COMPLEX(DP) :: e0(S%nat3, S%nat3)
    !
    TYPE(q_grid) :: qline
    !
    REAL(DP),PARAMETER :: dq = 0.0005_dp ! bohr^-1
    INTEGER :: i, iq0, nq0
    !
    nq0 = DSQRT(SUM((xq-self%xq_ref)**2))/(dq*S%tpiba)
    !print*, nq0

    CALL setup_path(self%xq_ref, 1, qline, S%at)
    CALL setup_path(xq, nq0, qline, S%at)
    !
    DO iq0 = 1, qline%nq
      CALL freq_phq(qline%xq(:,iq0), S, fc2, w0, e0)
      CALL set_idx(self, nat3, w0, e0, reset=(iq0==1))
    ENDDO
    !
    CALL q_grid_destroy(qline)
    !
  END SUBROUTINE set_idx_reference
  !
  SUBROUTINE set_xq_ref(self, nq, xq, xq_ref) !RESULT(order_out)
    USE constants, ONLY : RY_TO_CMM1
    !USE merge_degenerate, ONLY : merge_degen
    IMPLICIT NONE
    CLASS(order_type),INTENT(inout) :: self
    !
    INTEGER,INTENT(in)  :: nq
    REAL(DP),INTENT(in) :: xq(3,nq), xq_ref(3)
    INTEGER :: i
    !IF(input%sort_freq=="reference")THEN
    IF(ALL(xq_ref==0._dp))THEN
      ! use average point of path to sort frequency, to minimize distance
      DO i = 1,3
        self%xq_ref(i) = SUM(xq(i,:))/nq
      ENDDO
      ! fall-back to 1/4 1/4 1/4 if it was zero
      IF(SUM(ABS(self%xq_ref))<1.d-4) self%xq_ref = 0.25_dp
    ELSE
      self%xq_ref = xq_ref
    ENDIF
    ioWRITE(*,'(2x,a,3f12.6)') "Reference point for sorting phonon bands:", self%xq_ref
    !ENDIF
  END SUBROUTINE 

  !
  SUBROUTINE set_idx(self, nat3, w, e, keep_old_e, reset)
    USE constants, ONLY : RY_TO_CMM1
    !USE merge_degenerate, ONLY : merge_degen
    IMPLICIT NONE
    CLASS(order_type),INTENT(inout) :: self
    !
    INTEGER,INTENT(in) :: nat3
    REAL(DP),INTENT(in) :: w(nat3)
    COMPLEX(DP),INTENT(in) :: e(nat3,nat3)
    LOGICAL,OPTIONAL,INTENT(in) :: keep_old_e, reset
    LOGICAL :: update_e
    !INTEGER :: order_out(nat3)
    !
    COMPLEX(DP) :: e_new(nat3,nat3)
    !
    REAL(DP)   :: olap
    COMPLEX(DP) :: dolap
    REAL(DP)   :: maxx
    INTEGER    :: i,j,j_,jmax, orderx(nat3)
    LOGICAL    :: assigned(nat3)
    REAL(DP),PARAMETER :: olap_thr = 0._dp
    REAL(DP),PARAMETER :: w_thr = 1.e+10_dp/RY_TO_CMM1
    REAL(DP),PARAMETER :: eps = 0._dp
    !
    IF(present(keep_old_e))THEN
      update_e = .not. keep_old_e
    ELSE
      update_e = .true.
    ENDIF
    !
    ! On first call, order is just 1..3*nat, then quick return
    IF(.not.ALLOCATED(self%e_prev))THEN
      ALLOCATE(self%e_prev(nat3,nat3))
      self%e_prev = e
      ALLOCATE(self%idx(nat3))
      FORALL(i=1:nat3) self%idx(i) = i
      !order_out = order
      RETURN
    ELSEIF(present(reset))THEN
        IF(reset)THEN
          self%e_prev = e
          FORALL(i=1:nat3) self%idx(i) = i
          RETURN
        ENDIF
    ELSE
      IF(size(self%e_prev,1)/=nat3) CALL errore("overlap","inconsistent nat3",1)
    ENDIF
    !
    e_new = e
    !CALL merge_degen(nat3, nat3, e_new, w)
    !
    assigned = .false.
    !WRITE(*,'(99i2)') order
    DO i = nat3,1,-1
      maxx = -1._dp
      jmax = -1
      DO j_ = 1,nat3
        j = self%idx(j_)
        IF(assigned(j)) CYCLE
        ! Square modulus:
        dolap = SUM( self%e_prev(:,i)*DCONJG(e_new(:,j)) )
        olap = DBLE(dolap*DCONJG(dolap))
        ! Do not change the order if the overlap condition is weak
        ! i.e. where some modes are degenerate, or if points are too far apart
        IF(olap>(maxx+eps) .and.olap>olap_thr .and. ABS(w(i)-w(j))<w_thr )THEN
        !IF(olap>(maxx+eps) .and.olap>olap_thr )THEN
          jmax = j
          maxx = olap
        ENDIF
        !
      ENDDO
      !
      IF(jmax<0) THEN
        IF (assigned(i)) CALL errore("overlap_sort","dont know what to do",1)
        jmax = i
      ENDIF
      !
      !IF(i/=jmax) print*, i, "->", jmax, maxx
      ! Destroy polarizations already assigned:
      !e_new(:,jmax) = 0._dp
      IF(update_e) self%e_prev(:,i) = e_new(:,jmax)
      !
      assigned(jmax) = .true.
      orderx(i) = jmax !order(jmax)
    ENDDO
    !
    ! Prepare for next call
    !self%e_prev = e
    self%idx = orderx
    !
    ! Set output
    !order_out = orderx
    !
  END SUBROUTINE
  !  
END MODULE
