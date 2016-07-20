!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE ph_velocity

  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info
                       
  PRIVATE
  REAL(DP),PARAMETER :: h = 1.e-5_dp
         
  PUBLIC :: velocity_simple, velocity_proj, velocity_var
         
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute ph group velocity by simple straightforward finite difference
  ! This algorithm can fail at band crossings
  FUNCTION velocity_simple(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: velocity_simple(3,S%nat3)
    !
    REAL(DP) :: xvel(3,S%nat3)
    COMPLEX(DP),ALLOCATABLE :: D2(:,:)
    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:)
    !
    INTEGER :: ix, nu
    REAL(DP) :: xqp(3), xqm(3), dh
    !
    CALL errore("velocity_simple","This algorithm is wrong and should not be used",1)
    ALLOCATE(D2(S%nat3,S%nat3), w2p(S%nat3), w2m(S%nat3))
    !
    xvel = 0._dp
    !
! NOTE: not using OMP here because it is used in fftinterp_mat2
!-!$OMP PARALLEL DO DEFAULT(shared) &
!-!$OMP             PRIVATE(ix, nu, xqp, xqm, w2p, w2m, D2, dh) &
!-!$OMP             REDUCTION(+:xvel)
    DO ix = 1,3
      xqp = xq
      xqp(ix) = xq(ix)+h
      CALL fftinterp_mat2(xqp, S, fc, D2)
      CALL mat2_diag(S%nat3, D2, w2p)
      !
      xqm = xq
      xqm(ix) = xq(ix)-h
      CALL fftinterp_mat2(xqm, S, fc, D2)
      CALL mat2_diag(S%nat3, D2, w2m)
      !
      dh = ( xqp(ix)-xqm(ix) )*S%tpiba
      FORALL (nu = 1:S%nat3)
        xvel(ix, nu) = ( SQRT(w2p(nu))-SQRT(w2m(nu)) ) / dh
      END FORALL
      !
    ENDDO
!-!$OMP END PARALLEL DO
    !
    velocity_simple = xvel
    !
    DEALLOCATE(D2,w2p, w2m)

  END FUNCTION velocity_simple
  ! \/o\________\\\_________________________________________/^>
  ! Compute ph group velocity by diagonalizing D2 at q:
  !     D u_i = w^2_i u_i ; U = (u_1 u_2 .. u_2*nat)
  ! then at q+h and q-h instead of diagonalizing we just rotate:
  !     W(q+h) = U(q)^H D(q+h) U(q)
  !     w^2_i(q+h) = W(q+h)_i,i
  ! This algorithm does not fail at band crossings
  FUNCTION velocity_proj(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: velocity_proj(3,S%nat3)
    !
    REAL(DP),ALLOCATABLE :: xvel(:,:)
    COMPLEX(DP),ALLOCATABLE :: D2(:,:), U(:,:), W(:,:)
    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:)
    !
    INTEGER :: ix, nu
    REAL(DP) :: xqp(3), xqm(3), dh
    !
    ALLOCATE(D2(S%nat3,S%nat3), U(S%nat3,S%nat3), W(S%nat3,S%nat3), &
             w2p(S%nat3), w2m(S%nat3), xvel(3,S%nat3))
    !
    CALL fftinterp_mat2(xq, S, fc, U)
    CALL mat2_diag(S%nat3, U, w2p)
    !
    xvel = 0._dp
    !
! NOTE: the rotation U must be shared!
! NOTE2: not using OMP here because it is used in fftinterp_mat2
!-!$OMP PARALLEL DO DEFAULT(shared) &
!-!$OMP             PRIVATE(ix, nu, xqp, xqm, w2p, w2m, D2, W, dh) &
!-!$OMP             REDUCTION(+:xvel)
    DO ix = 1,3
      xqp = xq
      xqp(ix) = xq(ix)+h
      CALL fftinterp_mat2(xqp, S, fc, D2)
      W = rotate_d2(S%nat3, D2, U)
      FORALL(nu = 1:S%nat3) w2p(nu) = REAL(W(nu,nu),kind=DP)
      WHERE(w2p>=0._dp)
        w2p = DSQRT(w2p)
      ELSEWHERE
        w2p = -DSQRT(-w2p)
      END WHERE
      
!       DO nu = 1,S%nat3
!         w2p(nu) = REAL(W(nu,nu),kind=DP)
!         IF(w2p(nu)>=0._dp) THEN
!           w2p(nu) = SQRT(w2p(nu))
!         ELSE
!           w2p(nu) = -SQRT(-w2p(nu))
!         ENDIF
!       ENDDO
      !
      xqm = xq
      xqm(ix) = xq(ix)-h
      CALL fftinterp_mat2(xqm, S, fc, D2)
      W = rotate_d2(S%nat3, D2, U)
      FORALL(nu = 1:S%nat3) w2m(nu) = REAL(W(nu,nu),kind=DP)
      WHERE(w2m>=0._dp)
        w2m = DSQRT(w2m)
      ELSEWHERE
        w2m = -DSQRT(-w2m)
      END WHERE
      !       DO nu = 1,S%nat3
!         w2m(nu) = REAL(W(nu,nu),kind=DP)
!         IF(w2m(nu)>=0._dp) THEN
!           w2m(nu) = SQRT(w2m(nu))
!         ELSE
!           w2m(nu) = -SQRT(-w2m(nu))
!         ENDIF
!       ENDDO
      !
      dh = ( xqp(ix)-xqm(ix) ) *S%tpiba
      FORALL (nu = 1:S%nat3)
        xvel(ix, nu) = ( w2p(nu)-w2m(nu) ) / dh
      END FORALL
      !
    ENDDO
!-!$OMP END PARALLEL DO
    !
    velocity_proj = xvel
    DEALLOCATE(D2, U, W, w2p, w2m, xvel)
    !
  END FUNCTION velocity_proj
  ! \/o\________\\\_________________________________________/^>
  PURE FUNCTION rotate_d2(nat3, D, U)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    COMPLEX(DP),INTENT(in) :: D(nat3,nat3), U(nat3,nat3)
    COMPLEX(DP) :: rotate_d2(nat3,nat3)
    rotate_d2 = matmul(transpose(conjg(U)), matmul(D,U))
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  ! As in Fugallo et. al. PRB 
  SUBROUTINE velocity_var(S,fc, xq, xvel)
    USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    REAL(DP),INTENT(out) :: xvel(3, S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: D2(:,:)
    REAL(DP),ALLOCATABLE    :: w2(:)

    ALLOCATE(D2(S%nat3,S%nat3), w2(S%nat3))
    xvel = 0._dp
    CALL errore('velocity_var', 'velocity_var not implemented', 1)
    DEALLOCATE(D2,w2)

  END SUBROUTINE velocity_var
  
END MODULE ph_velocity

