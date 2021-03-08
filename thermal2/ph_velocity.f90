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
  REAL(DP),PARAMETER :: h = 1.e-7_dp
  
!  INTERFACE velocity
!    MODULE PROCEDURE velocity_fdiff
!  END INTERFACE  
  !
  PUBLIC :: velocity
  !_simple, velocity_fdiff, velocity_ft
         
  CONTAINS
  FUNCTION velocity(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, fftinterp_dmat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    REAL(DP) :: velocity(3,S%nat3)
    !
    ! If we have effective charges, use finite differences derivation
    IF(S%lrigid)THEN
      velocity = velocity_fdiff(S,fc, xq) 
    ELSE
    ! Otherwise, use the property of Fourier transform to get dD2/dq = \sum_R R e(iqR) F_R
      velocity = velocity_ft(S,fc, xq) 
    ENDIF
  END FUNCTION 
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute ph group velocity by diagonalizing D2 at q:
  !     D u_i = w^2_i u_i ; U = (u_1 u_2 .. u_2*nat)
  ! then at q+h and q-h instead of diagonalizing we just rotate:
  !     W(q+h) = U(q)^H D(q+h) U(q)
  !     w^2_i(q+h) = W(q+h)_i,i
  ! This algorithm does not fail at band crossings
  FUNCTION velocity_fdiff(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: velocity_fdiff(3,S%nat3)
    !
    REAL(DP),ALLOCATABLE :: xvel(:,:)
    COMPLEX(DP),ALLOCATABLE :: D2(:,:), U(:,:), W(:,:)
    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:), w2(:)
    !
    INTEGER :: ix, nu
    REAL(DP) :: xqp(3), xqm(3), dh
    !
    ALLOCATE(D2(S%nat3,S%nat3), U(S%nat3,S%nat3), W(S%nat3,S%nat3), &
             w2p(S%nat3), w2m(S%nat3), w2(S%nat3), xvel(3,S%nat3))
    !
    CALL fftinterp_mat2(xq, S, fc, U)
    CALL mat2_diag(S%nat3, U, w2)
    !
    xvel = 0._dp
    !
! NOTE: the rotation U must be shared!
!$OMP PARALLEL DO DEFAULT(shared) &
!$OMP             PRIVATE(ix, nu, xqp, xqm, w2p, w2m, D2, W, dh) &
!$OMP             REDUCTION(+:xvel)
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
      !
      dh = ( xqp(ix)-xqm(ix) ) *S%tpiba
      FORALL (nu = 1:S%nat3)
        xvel(ix, nu) = ( w2p(nu)-w2m(nu) ) / dh
      END FORALL
      !
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degenerate_velocity(S%nat3, xvel, w2)
    velocity_fdiff = xvel
    DEALLOCATE(D2, U, W, w2p, w2m, w2, xvel)
    !
  END FUNCTION velocity_fdiff
  !
  ! \/o\________\\\_________________________________________/^>
  PURE FUNCTION rotate_d2(nat3, D, U)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    COMPLEX(DP),INTENT(in) :: D(nat3,nat3), U(nat3,nat3)
    COMPLEX(DP) :: rotate_d2(nat3,nat3)
    rotate_d2 = MATMUL(TRANSPOSE(CONJG(U)), MATMUL(D,U))
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  ! Compute the derivative of the dynamical matrix using properties
  ! of fourier transform
  FUNCTION velocity_ft(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, fftinterp_dmat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: velocity_ft(3,S%nat3)
    !
    REAL(DP),ALLOCATABLE :: xvel(:,:)
    COMPLEX(DP),ALLOCATABLE :: U(:,:), rD(:,:), dD(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:)
    INTEGER :: ix, nu
    !
    ALLOCATE(U(S%nat3,S%nat3), rD(S%nat3,S%nat3), w2(S%nat3), &
             xvel(3,S%nat3), dD(S%nat3,S%nat3,3))
    !
    ! We need to get and diagonalize the dyn.mat. at q to get its eigenvectors
    ! (aka the rotation U that makes it diagonal)
    CALL fftinterp_mat2(xq, S, fc, U)
    CALL mat2_diag(S%nat3, U, w2)
    !
    ! The next call gives us the derivative of the dynamical matrix
    CALL fftinterp_dmat2(xq, S, fc, dD)
    !
    xvel = 0._dp
    !
    DO ix = 1,3
      ! Instead of diagonalizing dD, we rotate it with the same patterns as D
      rD = rotate_d2(S%nat3, dD(:,:,ix), U)
      ! The diagonal terms are the derivatives of the SQUARE of the frequencies
      ! to get the derivatives of the frequencies, we need to multiply by
      ! 1/(2 \omega)
      DO nu = 1, S%nat3
      IF(ABS(w2(nu)) > 0._dp)THEN
        xvel(ix,nu) = 0.5_dp*REAL(rD(nu,nu),kind=dp)/SIGN(DSQRT(ABS(w2(nu))),w2(nu))
      ELSE
        ! we set at zero velocity at Gamma for the acoustic branches
        ! (it is not well defined anyway)
        xvel(ix,nu) = 0._dp
      ENDIF
      ENDDO
    ENDDO
    !
    CALL merge_degenerate_velocity(S%nat3, xvel, w2)
    velocity_ft = xvel
    DEALLOCATE(U, xvel)

  END FUNCTION velocity_ft
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the derivative of the dynamical matrix on the basis
  ! of the phonon eigenvectors
  FUNCTION velocity_matrix_ft(S,fc, xq) &
  RESULT (rD)
    USE fc2_interpolate, ONLY : fftinterp_mat2, fftinterp_dmat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    REAL(DP) :: rD(3,S%nat3,S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:), dD(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:)
    INTEGER :: ix
    !
    ALLOCATE(U(S%nat3,S%nat3), w2(S%nat3), dD(S%nat3,S%nat3,3))
    !
    ! We need to get and diagonalize the dyn.mat. at q to get its eigenvectors
    ! (aka the rotation U that makes it diagonal)
    CALL fftinterp_mat2(xq, S, fc, U)
    CALL mat2_diag(S%nat3, U, w2)
    !
    ! The next call gives us the derivative of the dynamical matrix
    CALL fftinterp_dmat2(xq, S, fc, dD)
    !
    !
    DO ix = 1,3
      ! Instead of diagonalizing dD, we rotate it with the same patterns as D
      rD(ix,:,:) = rotate_d2(S%nat3, dD(:,:,ix), U)
    ENDDO
    !
    !CALL merge_degenerate_velocity(S%nat3, xvel, w2)
    DEALLOCATE(U,w2,dD)

  END FUNCTION velocity_matrix_ft
  !
  END MODULE ph_velocity














