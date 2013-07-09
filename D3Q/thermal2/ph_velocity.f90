!
! Copyright Lorenzo Paulatto 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE ph_velocity

  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       ph_system_info
                       
  REAL(DP),PARAMETER :: h = 1.e-5_dp
                       
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute ph group velocity by simple straightforward finite difference
  ! This algorithm can fail at band crossings
  SUBROUTINE velocity_simple(S,fc, xq, xvel)
    USE interp_fc, ONLY : fftinterp_mat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    REAL(DP),INTENT(out) :: xvel(3,S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: D2(:,:)
    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:)
    !
    INTEGER :: ix, nu
    REAL(DP) :: xqp(3), xqm(3), dh
    !
    ALLOCATE(D2(S%nat3,S%nat3), w2p(S%nat3), w2m(S%nat3))
    !
    DO ix = 1,3
      xqp = xq
      xqp(ix) = xq(ix)+h
      CALL fftinterp_mat2(xqp, S, fc, D2)
      CALL mat2_diag(S, D2, w2p)
      !
      xqm = xq
      xqm(ix) = xq(ix)-h
      CALL fftinterp_mat2(xqm, S, fc, D2)
      CALL mat2_diag(S, D2, w2m)
      !
      dh = ( xqp(ix)-xqm(ix) )
      FORALL (nu = 1:S%nat3)
        xvel(ix, nu) = ( SQRT(w2p(nu))-SQRT(w2m(nu)) ) / dh
      END FORALL
      !
    ENDDO
    !
    DEALLOCATE(D2,w2p, w2m)

  END SUBROUTINE velocity_simple
  ! \/o\________\\\_________________________________________/^>
  ! Compute ph group velocity by diagonalizing D2 at q:
  !     D u_i = w^2_i u_i ; U = (u_1 u_2 .. u_2*nat)
  ! then at q+h and q-h instead of diagonalizing we just rotate:
  !     W(q+h) = U(q)^H D(q+h) U(q)
  !     w^2_i(q+h) = W(q+h)_i,i
  ! This algorithm does not fail at band crossings
  SUBROUTINE velocity_proj(S,fc, xq, xvel)
    USE interp_fc, ONLY : fftinterp_mat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    REAL(DP),INTENT(out) :: xvel(3,S%nat3)
    !
    COMPLEX(DP),ALLOCATABLE :: D2(:,:), U(:,:), W(:,:)
    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:)
    !
    INTEGER :: ix, nu
    REAL(DP) :: xqp(3), xqm(3), dh
    !
    ALLOCATE(D2(S%nat3,S%nat3), U(S%nat3,S%nat3), W(S%nat3,S%nat3), &
             w2p(S%nat3), w2m(S%nat3))
    !
    CALL fftinterp_mat2(xq, S, fc, U)
    CALL mat2_diag(S, U, w2p)
    !
    DO ix = 1,3
      xqp = xq
      xqp(ix) = xq(ix)+h
      CALL fftinterp_mat2(xqp, S, fc, D2)
      W = rotate(S%nat3, D2, U)
      FORALL(nu = 1:S%nat3) w2p(nu) = REAL(W(nu,nu),kind=DP)
      !
      xqm = xq
      xqm(ix) = xq(ix)-h
      CALL fftinterp_mat2(xqm, S, fc, D2)
      W = rotate(S%nat3, D2, U)
      FORALL(nu = 1:S%nat3) w2m(nu) = REAL(W(nu,nu),kind=DP)
      !
      dh = ( xqp(ix)-xqm(ix) )
      FORALL (nu = 1:S%nat3)
        xvel(ix, nu) = ( SQRT(w2p(nu))-SQRT(w2m(nu)) ) / dh
      END FORALL
      !
    ENDDO
    !
    DEALLOCATE(D2, U, W, w2p, w2m)

  END SUBROUTINE velocity_proj
  ! \/o\________\\\_________________________________________/^>
  FUNCTION rotate(nat3, D, U)
    COMPLEX(DP),INTENT(in) :: D(nat3,nat3), U(nat3,nat3)
    COMPLEX(DP) :: rotate(nat3,nat3)
    rotate = matmul(transpose(conjg(U)), matmul(D,U))
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  ! As in Fugallo et. al. PRB 
  SUBROUTINE velocity_var(S,fc, xq, xvel)
    USE interp_fc, ONLY : fftinterp_mat2, mat2_diag
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

