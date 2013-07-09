!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE interp_fc
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid, forceconst3_grid

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat2(xq, S, fc, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:, :) = D(:, :) + phase * fc%fc(:, :, i)
    END DO
  END SUBROUTINE fftinterp_mat2
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat3(xq2,xq3, S, fc, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst3_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq2(3), xq3(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3, S%nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:,:,:) = D(:,:,:) + phase * fc%fc(:,:,:, i)
    END DO
  END SUBROUTINE fftinterp_mat3
  ! \/o\________\\\_________________________________________/^>
  ! IN PLACE diagonalization of D
  SUBROUTINE mat2_diag(S, D, w2)
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    COMPLEX(DP),INTENT(inout) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(out)      :: w2(S%nat)
    !
    INTEGER  :: nb    ! block size
    INTEGER,save :: nat3=-1, lwork=-1 ! aux. var

    INTEGER :: info      ! flag saying if the exec. of libr. routines was ok

    INTEGER,EXTERNAL ::ILAENV ! function which gives block size
    !
    REAL(DP), ALLOCATABLE    :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    IF ( nat3 /= S%nat3 .or. lwork < 0 ) THEN
      !     check for the block size
      nb = ILAENV( 1, 'ZHETRD', 'U', S%nat3, -1, -1, -1 )
      IF (nb<1) nb = MAX(1,S%nat3)
      IF (nb==1 .or. nb>=S%nat3) then
        lwork=2*S%nat3-1
      ELSE
        lwork = (nb+1)*S%nat3
      ENDIF
      !
      IF(nat3 > 0 ) PRINT*, "WARNING! Redoing ILAENV"
      nat3 = S%nat3
      !
    ENDIF
    !
    ALLOCATE(work (lwork))
    ALLOCATE(rwork (3*S%nat3-2))
    !
    CALL ZHEEV('V','U',S%nat3,D,S%nat3,w2,work,lwork,rwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(rwork)
    DEALLOCATE(work)
    !
  END SUBROUTINE mat2_diag

  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE calc_d3mm(d3mm, mat3q, zz1, zz2, zz3, nat3, nat32, nat33)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    REAL(DP),INTENT(out) :: d3mm(nat3, nat3, nat3)
    COMPLEX(DP),INTENT(in) :: mat3q(nat3, nat3, nat3)
    COMPLEX(DP),INTENT(in) :: zz3(nat3, nat3), zz2(nat3, nat3), zz1(nat3, nat3)
    INTEGER,INTENT(in) :: nat3, nat32, nat33

    INTEGER :: im1, im2, im3
    COMPLEX(DP),ALLOCATABLE :: aux1(:,:,:), aux2(:,:,:), out(:), in(:)
    ALLOCATE(aux1(nat3,nat3,nat3), aux2(nat3,nat3,nat3), out(nat3), in(nat3))

    DO im3 = 1, nat3
        DO im2 = 1, nat3
          in =  mat3q(im2, im3, :)
          CALL ZGEMV('T', nat3, nat3, (1._dp,0._dp), zz1, nat3,in, &
                          nat32, (0._dp,0._dp), out, nat32)
          aux1(im2,im3,:) = out
        END DO
    END DO
    !
    DO im1 = 1, nat3
        DO im3 = 1, nat3
          in = aux1(:, im3, im1)
          CALL ZGEMV('T', nat3, nat3, (1._dp,0._dp), zz2, nat3, in, &
                          1, (0._dp,0._dp), out, 1)
          aux2(:, im3, im1) = out
        END DO
    END DO
    !
    DO im1 = 1, nat3
        DO im2 = 1, nat3
          in = aux2(im2, 1, im1)
          CALL ZGEMV('T', nat3, nat3, (1._dp,0._dp), zz3, nat3, in, &
                          nat3, (0._dp,0._dp), out, nat3)
          aux1(im2, :, im1) = out
        END DO
    END DO
    !
    d3mm = REAL( CONJG(aux1)*aux1 , kind=DP)

    DEALLOCATE(aux1,aux2,out,in)
    !
  END SUBROUTINE calc_d3mm
  ! \/o\________\\\_________________________________________/^>
  ! Returns the cross section of a 3-phonon scattering ;
  ! this is | V^3(q1,q2,q3) |^2 in Fugallo et. al. PRB
  SUBROUTINE scatter_3q(S, fc2, fc3, xq1, xq2, xq3, V3sq)
    USE d3_basis,  ONLY : d3_cart2pat
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in)  :: xq1(3), xq2(3), xq3(3)
    REAL(DP),INTENT(out) :: V3sq(:,:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:)
    REAL(DP) :: xq(3,3)
    INTEGER :: i
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(w2(S%nat3))
    !
    xq(:,1) = xq1 ; xq(:,2) = xq2; xq(:,3) = xq3
    DO i = 1,3
      CALL fftinterp_mat2(xq(:,i), S, fc2, U(:,:,i))
      CALL mat2_diag(S, U(:,:,i), w2)
    ENDDO
    !
    DEALLOCATE(w2)
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
    !
    CALL d3_cart2pat(D3, S%nat, U(:,:,1), U(:,:,2), U(:,:,3))
    V3sq = REAL( CONJG(D3)*D3 , kind=DP)
    !
    DEALLOCATE(U, D3)
    !
  END SUBROUTINE scatter_3q

END MODULE interp_fc
! <<^V^\\=========================================//-//-//========//O\\//



