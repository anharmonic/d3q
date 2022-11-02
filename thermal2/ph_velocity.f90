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
  PUBLIC :: velocity, velocity_operator
  !_simple, velocity_fdiff, velocity_ft
         
  CONTAINS
  FUNCTION velocity(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, fftinterp_dmat2, mat2_diag
    USE merge_degenerate, ONLY : merge_degenerate_velocity
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP) :: xvel_operator(3,S%nat3,S%nat3)
    REAL(DP) :: velocity(3,S%nat3)
    !
    INTEGER :: ix, im1
    ! 
    xvel_operator = velocity_operator_fdiff(S,fc, xq) 
    ! group velocities are diagonal elements of velocity operator
    DO ix=1,3
      DO im1=1,S%nat3
        velocity(ix,im1)=REAL(xvel_operator(ix,im1,im1))
      ENDDO
    END DO
    !! If we have effective charges, use finite differences derivation
    !IF(S%lrigid)THEN
    !  velocity = velocity_fdiff(S,fc, xq) 
    !ELSE
    !! Otherwise, use the property of Fourier transform to get dD2/dq = \sum_R R e(iqR) F_R
    !  velocity = velocity_ft(S,fc, xq) 
    !ENDIF
  END FUNCTION 

  FUNCTION velocity_operator(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, fftinterp_dmat2, mat2_diag
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP) :: velocity_operator(3,S%nat3,S%nat3)
    INTEGER :: im1, dir
    !
    velocity_operator = velocity_operator_fdiff(S,fc, xq) 
    !
    ! the diagonal elements of the velocity operator are the group velocities
    !do dir=1,3
    !  do im1=1,S%nat3
    !      velocity(dir,im1) = REAL(velocity_operator(dir, im1,im1))
    !  enddo
    !enddo
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
  ! Compute the velocity operator derived in 
  ! Simoncelli, M., Marzari, N., & Mauri, F. Phys. Rev. X (2022). 
  FUNCTION velocity_operator_fdiff(S,fc, xq)
    USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
    USE constants, ONLY : tpi,RY_TO_CMM1
    USE more_constants, ONLY : complex_i
    USE kinds,    ONLY : DP

    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3)
    ! FUNCTION RESULT:
    complex(DP) :: velocity_operator_fdiff(3,S%nat3,S%nat3) 
    !
    COMPLEX(DP),ALLOCATABLE :: eig_Z(:,:),eig_plus_Z(:,:),eig_minus_Z(:,:)  ! eigenvectors in Ziman's phase
    COMPLEX(DP),ALLOCATABLE :: eig_W(:,:), eig_plus_W(:,:),eig_minus_W(:,:) ! eigenvectors in Wallace's phase
    COMPLEX(DP),ALLOCATABLE :: matrix_U(:),matrix_U_plus(:), matrix_U_minus(:)

    COMPLEX(DP),ALLOCATABLE :: xvel_operator(:,:,:), velph_matrix_rotated(:,:)
    COMPLEX(DP),ALLOCATABLE :: sqrt_D_plus_W(:,:), sqrt_D_minus_W(:,:)
    COMPLEX(DP),ALLOCATABLE :: diffM_W(:,:), diffM_W_sym(:,:),velph_D(:,:)

    REAL(DP),ALLOCATABLE    :: w2p(:), w2m(:), w2(:),freqcp(:),freqcm(:), sqrt_w2(:)
    COMPLEX(DP) :: average_complx
    LOGICAL :: gamma, flag_cycle, diagonalize_V_in_degenerate_subspace
    !
    INTEGER :: ix, nu, ip, jp, sp, im1, im2,i_alpha,ia, i_atom, lbl_deg, dim_dg, id_m1,id_m2, s2, id_e
    REAL(DP) :: xqp(3), xqm(3), dh, cq(3), full_infos_freq(S%nat3,3), arg
    !
    REAL(DP) :: threshold_f_degeneracy_cmm1, dim_deg_sub(S%nat3), Egval_D(S%nat3)
    REAL(DP),PARAMETER :: eps12 = 1.d-12
    !
    ALLOCATE(eig_Z(S%nat3,S%nat3),eig_plus_Z(S%nat3,S%nat3),eig_minus_Z(S%nat3,S%nat3),          &
             eig_W(S%nat3,S%nat3),eig_plus_W(S%nat3,S%nat3),eig_minus_W(S%nat3,S%nat3),          &
             matrix_U(S%nat3),matrix_U_plus(S%nat3), matrix_U_minus(S%nat3),                     &
             xvel_operator(3,S%nat3,S%nat3), velph_matrix_rotated(S%nat3,S%nat3),           &
             sqrt_D_plus_W(S%nat3,S%nat3),sqrt_D_minus_W(S%nat3,S%nat3),                         &
             diffM_W(S%nat3,S%nat3), diffM_W_sym(S%nat3,S%nat3),                                 &
             w2p(S%nat3),w2m(S%nat3),w2(S%nat3),sqrt_w2(S%nat3),freqcp(S%nat3),freqcm(S%nat3),  &
             velph_D(S%nat3,S%nat3))
    !
    threshold_f_degeneracy_cmm1=0.2
    diagonalize_V_in_degenerate_subspace=.true.
    ! compute frequencies and eigenvectors at point q
    CALL fftinterp_mat2(xq, S, fc, eig_Z)
    CALL mat2_diag(S%nat3, eig_Z, w2)
    !
    ! enforce acoustic modes to have zero frequency
    cq = xq
    CALL cryst_to_cart(1,cq,S%at,-1)  ! -1 means convert from Cartesian to Crystal
    gamma = ALL( ABS(cq-NINT(cq))<(h/10.0))
    IF( gamma )THEN
      w2(1:3) = 0._dp
      ! "bar" means eigevectors in the Ziman phase
      eig_Z(:,1:3) = (0._dp, 0._dp)
    ENDIF
    !
    sqrt_w2=SQRT(w2)
    xvel_operator(:,:,:)=(0.0d0,0.0d0)
    !
    ! DETECTION OF DEGENERACIES
    !use real to detect if the sqrt is not an integer 
    ! full_infos_freq(:,1) contains the frequencies at the q_point at which we are computing the velocity (i_q)
    ! full_infos_freq(:,2) contains a label to identify different degeneracies
    ! full_infos_freq(:,3) contains info on the degeneracy. If it contains a non-zero value at index "im",
    ! the frequency "im" is degenerate. Different non-zero values identifies the degeneracy subspace.
    ! 
    ! EXAMPLES:
    ! if  full_infos_freq(im0,2) = im0 
    ! and full_infos_freq(im0,3) = 0
    !then the frequency "im0" is not degenerate
    !
    ! if full_infos_freq(im1,2) = im1 and 
    !    full_infos_freq(im1,3) = 1
    !    full_infos_freq(im2,2) = im1 and 
    !    full_infos_freq(im2,3) = 1
    ! then the frequency "im1" is degenerate with the freq "im2"
    !The second index will be changed, all the degenerate frequencies belonging to
    !the same subspace will have the same second index
    !initialize the array containint infos on the degeneracies
    do im1=1,S%nat3
      full_infos_freq(im1,1)=sqrt_w2(im1)
      full_infos_freq(im1,2)=im1
      full_infos_freq(im1,3)=0.0d0
    end do
    !
    !
    im1=1
    im2=1
    do while ((im1<=S%nat3) .and. (im2<S%nat3))
      im2=im1+1
      flag_cycle=.true.
      do while (im2<=S%nat3 .and. flag_cycle)
        if (dabs(full_infos_freq(im1,1)-sqrt_w2(im2))<=(threshold_f_degeneracy_cmm1/RY_TO_CMM1)) then
          !thirs index set to 1 means that the freq is degenerate with some other freq.
          !This is needed to remove ambiguity in the case in which for example
          !freq(1) is degenerate with freq(2) and 1 labels also the degeneracy.
          full_infos_freq(im1,3)=1.0d0
          full_infos_freq(im2,3)=1.0d0
          ! ensure that you use the smallest index to label the degenerate subspace
          !if (full_infos_freq(im2,2)==im2) then  
          ! in this way we memorize that the frequency im2 is degenerate with the freq im1
          full_infos_freq(im2,2)=im1
          im2=im2+1
        else
          flag_cycle=.false.
        end if
      end do
      im1=im2
    end do
    !
    !
    matrix_U(:)=(0.0d0,0.0d0)
    !
    ! these are the eigevectors in the Wallace phase
    eig_W(:,:)=(0.0d0,0.0d0)
    eig_minus_W(:,:)=(0.0d0,0.0d0)
    eig_plus_W(:,:)=(0.0d0,0.0d0)
    !
    ! U = exp (-i*q*[2pi/alat] * tau*[alat] ) 
    do i_atom=1,S%nat
      do i_alpha=1,3
        ip = i_alpha + (i_atom-1)*3
        ! xq = q-point in Cartesian Coordinates [2pi/alat] units
        ! S%tau = atomic position in primitive cell, Cartesian Coordinates [alat] units
        arg = tpi * SUM(xq(:)*S%tau(:,i_atom)) 
        matrix_U(ip)=cmplx(cos(arg),-sin(arg),kind=DP)
      enddo
    enddo
    !
    ! compute eigevectors in Wallace phase
    DO sp=1,S%nat3
      DO ip = 1,S%nat3
        eig_W(ip,sp)=matrix_U(ip)*eig_Z(ip,sp)
      ENDDO
    ENDDO    
    !

    !
    DO ix = 1,3
      sqrt_D_plus_W = 0._dp
      sqrt_D_minus_W = 0._dp
      diffM_W =0._dp
      diffM_W_sym =0._dp
      xqp = xq  ! qpoint in cartesian Coordinates [2pi/alat] units
      xqp(ix) = xq(ix)+h
      CALL fftinterp_mat2(xqp, S, fc, eig_plus_Z)
      CALL mat2_diag(S%nat3, eig_plus_Z, w2p)
      !
      ! enforce acoustic modes to have zero frequency
      cq = xq
      CALL cryst_to_cart(1,cq,S%at,-1)
      gamma = ALL( ABS(cq-NINT(cq))<(h/10.0))
      IF( gamma )THEN
        w2p(1:3) = 0._dp
        eig_plus_Z(:,1:3) = (0._dp, 0._dp)
      ENDIF
      !
      freqcp = SQRT(w2p)
      !
      matrix_U_plus(:)=(0.0d0,0.0d0)
      do i_atom=1,S%nat
        do i_alpha=1,3
          ip = i_alpha + (i_atom-1)*3
          ! xqp = q-point in Cartesian Coordinates [2pi/alat] units
          ! S%tau = atomic position in primitive cell, Cartesian Coordinates [alat] units
          arg = tpi * SUM(xqp(:)*S%tau(:,i_atom))
          matrix_U_plus(ip)=cmplx(cos(arg),-sin(arg),kind=DP)
        enddo
      enddo
      DO sp=1,S%nat3
        DO ip = 1,S%nat3
          eig_plus_W(ip,sp)=matrix_U_plus(ip)*eig_plus_Z(ip,sp)
        ENDDO
      ENDDO
      !
      DO ip=1,S%nat3
        DO jp=1,S%nat3
          DO sp=1,S%nat3
            sqrt_D_plus_W(ip,jp) =sqrt_D_plus_W(ip,jp)+eig_plus_W(ip,sp)*dcmplx(freqcp(sp),0.0d0)*DCONJG(eig_plus_W(jp,sp))
          ENDDO
        ENDDO
      ENDDO
      !
      !
      xqm = xq
      xqm(ix) = xq(ix)-h
      CALL fftinterp_mat2(xqm, S, fc, eig_minus_Z)
      CALL mat2_diag(S%nat3, eig_minus_Z, w2m)
      !
      ! enforce acoustic modes to have zero frequency
      cq = xq
      CALL cryst_to_cart(1,cq,S%at,-1)
      gamma = ALL( ABS(cq-NINT(cq))<(h/10.0))
      IF( gamma )THEN
        w2m(1:3) = 0._dp
        eig_minus_Z(:,1:3) = (0._dp, 0._dp)
      ENDIF
      !
      freqcm = SQRT(w2m)
      !
      matrix_U_minus(:)=(0.0d0,0.0d0)
      do i_atom=1,S%nat
        do i_alpha=1,3
          ip = i_alpha + (i_atom-1)*3
          arg = tpi * SUM(xqm(:)*S%tau(:,i_atom))
          matrix_U_minus(ip)=cmplx(cos(arg),-sin(arg),kind=DP)
        enddo
      enddo
      DO sp=1,S%nat3
        DO ip = 1,S%nat3
          eig_minus_W(ip,sp)=matrix_U_minus(ip)*eig_minus_Z(ip,sp)
        ENDDO
      ENDDO      
      DO ip=1,S%nat3
        DO jp=1,S%nat3
          DO sp=1,S%nat3
            sqrt_D_minus_W(ip,jp) =sqrt_D_minus_W(ip,jp)+eig_minus_W(ip,sp)*dcmplx(freqcm(sp),0.0d0)*DCONJG(eig_minus_W(jp,sp))
          ENDDO
        ENDDO
      ENDDO
      !
      ! compute the derivative of the sqrt(dynamical matrix) in the Cartesian basis
      dh = ( xqp(ix)-xqm(ix) ) *S%tpiba ! S%tpiba = tpi/S%alat, where alat is in Bohr
      !
      DO ip=1,S%nat3
        DO jp=1,S%nat3
          diffM_W(ip,jp) =(sqrt_D_plus_W(ip,jp) - sqrt_D_minus_W(ip,jp)) /dcmplx(dh,0.0)
        ENDDO
      ENDDO
      ! enforce the matrix to be hermitian (in practice, it is already Hermitian)
      DO ip=1,S%nat3
        DO jp=1,S%nat3
          diffM_W_sym(ip,jp) = 0.5d0*(diffM_W(ip,jp)+DCONJG(diffM_W(jp,ip)))
        ENDDO
      ENDDO
      !
      ! now go back in the eigenstate basis
      !
      DO im1=1,S%nat3
        DO im2=1,S%nat3
          DO ip=1,S%nat3
            DO jp=1,S%nat3
              xvel_operator(ix,im1,im2)=xvel_operator(ix,im1,im2)+           &
                                        (DCONJG(eig_W(ip,im1)))*diffM_W_sym(ip,jp)*eig_W(jp,im2) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      !
    ENDDO
    !
    !_________________________________________________________________________________________________ 
    ! diagonalization of the velocity operator in the degenerate subspaces
    !
    if (diagonalize_V_in_degenerate_subspace) then
      do ix=1,3
        dim_deg_sub=0.0d0
        do sp =1, S%nat3
          !
          if ( (nint(full_infos_freq(sp,3))>0) .and.( nint(full_infos_freq(sp,2))==sp )) then
            !--------------------------------------------------------------
            ! Degeneracy, diagonalize proj. of V_ix in degenerate subspace 
            !--------------------------------------------------------------
            !this ensures that 's' is the label of the degeneracy
            lbl_deg=nint(full_infos_freq(sp,2))
            velph_D=(0.0d0,0.0d0)
            do im1 = lbl_deg,S%nat3
              do im2 = lbl_deg,S%nat3
                if (((nint(full_infos_freq(im1,2))==lbl_deg).and.(nint(full_infos_freq(im1,3))>0)) .and. & 
                    ((nint(full_infos_freq(im2,2))==lbl_deg).and.(nint(full_infos_freq(im2,3))>0)) ) then
                     dim_deg_sub(sp)=dim_deg_sub(sp)+1
                    !in this case we have the frequency im and im2 are degenerate
                    !and we are selecting all the eigenvectors which spans the 
                    !degenerate subspace
                    velph_D(im1, im2)=xvel_operator(ix,im1,im2)
                end if
              end do
            end do
            dim_deg_sub(sp)=sqrt(dim_deg_sub(sp))
            dim_dg=nint(dim_deg_sub(sp))
            !enforce to be hermitian
            do im1 =lbl_deg,lbl_deg+ dim_dg-1
              do im2 =lbl_deg,lbl_deg+ dim_dg-1
                average_complx=0.5d0*(velph_D(im1,im2)+dconjg(velph_D(im2,im1)))
                if (im1==im2) then
                  velph_D(im1,im2)=real( average_complx )
                else
                  velph_D(im1,im2)=average_complx
                  velph_D(im2,im1)=dconjg(average_complx)
                end if
              end do
            end do
            !
            call mat2_diag(dim_dg,velph_D(lbl_deg:lbl_deg+dim_dg-1,lbl_deg:lbl_deg+dim_dg-1),Egval_D(lbl_deg:lbl_deg+dim_dg-1))
            !
            do im1=1,S%nat3
              if ((im1<lbl_deg).or.(im1>=lbl_deg+dim_dg)) then
                ! identity outside the degenerate subspace
                velph_D(im1,im1)=1.0d0
              end if
            end do
            !
            !write the velocity operator with respect to the basis that diagonalizes V in degenerate subspace
            !
            velph_matrix_rotated=(0.0d0,0.0d0)
            do id_m1=lbl_deg,lbl_deg+dim_dg-1
              do id_m2=lbl_deg,lbl_deg+dim_dg-1
                do id_e = lbl_deg,lbl_deg+dim_dg-1
                  do s2 = lbl_deg,lbl_deg+dim_dg-1
                    velph_matrix_rotated(id_m1,id_m2)= velph_matrix_rotated(id_m1,id_m2)+ &
                    dconjg(velph_D(id_e,id_m1))*xvel_operator(ix,id_e,s2)*velph_D(s2,id_m2)
                  enddo
                enddo
              enddo
            enddo
            do id_m1 = lbl_deg,lbl_deg+dim_dg-1
              do id_m2 = lbl_deg,lbl_deg+dim_dg-1
                xvel_operator(ix,id_m1,id_m2)= velph_matrix_rotated(id_m1,id_m2)
              enddo
            enddo
            !         
          elseif ( (nint(full_infos_freq(sp,3)).eq.0) .and.( nint(full_infos_freq(sp,2))==sp )) then
            !--------------------------------------------------------------
            ! NO degeneracy, do nothing
            !--------------------------------------------------------------
            dim_deg_sub(sp)=dim_deg_sub(sp)+1
            Egval_D(sp)=(real(xvel_operator(ix,sp,sp)))
            !
          end if
          !
        end do
      end do
      !
      if ((sum(dim_deg_sub)-S%nat3>eps12 )) then 
        print*,'ERROR: sum of the dim of subsp. xq=',xq, sum(dim_deg_sub), S%nat3
      end if
    end if 
    !_________________________________________________________________________________________________ 
    velocity_operator_fdiff(:,:,:)=xvel_operator(:,:,:)

    DEALLOCATE(eig_Z,eig_plus_Z,eig_minus_Z,                        &
             eig_W,eig_plus_W,eig_minus_W,                          &
             matrix_U,matrix_U_plus, matrix_U_minus,                &
             xvel_operator,velph_matrix_rotated ,                   &
             sqrt_D_plus_W,sqrt_D_minus_W,                          &
             diffM_W, diffM_W_sym,w2p,w2m,w2,sqrt_w2,freqcp,freqcm)
    !
  END FUNCTION velocity_operator_fdiff
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














