!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
MODULE linewidth_program
  !
  USE kinds,    ONLY : DP
  USE input_fc, ONLY : forceconst2_grid, &
                       forceconst3_grid, &
                       ph_system_info
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE INPUT(s, fc2, fc3)
    USE iso_c_binding, ONLY : c_int
    USE input_fc,      ONLY : same_system, read_fc2, read_fc3, &
                              aux_system, div_mass_fc2, div_mass_fc3
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(inout) :: fc2
    TYPE(forceconst3_grid),INTENT(inout) :: fc3
    TYPE(ph_system_info),INTENT(inout)   :: s
    TYPE(ph_system_info) :: s3
    !
    INTEGER(kind=c_int) :: kb
    !
    CALL read_fc2("mat2R", s,  fc2)
    CALL read_fc3("mat3R", s3, fc3)
    IF(.not.same_system(s, s3)) THEN
      !PRINT*, "WARNING! FC2 and FC3 systems DO NOT MATCH !!!"
      CALL errore("INPUT", "FC2 and FC3 crystals DO NOT MATCH !!!", 1)
    ENDIF
    !
    CALL aux_system(s)
    !
    CALL memstat(kb)
    PRINT*, "Reading : done."
    PRINT*, "Memory used : ", kb/1000, "Mb"
    !
    CALL div_mass_fc2(S, fc2)
    CALL div_mass_fc3(S, fc3)
    !
  END SUBROUTINE INPUT


  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE QBZ_LINE(xq0, xq1, nq, S, fc)
    USE interp_fc, ONLY : fftinterp_mat2, mat2_diag
    USE constants, ONLY : RY_TO_CMM1
    USE ph_velocity 
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq0(3),xq1(3) ! line from xq0 to xq1
    INTEGER,INTENT(in)  :: nq            ! number of points
    !
    COMPLEX(DP),ALLOCATABLE :: D(:,:)
    REAL(DP),ALLOCATABLE    :: w2(:), xvels(:,:), xvelp(:,:)
    REAL(DP) :: xq(3), dxq, pl,dpl
    INTEGER :: i
    !
    ALLOCATE(D(S%nat3, S%nat3))
    ALLOCATE(w2(S%nat3), xvels(3,S%nat3),  xvelp(3,S%nat3))
    !
    dxq = 1._dp / REAL(nq-1,kind=DP)
    pl = 0._dp
    dpl = SQRT(SUM( (dxq*(xq1-xq0))**2 ))
    !
    DO i = 0,nq-1
      xq = (xq0*(nq-i) + xq1*i)*dxq
    
      CALL fftinterp_mat2(xq, S, fc, D)
      CALL mat2_diag(S, D, w2)
      WRITE(666, '(i4,f12.6,4x,3f12.6,5x,9f12.6)') i,pl,xq, SQRT(w2)*RY_TO_CMM1 

      CALL velocity_simple(S,fc, xq, xvels)
      WRITE(667, '(i4,f12.6,4x,3f12.6,5x,9(3e14.3,2x))') i,pl,xq, xvels
      
      CALL velocity_proj(S,fc, xq, xvelp)
      WRITE(668, '(i4,f12.6,4x,3f12.6,5x,9(3e14.3,2x))') i,pl,xq, xvelp

      pl = pl + dpl
    ENDDO
    !
    DEALLOCATE(D, w2, xvels, xvelp)
    !
  END SUBROUTINE QBZ_LINE
  
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE LW_QBZ_LINE(xq0, xq1, nq, S, fc2, fc3)
    USE interp_fc, ONLY : fftinterp_mat2, mat2_diag
    USE constants,  ONLY : RY_TO_CMM1
    USE ph_velocity 
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq0(3),xq1(3) ! line from xq0 to xq1
    INTEGER,INTENT(in)  :: nq            ! number of points
    !
    COMPLEX(DP),ALLOCATABLE :: D(:,:)
    REAL(DP),ALLOCATABLE    :: w2(:), lw(:)
    REAL(DP) :: xq(3), dxq, pl,dpl
    INTEGER :: i
    !
    ALLOCATE(D(S%nat3, S%nat3))
    ALLOCATE(w2(S%nat3), lw(S%nat3))
    !
    dxq = 1._dp / REAL(nq-1,kind=DP)
    pl = 0._dp
    dpl = SQRT(SUM( (dxq*(xq1-xq0))**2 ))
    !
    DO i = 0,nq-1
      xq = (xq0*(nq-1-i) + xq1*i)*dxq
      CALL fftinterp_mat2(xq, S, fc2, D)
      CALL mat2_diag(S, D, w2)
      
      lw = LW_TEST2(xq, 300._dp, S, fc2, fc3)
      WRITE(666, '(i4,f12.6,2x,3f12.6,2x,9f12.6,2x,9e15.4)') &
                   i,pl,xq, SQRT(w2)*RY_TO_CMM1, lw*RY_TO_CMM1

      pl = pl + dpl
    ENDDO
    !
    DEALLOCATE(D, w2, lw)
    !
  END SUBROUTINE LW_QBZ_LINE
 
   FUNCTION LW_TEST2(xq0,T, S, fc2, fc3) &
   RESULT(lw)
    !
    USE nanoclock
    !
    USE interp_fc,      ONLY : fftinterp_mat2, fftinterp_mat3, &
                               mat2_diag, scatter_3q, sum_modes, &
                               ip_cart2pat
    USE q_grid,         ONLY : q_grid_type, setup_simple_grid
    USE functions,      ONLY : f_bose
    
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: T ! Kelvin
    REAL(DP),INTENT(in) :: xq0(3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:), V3sq(:,:,:)
    INTEGER :: iq, jq, nu
    !
    INTEGER,PARAMETER :: n1=16,n2=16,n3=1
    TYPE(q_grid_type) :: grid
    !
    REAL(DP) :: freq(S%nat3,3)
    REAL(DP) :: bose(S%nat3,3)
    REAL(DP) :: xq(3,3)
    !
    REAL(DP) :: lw(S%nat3)
    !

    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(w2(S%nat3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    CALL setup_simple_grid(S, n1,n2,n3, grid)
    !
    lw = 0._dp
    !
    !CALL velocity_proj(S,fc2, xq, xvel)
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL fftinterp_mat2(xq(:,1), S, fc2, U(:,:,1))
    CALL mat2_diag(S, U(:,:,1), freq(:,1))
    freq(:,1) = SQRT(freq(:,1))
    bose(:,1) = f_bose(freq(:,1), T)
    U(:,:,1) = (CONJG(U(:,:,1)))

    !
    CALL start_nanoclock(lwtot)
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      CALL start_nanoclock(i_ph)
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(grid%xq(:,iq)+xq(:,1))
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL fftinterp_mat2(xq(:,jq), S, fc2, U(:,:,jq))
        CALL mat2_diag(S, U(:,:,jq), freq(:,jq))
        freq(:,jq) = SQRT(freq(:,jq))
        bose(:,jq) = f_bose(freq(:,jq), T)
        U(:,:,jq) = (CONJG(U(:,:,jq)))
      ENDDO
!$OMP END PARALLEL DO
      CALL stop_nanoclock(i_ph)
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      CALL start_nanoclock(d3time)
      CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)
      CALL stop_nanoclock(d3time)
      !
      CALL start_nanoclock(c2pat)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      CALL stop_nanoclock(c2pat)
      
      CALL start_nanoclock(tv3sq)
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      CALL stop_nanoclock(tv3sq)
      ! ------ end of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
      !
      CALL start_nanoclock(tsum)
      lw = lw + sum_modes( S, freq, bose, V3sq )
      CALL stop_nanoclock(tsum)
      !
    ENDDO
    !
    lw = lw/grid%nq
    !
    CALL stop_nanoclock(lwtot)
    CALL print_line()
    CALL print_nanoclock(lwtot)
    CALL print_nanoclock(d3time)
    !
    CALL print_nanoclock(i_ph)
    CALL print_nanoclock(c2pat)
    CALL print_nanoclock(tv3sq)
    CALL print_nanoclock(tsum)
    CALL print_memory()
    !
    DEALLOCATE(U, w2, V3sq)
    !
  END FUNCTION LW_TEST2  
  
END MODULE linewidth_program


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds, ONLY : DP
  USE linewidth_program
  USE environment, ONLY : environment_start, environment_end
  
  TYPE(forceconst2_grid) :: fc2
  TYPE(forceconst3_grid) :: fc3
  TYPE(ph_system_info)   :: s
  
  REAL(DP) :: xq(3,3)
  
  CALL environment_start('LW')
  
  CALL INPUT(s, fc2, fc3)
  
!   CALL QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/), &
!                 101, S, fc2)

  !CALL setup_symmetry(S)

  CALL LW_QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/),&
                   200, S, fc2, fc3)

!   xq(:,1) = (/0.5_dp,0.288675_dp,0._dp/)
!   xq(:,2) = - xq(:,1)
!   xq(:,3) = - (xq(:,1)+xq(:,2))
!   !CALL LW_TEST(xq, S, fc2, fc3)
! 
!   CALL LW_TEST2(xq(:,1),300._dp, S, fc2, fc3)
!   CALL LW_TEST2(xq(:,2),300._dp, S, fc2, fc3)
!   CALL LW_TEST2(xq(:,3),300._dp, S, fc2, fc3)

  CALL environment_end('LW')
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













