!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
MODULE more_constants
  USE kinds, ONLY : DP
  REAL(DP),PARAMETER :: RY_TO_JOULE =  0.5* 4.35974394e-18
  REAL(DP),PARAMETER :: RY_TO_SECOND = 2* 2.418884326505e-17
  REAL(DP),PARAMETER :: RY_TO_METER = 5.2917721092e-11
  REAL(DP),PARAMETER :: RY_TO_WATTMM1KM1 = RY_TO_JOULE / (RY_TO_SECOND * RY_TO_METER)
END MODULE more_constants

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
    USE iso_c_binding,  ONLY : c_int
    USE input_fc,       ONLY : same_system, read_fc2, read_fc3, &
                               aux_system, div_mass_fc2, div_mass_fc3
    USE io_global,	ONLY : stdout
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
!     CALL write_fc3("mat3bis", S, fc3)
    !
    CALL aux_system(s)
    !
    CALL memstat(kb)
    WRITE(stdout,*) "Reading : done."
    WRITE(stdout,*) "Memory used : ", kb/1000, "Mb"
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
      WRITE(766, '(i4,f12.6,4x,3f12.6,5x,9f12.6)') i,pl,xq, SQRT(w2)*RY_TO_CMM1 

      xvels = velocity_simple(S,fc, xq)
      WRITE(767, '(i4,f12.6,4x,3f12.6,5x,9(3e14.3,2x))') i,pl,xq, xvels
      
      xvelp = velocity_proj(S,fc, xq)
      WRITE(768, '(i4,f12.6,4x,3f12.6,5x,9(3e14.3,2x))') i,pl,xq, xvelp

      pl = pl + dpl
    ENDDO
    !
    DEALLOCATE(D, w2, xvels, xvelp)
    !
  END SUBROUTINE QBZ_LINE
  
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE LW_QBZ_LINE(xq0, xq1, nq, S, fc2, fc3)
    USE interp_fc,   ONLY : fftinterp_mat2, mat2_diag, linewidth_q
    USE constants,   ONLY : RY_TO_CMM1
    USE ph_velocity 
    USE q_grid,      ONLY : q_grid_type, setup_simple_grid
    USE nanoclock
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
    TYPE(q_grid_type) :: grid
    !
    ALLOCATE(D(S%nat3, S%nat3))
    ALLOCATE(w2(S%nat3), lw(S%nat3))
    !
    dxq = 1._dp / REAL(nq-1,kind=DP)
    pl = 0._dp
    dpl = SQRT(SUM( (dxq*(xq1-xq0))**2 ))
    !
    CALL setup_simple_grid(S, 16,16,1, grid)
    !
    CALL print_header()
    !
    DO i = 0,nq-1
      xq = (xq0*(nq-1-i) + xq1*i)*dxq
      CALL fftinterp_mat2(xq, S, fc2, D)
      CALL mat2_diag(S, D, w2)
      
!       lw = LW_TEST2(xq, 300._dp, S, fc2, fc3)
      lw = linewidth_q(xq, 300._dp, 10._dp/RY_TO_CMM1, S, grid, fc2, fc3)
      WRITE(666, '(i4,f12.6,2x,3f12.6,2x,9f12.6,2x,9e15.4)') &
                   i,pl,xq, SQRT(w2)*RY_TO_CMM1, lw*RY_TO_CMM1

      pl = pl + dpl
    ENDDO
    !
    DEALLOCATE(D, w2, lw)
    !
  END SUBROUTINE LW_QBZ_LINE
  
  !  
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE SMA_TRANSPORT(S, fc2, fc3, n1,n2,n3, T, sigma)
    USE interp_fc,      ONLY : fftinterp_mat2, mat2_diag, linewidth_q
    USE constants,      ONLY : RY_TO_CMM1, K_BOLTZMANN_RY
    USE more_constants, ONLY : RY_TO_WATTMM1KM1
    USE ph_velocity ,   ONLY : velocity_proj
    USE q_grid,         ONLY : q_grid_type, setup_simple_grid
    USE functions,      ONLY : f_bose
    USE mp_world,       ONLY : mpime, nproc, world_comm
    USE mp,             ONLY : mp_sum
    USE io_global,	ONLY : stdout
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    INTEGER,INTENT(in) :: n1,n2,n3
    REAL(DP),INTENT(in) :: T,sigma
    !
    COMPLEX(DP),ALLOCATABLE :: D(:,:)
    REAL(DP),ALLOCATABLE    :: w2(:), lw(:), xvel(:,:)
    REAL(DP) :: k(3,3,0:S%nat3), bose(S%nat3)
    INTEGER :: i, a,b, nu
    TYPE(q_grid_type) :: grid
    REAL(DP) :: dk, k0
    !
    ALLOCATE(D(S%nat3, S%nat3))
    ALLOCATE(w2(S%nat3), lw(S%nat3), xvel(3,S%nat3))
    !
    CALL setup_simple_grid(S, n1,n2,n3, grid)
    WRITE(stdout,*) "INTEGRATING:", grid%nq
    !
    k0 = 1/(S%omega*K_BOLTZMANN_RY*T**2)
    k = 0._dp
    !
    DO i = 1+mpime, grid%nq, nproc
      !
      WRITE(stdout,*) i
      CALL fftinterp_mat2(grid%xq(:,i), S, fc2, D)
      CALL mat2_diag(S, D, w2)
      xvel = velocity_proj(S,fc2, grid%xq(:,i))
      !
      lw = linewidth_q(grid%xq(:,i), T, sigma, S, grid, fc2, fc3)
      WRITE(1000+mpime, '(i4,3f12.6,2x,9f12.6,2x,9e15.4)') &
                   i,grid%xq(:,i), SQRT(w2)*RY_TO_CMM1, lw*RY_TO_CMM1      !
      IF(i<grid%nq)THEN
        IF(grid%xq(1,i)/=grid%xq(1,i+1)) WRITE(1000+mpime,*) i, "x"
      ENDIF
      
      bose = f_bose(SQRT(w2), T)
      
      DO nu = 1,S%nat3
        IF(lw(nu)/=0._dp)THEN
        DO a = 1,3
        DO b = 1,3
          dk =  k0*xvel(a,nu)*xvel(b,nu)*w2(nu)*bose(nu)*(bose(nu)+1)/lw(nu)/grid%nq
          k(a,b,nu) = k(a,b,nu) +dk
          k(a,b,0)  = k(a,b,0)  +dk
        ENDDO
        ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    DEALLOCATE(D, w2, lw)
    !
    CALL mp_sum(k, world_comm)
    !
    DO nu = 0,S%nat3
      write(667,*) "--------------", nu
      DO b = 1,3
        write(667, *) RY_TO_WATTMM1KM1*k(:,b,nu)
      ENDDO
    ENDDO
    !
  END SUBROUTINE SMA_TRANSPORT  
 
END MODULE linewidth_program


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE linewidth_program
  USE environment,      ONLY : environment_start, environment_end
  USE constants,        ONLY : RY_TO_CMM1
  USE mp_world,         ONLY : mp_global_start, mp_global_end
  USE input_fc
  
  TYPE(forceconst2_grid) :: fc2
  TYPE(forceconst3_grid) :: fc3
  TYPE(ph_system_info)   :: s
  
  REAL(DP) :: xq(3,3)
  
  CALL mp_global_start()
  CALL environment_start('LW')
  
  CALL INPUT(s, fc2, fc3)
 
  CALL QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/),&
                   200, S, fc2)
  
!   CALL LW_QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/),&
!                    200, S, fc2, fc3)

  CALL SMA_TRANSPORT(S, fc2, fc3, 3,3,1, 300._dp, 10._dp/RY_TO_CMM1)

  CALL environment_end('LW')
  CALL mp_global_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













