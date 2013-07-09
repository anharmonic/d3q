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
  SUBROUTINE LW_TEST(xq, S, fc2, fc3)
    USE interp_fc, ONLY : fftinterp_mat2, fftinterp_mat3, mat2_diag
    USE constants, ONLY : RY_TO_CMM1
    USE nanoclock, ONLY : f_nanosec
    USE d3_basis,  ONLY : d3_cart2pat
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(forceconst3_grid),INTENT(in) :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: xq(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: w2(:), d3mm(:,:,:)
    REAL(DP) :: t0,t1
    INTEGER :: i
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(w2(S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(d3mm(S%nat3, S%nat3, S%nat3))
    !
    DO i = 1,3
      CALL fftinterp_mat2(xq(:,i), S, fc2, U(:,:,i))
      CALL mat2_diag(S, U(:,:,i), w2)
    ENDDO
    
    CALL fftinterp_mat3(xq(:,2), xq(:,3), S, fc3, D3)

    CALL d3_cart2pat(D3, S%nat, U(:,:,1), U(:,:,2), U(:,:,3))
    d3mm = REAL( CONJG(D3)*D3 , kind=DP)
    !
    WRITE(*, '(3f12.6)') xq(:,:)
    DO i = 1,s%nat3
      WRITE(*, '(9f12.6)') d3mm(:,:,i)*1.d+15
      WRITE(*,*)
    ENDDO
    !
    DEALLOCATE(U, w2, D3, d3mm)
    !
  END SUBROUTINE LW_TEST  
  
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
  
  CALL QBZ_LINE((/0.5_dp,0.288675_dp,0._dp/), (/0.0_dp,0._dp,0._dp/), &
                101, S, fc2)


  xq(:,1) = (/0.5_dp,0.288675_dp,0._dp/)
  xq(:,2) = - xq(:,1)
  xq(:,3) = - (xq(:,1)+xq(:,2))
  CALL LW_TEST(xq, S, fc2, fc3)

  CALL environment_end('LW')
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!













