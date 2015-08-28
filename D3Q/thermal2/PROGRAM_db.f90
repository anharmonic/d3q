!
! Written by Lorenzo Paulatto (2014-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE add_bubble_program
  CONTAINS
  FUNCTION multiply_mass_dyn (S,dyn) RESULT(dyn_mass)
    USE kinds, only : DP
    USE input_fc,         ONLY : ph_system_info
    IMPLICIT NONE
    TYPE(ph_system_info)   :: S
    COMPLEX(DP),INTENT(in) :: dyn(S%nat3, S%nat3)
    COMPLEX(DP) :: dyn_mass(S%nat3, S%nat3)
    !
    INTEGER :: i, j
    !
    DO j = 1, S%nat3
    DO i = 1, S%nat3
      dyn_mass(i,j) = dyn(i,j)/(S%sqrtmm1(i)*S%sqrtmm1(j))
    ENDDO
    ENDDO
    !
  END FUNCTION multiply_mass_dyn
  
END MODULE add_bubble_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM add_bubble

  USE kinds,            ONLY : DP
  USE add_bubble_program
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, ph_system_info
  USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE q_grids,          ONLY : q_grid, setup_simple_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE dynbubble,        ONLY : dynbubble_q
  USE constants,        ONLY : RY_TO_CMM1
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: dbinput
  TYPE(q_grid)      :: qpts, grid
  !
  INTEGER :: iq, it, nu, mu
  COMPLEX(DP),ALLOCATABLE :: dyn(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dyn0(:,:), U(:,:), dynX(:,:)
  REAL(DP),ALLOCATABLE    :: freq(:), freq0(:)

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("DB", dbinput, qpts, S, fc2, fc3)
  CALL setup_simple_grid(S, dbinput%nk(1), dbinput%nk(2), dbinput%nk(3), grid)
  !
  ALLOCATE(dyn(S%nat3,S%nat3,dbinput%nconf))
  ALLOCATE(dyn0(S%nat3,S%nat3))
  ALLOCATE(dynX(S%nat3,S%nat3))
  ALLOCATE(U(S%nat3,S%nat3))
  ALLOCATE(freq(S%nat3), freq0(S%nat3))
  DO iq = 1, qpts%nq
    WRITE(*, *) "<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>"
    WRITE(*, '(i6,3f12.6)') iq, qpts%xq(:,iq)
    CALL fftinterp_mat2(qpts%xq(:,iq), S%nat3, fc2, dyn0)
    U = dyn0
    WRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, U)
    CALL mat2_diag(S%nat3, U, freq0)
    freq0 = SIGN(SQRT(ABS(freq0)), freq0)
    WRITE(*,"(6f20.12)") freq0*RY_TO_CMM1
    WRITE(*,*)

    ! Allocation is automatic in fortran 2003, of course it does not work with most compilers
    dyn = dynbubble_q(qpts%xq(:,iq), dbinput%nconf, dbinput%T, dbinput%sigma, S, grid, fc2, fc3)
    
    DO it = 1, dbinput%nconf
       WRITE(*, "(2f12.6)")  dbinput%T(it),  dbinput%sigma(it)
        dynX = dyn(:,:,it)
        
        !CALL dyn_cart2pat(dynX, S%nat3, U, +1)
!         DO mu = 1,S%nat3
!         DO nu = 1,S%nat3
!           dynX(nu,mu) = dynX(nu,mu)*DSQRT(ABS(freq0(nu))*ABS(freq0(mu)))
!         ENDDO
!         ENDDO

      !WRITE(*, "(6(2f12.6,4x))") RY_TO_CMM1*multiply_mass_dyn(S, dynX)
      WRITE(*, "(6(2f12.6,4x))") RY_TO_CMM1*dynX
        
        !dynX = multiply_mass_dyn(S, dynX) + dyn0
        dynX = dynX + dyn0
        CALL mat2_diag(S%nat3, dynX, freq)
        freq = SIGN(SQRT(ABS(freq)), freq)
        
        WRITE(*,"(6f20.12)") freq*RY_TO_CMM1
        WRITE(*,"(6e20.6)") (freq-freq0)*RY_TO_CMM1
        WRITE(*,*)
    ENDDO
  ENDDO
 
END PROGRAM add_bubble
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

