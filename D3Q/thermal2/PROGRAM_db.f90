!
! Written by Lorenzo Paulatto (2014-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE add_bubble_program
#include "mpi_thermal.h"
  
END MODULE add_bubble_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM add_bubble

  USE kinds,            ONLY : DP
  USE add_bubble_program
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, &
                               ph_system_info, multiply_mass_dyn, write_dyn
  USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE q_grids,          ONLY : q_grid, setup_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE dynbubble,        ONLY : dynbubble_q
  USE constants,        ONLY : RY_TO_CMM1
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: input
  TYPE(q_grid)      :: qpts, grid
  !
  INTEGER :: iq, it, nu, mu
  COMPLEX(DP),ALLOCATABLE :: dyn(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: dyn0(:,:), U(:,:), dynX(:,:), dynY(:,:)
  REAL(DP),ALLOCATABLE    :: freq(:), freq0(:), freqY(:)
  CHARACTER (LEN=6),  EXTERNAL :: int_to_char
  CHARACTER(len=512) :: filename

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL start_mpi()
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("DB", input, qpts, S, fc2, fc3)
  CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
  CALL grid%scatter()
  !
  ALLOCATE(dyn(S%nat3,S%nat3,input%nconf))
  ALLOCATE(dyn0(S%nat3,S%nat3))
  ALLOCATE(dynX(S%nat3,S%nat3))
  ALLOCATE(dynY(S%nat3,S%nat3))
  ALLOCATE(U(S%nat3,S%nat3))
  ALLOCATE(freq(S%nat3), freq0(S%nat3), freqY(S%nat3))
  DO iq = 1, qpts%nq
    ioWRITE(*, *) "<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>"
    ioWRITE(*, '(i6,3f12.6)') iq, qpts%xq(:,iq)
    CALL fftinterp_mat2(qpts%xq(:,iq), S%nat3, fc2, dyn0)
    U = dyn0
    ioWRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, U)
    CALL mat2_diag(S%nat3, U, freq0)
    freq0 = SIGN(SQRT(ABS(freq0)), freq0)
    ioWRITE(*,"(6f20.12)") freq0*RY_TO_CMM1
    ioWRITE(*,*)

    ! Allocation is automatic in fortran 2003, of course it does not work with most compilers
    dyn = dynbubble_q(qpts%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, S, grid, fc2, fc3)
    
    DO it = 1, input%nconf
      ioWRITE(*, "(2f12.6)")  input%T(it),  input%sigma(it)
      dynX = dyn(:,:,it)
       
      !CALL dyn_cart2pat(dynX, S%nat3, U, +1)
      !DO mu = 1,S%nat3
      !DO nu = 1,S%nat3
      !  dynX(nu,mu) = dynX(nu,mu)*DSQRT(ABS(freq0(nu))*ABS(freq0(mu)))
      !ENDDO
      !ENDDO
      !CALL dyn_cart2pat(dynX, S%nat3, U, -1)

      ioWRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, dynX)
     !ioWRITE(*, "(6(2f12.6,4x))") RY_TO_CMM1*dynX
       
      !dynX = multiply_mass_dyn(S, dynX) + dyn0
      dynX = dynX + dyn0
      !
      dynY = multiply_mass_dyn(S, dynX)
      filename = "dyn_conf"//TRIM(int_to_char(it))//"_"//TRIM(int_to_char(iq))
      CALL write_dyn(filename, qpts%xq(:,iq), dynY, S)
      !
      dynY = dynX
      CALL mat2_diag(S%nat3, dynX, freq)
      freq = SIGN(SQRT(ABS(freq)), freq)
      
      CALL dyn_cart2pat(dynY, S%nat3, U, +1)
      FORALL(nu=1:S%nat3) freqY(nu) = dynY(nu,nu)
      freqY = SIGN(SQRT(ABS(freqY)), freqY)
             
      ioWRITE(*,"('freqX ',6f20.12)") freq*RY_TO_CMM1
      ioWRITE(*,"('freqY ',6f20.12)") freqY*RY_TO_CMM1
      ioWRITE(*,"('shiftX',6e20.6)") (freq-freq0)*RY_TO_CMM1
      ioWRITE(*,"('shiftY',6e20.6)") (freqY-freq0)*RY_TO_CMM1
      ioWRITE(*,*)
      !
      !
    ENDDO
  ENDDO

  CALL stop_mpi()
 
END PROGRAM add_bubble
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

