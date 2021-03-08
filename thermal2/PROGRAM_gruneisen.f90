!
! Written by Lorenzo Paulatto (2020) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE add_gruneisen_program
#include "mpi_thermal.h"
  USE kinds, ONLY : DP

  CONTAINS 


  FUNCTION gruneisen_q(xq0, nconf, T, sigma, S, grid, fc2, fc3)
    USE q_grids,          ONLY : q_grid
    USE input_fc,         ONLY : ph_system_info, div_mass_dyn
    USE fc2_interpolate,  ONLY : forceconst2_grid, freq_phq_safe, bose_phq,&
                                 set_nu0, ip_cart2pat, dyn_cart2pat
    USE fc3_interpolate,  ONLY : forceconst3
    USE constants,        ONLY : RY_TO_CMM1
    USE mpi_thermal,      ONLY : mpi_bsum
    USE ph_velocity,      ONLY : velocity
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    ! FUNCTION RESULT:
    REAL(DP) :: gruneisen_q(S%nat3,nconf)
    REAL(DP) :: grun(S%nat3,nconf)
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    INTEGER :: iq, jq, nu, it
    INTEGER :: nu0(3)
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3), vel(3,S%nat3)
    LOGICAL :: only_diag, scale_D3
    REAL(DP) :: scale_D3_fac, freqm1(S%nat3)
    INTEGER :: mu
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    grun = (0._dp, 0._dp)
    !
    xq(:,1) = -xq0
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
      vel = velocity(S, fc2, xq(:,2))
      !
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      !
      DO it = 1,nconf
        DO jq = 1,3
          CALL bose_phq(T(it),s%nat3, freq(:,jq), bose(:,jq))
        ENDDO
        !
      ENDDO
      !
      DO nu = 1,S%nat3
        
      ENDDO
      !
    ENDDO
    !
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf, grun)
    !
    gruneisen_q = -0.5_dp * grun
    !
    DEALLOCATE(U)
    !
  END FUNCTION gruneisen_q  

  SUBROUTINE GRUNEISEN_PATH(input,qpath, S,fc2,fc2b,fc3)
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info, multiply_mass_dyn, &
                                 write_dyn, read_fc2, aux_system, div_mass_fc2
    USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat, set_nu0, freq_phq_safe
    USE q_grids,          ONLY : q_grid, setup_grid
    USE fc3_interpolate,  ONLY : forceconst3
    USE code_input,       ONLY : READ_INPUT, code_input_type
    USE constants,        ONLY : RY_TO_CMM1
    USE more_constants,   ONLY : write_conf
    IMPLICIT NONE
    TYPE(forceconst2_grid) :: fc2, fc2b
    CLASS(forceconst3),INTENT(in) :: fc3
    TYPE(ph_system_info)   :: S
    TYPE(code_input_type)     :: input
    TYPE(q_grid)      :: qpath
    
    TYPE(q_grid) :: grid

    INTEGER :: iq, it
    COMPLEX(DP),ALLOCATABLE :: D(:,:)
    REAL(DP),ALLOCATABLE    :: w(:), grun(:,:)
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    CHARACTER(len=512) :: filename

    ALLOCATE(grun(S%nat3,input%nconf))
    ALLOCATE(D(S%nat3,S%nat3))
    ALLOCATE(w(S%nat3))

    ioWRITE(*,*) 'PROGRAM GRUNEISEN'

    DO it = 1,input%nconf
      filename=TRIM(input%outdir)//"/"//&
              TRIM(input%prefix)//"_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out"
      OPEN(unit=1000+it, file=filename)
      ioWRITE(1000+it, *) "# calculation of gruneisen (gamma_n)"
      ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
    ENDDO

    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
    CALL grid%scatter()
 
    PATH : & 
    DO iq = 1, qpath%nq
      ioWRITE(stdout, *) "<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>"
      ioWRITE(stdout, '(i6,3f12.6)') iq, qpath%xq(:,iq)
      CALL freq_phq_safe(qpath%xq(:,iq), S,fc2,w,D)
      ioWRITE(stdout,"(6f20.12)") w*RY_TO_CMM1
      ioWRITE(stdout,*)
      IF(ionode) FLUSH(stdout)

      grun = gruneisen_q(qpath%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, S, grid, fc2, fc3)
      !
      DO it = 1, input%nconf
       ioWRITE(1000+it, *) iq, w, grun(:,it) 
       !       
      ENDDO
    ENDDO &
    PATH ! <--- path end

    DO it = 1,input%nconf
      CLOSE(1000+it)
    ENDDO
  END SUBROUTINE GRUNEISEN_PATH

END MODULE add_gruneisen_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM add_gruneisen

  USE kinds,            ONLY : DP
  USE add_gruneisen_program
  USE fc3_interpolate,  ONLY : forceconst3
  USE input_fc,         ONLY : forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE more_constants,   ONLY : print_citations_linewidth

  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2, fc2b
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S, Sb
  TYPE(code_input_type)     :: input
  TYPE(q_grid)      :: qpts

  CALL start_mpi()
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("GRUN", input, qpts, S, fc2, fc3)
  !
  IF(input%calculation == "grun")THEN
    CALL GRUNEISEN_PATH(input,qpts, S,fc2,fc2b,fc3)
  ELSE
    CALL errore("gruneisen", "no such calculation: "//TRIM(input%calculation),1)
  ENDIF
  !
  CALL stop_mpi()
  CALL print_citations_linewidth()
 
END PROGRAM add_gruneisen
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

