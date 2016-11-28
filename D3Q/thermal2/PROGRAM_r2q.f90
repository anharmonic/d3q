!
! Written by Lorenzo Paulatto (2014-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE r2q_program

 ! this space unintentionally left non-blank
#include "mpi_thermal.h"
  CONTAINS

  SUBROUTINE joint_dos(input, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, q_basis, setup_grid, prepare_q_basis
    USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    USE random_numbers,   ONLY : randy
    IMPLICIT NONE
    TYPE(code_input_type) :: input
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: qgrid, sgrid
    TYPE(q_basis) :: qbasis, sbasis
    !
    REAL(DP) :: freqj(S%nat3), freqk(S%nat3)
    REAL(DP) :: bosej(S%nat3), bosek(S%nat3)
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    !
    REAL(DP) :: xq0(3) = (/ 0._dp, 0._dp, 0._dp /)
    REAL(DP) :: xq_random(3)
    !
    REAL(DP) :: nrg(input%ne), jd_C(input%ne), jd_X(input%ne), xq_j(3), xq_k(3)
    REAL(DP) :: sigma_ry
    REAL(DP) :: dom_C(input%ne), dom_X(input%ne)
    REAL(DP) ::  ctm_C(input%ne), ctm_X(input%ne), bose_C, bose_X
    INTEGER :: jq, k,j,i
    !
    !
    FORALL(i=1:input%ne) nrg(i) = input%de * (i-1) + input%e0
    nrg = nrg/RY_TO_CMM1
    !
    sigma_ry = input%sigma(1)/RY_TO_CMM1
    
    jd_C = 0._dp
    jd_X = 0._dp

    xq_random  = (/ randy(), randy(), randy() /)
    CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
                qgrid, scatter=.false.)
    
    DO jq = 1, qgrid%nq
      xq_j = qgrid%xq(:,jq)
      CALL freq_phq(xq_j, S, fc, freqj, U)
      bosej(:) = f_bose(freqj, input%T(1))
      !WRITE(20001, '(6f14.6)') freqj*RY_TO_CMM1
      
      xq_k = -(xq0 + xq_j)
      CALL freq_phq(xq_k, S, fc, freqk, U)
      bosek(:) = f_bose(freqk, input%T(1))
      !
      DO k = 1,S%nat3
        DO j = 1,S%nat3
          !
          bose_C = 2*(bosej(j) - bosek(k))
          dom_C(:) =nrg(:)+freqj(j)-freqk(k) ! cohalescence
          ctm_C = bose_C * f_gauss(dom_C, sigma_ry) !delta 
          !
          bose_X = bosej(j) + bosek(k) + 1
          dom_X(:) =nrg(:)-freqj(j)-freqk(k) ! scattering/decay
          ctm_X = bose_X * f_gauss(dom_X, sigma_ry) !delta
          !
          jd_C(:) = jd_C(:) + qgrid%w(jq)*ctm_C(:)
          jd_X(:) = jd_X(:) + qgrid%w(jq)*ctm_X(:)
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
    DO i = 1,input%ne
      WRITE(10000,'(4e24.15)') RY_TO_CMM1*nrg(i),jd_C(i)+jd_X(i),jd_C(i),jd_X(i)
    ENDDO
    !
  END SUBROUTINE

END MODULE r2q_program

PROGRAM r2q 

  USE kinds,            ONLY : DP
  USE r2q_program
  USE input_fc,         ONLY : read_fc2, aux_system, div_mass_fc2, &
                              forceconst2_grid, ph_system_info, &
                              multiply_mass_dyn, write_dyn
  USE asr2_module,      ONLY : impose_asr2
  !USE q_grids,          ONLY : q_grid
  USE constants,        ONLY : RY_TO_CMM1
  USE fc2_interpolate,  ONLY : freq_phq
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : code_input_type, READ_INPUT
  USE ph_velocity,      ONLY : velocity
  USE more_constants,   ONLY : print_citations_linewidth
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(ph_system_info)   :: S
  TYPE(q_grid)           :: qpath
  TYPE(code_input_type)  :: input
  !
  CHARACTER(len=512) :: filename
  !
  REAL(DP) :: xq(3)
  REAL(DP),ALLOCATABLE :: freq(:), vel(:,:)
  COMPLEX(DP),ALLOCATABLE :: U(:,:)
  INTEGER :: i, output_unit=10000
  CHARACTER (LEN=6),  EXTERNAL :: int_to_char
  !
  CALL print_citations_linewidth()
  !  
  CALL READ_INPUT("R2Q", input, qpath, S, fc2)
!   CALL read_fc2(file_mat2, S,  fc2)
!   CALL aux_system(S)
!   IF(asr2) CALL impose_asr2("simple",S%nat, fc2)
!   CALL div_mass_fc2(S, fc2)

  IF( input%calculation=="jdos") THEN
    CALL joint_dos(input,S,fc2)
  ELSE
    ALLOCATE(freq(S%nat3))
    ALLOCATE(U(S%nat3,S%nat3))

    filename=TRIM(input%outdir)//"/"//TRIM(input%prefix)//".out"
    OPEN(unit=output_unit, file=filename)

    IF(input%print_velocity) THEN
      filename=TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_vel.out"
      OPEN(unit=output_unit+1, file=filename)
      ALLOCATE(vel(3,S%nat3))
    ENDIF
    !
    DO i = 1,qpath%nq
      CALL freq_phq(qpath%xq(:,i), S, fc2, freq, U)
      !WRITE(*, '(a16,999f12.4)') "freq", freq*RY_TO_CMM1
      !WRITE(*, '(999f12.4)') freq*RY_TO_CMM1
      ioWRITE(output_unit, '(i6,f12.6,3x,3f12.6,999e16.6)') &
        i, qpath%w(i), qpath%xq(:,i), freq*RY_TO_CMM1
      ioFLUSH(output_unit)
      
      IF(input%print_dynmat) THEN
        U = multiply_mass_dyn(S, U)
        filename = TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_dyn"//TRIM(int_to_char(i))
        CALL write_dyn(filename, qpath%xq(:,i), U, S)
      ENDIF

      IF(input%print_velocity) THEN
        vel = velocity(S, fc2, qpath%xq(:,i))
        ioWRITE(output_unit+1, '(i6,f12.6,3x,3f12.6,999(3e16.8,3x))') &
          i, qpath%w(i), qpath%xq(:,i), vel*RY_TO_CMM1
        ioFLUSH(output_unit+1)
      ENDIF

    ENDDO
    !
    CLOSE(output_unit)
    DEALLOCATE(freq, U)
    IF(input%print_velocity) THEN
      CLOSE(output_unit+1)
      DEALLOCATE(vel)
    ENDIF
    !
  ENDIF
  
END PROGRAM r2q
