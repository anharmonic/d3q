!
! Written by Lorenzo Paulatto (2014-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE r2q_program

 ! this space unintentionally left non-blank
  
END MODULE r2q_program

PROGRAM r2q 

  USE kinds,            ONLY : DP
  USE r2q_program
  USE input_fc,         ONLY : read_fc2, aux_system, div_mass_fc2, &
                              forceconst2_grid, ph_system_info, &
                              print_citations_linewidth, &
                              multiply_mass_dyn, write_dyn
  USE asr2_module,      ONLY : impose_asr2
  !USE q_grids,          ONLY : q_grid
  USE constants,        ONLY : RY_TO_CMM1
  USE fc2_interpolate,  ONLY : freq_phq
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : code_input_type, READ_INPUT
  USE ph_velocity,      ONLY : velocity_proj
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
    WRITE(output_unit, '(i6,f12.6,3x,3f12.6,999e16.6)') &
      i, qpath%w(i), qpath%xq(:,i), freq*RY_TO_CMM1
    FLUSH(output_unit)
    
    IF(input%print_dynmat) THEN
      U = multiply_mass_dyn(S, U)
      filename = TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_dyn"//TRIM(int_to_char(i))
      CALL write_dyn(filename, qpath%xq(:,i), U, S)
    ENDIF

    IF(input%print_velocity) THEN
      vel = velocity_proj(S, fc2, qpath%xq(:,i))
      WRITE(output_unit+1, '(i6,f12.6,3x,3f12.6,999(3e16.8,3x))') &
        i, qpath%w(i), qpath%xq(:,i), vel*RY_TO_CMM1
      FLUSH(output_unit+1)
    ENDIF

  ENDDO
  !
  CLOSE(output_unit)
  DEALLOCATE(freq, U)
  IF(input%print_velocity) THEN
    CLOSE(output_unit+1)
    DEALLOCATE(vel)
  ENDIF
  
END PROGRAM r2q
