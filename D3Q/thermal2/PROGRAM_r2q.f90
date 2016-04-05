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
                              print_citations_linewidth
  USE asr2_module,      ONLY : impose_asr2
  !USE q_grids,          ONLY : q_grid
  USE constants,        ONLY : RY_TO_CMM1
  USE fc2_interpolate,       ONLY : freq_phq_safe
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(ph_system_info)   :: S
  !TYPE(q_grid)      :: qpath
  !
  LOGICAL :: asr2 = .true.
  CHARACTER(len=512) :: file_mat2 = "quter.fc"
  !
  REAL(DP) :: xq(3)
  REAL(DP),ALLOCATABLE :: freq(:)
  COMPLEX(DP),ALLOCATABLE :: U(:,:)
  INTEGER :: i
  !
  CALL print_citations_linewidth()
  !  
  CALL read_fc2(file_mat2, S,  fc2)
  CALL aux_system(S)
  IF(asr2) CALL impose_asr2("simple",S%nat, fc2)
  CALL div_mass_fc2(S, fc2)


  ALLOCATE(freq(S%nat3))
  ALLOCATE(U(S%nat3,S%nat3))
  !xq = (/0.333333333333_dp,  0.1924500_dp,  0.0_dp /)
  xq = 0._dp
  DO i = 1,26
    xq(3) = DBLE(i-1)/DBLE(25)
    CALL freq_phq_safe(xq, S, fc2, freq, U)
    !WRITE(*, '(a16,999f12.4)') "freq", freq*RY_TO_CMM1
    WRITE(*, '(999f12.4)') freq*RY_TO_CMM1
    WRITE(998, '(999f12.4)') xq
    WRITE(999, '(999f12.4)') freq*RY_TO_CMM1
  ENDDO

END PROGRAM r2q
