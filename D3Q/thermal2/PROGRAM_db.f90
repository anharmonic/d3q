!
! Written by Lorenzo Paulatto (2014-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE add_bubble_program
#include "mpi_thermal.h"
  USE kinds, ONLY : DP

  CONTAINS 
  SUBROUTINE DYNBUBBLE_PATH(input,qpath, S,fc2,fc2b,fc3)
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info, multiply_mass_dyn, &
                                 write_dyn, read_fc2, aux_system, div_mass_fc2
    USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
    USE q_grids,          ONLY : q_grid, setup_grid
    USE fc3_interpolate,  ONLY : forceconst3
    USE code_input,       ONLY : READ_INPUT, code_input_type
    USE dynbubble,        ONLY : dynbubble_q
    USE constants,        ONLY : RY_TO_CMM1
    USE more_constants,   ONLY : write_conf
    IMPLICIT NONE
    TYPE(forceconst2_grid) :: fc2, fc2b
    CLASS(forceconst3),INTENT(in) :: fc3
    TYPE(ph_system_info)   :: S
    TYPE(code_input_type)     :: input
    TYPE(q_grid)      :: qpath
    
    TYPE(q_grid) :: grid

    INTEGER :: iq, it, nu, mu
    COMPLEX(DP),ALLOCATABLE :: dyn(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: dyn0(:,:), U(:,:), dynX(:,:), dynY(:,:)
    REAL(DP),ALLOCATABLE    :: freq(:), freq0(:), freqY(:)
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    CHARACTER(len=512) :: filename
  
    ALLOCATE(dyn(S%nat3,S%nat3,input%nconf))
    ALLOCATE(dyn0(S%nat3,S%nat3))
    ALLOCATE(dynX(S%nat3,S%nat3))
    ALLOCATE(dynY(S%nat3,S%nat3))
    ALLOCATE(U(S%nat3,S%nat3))
    ALLOCATE(freq(S%nat3), freq0(S%nat3), freqY(S%nat3))

    DO it = 1,input%nconf
      filename=TRIM(input%outdir)//"/"//&
              TRIM(input%prefix)//"_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out"
      !ioWRITE(*,*) "filename", TRIM(filename)
      OPEN(unit=1000+it, file=filename)
      IF (TRIM(input%mode) == "full") THEN
        ioWRITE(1000+it, *) "# calculation of linewidth (gamma_n) [and lineshift (delta_n)]"
      ELSE
        ioWRITE(1000+it, *) "# calculation of linewidth (gamma_n)"
      ENDIF
      ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      FLUSH(1000+it)
    ENDDO

    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
    CALL grid%scatter()
  
    DO iq = 1, qpath%nq
      ioWRITE(*, *) "<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>"
      ioWRITE(*, '(i6,3f12.6)') iq, qpath%xq(:,iq)
      CALL fftinterp_mat2(qpath%xq(:,iq), S%nat3, fc2b, dyn0)
      U = dyn0
      ioWRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, U)
      CALL mat2_diag(S%nat3, U, freq0)
      freq0 = SIGN(SQRT(ABS(freq0)), freq0)
      ioWRITE(*,"(6f20.12)") freq0*RY_TO_CMM1
      ioWRITE(*,*)

      ! Allocation is automatic in fortran 2003, but of course it does not work with most compilers
      dyn = dynbubble_q(qpath%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, S, grid, fc2, fc3)
      
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
        !filename = "dyn_conf"//TRIM(int_to_char(it))//"_"//TRIM(int_to_char(iq))
        IF(input%print_dynmat) THEN
          filename = "dyn_"//TRIM(input%prefix)//&
                    "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                    "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//&
                    "_q"//TRIM(int_to_char(iq))
          CALL write_dyn(filename, qpath%xq(:,iq), dynY, S)
        ENDIF
        !
        dynY = dynX
        CALL mat2_diag(S%nat3, dynX, freq)
        freq = SIGN(SQRT(ABS(freq)), freq)
        ioWRITE(1000+it,"(i6,4f12.6,2(12e20.6,5x))") iq,qpath%w(iq), qpath%xq(:,iq), freq0*RY_TO_CMM1,freq*RY_TO_CMM1
        FLUSH(1000+it)
        
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

    DO it = 1,input%nconf
      CLOSE(1000+it)
    ENDDO
  END SUBROUTINE DYNBUBBLE_PATH


  SUBROUTINE DYNSPECTRE_PATH(input, qpath, S, fc2, fc2b, fc3)
    USE fc2_interpolate,     ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,      ONLY : add_exp_t_factor
    USE dynbubble,      ONLY : tracebubble_q
    USE constants,      ONLY : RY_TO_CMM1
    USE q_grids,        ONLY : q_grid, setup_grid
    USE more_constants, ONLY : write_conf
    USE input_fc,       ONLY : forceconst2_grid, ph_system_info
    USE fc3_interpolate,ONLY : forceconst3
    USE code_input,     ONLY : code_input_type
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2, fc2b
    CLASS(forceconst3),INTENT(inout)  :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: qpath
    !
    INTEGER :: iq, it, ie
    TYPE(q_grid) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP)   :: sigma_ry(input%nconf)
    !
    REAL(DP),ALLOCATABLE :: ener(:), spectralf(:,:,:)
    !
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)
    ALLOCATE(spectralf(input%ne,S%nat3,input%nconf))
    !
    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
    CALL grid%scatter()
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    ener = ener/RY_TO_CMM1
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    IF(ionode)THEN
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//"_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                              "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
      ioWRITE(1000+it, *) "# spectral function mode: ", input%mode
      ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      ioWRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      FLUSH(1000+it)
    ENDDO
    ENDIF
    !
    ioWRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points (2)"
    
    DO iq = 1,qpath%nq
      CALL freq_phq(qpath%xq(:,iq), S, fc2, w2, D)
      ioWRITE(*,'(i6,3f12.6,5x,6f12.6,100(/,47x,6f12.6))') iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      !
      DO it = 1,input%nconf
        ioWRITE(1000+it, *)
        ioWRITE(1000+it, '(a,i6,3f15.8,5x,100(/,"#",6f12.6))') "#  xq",  iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      ENDDO
      !
      spectralf = tracebubble_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                                S, grid, fc2, fc3, input%ne, ener)
      !
      IF(input%exp_t_factor) CALL add_exp_t_factor(input%nconf, input%T, input%ne, S%nat3, ener, spectralf)
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          ioWRITE(1000+it, '(2f14.8,100e14.6)') &
                qpath%w(iq), ener(ie)*RY_TO_CMM1, SUM(spectralf(ie,:,it))/RY_TO_CMM1, &
                spectralf(ie,:,it)/RY_TO_CMM1
          FLUSH(1000+it)
        ENDDO
      ENDDO
      !
    ENDDO
    !
    !
    IF(ionode)THEN
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    ENDIF
    !
    DEALLOCATE(spectralf, ener)
    !
  END SUBROUTINE DYNSPECTRE_PATH

END MODULE add_bubble_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM add_bubble

  USE kinds,            ONLY : DP
  USE add_bubble_program
  USE fc3_interpolate,  ONLY : forceconst3
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, &
                               ph_system_info, same_system, read_fc2, aux_system, div_mass_fc2
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE asr2_module,      ONLY : impose_asr2

  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2, fc2b
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S, Sb
  TYPE(code_input_type)     :: input
  TYPE(q_grid)      :: qpts
  !

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL start_mpi()
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("DB", input, qpts, S, fc2, fc3)
  !
  ! Do some extra input, only for dynbubble
  CALL read_fc2(input%file_mat2_final, Sb,  fc2b)
  CALL aux_system(Sb)
  IF(.not. same_system(S,Sb)) CALL errore("PROGRAM_db", "S and Sb do not match", 1)
  CALL impose_asr2(input%asr2, Sb%nat, fc2b)
  CALL div_mass_fc2(Sb, fc2b)
  !
  IF(input%calculation == "db")THEN
    CALL DYNBUBBLE_PATH(input,qpts, S,fc2,fc2b,fc3)
  ELSEIF(input%calculation == "dbs")THEN
    CALL DYNSPECTRE_PATH(input,qpts, S,fc2,fc2b,fc3)
  
  ELSE
    CALL errore("db", "no such calculation: "//TRIM(input%calculation),1)
  ENDIF
  !

  CALL stop_mpi()
 
END PROGRAM add_bubble
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

