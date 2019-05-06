!
! Written by Lorenzo Paulatto (2014-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE add_bubble_program
#include "mpi_thermal.h"
  USE kinds, ONLY : DP

  CONTAINS 
  SUBROUTINE DYNBUBBLE_PATH(input,qpath, S,fc2,fc2b,fc3)
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info, multiply_mass_dyn, &
                                 write_dyn, read_fc2, aux_system, div_mass_fc2
    USE fc2_interpolate,  ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat, set_nu0, freq_phq_safe
    USE q_grids,          ONLY : q_grid, setup_grid
    USE fc3_interpolate,  ONLY : forceconst3
    USE code_input,       ONLY : READ_INPUT, code_input_type
    USE dynbubble,        ONLY : dynbubble_q
    USE constants,        ONLY : RY_TO_CMM1
    USE more_constants,   ONLY : write_conf
    USE functions,        ONLY : quicksort
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
    REAL(DP),ALLOCATABLE    :: freq(:), freq0(:), freqY(:), freqZ(:)
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    CHARACTER(len=512) :: filename
    ! RAFTEST
    INTEGER :: nu0
    REAL(DP) :: w(S%nat3),wm1(S%nat3)
    COMPLEX(DP) :: D(S%nat3,S%nat3)

    ALLOCATE(dyn(S%nat3,S%nat3,input%nconf))
    ALLOCATE(dyn0(S%nat3,S%nat3))
    ALLOCATE(dynX(S%nat3,S%nat3))
    ALLOCATE(dynY(S%nat3,S%nat3))
    ALLOCATE(U(S%nat3,S%nat3))
    ALLOCATE(freq(S%nat3), freq0(S%nat3), freqY(S%nat3), freqZ(S%nat3) )

    ! RAFTEST
    !write(*,*)  RY_TO_CMM1
    write(*,*) 'PROGRAM STATIC'

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
      ioFLUSH(1000+it)
    ENDDO

    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
    CALL grid%scatter()
  
    DO iq = 1, qpath%nq
      ioWRITE(stdout, *) "<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>"
      ioWRITE(stdout, '(i6,3f12.6)') iq, qpath%xq(:,iq)
      CALL fftinterp_mat2(qpath%xq(:,iq), S, fc2b, dyn0)
      U = dyn0
      !ioWRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, U)
      CALL mat2_diag(S%nat3, U, freq0)
      freq0 = SIGN(DSQRT(ABS(freq0)), freq0)
      ioWRITE(stdout,"(6f20.12)") freq0*RY_TO_CMM1
      ioWRITE(stdout,*)
      IF(ionode) FLUSH(stdout)

     
      ! RAFTEST
      !nu0 = set_nu0(qpath%xq(:,iq), S%at)
      !freqm1 = 0._dp
      !DO nu= 1,S%nat3
      ! IF (nu >=nu0) freqm1(nu)=0.5_dp/freq0(nu) 
      !END DO

      ! Compute 3rd order correction to dynamical matrix
      dyn = dynbubble_q(qpath%xq(:,iq), input%nconf, input%T, input%sigma/RY_TO_CMM1, S, grid, fc2, fc3)
     
      !
      nu0 = set_nu0(qpath%xq(:,iq), S%at)
      CALL freq_phq_safe(qpath%xq(:,iq), S,fc2,w,D)
      wm1 = 0._dp
      DO nu= 1,S%nat3
       IF (nu >=nu0) wm1(nu)=0.5_dp/w(nu)
      END DO
      
      !do nu=1,S%nat3      
      ! write(*,*) 'AAA', U(nu,nu), D(nu,nu)
      !end do

      DO it = 1, input%nconf
        ioWRITE(stdout, "(2f12.6)")  input%T(it),  input%sigma(it)
        dynX = dyn(:,:,it)
        
        !CALL dyn_cart2pat(dynX, S%nat3, U, +1)
        !DO mu = 1,S%nat3
        !DO nu = 1,S%nat3
        !  dynX(nu,mu) = dynX(nu,mu)*DSQRT(ABS(freq0(nu))*ABS(freq0(mu)))
        !ENDDO
        !ENDDO
        !CALL dyn_cart2pat(dynX, S%nat3, U, -1)

        !ioWRITE(*, "(6(2f12.6,4x))") multiply_mass_dyn(S, dynX)
        !ioWRITE(*, "(6(2f12.6,4x))") RY_TO_CMM1*dynX
        
        !dynX = multiply_mass_dyn(S, dynX) + dyn0
        dynX = dynX + dyn0
        !
        !
        IF(input%print_dynmat) THEN
          dynY = multiply_mass_dyn(S, dynX)
          filename = "dyn_"//TRIM(input%prefix)//&
                    "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                    "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//&
                    "_q"//TRIM(int_to_char(iq))
          CALL write_dyn(filename, qpath%xq(:,iq), dynY, S)
        ENDIF
        !
        dynY = dynX
        CALL mat2_diag(S%nat3, dynX, freq)
        freq = SIGN(DSQRT(ABS(freq)), freq)
        !
        ioWRITE(1000+it,"(i6,4f12.6,2(12e20.6,5x))") iq,qpath%w(iq), qpath%xq(:,iq),& 
         freq0*RY_TO_CMM1,freq*RY_TO_CMM1
        ! 
        ioFLUSH(1000+it)
        
        CALL dyn_cart2pat(dynY, S%nat3, U, +1)
        ! RAFTEST
        !DO nu=1,S%nat3
        ! DO mu=1,S%nat3
          !IF (nu/=mu) 
        !  write(*,*) nu,mu,dynY(nu,mu)
        ! END DO
        !END DO
        !       
        FORALL(nu=1:S%nat3) freqY(nu) = dynY(nu,nu)
        freqY = SIGN(DSQRT(ABS(freqY)), freqY)
        ! RAF
        call quicksort(freqY,1,S%nat3)
        dynX = dyn(:,:,it)
        CALL dyn_cart2pat(dynX, S%nat3, D, +1)
        DO nu=1, S%nat3
         freqZ(nu) = freq0(nu)+dynX(nu,nu)*wm1(nu)
        END DO
        call quicksort(freqZ,1,S%nat3)
        ! 
        ioWRITE(stdout,"('freq full correction')")
        ioWRITE(stdout,"('freqX ',6f20.12)") freq*RY_TO_CMM1
        ioWRITE(stdout,"('shiftX',6e20.6)") (freq-freq0)*RY_TO_CMM1
        ioWRITE(stdout,"('freq diag correction - from squared freq corr')")
        ioWRITE(stdout,"('freqY ',6f20.12)") freqY*RY_TO_CMM1
        ioWRITE(stdout,"('shiftY',6e20.6)") (freqY-freq0)*RY_TO_CMM1
        ioWRITE(stdout,"('freq diag correction')")
        ioWRITE(stdout,"('freqZ ',6f20.12)") freqZ*RY_TO_CMM1
        ioWRITE(stdout,"('shiftZ',6e20.6)") (freqZ-freq0)*RY_TO_CMM1
        ioWRITE(stdout,*)
        IF(ionode) FLUSH(stdout)
        !
        !
        ! RAFTEST
        !dynX=dynY
        !CALL mat2_diag(S%nat3, dynY, freq)
        !DO nu=1,S%nat3
        ! DO mu=1,S%nat3
          !IF (nu/=mu) write(*,*) dynX(nu,mu),dynY(nu,mu)
        !  write(*,*) freq(mu),dynX(mu,mu)
        ! END DO
        !END DO
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
    !
    ALLOCATE(spectralf(input%ne,S%nat3,input%nconf))
    !
    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), grid)
    CALL grid%scatter()
    !
    ALLOCATE(ener(input%ne))

    ! RAFTEST
    write(*,*) 'PROGRAM DYNAMIC BUBBLE'

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
      ioFLUSH(1000+it)
    ENDDO
    ENDIF
    !
    ioWRITE(*,'(2x,a,i6,a)') "Going to compute", qpath%nq, " points (2)"
    
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
      !IF(input%exp_t_factor) CALL add_exp_t_factor(input%nconf, input%T, input%ne, S%nat3, ener, spectralf)
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          ioWRITE(1000+it, '(2f14.8,100e14.6)') &
                qpath%w(iq), ener(ie)*RY_TO_CMM1, SUM(spectralf(ie,:,it))/RY_TO_CMM1, &
                spectralf(ie,:,it)/RY_TO_CMM1
          ioFLUSH(1000+it)
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

! RAF
END MODULE add_bubble_program

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM add_bubble

  USE kinds,            ONLY : DP
  USE add_bubble_program
  USE fc3_interpolate,  ONLY : forceconst3
  USE input_fc,         ONLY : forceconst2_grid, read_fc2, div_mass_fc2, &
                               ph_system_info, same_system, aux_system
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE asr2_module,      ONLY : impose_asr2
  USE more_constants,   ONLY : print_citations_linewidth

  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2, fc2b
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S, Sb
  TYPE(code_input_type)     :: input
  TYPE(q_grid)      :: qpts
  ! RAFTEST
  !LOGICAL  :: scale_fc3
  !REAL(DP) :: scale_fc3_fac
  !scale_fc3 = .true.
  !scale_fc3_fac = 1.0d-3
  !
  !

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL start_mpi()
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("DB", input, qpts, S, fc2, fc3)


 
  ! RAFTEST
  !print*, fc3%n_R
  !print*, fc3%nq
  !print*, fc3%xR3


  !IF (scale_fc3) THEN
  !  write(*,*) 'scaling of fc3 with factor', scale_fc3_fac
  !  fc3%sFC=fc3%sFC*scale_fc3_fac
  !END IF
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
  CALL print_citations_linewidth()
 
END PROGRAM add_bubble
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

