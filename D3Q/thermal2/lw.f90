!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! CONVENTIONS :
! xR, xq --> cartesian coordinates
! yR, yq --> crystalline coordinates
!
#define timer_CALL CALL
!
MODULE linewidth_program
  !
  USE kinds,    ONLY : DP
  !
  CONTAINS
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE LW_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE fc2_interpolate,         ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q, selfnrg_q, spectre_q
    USE constants,          ONLY : RY_TO_CMM1
    USE q_grids,            ONLY : q_grid, setup_simple_grid
    USE more_constants,     ONLY : write_conf
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: qpath
    !
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)
    REAL(DP) :: pl,dpl
    INTEGER :: iq, it
    TYPE(q_grid) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP)   :: lw(S%nat3,input%nconf)
    REAL(DP)   :: lw_casimir(S%nat3)
    REAL(DP)   :: sigma_ry(input%nconf)
    CHARACTER(len=15) :: f1, f2
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
      IF (TRIM(input%mode) == "full") THEN
        WRITE(1000+it, *) "# calculation of linewidth (gamma_n) [and lineshift (delta_n)]"
      ELSE
        WRITE(1000+it, *) "# calculation of linewidth (gamma_n)"
      ENDIF
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
    ENDDO
    ! Prepare formats to write out data
    WRITE(f1,'(i6,a)') S%nat3, "f12.6,6x,"
    WRITE(f2,'(i6,a)') S%nat3, "e15.5,6x,"
    !
    ! Gaussian: exp(x^2/(2s^2)) => FWHM = 2sqrt(2log(2)) s
    ! Wrong Gaussian exp(x^2/c^2) => FWHM = 2 sqrt(log(2)) c
    ! Lorentzian: (g/2)/(x^2 + (g/2)^2) => FWHM = g
    ! Wrong Lorentzian: d/(x^2+d^2) => FWHM = 2d
    !  => 2d = 2 sqrt(log(2) c
    !      d = sqrt(log(2)) c
    !      d = 0.83255 c
    ! IMHO: you need to use a sigma that is 0.6 (=0.5/0.83255) times smaller when using
    ! linewidth_q than when using selfnrg_q in order to get the same values
    ! THIS FACTOR IS NOW INCLUDED DIRECTLY IN SUM_LINEWIDTH_MODES!
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    !
    DO iq = 1,qpath%nq
      WRITE(*,'(i6,3f15.8)') iq, qpath%xq(:,iq)
      !
      CALL freq_phq(qpath%xq(:,iq), S, fc2, w2, D)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      IF (TRIM(input%mode) == "full") THEN
        ls = selfnrg_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                        S, grid, fc2, fc3)
        DO it = 1,input%nconf
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,'//f1//f2//f2//'x)') &
                iq,pl,qpath%xq(:,iq), w2*RY_TO_CMM1, -DIMAG(ls(:,it))*RY_TO_CMM1, DBLE(ls(:,it))*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE IF (TRIM(input%mode) == "real") THEN
          timer_CALL t_lwphph%start()
        lw = linewidth_q(qpath%xq(:,iq), input%nconf, input%T,  sigma_ry, &
                        S, grid, fc2, fc3)
          timer_CALL t_lwphph%stop()
        IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw = lw + isotopic_linewidth_q(qpath%xq(:,iq), input%nconf, input%T,  sigma_ry, &
                                      S, grid, fc2)
          timer_CALL t_lwisot%stop()
        ENDIF
        IF(input%casimir_scattering)THEN
          timer_CALL t_lwcasi%start()
         lw_casimir = casimir_linewidth_q(qpath%xq(:,iq), input%casimir_length, &
                                      input%casimir_dir, S, fc2)
         DO it = 1,input%nconf
          lw(:,it) = lw(:,it) + lw_casimir
         ENDDO
          timer_CALL t_lwcasi%stop()
        ENDIF
        
        DO it = 1,input%nconf
          WRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,'//f1//f2//'x)') &
                iq,pl,qpath%xq(:,iq), w2*RY_TO_CMM1, lw(:,it)*RY_TO_CMM1
          CALL flush_unit(1000+it)
        ENDDO
      ELSE
        CALL errore('LW_QBZ_LINE', 'wrong mode (real or full)', 1)
      ENDIF
      !
    ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
#ifdef timer_CALL
      timer_CALL t_lwisot%print()
      timer_CALL t_lwcasi%print()
      timer_CALL t_lwphph%print()
      WRITE(*,'(a)') "*** * Contributions to ph-ph linewidth time:"
      timer_CALL t_freq%print() 
      timer_CALL t_bose%print() 
      timer_CALL t_sum%print() 
      timer_CALL t_fc3int%print() 
      timer_CALL t_fc3m2%print() 
      timer_CALL t_fc3rot%print() 
#endif
    !
  END SUBROUTINE LW_QBZ_LINE
  !   
  !  
  SUBROUTINE SPECTR_QBZ_LINE(input, qpath, S, fc2, fc3)
    USE fc2_interpolate,     ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,      ONLY : spectre_q, simple_spectre_q, add_exp_t_factor
    USE constants,      ONLY : RY_TO_CMM1
    USE q_grids,        ONLY : q_grid, setup_simple_grid
    USE more_constants, ONLY : write_conf
    USE input_fc,       ONLY : forceconst2_grid, ph_system_info
    USE fc3_interpolate,ONLY : forceconst3
    USE code_input,     ONLY : code_input_type
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(inout)  :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: qpath
    !
    REAL(DP) :: pl,dpl
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
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    ener = ener/RY_TO_CMM1
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                              "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "# spectral function mode: ", input%mode
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      WRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      CALL flush_unit(1000+it)
    ENDDO
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
    DO iq = 1,qpath%nq
      CALL freq_phq(qpath%xq(:,iq), S, fc2, w2, D)
      WRITE(*,'(i6,3f15.8,5x,6f12.6)') iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      !
      DO it = 1,input%nconf
        WRITE(1000+it, *)
        WRITE(1000+it, '(a,i6,3f15.8,5x,6f12.6)') "#  xq",  iq, qpath%xq(:,iq), w2*RY_TO_CMM1
      ENDDO
      !
      IF (TRIM(input%mode) == "full") THEN
        spectralf = spectre_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                                  S, grid, fc2, fc3, input%ne, ener)
      ELSE IF (TRIM(input%mode) == "simple") THEN
        spectralf = simple_spectre_q(qpath%xq(:,iq), input%nconf, input%T, sigma_ry, &
                                  S, grid, fc2, fc3, input%ne, ener)
      ELSE
        CALL errore("SPECTR_QBZ_LINE", 'unknown mode "'//TRIM(input%mode)//'"', 1)
      ENDIF
      !
      IF(input%exp_t_factor) CALL add_exp_t_factor(input%nconf, input%T, input%ne, S%nat3, ener, spectralf)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          WRITE(1000+it, '(2f14.8,100e14.6)') &
                pl, ener(ie)*RY_TO_CMM1, SUM(spectralf(ie,:,it))/RY_TO_CMM1**2, spectralf(ie,:,it)/RY_TO_CMM1**2
          CALL flush_unit(1000+it)
        ENDDO
      ENDDO
      !
    ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    DEALLOCATE(spectralf, ener)
    !
  END SUBROUTINE SPECTR_QBZ_LINE
  !   
  !  
  SUBROUTINE FINAL_STATE_LINE(input, qpath, S, fc2, fc3)
    USE fc2_interpolate,     ONLY : fftinterp_mat2, mat2_diag
    USE  final_state,   ONLY : final_state_q
    USE constants,      ONLY : RY_TO_CMM1
    USE q_grids,        ONLY : q_grid, setup_simple_grid
    USE more_constants, ONLY : write_conf
    USE input_fc,       ONLY : forceconst2_grid, ph_system_info
    USE fc3_interpolate,ONLY : forceconst3
    USE code_input,     ONLY : code_input_type
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(inout)  :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: qpath
    !
    REAL(DP) :: pl,dpl, e_inital_ry 
    INTEGER :: iq, it, ie
    TYPE(q_grid) :: grid
    COMPLEX(DP):: ls(S%nat3,input%nconf)
    REAL(DP) :: sigma_ry(input%nconf)
    !
    REAL(DP),PARAMETER :: inv2_RY_TO_CMM1 = 1/(RY_TO_CMM1)**2
    !
    REAL(DP),ALLOCATABLE :: ener(:), fstate(:,:,:)
    CHARACTER(len=256) :: filename
    !
    ALLOCATE(fstate(input%ne,S%nat3,input%nconf))
    !
    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), grid)
    !
    ALLOCATE(ener(input%ne))
    FORALL(ie = 1:input%ne) ener(ie) = (ie-1)*input%de+input%e0
    ! Convert to Rydberg:
    ener = ener/RY_TO_CMM1
    e_inital_ry = input%e_initial/RY_TO_CMM1
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    DO it = 1,input%nconf
      filename = TRIM(input%outdir)//"/"//&
                              TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                              "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out"
      OPEN(unit=1000+it, file=filename)
      WRITE(*,*) "opening ", TRIM(filename)
      WRITE(1000+it, *) "# final state decompositions, mode: ", input%mode
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "#", it, "T=",input%T(it), "sigma=", input%sigma(it)
      WRITE(1000+it, *) "#   q-path     energy (cm^-1)         total      band1      band2    ....     "
      CALL flush_unit(1000+it)
    ENDDO
    !
    dpl = 0._dp; pl = 0._dp
    !
    WRITE(*,'(1x,a,i6,a)') "Going to compute", qpath%nq, " points"
    
!     DO iq = 1,qpath%nq
!       WRITE(*,'(i6,3f15.8)') iq, qpath%xq(:,iq)
!       !
!       DO it = 1,input%nconf
!         WRITE(1000+it, *)
!         WRITE(1000+it, '(a,i6,3f15.8)') "#  xq",  iq, qpath%xq(:,iq)
!       ENDDO
      !
      fstate = final_state_q(input%q_initial, qpath, input%nconf, input%T, sigma_ry, &
                             S, grid, fc2, fc3, e_inital_ry, input%ne, ener, &
                             input%q_resolved, input%sigmaq, input%outdir, input%prefix)
      !
      IF(iq>1) dpl = SQRT(SUM( (qpath%xq(:,iq-1)-qpath%xq(:,iq))**2 ))
      pl = pl + dpl
      !
      DO it = 1,input%nconf
        DO ie = 1,input%ne
          WRITE(1000+it, '(2f14.8,100e18.6e4)') &
                pl, ener(ie), SUM(fstate(ie,:,it)), fstate(ie,:,it)
          CALL flush_unit(1000+it)
        ENDDO
      ENDDO
      !
!     ENDDO
    !
    !
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    !
    DEALLOCATE(fstate, ener)
    !
  END SUBROUTINE FINAL_STATE_LINE
  !   
  END MODULE linewidth_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!               !               !               !               !               !               !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE linewidth_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid !, setup_simple_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: lwinput
  TYPE(q_grid)      :: qpath

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("LW", lwinput, qpath, S, fc2, fc3)
  !
  IF(    TRIM(lwinput%calculation) == "lw"   &
    .or. TRIM(lwinput%calculation) == "grid" &
    ) THEN
    !
    CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE &
  IF(TRIM(lwinput%calculation) == "spf") THEN
    !
    CALL SPECTR_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
    !
!   ELSE &
!   IF(TRIM(lwinput%calculation) == "grid") THEN
!     !
!     ! Compute the linewidth over the grid, will be reused later to do self-consistent linewidth
!     CALL LW_QBZ_LINE(lwinput, qpath, S, fc2, fc3)
  ELSE &
  IF(TRIM(lwinput%calculation) == "final") THEN
    !
    CALL FINAL_STATE_LINE(lwinput, qpath, S, fc2, fc3)
    !
  ELSE
    CALL errore("lw", "what else to do?", 1)
  ENDIF
  !
!   CALL environment_end('LW')
!   CALL mp_world_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

