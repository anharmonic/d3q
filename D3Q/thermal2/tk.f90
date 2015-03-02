!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
#define timer_CALL CALL
!
MODULE thermalk_program
  !
  USE kinds,    ONLY : DP
  !
  CONTAINS
  ! Test subroutine: compute phonon frequencies along a line and save them to unit 666  
  SUBROUTINE TK_SMA(input, qgrid, S, fc2, fc3)
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf
    USE q_grid,             ONLY : q_grid_type, setup_simple_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE fc2_interpolate,    ONLY : freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity_proj
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qgrid
    !
!    TYPE(q_grid_type) :: kgrid
    REAL(DP) :: sigma_ry(input%nconf)
    REAL(DP) :: lw(S%nat3,input%nconf)
    REAL(DP) :: lw_isotopic(S%nat3,input%nconf)
    REAL(DP) :: lw_phph(S%nat3,input%nconf)
    REAL(DP) :: lw_casimir(S%nat3)
    REAL(DP) :: freq(S%nat3)
    REAL(DP) :: bose(S%nat3,input%nconf)
    COMPLEX(DP) :: U(S%nat3,S%nat3)
    !
    REAL(DP) :: tk(3,3,input%nconf)
    REAL(DP) :: vel(3,S%nat3)
    REAL(DP) :: dq, pref
    INTEGER  :: iq, it, a, b, nu
    !
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    !
    sigma_ry = input%sigma/RY_TO_CMM1
!    CALL setup_simple_grid(S, input%nk(1), input%nk(2), input%nk(3), kgrid)
    !
    WRITE(*,'(1x,a,i10,a)') "Integrating over a grid of", qgrid%nq, " points"
    !
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              "lw."//TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
      WRITE(1000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
      WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
      !
      IF(input%isotopic_disorder) THEN
        OPEN(unit=2000+it, file=TRIM(input%outdir)//"/"//&
                                "lwiso."//TRIM(input%prefix)//".T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                  "s"//TRIM(write_conf(it,input%nconf,input%sigma))//"out")
        WRITE(1000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
        WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        CALL flush_unit(2000+it)
      ENDIF
    ENDDO
      IF(input%casimir_scattering) THEN
        OPEN(unit=3000, file=TRIM(input%outdir)//"/"//&
                                "lwcas."//TRIM(input%prefix)//".out")
        WRITE(3000, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
!         WRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        CALL flush_unit(3000)
      ENDIF
    !
    tk = 0._dp
    dq = S%Omega/qgrid%nq
    !
      timer_CALL t_tksma%start()
    QPOINT_LOOP : &
    DO iq = 1,qgrid%nq
      WRITE(*,'(i6,3f15.8)') iq, qgrid%xq(:,iq)
      !
        timer_CALL t_lwphph%start()
      lw_phph = linewidth_q(qgrid%xq(:,iq), input%nconf, input%T,&
                       sigma_ry, S, qgrid, fc2, fc3)
      IF(ANY(lw<0._dp)) WRITE(*,'(a,99e12.4)') "Negative LW!", lw
        timer_CALL t_lwphph%stop()
      !
      ! Compute contribution of isotopic disorder
      IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw_isotopic = isotopic_linewidth_q(qgrid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, qgrid, fc2)
        IF(ANY(lw_isotopic<0._dp)) WRITE(*,'(a,99e12.4)') "Negative LW isotopic!", lw
          timer_CALL t_lwisot%stop()
      ELSE
        lw_isotopic = 0._dp
      ENDIF
      !
      !
        timer_CALL t_freq%start() 
      ! Velocity
      vel = velocity_proj(S, fc2, qgrid%xq(:,iq))
      CALL  freq_phq_safe(qgrid%xq(:,iq), S, fc2, freq, U)
        timer_CALL t_freq%stop() 
      !
      ! Compute anisotropic Casimir linewidth
      IF(input%casimir_scattering) THEN
          timer_CALL t_lwcasi%start() 
        lw_casimir = casimir_linewidth_vel(vel, input%casimir_length, input%casimir_dir, S%nat3)
          timer_CALL t_lwcasi%stop() 
      ELSE
        lw_casimir = 0._dp
      ENDIF
      !
      DO it = 1, input%nconf
        WRITE(1000+it,'(3f12.6,99e20.10)') qgrid%xq(:,iq), lw(:,it)*RY_TO_CMM1
        IF(input%isotopic_disorder) &
          WRITE(2000+it,'(3f12.6,99e20.10)') qgrid%xq(:,iq), lw_isotopic(:,it)*RY_TO_CMM1
      ENDDO
      IF(input%casimir_scattering) &
        WRITE(3000,'(3f12.6,99e20.10)') qgrid%xq(:,iq), lw_casimir(:)*RY_TO_CMM1
      !
        timer_CALL t_tksum%start()
      !
      DO it = 1,input%nconf
        lw(:,it) = lw_phph(:,it) + lw_isotopic(:,it) + lw_casimir
      ENDDO
      !
      CONF_LOOP : &
      DO it = 1, input%nconf
        !
        IF(input%T(it)==0._dp) CYCLE CONF_LOOP
        !
        CALL  bose_phq(input%T(it), S%nat3, freq, bose(:,it))
        !
        MODE_LOOP : &
        DO nu = 1, S%nat3
          ! Check if we have zero linewidth and non-zero velocity it is a problem
          ! lw can be NaN when T=0 and xq=0, check for lw>0 insteand, because NaN/=0 is true
          IF(lw(nu,it)<0._dp)THEN
            WRITE(*,"(3x,a,e12.4,3i6)") "warning: negative lw (iq,nu,it):", lw(nu,it), iq, nu, it 
            lw(nu,it) = - lw(nu,it)
          ENDIF
          IF(.not. lw(nu,it)>0._dp)THEN
            IF(ANY(ABS(vel(:,nu))>eps_vel ))THEN
              WRITE(*,*) iq, it, lw(nu,it), vel(:,nu)
              CALL errore("TK_SMA", "cannot threat this case", 1)
            ELSE
              WRITE(*,"(3x,a,3i6)") "skip (iq,nu,it):", iq, nu, it
              CYCLE MODE_LOOP 
            ENDIF
          ENDIF
          !
          pref = freq(nu)**2 *bose(nu,it)*(1+bose(nu,it))/input%T(it)**2 *dq /lw(nu,it)
          !WRITE(*,"(3x,a,3i6,4e15.6)") "do:", iq, nu, it, pref, freq(nu), bose(nu,it), lw(nu,it)
          DO a = 1,3
          DO b = 1,3
            tk(a,b,it) = tk(a,b,it) + pref*vel(a,nu)*vel(b,nu)
          ENDDO
          ENDDO
        ENDDO MODE_LOOP 
        !
      ENDDO CONF_LOOP 
        timer_CALL t_tksum%stop()
      !
    ENDDO QPOINT_LOOP
      timer_CALL t_tksma%stop()
      
    DO it = 1, input%nconf
      CLOSE(1000+it)
      IF(input%isotopic_disorder) CLOSE(2000+it)
    ENDDO
    IF(input%casimir_scattering) CLOSE(3000)
    !
    DO it = 1,input%nconf
      WRITE(*,"(a,i3,2f12.6)") "conf:", it, input%sigma(it), input%T(it)
      WRITE(*,"(3x,3e20.6)") tk(:,1,it)*RY_TO_WATTMM1KM1
      WRITE(*,"(3x,3e20.6)") tk(:,2,it)*RY_TO_WATTMM1KM1
      WRITE(*,"(3x,3e20.6)") tk(:,3,it)*RY_TO_WATTMM1KM1
      WRITE(*,"(3x,a)") "------------"
    ENDDO
    
#ifdef timer_CALL
      CALL t_tksma%print()
      WRITE(*,'(a)') "*** * Contributions to SMA conductivity:"
      CALL t_tksum%print()
      CALL t_lwisot%print()
      CALL t_lwcasi%print()
      CALL t_lwphph%print()
      WRITE(*,'(a)') "*** * Contributions to ph-ph linewidth time:"
      CALL t_freq%print() 
      CALL t_bose%print() 
      CALL t_sum%print() 
      CALL t_fc3int%print() 
      CALL t_fc3m2%print() 
      CALL t_fc3rot%print() 
#endif

    !
    !
  END SUBROUTINE TK_SMA
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  SUBROUTINE TK_ITERATIVE(input, qgrid, S, fc2, fc3)
    USE fc2_interpolate,         ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE q_grid,             ONLY : q_grid_type, setup_simple_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: qgrid
    !
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)  
    INTEGER :: iq, it
!    TYPE(q_grid_type) :: kgrid
    REAL(DP)   :: sigma_ry(input%nconf)
  END SUBROUTINE TK_ITERATIVE
  !   
  END MODULE thermalk_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!               !               !               !               !               !               !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM linewidth

  USE kinds,            ONLY : DP
  USE thermalk_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, ph_system_info
  USE q_grid,           ONLY : q_grid_type !, setup_simple_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)   :: S
  TYPE(code_input_type)     :: tkinput
  TYPE(q_grid_type)      :: qgrid

!   CALL mp_world_start(world_comm)
!   CALL environment_start('LW')
  CALL print_citations_linewidth()

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("TK", tkinput, qgrid, S, fc2, fc3)
  !
  IF(TRIM(tkinput%calculation) == "sma") THEN
    !
    CALL TK_SMA(tkinput, qgrid, S, fc2, fc3)
    !
  ELSE &
  IF(TRIM(tkinput%calculation) == "iterative") THEN
    !
    CALL TK_ITERATIVE(tkinput, qgrid, S, fc2, fc3)
    !
  ELSE
    CALL errore("lw", "what else to do?", 1)
  ENDIF
  !
!   CALL environment_end('LW')
!   CALL mp_world_end()
 
END PROGRAM
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

