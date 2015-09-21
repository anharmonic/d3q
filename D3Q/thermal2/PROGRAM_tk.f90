!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!
MODULE thermalk_program
  !
#include "mpi_thermal.h"
  USE kinds,       ONLY : DP
  USE mpi_thermal, ONLY : ionode
  USE posix_signal,  ONLY : check_graceful_termination
  USE timers
  !
  CONTAINS
  ! 
  ! This subroutine computes the SMA thermal conducivity, it is mainly just a driver
  ! that uses other subroutines to compute th intrinsic, isotopic and casimir linewidths,
  ! than it sums everything up and takes care of input/output.
  !
  ! This subroutine is obsoleted by the first iteration of the variational method, 
  ! but we keep it for didactical purposes
  SUBROUTINE TK_SMA(input, out_grid, S, fc2, fc3)
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf
    USE q_grids,            ONLY : q_grid, setup_grid, setup_bz_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : ph_system_info
    USE code_input,         ONLY : code_input_type
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity_proj
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: out_grid
    !
    TYPE(q_grid) :: in_grid
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
    !
    ! We are using the same grid size for the inner and outer grid    
    ! this is not relly necessary in SMA, will fix later
    ! In any case we need to use to separate grids, because the 
    ! inner one (in_grid) is scatterd over MPI
    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), in_grid)
    CALL in_grid%scatter()
    !
    ioWRITE(stdout,'(1x,a,i10,a)') "Integrating over an inner grid of", in_grid%nq, " points"
    ioWRITE(stdout,'(1x,a,i10,a)') "Integrating over an outer grid of", out_grid%nq, " points"
    !
    IF(ionode)THEN
    DO it = 1,input%nconf
      OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                              "lw."//TRIM(input%prefix)//&
                               "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                               "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
      ioWRITE(1000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
      ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
      CALL flush_unit(1000+it)
      !
      IF(input%isotopic_disorder) THEN
        OPEN(unit=2000+it, file=TRIM(input%outdir)//"/"//&
                                "lwiso."//TRIM(input%prefix)//&
                                    "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                    "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
        ioWRITE(2000+it, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
        ioWRITE(2000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        CALL flush_unit(2000+it)
      ENDIF
    ENDDO
      IF(input%casimir_scattering) THEN
        OPEN(unit=3000, file=TRIM(input%outdir)//"/"//&
                                "lwcas."//TRIM(input%prefix)//".out")
        ioWRITE(3000, *) "# qpoint [2pi/alat], linewidth [cm^-1]"
!         ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it), "    sigma=", input%sigma(it)
        CALL flush_unit(3000)
      ENDIF
    ENDIF
    !
    tk = 0._dp
    dq = S%Omega/out_grid%nq
    !
      timer_CALL t_tksma%start()
    QPOINT_LOOP : &
    DO iq = 1,out_grid%nq
      ioWRITE(stdout,'(i6,3f15.8)') iq, out_grid%xq(:,iq)
      !
        timer_CALL t_lwphph%start()
      lw_phph = linewidth_q(out_grid%xq(:,iq), input%nconf, input%T,&
                       sigma_ry, S, in_grid, fc2, fc3)
      IF(ANY(lw<0._dp)) WRITE(stdout,'(a,99e12.4)') "Negative LW!", lw
        timer_CALL t_lwphph%stop()
      !
      ! Compute contribution of isotopic disorder
      IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw_isotopic = isotopic_linewidth_q(out_grid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, out_grid, fc2)
        IF(ANY(lw_isotopic<0._dp)) WRITE(stdout,'(a,99e12.4)') "Negative LW isotopic!", lw
          timer_CALL t_lwisot%stop()
      ELSE
        lw_isotopic = 0._dp
      ENDIF
      !
      !
        timer_CALL t_velcty%start() 
      ! Velocity
      vel = velocity_proj(S, fc2, out_grid%xq(:,iq))
      CALL  freq_phq_safe(out_grid%xq(:,iq), S, fc2, freq, U)
        timer_CALL t_velcty%stop() 
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
        timer_CALL t_lwinout%start()
      DO it = 1, input%nconf
        ioWRITE(1000+it,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw(:,it)*RY_TO_CMM1
        IF(input%isotopic_disorder) THEN
          ioWRITE(2000+it,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw_isotopic(:,it)*RY_TO_CMM1
        ENDIF
      ENDDO
      IF(input%casimir_scattering) THEN
        ioWRITE(3000,'(3f12.6,99e20.10)') out_grid%xq(:,iq), lw_casimir(:)*RY_TO_CMM1
      ENDIF
        timer_CALL t_lwinout%stop()
      !
        timer_CALL t_tksum%start()
      ! Casimir linewidth is temperature/smearing-independent, sum it to all configurations
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
          IF(lw(nu,it)<0._dp)THEN ! true for NaN
            WRITE(stdout,"(3x,a,e12.4,3i6)") "WARNING! Negative lw (idx q, mode, conf):", lw(nu,it), iq, nu, it 
            lw(nu,it) = - lw(nu,it)
          ENDIF
          IF(.not. lw(nu,it)>0._dp)THEN ! false for NaN
            IF(ANY(ABS(vel(:,nu))>eps_vel ))THEN
              WRITE(stdout,'(3i6,1e20.10,5x,3e20.10)') iq, nu, it, lw(nu,it), vel(:,nu)
              CALL errore("TK_SMA", "cannot threat this case", 1)
            ELSE
              !ioWRITE(stdout,"(3x,a,3i6)") "skip (iq,nu,it):", iq, nu, it
              CYCLE MODE_LOOP 
            ENDIF
          ENDIF
          !
          pref = freq(nu)**2 *bose(nu,it)*(1+bose(nu,it))/input%T(it)**2 *dq /lw(nu,it)
          !ioWRITE(stdout,"(3x,a,3i6,4e15.6)") "do:", iq, nu, it, pref, freq(nu), bose(nu,it), lw(nu,it)
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

    IF(ionode)THEN      
    DO it = 1, input%nconf
      CLOSE(1000+it)
      IF(input%isotopic_disorder) CLOSE(2000+it)
    ENDDO
    IF(input%casimir_scattering) CLOSE(3000)
    ENDIF
    !
    ! Write to disk
    IF(ionode) OPEN(unit=10000, file=TRIM(input%outdir)//"/"//&
                                     TRIM(input%prefix)//"."//"out")
    ioWRITE(10000,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                        "    K_x            K_y            K_z             ",&
                        "    K_xy           K_xz           K_yz            ",&
                        "    K_yx           K_zx           K_zy       "
    DO it = 1,input%nconf
      ioWRITE(10000,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
      tk(1,1,it)*RY_TO_WATTMM1KM1,tk(2,2,it)*RY_TO_WATTMM1KM1,tk(3,3,it)*RY_TO_WATTMM1KM1, &
      tk(1,2,it)*RY_TO_WATTMM1KM1,tk(1,3,it)*RY_TO_WATTMM1KM1,tk(2,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,1,it)*RY_TO_WATTMM1KM1,tk(3,1,it)*RY_TO_WATTMM1KM1,tk(3,2,it)*RY_TO_WATTMM1KM1
    ENDDO
    IF(ionode) CLOSE(10000)
    !
    ! Write to screen
    ioWRITE(stdout,"(3x,a,/,3x,a)") "************", "SMA thermal conductivity, also stored in file:"
    ioWRITE(stdout,'(5x,a)') TRIM(input%outdir)//"/"//TRIM(input%prefix)//"."//"out"
    DO it = 1,input%nconf
      ioWRITE(stdout,"(3x,a)") "**"
      ioWRITE(stdout,"(a,i3,2f12.6)") "conf:", it, input%sigma(it), input%T(it)
      ioWRITE(stdout,"(3x,3e20.6)") tk(:,1,it)*RY_TO_WATTMM1KM1
      ioWRITE(stdout,"(3x,3e20.6)") tk(:,2,it)*RY_TO_WATTMM1KM1
      ioWRITE(stdout,"(3x,3e20.6)") tk(:,3,it)*RY_TO_WATTMM1KM1
    ENDDO
    !
#ifdef timer_CALL
    ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
    CALL print_timers_header()
    CALL t_tksma%print()
    ioWRITE(stdout,'(a)') "*** * Contributions to SMA conductivity:"
    CALL t_tksum%print()
    CALL t_lwisot%print()
    CALL t_lwcasi%print()
    CALL t_lwphph%print()
    CALL t_velcty%print()
    CALL t_lwinout%print()
    CALL t_mpicom%print()
    CALL t_readdt%print()
    ioWRITE(stdout,'(a)') "*** * Contributions to ph-ph linewidth time:"
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
  SUBROUTINE TK_CG(input, out_grid, S, fc2, fc3)
    USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE q_grids,            ONLY : q_grid, q_basis, setup_grid, setup_bz_grid, &
                                   prepare_q_basis, qbasis_dot, qbasis_ax, &
                                   qbasis_a_over_b
    USE variational_tk
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: out_grid
    TYPE(q_grid)                 :: in_grid ! inner grid is MPI-scattered it is used 
                                            ! for integrating the ph-ph scattering terms and linewidth
    !
    INTEGER :: ix, nu, iq, it, nu0, iter
    !
    TYPE(q_basis) :: qbasis
    REAL(DP),ALLOCATABLE :: A_out(:,:,:), inv_sqrt_A_out(:,:,:), inv_A_out(:,:,:)
    REAL(DP),ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), h(:,:,:,:), t(:,:,:,:)
    REAL(DP),ALLOCATABLE :: Af(:,:,:,:)
    REAL(DP),ALLOCATABLE :: g_dot_h(:,:), h_dot_t(:,:), &
                            g_mod2(:,:), g_mod2_old(:,:), &
                            pref(:,:)
    INTEGER,PARAMETER :: niter_max = 10000
    REAL(DP) :: tk(3,3,input%nconf)
    !
    ! For code readability:
    INTEGER :: nconf, nat3, nq
    nconf = input%nconf
    nat3  = S%nat3
    nq    = out_grid%nq
    
    ! make the inner grid on top of the outer one
    CALL setup_grid(input%grid_type, S%bg, out_grid%n(1),out_grid%n(2),out_grid%n(3), in_grid)
    CALL in_grid%scatter()

    CALL prepare_q_basis(out_grid, qbasis, nconf, input%T, S, fc2)
    ! Compute A_out diagonal matrix
    ALLOCATE(A_out(nconf, nat3, nq))
    CALL compute_A_out(A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
    ! Compute 1/sqrt(A_out) from A_out
!     ALLOCATE(inv_sqrt_A_out(nconf, nat3, nq))
!     CALL compute_inv_sqrt_A_out(A_out, inv_sqrt_A_out, nconf, nat3, nq)
!     ! Compute 1/A_out from A_out
!     ALLOCATE(inv_A_out(nconf, nat3, nq))
!     CALL compute_inv_A_out(A_out, inv_A_out, nconf, nat3, nq)
    !
    !
    ! f0 = f_SMA = 1/(A_out) b
    ALLOCATE(f(3, nconf, nat3, nq))
    ALLOCATE(g(3, nconf, nat3, nq))
    ALLOCATE(h(3, nconf, nat3, nq))
    ALLOCATE(t(3, nconf, nat3, nq))
    ALLOCATE(   g_dot_h(3, nconf) )
    ALLOCATE(   h_dot_t(3, nconf) )
    ALLOCATE(    g_mod2(3, nconf) )
    ALLOCATE(g_mod2_old(3, nconf) )
    ALLOCATE(      pref(3, nconf) )
    !
    ioWRITE(stdout,"(3x,'\>\^\~',80('-'),'^v^v',40('-'),'=/~/o>',/,4x,a,i4)") "iter ", 0
      timer_CALL t_tkaout%start()
    ! f0 = f_SMA = 1/(A_out) b
    !f = A_diag_f(inv_A_out, qbasis%b, nconf, nat3, nq)
    f = A_diagm1_f(A_out, qbasis%b, nconf, nat3, nq)
    tk = calc_tk_simple(f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
    CALL print_tk(tk, input%sigma, input%T, nconf, "SMA tk")
    ! g0 = Af0 - b
    CALL A_times_f(f, g, A_out, input, qbasis, out_grid, S, fc2, fc3)
!     g = A_diag_f(A_out, f, nconf, nat3, nq)
    g = qbasis%b-g
    !g = g/100._dp
    g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
    !
    h = g
    tk = calc_tk_gf(g, f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
    CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb) - initial")
    CALL print_gradmod2_tk(g_mod2, "TK gradient mod", input%T, S%omega, nconf, nat3, nq)
      timer_CALL t_tkaout%stop()
    !
    DO iter = 1,niter_max
      ioWRITE(stdout,"(3x,'\>\^\~',80('-'),'^v^v',40('-'),'=/~/o>',/,4x,a,i4)") "iter ", iter
      ! t = (A_in+A_out)h
        timer_CALL t_tkain%start()
      CALL A_times_f(h, t, A_out, input, qbasis, out_grid, S, fc2, fc3)
        timer_CALL t_tkain%stop()
      !
        timer_CALL t_tkcg%start()
      ! f_(i+1) = f_i - (g_i.h_i) / (h_i.t_i) h_i
      ! g_(i+1) = g_i - (g_i.h_i) / (h_i.t_i) t_i
      g_dot_h = qbasis_dot(g, h,  nconf, nat3, nq )
      h_dot_t = qbasis_dot(h, t, nconf, nat3, nq )
      ! qbasis_b_over_a check for zero/zero, preventing NaN
      pref = qbasis_a_over_b(g_dot_h, h_dot_t, nconf)
      f = f - qbasis_ax(pref, h, nconf, nat3, nq)
      g = g - qbasis_ax(pref, t, nconf, nat3, nq)
      ! compute gradient explicitly:
      !
      !tk = -\lambda 1/2(f.g-f.b)
      tk = calc_tk_gf(g, f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
      CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb)")
      !
      ! h_(i+1) = -g_(i+1) + (g_(i+1).g_(i+1)) / (g_i.g_i) h_i
      g_mod2_old = g_mod2
      g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
      CALL print_gradmod2_tk(g_mod2, "TK gradient mod", input%T, S%omega, nconf, nat3, nq)
      !pref = g_mod2 / g_mod2_old
      pref = qbasis_a_over_b(g_mod2, g_mod2_old, nconf)
      h = qbasis_ax(pref, h, nconf, nat3, nq) - g
        timer_CALL t_tkcg%stop()
    ENDDO
    !
  END SUBROUTINE TK_CG
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  SUBROUTINE TK_CG_prec(input, out_grid, S, fc2, fc3)
    USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE q_grids,            ONLY : q_grid, q_basis, setup_grid, setup_bz_grid, &
                                   prepare_q_basis, qbasis_dot, qbasis_ax, &
                                   qbasis_a_over_b
    USE variational_tk
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: out_grid
    !
    INTEGER :: ix, nu, iq, it, nu0, iter
    !
    TYPE(q_grid)         :: in_grid ! inner grid is MPI-scattered it is used 
                                    ! for integrating the ph-ph scattering terms and linewidth
    !
    TYPE(q_basis) :: qbasis
    REAL(DP),ALLOCATABLE :: A_out(:,:,:), inv_sqrt_A_out(:,:,:), inv_A_out(:,:,:)
    REAL(DP),ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), h(:,:,:,:), t(:,:,:,:)
    REAL(DP),ALLOCATABLE :: Af(:,:,:,:)
    REAL(DP),ALLOCATABLE :: g_dot_h(:,:), h_dot_t(:,:), &
                            g_mod2(:,:), g_mod2_old(:,:), &
                            pref(:,:)
    INTEGER,PARAMETER :: niter_max = 10000
    REAL(DP) :: tk(3,3,input%nconf)
    !
    ! For code readability:
    INTEGER :: nconf, nat3, nq
    nconf = input%nconf
    nat3  = S%nat3
    nq    = out_grid%nq
    
    ! make the inner grid on top of the outer one
    !CALL setup_grid(input%grid_type, S%bg, out_grid%n(1),out_grid%n(2),out_grid%n(3), in_grid)
    CALL setup_bz_grid(S%bg, out_grid%n(1),out_grid%n(2),out_grid%n(3), in_grid)
    CALL in_grid%scatter()
    
    CALL prepare_q_basis(out_grid, qbasis, nconf, input%T, S, fc2)
    ! Compute A_out diagonal matrix
    ALLOCATE(A_out(nconf, nat3, nq))
    CALL compute_A_out(A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
    ! Compute 1/sqrt(A_out) from A_out
    ALLOCATE(inv_sqrt_A_out(nconf, nat3, nq))
    CALL compute_inv_sqrt_A_out(A_out, inv_sqrt_A_out, nconf, nat3, nq)
    ! Compute 1/A_out from A_out
!     ALLOCATE(inv_A_out(nconf, nat3, nq))
!     CALL compute_inv_A_out(A_out, inv_A_out, nconf, nat3, nq)
    !
    !
    ALLOCATE(f(3, nconf, nat3, nq))
    ALLOCATE(g(3, nconf, nat3, nq))
    ALLOCATE(h(3, nconf, nat3, nq))
    ALLOCATE(t(3, nconf, nat3, nq))
    ALLOCATE(   g_dot_h(3, nconf) )
    ALLOCATE(   h_dot_t(3, nconf) )
    ALLOCATE(    g_mod2(3, nconf) )
    ALLOCATE(g_mod2_old(3, nconf) )
    ALLOCATE(      pref(3, nconf) )
    !
    ioWRITE(stdout,"(3x,'\>\^\~',80('-'),'^v^v',40('-'),'=/~/o>',/,4x,a,i4)") "iter ", 0
        timer_CALL t_tkaout%start()
    ! \tilde{f0} = A_out^(-1/2) b
    ! \tilde{b} = \tilde{f0}
    f = A_diag_f(inv_sqrt_A_out, qbasis%b, nconf, nat3, nq)
    qbasis%b = f
    !
    tk = calc_tk_simple(f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
    CALL print_tk(tk, input%sigma, input%T, nconf, "SMA tk")
    !
    ! f0 = f_SMA = 1/(A_out) b
    ! g0 = Af0 - b
    CALL tilde_A_times_f(f, g, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
!     g = A_diag_f(A_out, f, nconf, nat3, nq)
    g = qbasis%b-g
    !g = g/100._dp
    g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
    !
    h = g
    tk = calc_tk_gf(g, f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
    CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb) - initial")
    CALL print_gradmod2_tk(g_mod2, "TK gradient mod", input%T, S%omega, nconf, nat3, nq)
        timer_CALL t_tkaout%stop()
    !
    DO iter = 1,niter_max
      CALL check_graceful_termination()
      ioWRITE(stdout,"(3x,'\>\^\~',80('-'),'^v^v',40('-'),'=/~/o>',/,4x,a,i4)") "iter ", iter
      ! t = (A_in+A_out)h
        timer_CALL t_tkain%start()
      CALL tilde_A_times_f(h, t, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
        timer_CALL t_tkain%stop()
      !
        timer_CALL t_tkcg%start()
      ! f_(i+1) = f_i - (g_i.h_i) / (h_i.t_i) h_i
      ! g_(i+1) = g_i - (g_i.h_i) / (h_i.t_i) t_i
      g_dot_h = qbasis_dot(g, h,  nconf, nat3, nq )
      h_dot_t = qbasis_dot(h, t, nconf, nat3, nq )
      ! qbasis_a_over_b check for zero/zero, preventing NaN
      pref = qbasis_a_over_b(g_dot_h, h_dot_t, nconf)
      f = f - qbasis_ax(pref, h, nconf, nat3, nq)
      g = g - qbasis_ax(pref, t, nconf, nat3, nq)
      !
      !tk = -\lambda 1/2(f.g-f.b)
      tk = calc_tk_gf(g, f, qbasis%b, input%T, S%omega, nconf, nat3, nq)
      CALL print_tk(tk, input%sigma, input%T, nconf, "TK from 1/2(fg-fb)")
      !
      ! h_(i+1) = -g_(i+1) + (g_(i+1).g_(i+1)) / (g_i.g_i) h_i
      g_mod2_old = g_mod2
      g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
      CALL print_gradmod2_tk(g_mod2, "TK gradient mod", input%T, S%omega, nconf, nat3, nq)
      !pref = g_mod2 / g_mod2_old
      pref = qbasis_a_over_b(g_mod2, g_mod2_old, nconf)
      h = qbasis_ax(pref, h, nconf, nat3, nq) - g
        timer_CALL t_tkcg%stop()
    ENDDO
    !
  END SUBROUTINE TK_CG_prec
  END MODULE thermalk_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!               !               !               !               !               !               !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM thermalk

  USE kinds,            ONLY : DP
  USE thermalk_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : print_citations_linewidth, forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid !, setup_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE nanoclock,        ONLY : init_nanoclock
  !
  USE posix_signal,       ONLY : set_TERMINATE_GRACEFULLY
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid)     :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)       :: S
  TYPE(code_input_type)      :: tkinput
  TYPE(q_grid)               :: out_grid

!   CALL mp_world_start(world_comm)
!   CALL environment_start('TK')
  CALL init_nanoclock()
  CALL start_mpi()
  CALL print_citations_linewidth()
  CALL set_TERMINATE_GRACEFULLY() !print_timers_and_die)

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  CALL READ_INPUT("TK", tkinput, out_grid, S, fc2, fc3)
  !
  IF(TRIM(tkinput%calculation) == "sma") THEN
    !
    CALL TK_SMA(tkinput, out_grid, S, fc2, fc3)
    !
  ELSEIF(TRIM(tkinput%calculation) == "cg") THEN
    !
    CALL TK_CG(tkinput, out_grid, S, fc2, fc3)
    !
  ELSEIF(TRIM(tkinput%calculation) == "cgp") THEN
    !
    CALL TK_CG_prec(tkinput, out_grid, S, fc2, fc3)
    !
  ELSE
    CALL errore("lw", "what else to do?", 1)
  ENDIF
  !
  CALL stop_mpi()
!   CALL environment_end('TK')
!   CALL mp_world_end()
 
END PROGRAM thermalk
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

