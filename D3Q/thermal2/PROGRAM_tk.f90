!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE thermalk_program
  !
#include "mpi_thermal.h"
  USE kinds,       ONLY : DP
  USE mpi_thermal, ONLY : ionode
  USE posix_signal,ONLY : check_graceful_termination
  USE timers
  !
  CONTAINS
  !
  SUBROUTINE check_negative_lw(lw, nat3, nconf, name)
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: lw(nat3, nconf)
    INTEGER,INTENT(in)  :: nat3, nconf
    CHARACTER(len=*),INTENT(in) :: name
    !
    INTEGER :: i,j
    IF(ANY(lw<0._dp))THEN
      DO i = 1,nconf
        DO j = 1,nat3
          IF(lw(j,i) < 0._dp) THEN
            WRITE(*,*) i, j, lw(j,i)
          ENDIF
        ENDDO
      ENDDO
      CALL errore(name, "negative linewidth", 1)
    ENDIF
    RETURN
  END SUBROUTINE
  !
  ! This subroutine computes the SMA thermal conducivity, it is mainly just a driver
  ! that uses other subroutines to compute th intrinsic, isotopic and casimir linewidths,
  ! than it sums everything up and takes care of input/output.
  !
  ! This subroutine is obsoleted by the first iteration of the variational method,
  ! but we keep it for didactical purposes. This is also faster, but it does not respect
  ! the detailed balance exactly (negative tk is possible at low temperature)
  !
  ! NOTE: all the *linewidth* functions return the HALF width half maximum, here 
  ! we multiply by 2 in order to get the full width.
  SUBROUTINE TK_SMA(input, out_grid, S, fc2, fc3)
    USE linewidth,          ONLY : linewidth_q
    USE constants,          ONLY : RY_TO_CMM1, K_BOLTZMANN_RY
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf
    USE q_grids,            ONLY : q_grid, setup_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel, mfp_scatter_vel
    USE input_fc,           ONLY : ph_system_info
    USE code_input,         ONLY : code_input_type
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity
    !USE overlap,            ONLY : order_type
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
    REAL(DP) :: pref
    INTEGER  :: iq, it, a, b, nu
    !
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    LOGICAL :: gamma
    !TYPE(order_type) :: order
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    ! the inner grid (in_grid) is scatterd over MPI
    ioWRITE(*,*) "--> Setting up inner grid"
    CALL setup_grid(input%grid_type_in, S%bg, input%nk_in(1), input%nk_in(2), input%nk_in(3),&
                    in_grid, scatter=.true., xq0=input%xk0_in)
    !
    ! Open files to store the linewidth
    IF(ionode.and.input%store_lw)THEN
      DO it = 1,input%nconf
        IF (input%intrinsic_scattering) THEN
          OPEN(unit=1000+it, file=TRIM(input%outdir)//"/"//&
                                  "lw."//TRIM(input%prefix)//&
                                  "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                  "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
          ioWRITE(1000+it, *) "# linewidth [cm^-1]"
          ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it),&
                                                      "    sigma=", input%sigma(it)
          ioFLUSH(1000+it)
        ENDIF
        !
        !
        IF(input%isotopic_disorder) THEN
          OPEN(unit=2000+it, file=TRIM(input%outdir)//"/"//&
                                  "lwiso."//TRIM(input%prefix)//&
                                      "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                      "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
          ioWRITE(2000+it, *) "# isotopic linewidth [cm^-1]"
          ioWRITE(2000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it),&
                                                      "    sigma=", input%sigma(it)
          ioFLUSH(2000+it)
        ENDIF
      ENDDO
      IF(input%casimir_scattering) THEN
        OPEN(unit=3000, file=TRIM(input%outdir)//"/"//&
                                "lwcas."//TRIM(input%prefix)//".out")
        ioWRITE(3000, *) "# casimir linewidth [cm^-1]"
        ioFLUSH(3000)
      ENDIF
      ! Velocity:
        OPEN(unit=5000, file=TRIM(input%outdir)//"/"//&
                                "vel."//TRIM(input%prefix)//".out")
        ioWRITE(5000, *) "# ph group velocity x1, y1, z1, x2, y2, z2, ..."
        ioFLUSH(5000)
        OPEN(unit=6000, file=TRIM(input%outdir)//"/"//&
                                "freq."//TRIM(input%prefix)//".out")
        ioWRITE(6000, *) "# phonon frequencies"
        ioFLUSH(6000)
        OPEN(unit=7000, file=TRIM(input%outdir)//"/"//&
                                "q."//TRIM(input%prefix)//".out")
        ioWRITE(7000, *) "# qpoint [2pi/alat], weight"
        ioFLUSH(7000)
    ENDIF
    !
    tk = 0._dp
    !
      timer_CALL t_tksma%start()
    QPOINT_LOOP : &
    DO iq = 1,out_grid%nq
      !ioWRITE(stdout,'(i6,3f15.8)') iq, out_grid%xq(:,iq)
      CALL print_percent_wall(10._dp, 300._dp, iq, out_grid%nq, (iq==1))
      !
      IF (input%intrinsic_scattering) THEN
          timer_CALL t_lwphph%start()
        lw_phph = linewidth_q(out_grid%xq(:,iq), input%nconf, input%T,&
                         sigma_ry, S, in_grid, fc2, fc3)
        CALL check_negative_lw(lw_phph, S%nat3, input%nconf, "SMA:phph")
          timer_CALL t_lwphph%stop()
      ELSE
        lw_phph = 0._dp
      ENDIF
      !
      ! Compute contribution of isotopic disorder
      IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw_isotopic = isotopic_linewidth_q(out_grid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, in_grid, fc2)
        CALL check_negative_lw(lw_isotopic, S%nat3, input%nconf, "SMA:isotopic")
          timer_CALL t_lwisot%stop()
      ELSE
        lw_isotopic = 0._dp
      ENDIF
      !
      !
        timer_CALL t_velcty%start()
      ! Velocity
      vel = velocity(S, fc2, out_grid%xq(:,iq))
        timer_CALL t_velcty%stop()
      !
      CALL  freq_phq_safe(out_grid%xq(:,iq), S, fc2, freq, U)
      !
      !CALL order%set(S%nat3, freq, U, keep_old_e=.true.)
      !
      ! Compute anisotropic Casimir linewidth
      IF(input%casimir_scattering) THEN
          timer_CALL t_lwcasi%start()
        lw_casimir = casimir_linewidth_vel(vel, input%sample_length, input%sample_dir, S%nat3)
          timer_CALL t_lwcasi%stop()
      ELSE
        lw_casimir = 0._dp
      ENDIF
      !
      IF(input%store_lw)THEN
        timer_CALL t_lwinout%start()
        DO it = 1, input%nconf
          IF(input%intrinsic_scattering.and.ionode) WRITE(1000+it,'(99e20.10)') lw_phph(:,it)*RY_TO_CMM1
          IF(input%isotopic_disorder.and.ionode)    WRITE(2000+it,'(99e20.10)') lw_isotopic(:,it)*RY_TO_CMM1
        ENDDO
        IF(input%casimir_scattering) THEN
          ioWRITE(3000,'(99e20.10)') lw_casimir(:)*RY_TO_CMM1
        ENDIF
        ioWRITE(5000,'(3(99e20.10,3x))') vel(:,:)
        ioWRITE(6000,'(99e20.10)') freq(:)*RY_TO_CMM1
        ioWRITE(7000,'(4e20.10)') out_grid%xq(:,iq), out_grid%w(iq)
        timer_CALL t_lwinout%stop()
      ENDIF
      !
        timer_CALL t_tksum%start()
      ! Casimir linewidth is temperature/smearing-independent, sum it to all configurations
      ! Also multiply by 2 in order to get the FULL width from the HALF width
      ! thanks to Francesco Macheda for reporting this
      DO it = 1,input%nconf
        lw(:,it) = 2*(lw_phph(:,it) + lw_isotopic(:,it) + lw_casimir)
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
          ! Check: if we have zero linewidth and non-zero velocity it is a problem
          ! lw can be NaN when T=0 and xq=0, check for lw>0 instead, because NaN/=0 is true
          IF(lw(nu,it)<0._dp)THEN ! true for NaN
            WRITE(stdout,"(3x,a,e12.4,3i6)") &
               "WARNING! Negative lw (idx q, mode, conf):", lw(nu,it), iq, nu, it
            lw(nu,it) = - lw(nu,it)
          ENDIF
          IF( (.not. lw(nu,it)>0._dp) .and. (input%intrinsic_scattering) )THEN ! false for NaN
            gamma=(ALL(ABS(out_grid%xq(:,iq))<1.d-12)) 
            IF((ANY(ABS(vel(:,nu))>eps_vel )) .and.(.not. gamma))THEN
              WRITE(stdout,'(3i6,1e20.10,5x,3e20.10)') iq, nu, it, lw(nu,it), vel(:,nu)
              CALL errore("TK_SMA", "Zero or NaN linewidth, with non-zero velocity", 1)
            ELSE
              !ioWRITE(stdout,"(3x,a,3i6)") "skip (iq,nu,it):", iq, nu, it
              CYCLE MODE_LOOP
            ENDIF
          ENDIF
          !
          IF(input%mfp_cutoff)THEN
            ! Note that lw is already multiplied by 2 here
            IF( mfp_scatter_vel(vel(:,nu), lw(nu,it), &
                         input%sample_length, input%sample_dir) ) &
               CYCLE MODE_LOOP
          ENDIF
          !
          pref = freq(nu)**2 *bose(nu,it)*(1+bose(nu,it))&
                             /(input%T(it)**2 *S%Omega*K_BOLTZMANN_RY )&
                             *out_grid%w(iq) /lw(nu,it)
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

    IF(ionode.and.input%store_lw)THEN
      DO it = 1, input%nconf
        IF(input%intrinsic_scattering) CLOSE(1000+it)
        IF(input%isotopic_disorder)    CLOSE(2000+it)
      ENDDO
      IF(input%casimir_scattering) CLOSE(3000)
      CLOSE(5000) !velocity
      CLOSE(6000) ! freq
      CLOSE(7000) ! q-points
    ENDIF
    !
    IF(input%intrinsic_scattering) THEN
      ! Write to disk
      IF(ionode) OPEN(unit=10000, file=TRIM(input%outdir)//"/"//&
                                       TRIM(input%prefix)//"."//"out")
      ioWRITE(10000,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                          "    K_x            K_y            K_z             ",&
                          "    K_xy           K_xz           K_yz            ",&
                          "    K_yx           K_zx           K_zy       "
      tk = tk*RY_TO_WATTMM1KM1
      DO it = 1,input%nconf
        ioWRITE(10000,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
        tk(1,1,it),tk(2,2,it),tk(3,3,it), &
        tk(1,2,it),tk(1,3,it),tk(2,3,it), &
        tk(2,1,it),tk(3,1,it),tk(3,2,it)
      ENDDO
      IF(ionode) CLOSE(10000)
      !
      ! Write to screen
      ioWRITE(stdout,"(3x,a,/,3x,a)") "************", "SMA thermal conductivity, stored to file:"
      ioWRITE(stdout,'(5x,a)') TRIM(input%outdir)//"/"//TRIM(input%prefix)//"."//"out"
      ioWRITE(stdout,"(3x,a)") "Diagonal components (conf, sigma, T, K_x, K_y, K_z):"
      DO it = 1,input%nconf
        ioWRITE(stdout,"(i3,2f12.6,3e16.8)")  it, input%sigma(it), input%T(it),&
                                        tk(1,1,it), tk(2,2,it), tk(3,3,it)
      ENDDO
    ELSE
      ioWRITE(stdout,"(3x,a,/,3x,a)") "************", "SMA thermal conductivity not computed (intrinsic_scattering is false)"
    ENDIF
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
    CALL t_merged%print()
#endif
    !
    !
  END SUBROUTINE TK_SMA
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  SUBROUTINE TK_CG_prec(input, out_grid, S, fc2, fc3)
    USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, freq_phq
    USE linewidth,          ONLY : linewidth_q
    !USE constants,          ONLY : RY_TO_CMM1
     USE more_constants,     ONLY : RY_TO_WATTMM1KM1!, write_conf
    USE fc3_interpolate,    ONLY : forceconst3
    !USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    !USE casimir_linewidth,  ONLY : casimir_linewidth_q
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
    INTEGER :: ix, nu, iq, it, nu0, iter, iter0
    !
    TYPE(q_grid)         :: in_grid ! inner grid is MPI-scattered it is used
                                    ! for integrating the ph-ph scattering terms and linewidth
    !
    TYPE(q_basis) :: qbasis
    REAL(DP),ALLOCATABLE :: A_out(:,:,:), inv_sqrt_A_out(:,:,:), inv_A_out(:,:,:)
    REAL(DP),ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), &
                            h(:,:,:,:), t(:,:,:,:)!, j(:,:,:,:)
    REAL(DP),ALLOCATABLE :: Af(:,:,:,:)
    REAL(DP),ALLOCATABLE :: g_dot_h(:,:), h_dot_t(:,:), &
                            g_mod2(:,:), g_mod2_old(:,:), &
                            pref(:,:)
    REAL(DP) :: tk(3,3,input%nconf), delta_tk(3,3,input%nconf), &
                tk_old(3,3,input%nconf)
    LOGICAL :: conv
!     CHARACTER(len=6) :: what
!     CHARACTER(len=256) :: filename
    !
    ! For code readability:
    LOGICAL :: restart_ok
    INTEGER :: nconf, nat3, nq
    nconf = input%nconf
    nat3  = S%nat3
    nq    = out_grid%nq

        timer_CALL t_tkprec%start()
    ! make the inner grid on top of the outer one, they are actually identical
    ! but the inner one is MPI-scattered, so we need two separate objects to store them
    ioWRITE(*,*) "--> Setting up inner grid from outer grid"
    ! Grid type here is always "simple" because we want to reuse the eventual random
    ! shift from the outer grid (using two different shifts would not work)
    CALL setup_grid("simple", S%bg, out_grid%n(1),out_grid%n(2),out_grid%n(3), &
                    in_grid, xq0=out_grid%xq0, scatter=.true.)
    CALL prepare_q_basis(out_grid, qbasis, nconf, input%T, S, fc2)

    !
    ALLOCATE(A_out(nconf, nat3, nq))
    ALLOCATE(inv_sqrt_A_out(nconf, nat3, nq))
    ALLOCATE(f(3, nconf, nat3, nq))
    ALLOCATE(g(3, nconf, nat3, nq))
    ALLOCATE(h(3, nconf, nat3, nq))
    ALLOCATE(t(3, nconf, nat3, nq))
    ALLOCATE(   g_dot_h(3, nconf) )
    ALLOCATE(   h_dot_t(3, nconf) )
    ALLOCATE(    g_mod2(3, nconf) )
    ALLOCATE(g_mod2_old(3, nconf) )
    ALLOCATE(      pref(3, nconf) )
      timer_CALL t_tkprec%stop()
    !
    ! Try if restart file can be read (will return false if input%restart is false):
      timer_CALL t_restart%start()
    restart_ok = read_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq, iter0, tk)
      timer_CALL t_restart%stop()
    !
    IF(.not. restart_ok) THEN
      ! Compute A_out diagonal matrix
    ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
            "Computing A_out"
        timer_CALL t_tkaout%start()
      CALL compute_A_out(A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
          timer_CALL t_tkaout%stop()
      !
      ! Save A_out, so we do not have to repeat it
        timer_CALL t_restart%start()
      CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq, iter0, tk)
        timer_CALL t_restart%stop()
      !
    ENDIF
    !
    ! Compute A_out^(-1/2) and A_out^-1 from A_out
          timer_CALL t_tkprec%start()
    CALL compute_inv_sqrt_A_out(A_out, inv_sqrt_A_out, nconf, nat3, nq)
          timer_CALL t_tkprec%stop()
    ! Go to the reduced variables:
    ! \tilde{b} = A_out^(-1/2) b
    qbasis%b = A_diag_f(inv_sqrt_A_out, qbasis%b, nconf, nat3, nq)
    !
    !
    ! Bootstrap the CG algorithm, not needed when restarting except if restarting from just after A_out
    IF(iter0 == 0) THEN
      ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
              "iter ", 0
      ! \tilde{f0} = A_out^(1/2) f0 = A_out^(-1/2) b = \tilde{b}
      f = qbasis%b
      ! "tilde" label dropped from here on
      !
      ! Compute SMA tk = lambda f0.b
          timer_CALL t_tkprec%start()
      g = 0._dp
      !
      tk = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
      ! Open files to output the SMA TK values separately
      CALL print_tk(input, tk, "SMA TK", 10000, -1)
          timer_CALL t_tkprec%stop()
      !
      ! Compute g0 = Af0 - b
        timer_CALL t_tkain%start()
      CALL tilde_A_times_f(f, g, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
        timer_CALL t_tkain%stop()
      !
        timer_CALL t_tkprec%start()
      g =  g - qbasis%b
      !
      ! Trivially set h0 = -g0
      h = -g
      !
      !
      tk = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
      CALL print_tk(input, tk, "Initial TK from 1/2(fg-fb)", 10000, 0)
        timer_CALL t_tkprec%stop()
      !
        timer_CALL t_restart%start()
      iter0 = 1
      CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq, iter0, tk)
        timer_CALL t_restart%stop()
    ENDIF
    !
    ! Compute initial variational tk0 =-\lambda (f0.g0-f0.b)
    !ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") "iter ", 0
    !
    g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
    !
    CGP_ITERATION : &
    DO iter = iter0,input%niter_max
      CALL check_graceful_termination()
      ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") "iter ", iter
      ! t = Ah = [A_out^(-1/2) (1+A_in) A_out^(-1/2)] h
        timer_CALL t_tkain%start()
      CALL tilde_A_times_f(h, t, inv_sqrt_A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
        timer_CALL t_tkain%stop()
      !
        timer_CALL t_tkcg%start()
      ! Do a CG step:
      ! f_(i+1) = f_i - (g_i.h_i / h_i.t_i) h_i
      ! g_(i+1) = g_i - (g_i.h_i / h_i.t_i) t_i
      g_dot_h = qbasis_dot(g, h,  nconf, nat3, nq )   ! g.h
      h_dot_t = qbasis_dot(h, t, nconf, nat3, nq )    ! h.t
      pref = qbasis_a_over_b(g_dot_h, h_dot_t, nconf) ! g.h/h.t
      !f_old = f
      f = f - qbasis_ax(pref, h, nconf, nat3, nq)

      g = g - qbasis_ax(pref, t, nconf, nat3, nq)
      !
      !In case you want to compute explicitly the gradient (i.e. for testing):
!       ALLOCATE(j(3, nconf, nat3, nq))
!       CALL tilde_A_times_f(f, j, inv_sqrt_A_out, input, qbasis, &
!                            out_grid, in_grid, S, fc2, fc3)
!       j = j-qbasis%b
      !
      tk_old = tk
      !tk = -\lambda (f.g-f.b)
      tk = calc_tk_gf(g, f, qbasis%b, input%T, out_grid%w, S%omega, nconf, nat3, nq)
      CALL print_tk(input, tk, "TK from 1/2(fg-fb)", 10000, iter)
      ! also compute the variation of tk and check for convergence
      WHERE(tk/=0._dp); delta_tk = (tk-tk_old)/ABS(tk)
      ELSEWHERE;        delta_tk = 0._dp
      END WHERE
      CALL print_deltatk(delta_tk, input%sigma, input%T, nconf, "Delta TK (relative)")
      conv = check_conv_tk(input%thr_tk, nconf, delta_tk)
      !
      IF(conv) THEN
        ioWRITE(stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
                      "Convergence achieved"
        CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq, iter+1, tk)
        EXIT CGP_ITERATION
      ENDIF
      !
      ! h_(i+1) = -g_(i+1) + (g_(i+1).g_(i+1)) / (g_i.g_i) h_i
      g_mod2_old = g_mod2
      g_mod2 = qbasis_dot(g, g, nconf, nat3, nq )
      !
      ! Compute the new conjugate gradient:
      ! h_(i+1) = (g_(i+1).g_(i+1) / g_i.g_i) h_i - g_(i+1)
      pref = qbasis_a_over_b(g_mod2, g_mod2_old, nconf)
      h = -g + qbasis_ax(pref, h, nconf, nat3, nq)
        timer_CALL t_tkcg%stop()
      !
      ! Store to file for restart
        timer_CALL t_restart%start()
      CALL save_cg_step(input, S, A_out, f, g, h, nconf, nat3, nq, iter+1, tk)
        timer_CALL t_restart%stop()
      !
    ENDDO &
    CGP_ITERATION
    !
    IF(iter>input%niter_max) THEN
      ioWRITE (stdout,"(3x,'\>\^\~',40('-'),'^v^v',20('-'),'=/~/o>',/,4x,a,i4)") &
                    "Maximum number of iterations reached."
    ENDIF
    !
#ifdef timer_CALL
      ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
      CALL print_timers_header()
      CALL t_tkprec%print()
      CALL t_tkaout%print()
      CALL t_tkain%print()
      CALL t_tkcg%print()
      CALL t_restart%print()
      ioWRITE(*,'(a)') "*** * High level components:"
      CALL t_tktld%print()
      CALL t_lwisot%print()
      CALL t_lwcasi%print()
      CALL t_lwphph%print()
      CALL t_lwchk%print()
      ioWRITE(*,'(a)') "*** * Low level subroutines: "
      CALL t_freq%print()
      CALL t_bose%print()
      CALL t_sum%print()
      CALL t_fc3int%print()
      CALL t_fc3dint%print()
      CALL t_fc3m2%print()
      CALL t_fc3rot%print()
      CALL t_mpicom%print()
      CALL t_merged%print()
#endif
    END SUBROUTINE TK_CG_prec
  END MODULE thermalk_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM thermalk

  USE kinds,            ONLY : DP
  USE thermalk_program
!   USE environment,      ONLY : environment_start, environment_end
!   USE mp_world,         ONLY : mp_world_start, mp_world_end, world_comm
  USE input_fc,         ONLY : forceconst2_grid, ph_system_info
  USE q_grids,          ONLY : q_grid
  USE fc3_interpolate,  ONLY : forceconst3
  USE code_input,       ONLY : READ_INPUT, code_input_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE nanoclock,        ONLY : init_nanoclock
  USE more_constants,   ONLY : print_citations_linewidth
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
  ELSEIF(TRIM(tkinput%calculation) == "cgp" .or. TRIM(tkinput%calculation) == "exact") THEN
    !
    CALL TK_CG_prec(tkinput, out_grid, S, fc2, fc3)
    !
  ELSE
    CALL errore("tk", "Unknown calculation type: "//TRIM(tkinput%calculation), 1)
  ENDIF
  !
  IF(ionode) CALL print_citations_linewidth()
  CALL stop_mpi()

END PROGRAM thermalk
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

