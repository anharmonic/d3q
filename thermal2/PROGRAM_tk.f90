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
    USE constants,          ONLY : RY_TO_CMM1, K_BOLTZMANN_RY, tpi
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf, ryvel_si
    USE q_grids,            ONLY : q_grid, setup_grid
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel, mfp_scatter_vel
    USE input_fc,           ONLY : ph_system_info
    USE code_input,         ONLY : code_input_type
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity, velocity_operator
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
    REAL(DP) :: lw_un(S%nat3,2,input%nconf)
    REAL(DP) :: freq(S%nat3)
    REAL(DP) :: bose(S%nat3,input%nconf)
    COMPLEX(DP) :: U(S%nat3,S%nat3)
    !
    REAL(DP) :: tk(3,3,input%nconf)
    REAL(DP) :: vel(3,S%nat3)
    !
    COMPLEX(DP) :: tk_P(3,3,input%nconf)
    COMPLEX(DP) :: tk_C(3,3,input%nconf)
    COMPLEX(DP) :: vel_operator(3,S%nat3,S%nat3)
    COMPLEX(DP) :: prefactor
    REAL(DP) :: Lorentzian, boseplus, boseplus_p, debug_max_imag_k_P, debug_max_imag_k_C
    !
    !for debug
    REAL(DP) :: vel_diag(3,S%nat3)
    !
    REAL(DP) :: pref
    INTEGER  :: iq, it, a, b, nu, id_dir,im, im2, nup, iu
    character(len=29) :: filename
    LOGICAL :: condition_print
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    LOGICAL :: gamma
    !TYPE(order_type) :: order
    !
    character*1 tab
    tab = char(9)
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
          !
          iu = 1000+input%nconf+it
          OPEN(unit=iu, file=TRIM(input%outdir)//"/"//&
                                  "lw_N."//TRIM(input%prefix)//&
                                  "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                  "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
          ioWRITE(iu, *) "# linewidth NORMAL [cm^-1]"
          ioWRITE(iu, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it),&
                                                      "    sigma=", input%sigma(it)
          ioFLUSH(iu)
          !
          iu = 1000+2*input%nconf+it
          OPEN(unit=iu, file=TRIM(input%outdir)//"/"//&
                                  "lw_U."//TRIM(input%prefix)//&
                                  "_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                                  "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out")
          ioWRITE(iu, *) "# linewidth UMKLAPP [cm^-1]"
          ioWRITE(iu, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, "     T=",input%T(it),&
                                                      "    sigma=", input%sigma(it)
          ioFLUSH(iu)

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
        OPEN(unit=5001, file=TRIM(input%outdir)//"/"//&
                                "vel_diag."//TRIM(input%prefix)//".out")
        ioWRITE(5001, *) "# ph group velocity x1, y1, z1, x2, y2, z2, ..."
        ioFLUSH(5001)

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
    tk_P = 0._dp
    tk_C = 0._dp
    !
    IF(input%print_all.and.ionode) THEN
      open(unit=7022,file=TRIM(input%outdir)//"/"//'phonon_velocity_operator.dat',status='replace')
      WRITE(7022,'(A)') '#i_q, im1, im2, f_1 (cm^-1), f_2(cm^-1), |V^x|^2, |V^y|^2, |V^z|^2, units (m/s)^2'
      DO it = 1, input%nconf
        write(filename, '(A22,i3.3,A4)') TRIM(input%outdir)//"/"//'phonon_properties_raw_',it,'.dat'
        open(unit=7023+it, file=filename, status='replace')
        write(7023+it,'(A)') "#  q, mod1, mod2, omega1(Ry),  omega2(Ry),gamma1(Ry),  gamma2(Ry)"
      END DO
    ENDIF    

    timer_CALL t_tksma%start()
    QPOINT_LOOP : &
    DO iq = 1,out_grid%nq
      !ioWRITE(stdout,'(i6,3f15.8)') iq, out_grid%xq(:,iq)
      CALL print_percent_wall(10._dp, 300._dp, iq, out_grid%nq, (iq==1))
      !
      IF (input%intrinsic_scattering) THEN
          timer_CALL t_lwphph%start()
        lw_phph = linewidth_q(out_grid%xq(:,iq), input%nconf, input%T,&
                         sigma_ry, S, in_grid, fc2, fc3, lw_un=lw_un)
        CALL check_negative_lw(lw_phph, S%nat3, input%nconf, "SMA:phph")
          timer_CALL t_lwphph%stop()
      ELSE
        lw_phph = 0._dp
        IF (input%debug_linewidth) THEN
          lw_phph = 15._dp/RY_TO_CMM1
        ENDIF
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
      ! Velocity operator
      !ioWRITE(stdout,"(3x,a)") "velocity_operator start"      
      vel_operator = velocity_operator( S, fc2, out_grid%xq(:,iq))
      !
      !ioWRITE(stdout,"(3x,a)") "velocity_operator end"      
      ! old subroutine, for debug
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
          IF(input%intrinsic_scattering.and.ionode) THEN 
            WRITE(1000+it,'(99e20.10)') lw_phph(:,it)*RY_TO_CMM1
            WRITE(1000+input%nconf+it,'(99e20.10)') lw_un(:,1,it)*RY_TO_CMM1
            WRITE(1000+2*input%nconf+it,'(99e20.10)') lw_un(:,2,it)*RY_TO_CMM1
          ENDIF
          IF(input%isotopic_disorder.and.ionode)    WRITE(2000+it,'(99e20.10)') lw_isotopic(:,it)*RY_TO_CMM1
        ENDDO
        IF(input%casimir_scattering) THEN
          ioWRITE(3000,'(99e20.10)') lw_casimir(:)*RY_TO_CMM1
        ENDIF
        ioWRITE(5000,'(3(99e20.10,3x))') vel(:,:)
        ioWRITE(5001,'(3(99e20.10,3x))') vel_diag(:,:)
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
        IF(input%print_all.and.ionode) THEN
          ! print velocity operator
          do im=1,S%nat3
             do im2=1,S%nat3 
                condition_print=.false.      
                if (input%workaround_print_v) then     
                  if (iq==1) then
                    if (im>3 .and. im2>3) then
                      condition_print=.true.
                    end if
                  else
                      condition_print=.true.
                  end if
                else
                  condition_print=((lw(im,it)>0.0d0).and.(lw(im2,it)>0.0d0))
                end if            
               if ( condition_print .and. (im<=im2)) then
                if (it==1) then
                  WRITE(7022 ,'(I6,I4,I4,2ES11.4, 3ES14.7)')       &
                   iq,im,im2,                                  & 
                   freq(im)*RY_TO_CMM1,                   &
                   freq(im2)*RY_TO_CMM1,                  &
                  ((real(vel_operator(1,im,im2)))**2+     &
                  (aimag(vel_operator(1,im,im2))**2 ) )   &
                                         *((tpi*ryvel_si)**2 ), &
                  ((real(vel_operator(2,im,im2)))**2+     &
                  (aimag(vel_operator(2,im,im2))**2 ) )   &
                                         *((tpi*ryvel_si)**2 ), &
                  ((real(vel_operator(3,im,im2)))**2+     &
                  (aimag(vel_operator(3,im,im2))**2 ) )   &
                                         *((tpi*ryvel_si)**2 )
               endif
               ! print phonon phonon_properties_raw
               write(7023+it,'(I6,I4,I4,A1,ES17.10,A1,ES17.10,A1,ES17.10,A1,ES17.10)') &
                  iq, im,im2,tab,freq(im),tab,freq(im2),tab,lw(im,it),tab,lw(im2,it)
              end if
           end do
         end do
        END IF
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
          ! end of the safety checks loop 1 
          ! coherence term
          boseplus=bose(nu,it)*(1.0+bose(nu,it))
          DO nup = 1, S%nat3
            boseplus_p=bose(nup,it)*(1.0+bose(nup,it))
            !
            prefactor=0.25d0*(freq(nu)+freq(nup))*                                     &
                          (boseplus*freq(nu)+boseplus_p*freq(nup))*                    &
                          out_grid%w(iq)/(input%T(it)**2 *S%Omega*K_BOLTZMANN_RY )
            !
            Lorentzian=  (0.5d0*(lw(nu,it)+lw(nup,it)))/            &
                         (  (freq(nu)-freq(nup))**2 + (0.5d0*(lw(nu,it)+lw(nup,it)))**2 )
            !
            IF ((lw(nu,it)>0.0d0) .and. (lw(nup,it)>0.0d0)) THEN
              DO a = 1,3
                DO b = 1,3
                  IF ( nu == nup) THEN
                    tk_P(a,b,it) = tk_P(a,b,it) + prefactor*Lorentzian*    &
                        (vel_operator(a,nu,nup)*vel_operator(b,nup,nu))
                  ELSE
                    tk_C(a,b,it) = tk_C(a,b,it) + prefactor*Lorentzian*    &
                        (vel_operator(a,nu,nup)*vel_operator(b,nup,nu))
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO ! end second mode loop
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
      CLOSE(5001) !diagonal elements velocity operator
      CLOSE(6000) ! freq
      CLOSE(7000) ! q-points
    ENDIF
    IF(ionode.and.input%print_all)THEN
      CLOSE(7022)
      DO it = 1, input%nconf
        close(7023+it)
      END DO
    ENDIF
    !
    IF(input%intrinsic_scattering .or. input%debug_linewidth ) THEN
      ! Write to disk
      !      
      IF(ionode) OPEN(unit=10001, file=TRIM(input%outdir)//"/"//&
                                       TRIM(input%prefix)//"_P_sma."//"out")
      IF(ionode) OPEN(unit=10002, file=TRIM(input%outdir)//"/"//&
                                       TRIM(input%prefix)//"_C."//"out")  
      IF(ionode) OPEN(unit=10003, file=TRIM(input%outdir)//"/"//&
                                       TRIM(input%prefix)//"_TOT."//"out")                                             

      ioWRITE(10001,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                          "    K_P^x            K_P^y            K_P^z             ",&
                          "    K_P^xy           K_P^xz           K_P^yz            ",&
                          "    K_P^yx           K_P^zx           K_P^zy       "

      ioWRITE(10002,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                          "    K_C^x            K_C^y            K_C^z             ",&
                          "    K_C^xy           K_C^xz           K_C^yz            ",&
                          "    K_C^yx           K_C^zx           K_C^zy       "

      ioWRITE(10003,'(4a)') "#conf  sigma[cmm1]   T[K]  ",&
                          "    K_TOT^x            K_TOT^y            K_TOT^z             ",&
                          "    K_TOT^xy           K_TOT^xz           K_TOT^yz            ",&
                          "    K_TOT^yx           K_TOT^zx           K_TOT^zy       "
      tk = tk*RY_TO_WATTMM1KM1

      tk_P = tk_P*RY_TO_WATTMM1KM1
      tk_C = tk_C*RY_TO_WATTMM1KM1

      debug_max_imag_k_P=maxval(abs(aimag(tk_P)))
      debug_max_imag_k_C=maxval(abs(aimag(tk_C)))

      DO it = 1,input%nconf
        !
        ioWRITE(10001,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
        REAL(tk_P(1,1,it)),REAL(tk_P(2,2,it)),REAL(tk_P(3,3,it)), &
        REAL(tk_P(1,2,it)),REAL(tk_P(1,3,it)),REAL(tk_P(2,3,it)), &
        REAL(tk_P(2,1,it)),REAL(tk_P(3,1,it)),REAL(tk_P(3,2,it))
        !
        ioWRITE(10002,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
        REAL(tk_C(1,1,it)),REAL(tk_C(2,2,it)),REAL(tk_C(3,3,it)), &
        REAL(tk_C(1,2,it)),REAL(tk_C(1,3,it)),REAL(tk_C(2,3,it)), &
        REAL(tk_C(2,1,it)),REAL(tk_C(3,1,it)),REAL(tk_C(3,2,it))
        !
        ioWRITE(10003,"(i3,2f12.6,3(3e15.6,5x))") it, input%sigma(it), input%T(it), &
        REAL(tk_P(1,1,it)+tk_C(1,1,it)),REAL(tk_P(2,2,it)+tk_C(2,2,it)),REAL(tk_P(3,3,it)+tk_C(3,3,it)), &
        REAL(tk_P(1,2,it)+tk_C(1,2,it)),REAL(tk_P(1,3,it)+tk_C(1,3,it)),REAL(tk_P(2,3,it)+tk_C(2,3,it)), &
        REAL(tk_P(2,1,it)+tk_C(2,1,it)),REAL(tk_P(3,1,it)+tk_C(3,1,it)),REAL(tk_P(3,2,it)+tk_C(3,2,it))
      ENDDO
      IF(ionode) CLOSE(10001)
      IF(ionode) CLOSE(10002)
      IF(ionode) CLOSE(10003)

      !
      ! Write to screen
      ioWRITE(stdout,"(3x,a,/,3x,a)") "************", "SMA thermal conductivity, stored to file:"
      ioWRITE(stdout,'(5x,a)') TRIM(input%outdir)//"/"//TRIM(input%prefix)//"."//"out"
      !
      ioWRITE(stdout,"(3x,a)") "Thermal conductivity  (conf, sigma, T, K_x, K_y, K_z):"
      DO it = 1,input%nconf
        ioWRITE(stdout,"(a7,i3,2f12.6,3e16.8)") 'K_P   =', it, input%sigma(it), input%T(it),&
                                        REAL(tk_P(1,1,it)), REAL(tk_P(2,2,it)), REAL(tk_P(3,3,it))
        ioWRITE(stdout,"(a7,i3,2f12.6,3e16.8)") 'K_C   =', it, input%sigma(it), input%T(it),&
                                        REAL(tk_C(1,1,it)), REAL(tk_C(2,2,it)), REAL(tk_C(3,3,it))
        ioWRITE(stdout,"(a7,i3,2f12.6,3e16.8)") 'K_TOT =', it, input%sigma(it), input%T(it),&
                           REAL(tk_P(1,1,it)+tk_C(1,1,it)), REAL(tk_P(2,2,it)+tk_C(2,2,it)),&
                           REAL(tk_P(3,3,it)+tk_C(3,3,it))                                                                
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
      !CALL compute_A_out_symq(A_out, input, qbasis, out_grid, in_grid, S, fc2, fc3)
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
  CALL remove_stack_limit()
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
    !CALL TK_SMA(tkinput, out_grid, S, fc2, fc3)    
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

