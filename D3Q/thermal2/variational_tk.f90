!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE variational_tk
  USE kinds, ONLY : DP
  ! 
  CONTAINS
  ! This subroutine computes the SMA thermal conducivity, it is mainly just a driver
  ! that uses other subroutines to compute th intrinsic, isotopic and casimir linewidths,
  ! than it sums everything up and takes care of input/output.
  SUBROUTINE gen_a_out(A_out, input, qbasis, S, fc2, fc3)
    USE constants,          ONLY : RY_TO_CMM1
    USE linewidth,          ONLY : linewidth_q
    USE q_grids,            ONLY : q_basis
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_vel
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)     :: qbasis
    REAL(DP),INTENT(out) :: A_out(3,input%nconf,S%nat3, qbasis%grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    REAL(DP) :: lw(S%nat3, input%nconf)
    REAL(DP) :: lw_isotopic(S%nat3, input%nconf)
    REAL(DP) :: lw_phph(S%nat3, input%nconf)
    REAL(DP) :: lw_casimir(S%nat3)
    REAL(DP) :: bose(S%nat3, input%nconf)
    !
    REAL(DP) :: vel(3,S%nat3)
    INTEGER  :: iq, it, ix, nu
    !
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    lw_isotopic = 0._dp
    lw_casimir  = 0._dp
    
    QPOINT_LOOP : &
    DO iq = 1,qbasis%grid%nq
      ! bang!
      lw_phph = linewidth_q(qbasis%grid%xq(:,iq), input%nconf, input%T,sigma_ry, S, qbasis%grid, fc2, fc3)
      !
      ! Compute contribution of isotopic disorder
      IF(input%isotopic_disorder) &
        lw_isotopic = isotopic_linewidth_q(qbasis%grid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, qbasis%grid, fc2)
      !
      vel = qbasis%c(:,:,iq)
      !
      ! Compute potentially anisotropic Casimir linewidth
      IF(input%casimir_scattering) &
        lw_casimir = casimir_linewidth_vel( qbasis%c(:,:,iq), input%casimir_length, input%casimir_dir, S%nat3)
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
        bose(:,it) = qbasis%be(:,it,iq)
        !
        MODE_LOOP : &
        DO nu = 1, S%nat3
          ! Check if we have zero linewidth and non-zero velocity it is a problem
          ! lw can be NaN when T=0 and xq=0, check for lw>0 insteand, because NaN/=0 is true
          IF(lw(nu,it)<0._dp)THEN ! true for NaN
            WRITE(*,"(3x,a,e12.4,3i6)") "WARNING! Negative lw (idx q, mode, conf):", lw(nu,it), iq, nu, it 
            lw(nu,it) = - lw(nu,it)
          ENDIF
          IF(.not. lw(nu,it)>0._dp)THEN ! false for NaN
            IF(ANY(ABS(qbasis%c(:,nu,iq))>eps_vel ))THEN
              WRITE(*,'(3i6,1e20.10,5x,3e20.10)') iq, nu, it, lw(nu,it), qbasis%c(:,nu,iq)
              CALL errore("TK_SMA", "cannot threat this case", 1)
            ELSE
              WRITE(*,"(3x,a,3i6)") "skip (iq,nu,it):", iq, nu, it
              CYCLE MODE_LOOP 
            ENDIF
          ENDIF
          !
          DO ix = 1,3
            A_out(ix,it,nu,iq) = bose(nu,it)*(1+bose(nu,it)) /lw(nu,it)
          ENDDO
        ENDDO MODE_LOOP 
      ENDDO CONF_LOOP 
    ENDDO QPOINT_LOOP
    !
  END SUBROUTINE gen_a_out
  !
END MODULE variational_tk
