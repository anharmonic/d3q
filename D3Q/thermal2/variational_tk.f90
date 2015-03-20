!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE variational_tk
  USE kinds, ONLY : DP
  ! 
  ! <<^V^\\=========================================//-//-//========//O\\//
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute the diagonal matrix A_out
  SUBROUTINE gen_A_out(A_out, input, qbasis, S, fc2, fc3)
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
    REAL(DP),INTENT(out) :: A_out(input%nconf,S%nat3, qbasis%grid%nq)
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
          !DO ix = 1,3
          A_out(it,nu,iq) = bose(nu,it)*(1+bose(nu,it)) * lw(nu,it)
          !ENDDO
        ENDDO MODE_LOOP 
      ENDDO CONF_LOOP 
    ENDDO QPOINT_LOOP
    !
  END SUBROUTINE gen_A_out
  ! \/o\________\\\_________________________________________/^>
  ! Apply the A_in matrix to a vector f
  SUBROUTINE A_times_f(f, Af, input, basis, S, fc2, fc3)
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
    TYPE(q_basis),INTENT(in)     :: basis
    REAL(DP),INTENT(in)  ::  f(input%nconf,S%nat3, basis%grid%nq)
    REAL(DP),INTENT(out) :: Af(3,input%nconf,S%nat3, basis%grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    REAL(DP) :: bose(S%nat3, input%nconf)
    !
    !REAL(DP) :: vel(3,S%nat3)
    REAL(DP) :: aux(input%nconf,S%nat3)
    INTEGER  :: iq, it, ix, nu
    !
    REAL(DP),PARAMETER :: eps_vel = 1.e-12_dp
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    QPOINT_OUTER_LOOP : & ! the inner loop in in A_times_f_q
    DO iq = 1,basis%grid%nq
      ! bang!
       aux = A_times_f_q(f, iq, input%nconf, input%T, sigma_ry, S, basis, fc2, fc3)
       DO nu = 1,S%nat3
       DO it = 1,input%nconf
       DO ix = 1,3
        Af(ix,it,nu,iq) = aux(it,nu)
       ENDDO
       ENDDO
       ENDDO
      !
!       ! Compute contribution of isotopic disorder
!       IF(input%isotopic_disorder) &
!         lw_isotopic = isotopic_linewidth_q(basis%grid%xq(:,iq), input%nconf, input%T, &
!                                            sigma_ry, S, basis%grid, fc2)
!       !
!       vel = basis%c(:,:,iq)
!       !
!       ! Compute potentially anisotropic Casimir linewidth
!       IF(input%casimir_scattering) &
!         lw_casimir = casimir_linewidth_vel( basis%c(:,:,iq), input%casimir_length, input%casimir_dir, S%nat3)
!       ! Casimir linewidth is temperature/smearing-independent, sum it to all configurations
!       DO it = 1,input%nconf
!         lw(:,it) = lw_phph(:,it) + lw_isotopic(:,it) + lw_casimir
!       ENDDO
      !
    ENDDO QPOINT_OUTER_LOOP
    !
    !DEALLOCATE(A_iq)
    !
  END SUBROUTINE A_times_f
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute 3*nat lines of the matrix A (those corresponding to the
  ! q-point iq0) and apply them to the vector f. In output in gives 
  ! 3*nat elements of Af
  FUNCTION A_times_f_q(f, iq0, nconf, T, sigma, S, basis, fc2, fc3)
    USE q_grids,          ONLY : q_basis
    USE fc2_interpolate,  ONLY : forceconst2_grid, bose_phq, set_nu0, &
                                 freq_phq_safe, bose_phq, ip_cart2pat
    USE fc3_interpolate,  ONLY : forceconst3
    USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)      :: basis
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    REAL(DP),INTENT(in) :: f(nconf,S%nat3,basis%grid%nq)
    INTEGER,INTENT(in)  :: iq0
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    ! FUNCTION RESULT:
    REAL(DP) :: A_times_f_q(nconf,S%nat3)
    REAL(DP) :: P3(S%nat3,S%nat3) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:), V3Bsq(:,:,:)
    INTEGER :: iq, jq, mu, nu, it, nu0(4)
    !
    REAL(DP) :: freq(S%nat3,4), bose(S%nat3,4), xq(3,4), dq
    REAL(DP),PARAMETER :: epsq = 1.e-6_dp
    !

    ALLOCATE(U(S%nat3, S%nat3,4))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(V3Bsq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    A_times_f_q = 0._dp
    !
    dq = 1._dp/basis%grid%nq
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = basis%grid%xq(:,iq0)
    nu0(1) = set_nu0(xq(:,1), S%at)
    !freq(:,1) = basis%w(:,iq0)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    QPOINT_INNER_LOOP : &
    DO iq = 1, basis%grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = basis%grid%xq(:,iq)
      xq(:,3) = -xq(:,2)-xq(:,1)
      xq(:,4) =  xq(:,2)-xq(:,1)
      !write(*,'(4i3,4(3f8.3,3x))') nu0, xq
      DO jq = 2,4
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
        !write(*,'(99f12.4)') freq(:,jq)*RY_TO_CMM1
      ENDDO
      !
      ! Interpolate D3(q1,q2,-q1-q2)
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      !
      ! Interpolate D3(q1,-q2, q2-q1)
      ! For this process, we send q2 -> -q2,
      ! i.e. D2(-q2,q2) -> D2(q2,-q2) = D2(-q2,q2)*
      IF( ALL(ABS(xq(:,2))<epsq) ) THEN
        ! When q2 == 0, just copy over
        V3Bsq = V3sq
      ELSE
        xq(:,2) = -xq(:,2)
        U(:,:,2) = CONJG(U(:,:,2))
        CALL fc3%interpolate(xq(:,2), xq(:,4), S%nat3, D3)
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,4))
        V3Bsq = REAL( CONJG(D3)*D3 , kind=DP)
      ENDIF
      !
      DO it = 1,nconf
        bose(:,1) = basis%be(:,it,iq0)
        bose(:,2) = basis%be(:,it,iq)
        DO jq = 3,4
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
        !
        ! P3 is a 3*nat x 3*nat minor of the A matrix, the implicit indexes are 
        ! iq0 and iq, the matrix A has dimension (3*nat*nq x 3*nat*nq)
        P3 = sum_scattering_modes( S%nat3, sigma(it), freq, bose, V3sq, V3Bsq, nu0 )
        !
        ! 3*nat lines of the A matrix are applied now to f to produce 3*nat elements of Af
        DO mu = 1,S%nat3
        DO nu = 1,S%nat3
          A_times_f_q(it,nu) = A_times_f_q(it,nu)+ P3(nu,mu)*f(it,mu,iq)!*dq
        ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDDO &
    QPOINT_INNER_LOOP
    !
    !A_times_f_q = p/basis%grid%nq
    !
    DEALLOCATE(U, V3sq, V3Bsq, D3)
    !
  END FUNCTION A_times_f_q
  !  
  ! \/o\________\\\_________________________________________/^>
  ! Compute the scattering probability required for the matrix A.
  ! The form is a bit more complex than in  PRB 88, 045430 (2013)
  ! because we have rearranged the indexes to only have cohalescence
  ! processes instead of mixed scattering-cohalescence. This trick ensures
  ! that everything is positive definite.
  FUNCTION sum_scattering_modes(nat3, sigma, freq, bose, V3sq, V3Bsq, nu0)
    USE functions, ONLY : f_gauss => f_gauss
    USE constants, ONLY : pi, RY_TO_CMM1
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: freq(nat3,4)
    REAL(DP),INTENT(in) :: bose(nat3,4)
    REAL(DP),INTENT(in) :: V3sq(nat3,nat3,nat3)
    REAL(DP),INTENT(in) :: V3Bsq(nat3,nat3,nat3)
    INTEGER,INTENT(in)  :: nu0(4)
    !
    REAL(DP) :: sum_scattering_modes(nat3,nat3)
    !
    ! _C -> scattering, _X -> cohalescence
    REAL(DP) :: bose_a, bose_b, bose_c ! final/initial state populations 
    REAL(DP) :: dom_a, dom_b, dom_c   ! \delta\omega
    REAL(DP) :: ctm_a, ctm_b, ctm_c   !
    REAL(DP) :: norm_a, norm_bc
    REAL(DP) :: sum_a, sum_bc
    !
    REAL(DP),PARAMETER :: norm = pi/8
    !
    INTEGER :: i,j,k
    REAL(DP) :: p(nat3,nat3)
    p(:,:) = 0._dp
    
!     DO j = 1,4
! !     write(*,'("--",99f12.4)') freq(:,j)*RY_TO_CMM1
!     ENDDO

    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(i,j,k,bose_a,bose_b,bose_c,dom_a,dom_b,dom_c,ctm_a,ctm_b,ctm_c,sum_a,sum_bc,norm_a,norm_bc) &
!$OMP             REDUCTION(+: p) COLLAPSE(3)
    DO k = 1,nat3
      DO j = nu0(2),nat3
        DO i = nu0(1),nat3
          !
          bose_a = bose(i,1) * bose(j,2) * (bose(k,3)+1)
          bose_b = bose(i,1) * (bose(j,2)+1) * bose(k,4)
          bose_c = (bose(i,1)+1) * bose(j,2) * bose(k,4)
          !
          dom_a =  freq(i,1) + freq(j,2) - freq(k,3)
          dom_b =  freq(i,1) - freq(j,2) + freq(k,4)
          dom_c = -freq(i,1) + freq(j,2) + freq(k,4)
          !
          ctm_a = bose_a *  f_gauss(dom_a, sigma)
          ctm_b = bose_b *  f_gauss(dom_b, sigma)
          ctm_c = bose_c *  f_gauss(dom_c, sigma)
          !
          IF(k>=nu0(3))THEN
!             IF((freq(i,1)*freq(j,2)*freq(k,3))/=0._dp) THEN
            norm_a  = norm/(freq(i,1)*freq(j,2)*freq(k,3))
!             ELSE
!               WRITE(*,"(a,6i3,3f12.4)") "problem A", i,j,k,nu0(1),nu0(2),nu0(3),freq(i,1)*RY_TO_CMM1,freq(j,2)*RY_TO_CMM1,freq(k,3)*RY_TO_CMM1
!             ENDIF
          ELSE
            norm_a = 0._dp
          ENDIF
          !
          IF(k>=nu0(4))THEN
!             IF((freq(i,1)*freq(j,2)*freq(k,4))/=0._dp) THEN
            norm_bc = norm/(freq(i,1)*freq(j,2)*freq(k,4))
!             ELSE
!               WRITE(*,"(a,6i3,3f12.4)") "problem B", i,j,k,nu0(1),nu0(2),nu0(4),freq(i,1)*RY_TO_CMM1,freq(j,2)*RY_TO_CMM1,freq(k,4)*RY_TO_CMM1
!             ENDIF
          ELSE
            norm_bc = 0._dp
          ENDIF
          !
          sum_a  = norm_a * ctm_a * V3sq(i,j,k)
          !IF(ISNAN(sum_a)) STOP 1
          sum_bc = norm_bc * (ctm_b+ctm_c) * V3Bsq(i,j,k)
          !IF(ISNAN(sum_bc)) STOP 2
          !
          p(i,j) = p(i,j) -sum_a + sum_bc
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    sum_scattering_modes = p
    !
  END FUNCTION sum_scattering_modes
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute thermal conductivity as \lambda \over { N T^2 } f \dot b
  SUBROUTINE calc_tk_simple(f, b, T, Omega, nconf, nat3, nq)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE timers
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: Omega ! cell volume (bohr^3)
    !
    REAL(DP) :: tk(3,3,nconf)
    INTEGER  :: iq, it, ix, jx, nu
    !
    !
    tk = 0._dp
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  f(ix,it,nu,iq)*b(jx,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
!!$OMP ENDDO      
    ENDDO
    ENDDO
!!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
    WRITE(*,*) "tk_sma"
    DO it = 1,nconf
      WRITE(*,'(i6,f12.6,9e20.8)') it, T(it), &
      tk(1,1,it)*RY_TO_WATTMM1KM1, &
      tk(2,2,it)*RY_TO_WATTMM1KM1, &
      tk(3,3,it)*RY_TO_WATTMM1KM1, &
      tk(1,2,it)*RY_TO_WATTMM1KM1, &
      tk(1,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,3,it)*RY_TO_WATTMM1KM1
    ENDDO
    !
  END SUBROUTINE calc_tk_simple
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute thermal conductivity as -2\lambda \over { N T^2 } (1/2 f \dot Af - f \dot b)
  SUBROUTINE calc_tk_variational(f, Af, b, T, Omega, nconf, nat3, nq)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE timers
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: Af(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: Omega ! cell volume (bohr^3)
    !
    REAL(DP) :: tk(3,3,nconf)
    INTEGER  :: iq, it, ix, jx, nu
    !
    !
    tk = 0._dp
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  0.5_dp*f(ix,it,nu,iq)*Af(jx,it,nu,iq) - f(ix,it,nu,iq)*b(jx,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
!!$OMP ENDDO      
    ENDDO
    ENDDO
!!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = -2*Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
    WRITE(*,*) "tk_var"
    DO it = 1,nconf
      WRITE(*,'(i6,f12.6,9e20.8)') it, T(it), &
      tk(1,1,it)*RY_TO_WATTMM1KM1, &
      tk(2,2,it)*RY_TO_WATTMM1KM1, &
      tk(3,3,it)*RY_TO_WATTMM1KM1, &
      tk(1,2,it)*RY_TO_WATTMM1KM1, &
      tk(1,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,3,it)*RY_TO_WATTMM1KM1
    ENDDO
    !
  END SUBROUTINE calc_tk_variational
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute thermal conductivity as -2\lambda \over { N T^2 } (1/2 f \dot Af - f \dot b)
  SUBROUTINE calc_tk_gf(g, f, b, T, Omega, nconf, nat3, nq)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE timers
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: g(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: Omega ! cell volume (bohr^3)
    !
    REAL(DP) :: tk(3,3,nconf)
    INTEGER  :: iq, it, ix, jx, nu
    !
    !
    tk = 0._dp
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  0.5_dp*( f(ix,it,nu,iq)*g(jx,it,nu,iq) &
                                                 -f(ix,it,nu,iq)*b(jx,it,nu,iq) )
        ENDDO
        ENDDO
      ENDDO
!!$OMP ENDDO      
    ENDDO
    ENDDO
!!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = -2*Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
    WRITE(*,*) "tk_gf"
    DO it = 1,nconf
      WRITE(*,'(i6,f12.6,9e20.8)') it, T(it), &
      tk(1,1,it)*RY_TO_WATTMM1KM1, &
      tk(2,2,it)*RY_TO_WATTMM1KM1, &
      tk(3,3,it)*RY_TO_WATTMM1KM1, &
      tk(1,2,it)*RY_TO_WATTMM1KM1, &
      tk(1,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,3,it)*RY_TO_WATTMM1KM1
    ENDDO
    !
  END SUBROUTINE calc_tk_gf
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION b_over_A(b, A, nconf, nat3, nq) RESULT(f)
    IMPLICIT NONE
    !
    REAL(DP):: f(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, ix, nu
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq,nu,it,ix) COLLAPSE(4)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
      DO ix = 1,3
        IF(A(it,nu,iq)/=0._dp)THEN
          f(ix,it,nu,iq) = b(ix,it,nu,iq)/A(it,nu,iq)
        ELSE
          f(ix,it,nu,iq) = 0._dp
        ENDIF
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO    
    !
  END FUNCTION b_over_A
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION A_out_times_f(A_out, f, nconf, nat3, nq) RESULT(Af)
    IMPLICIT NONE
    !
    REAL(DP):: Af(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A_out(nconf, nat3, nq)
    REAL(DP),INTENT(in) :: f(3, nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, ix, nu
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq,nu,it,ix) COLLAPSE(4)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
      DO ix = 1,3
        Af(ix,it,nu,iq) = A_out(it,nu,iq)*f(ix,it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO    
    !
  END FUNCTION A_out_times_f
  
END MODULE variational_tk
