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
  ! Multiply matrix in diagonal form, like A_out, A_out^-1, A_out^-1/2 
  ! with a vector (f, b, etc)
  FUNCTION A_diag_f(A, f, nconf, nat3, nq) RESULT(Af)
    IMPLICIT NONE
    !
    REAL(DP):: Af(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
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
        Af(ix,it,nu,iq) = f(ix,it,nu,iq)*A(it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO    
    !
  END FUNCTION A_diag_f
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the diagonal matrix A_out
  SUBROUTINE compute_A_out(A_out, input, qbasis, S, fc2, fc3)
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
      IF(input%isotopic_disorder) THEN
        lw_isotopic = isotopic_linewidth_q(qbasis%grid%xq(:,iq), input%nconf, input%T, &
                                           sigma_ry, S, qbasis%grid, fc2)
      ENDIF
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
  END SUBROUTINE compute_A_out
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE compute_inv_sqrt_A_out(A, inv_sqrt_A, nconf, nat3, nq)
    IMPLICIT NONE
    !
    REAL(DP):: f(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(out) :: inv_sqrt_A(nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, nu
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq,nu,it) COLLAPSE(3)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
        IF(A(it,nu,iq)>0._dp)THEN
          inv_sqrt_A(it,nu,iq) = 1/DSQRT(A(it,nu,iq))
        ELSE IF (A(it,nu,iq)<0._dp) THEN
          inv_sqrt_A(it,nu,iq) = -1/SQRT( -A(it,nu,iq) )
          IF(iq/=1) WRITE(*,*) "Warning: negative A_out", iq, nu, it
        ELSE
          ! This hould be infinity, but it should only happen at Gamma for
          ! acoustic bands, where we can ignore it because everything is zero
          inv_sqrt_A(it,nu,iq) = 0._dp
          IF(iq/=1) WRITE(*,*) "Suspicious null A_out", iq, nu, it
        ENDIF
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO    
    !
  END SUBROUTINE compute_inv_sqrt_A_out
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE compute_inv_A_out(A, inv_A, nconf, nat3, nq)
    IMPLICIT NONE
    !
    REAL(DP):: f(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
    REAL(DP),INTENT(out) :: inv_A(nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, nu
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq,nu,it) COLLAPSE(3)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
        IF(A(it,nu,iq)/=0._dp)THEN
          inv_A(it,nu,iq) = 1/(A(it,nu,iq))
          ! This hould be infinity, but it should only happen at Gamma for
          ! acoustic bands, where we can ignore it because everything is zero
        ELSE
          inv_A(it,nu,iq) = 0._dp
          IF(iq/=1) WRITE(*,*) "Suspicious null A_out", iq, nu, it
        ENDIF
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO    
    !
  END SUBROUTINE compute_inv_A_out
  !
  ! \/o\________\\\_________________________________________/^>
  ! Apply the A_in matrix to a vector f
  ! this is left as a subroutine to underline that it is the most time consuming step.
  SUBROUTINE A_in_times_f(f, Af, input, basis, S, fc2, fc3)
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
    REAL(DP),INTENT(in)  ::  f(3,input%nconf,S%nat3, basis%grid%nq)
    REAL(DP),INTENT(out) :: Af(3,input%nconf,S%nat3, basis%grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    !
    INTEGER  :: iq, it, ix, nu
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    QPOINT_OUTER_LOOP : & ! the inner loop in in A_in_times_f_q
    DO iq = 1,basis%grid%nq
      ! bang!
       Af(:,:,:,iq) = A_in_times_f_q(f, iq, input%nconf, input%T, sigma_ry, S, basis, &
                                     fc2, fc3, input%isotopic_disorder)
      !
      !
    ENDDO QPOINT_OUTER_LOOP
    !
  END SUBROUTINE A_in_times_f
  !
  ! \/o\________\\\_________________________________________/^>
  ! Apply the A = (A_out+A_in) matrix to a vector f, A_out must be already computed
  ! this is left as a subroutine to underline that it is the most time consuming step.
  SUBROUTINE A_times_f(f, Af, A_out, input, basis, S, fc2, fc3)
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
    REAL(DP),INTENT(in)  ::     f(3,input%nconf,S%nat3, basis%grid%nq)
    REAL(DP),INTENT(in)  :: A_out(input%nconf,S%nat3, basis%grid%nq)
    !
    REAL(DP),INTENT(out) ::    Af(3,input%nconf,S%nat3, basis%grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    !
    INTEGER  :: iq, it, ix, nu
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    QPOINT_OUTER_LOOP : & ! the inner loop is in A_in_times_f_q
    DO iq = 1,basis%grid%nq
      ! bang!
      Af(:,:,:,iq) = A_in_times_f_q(f, iq, input%nconf, input%T, sigma_ry, S, basis, &
                                    fc2, fc3, input%isotopic_disorder)
      DO nu = 1,S%nat3
        DO it = 1,input%nconf
        DO ix = 1,3
          Af(ix,it,nu,iq) = Af(ix,it,nu,iq) + A_out(it,nu,iq)*f(ix,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
      !
    ENDDO QPOINT_OUTER_LOOP
    !
  END SUBROUTINE A_times_f
  !
  ! \/o\________\\\_________________________________________/^>
  ! Apply the \tilde{A} = (1+A_out^{-1/2} A_in A_out^{-1/2}) matrix to a vector f,
  ! A_out^{-1/2} must be already computed and provided in input
  ! this is left as a subroutine to underline that it is the most time consuming step.
  SUBROUTINE tilde_A_times_f(f, Af, inv_sqrt_A_out, input, basis, S, fc2, fc3)
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
    REAL(DP),INTENT(in)  ::     f(3,input%nconf,S%nat3,basis%grid%nq)
    REAL(DP),INTENT(in)  :: inv_sqrt_A_out(input%nconf,S%nat3,basis%grid%nq)
    !
    REAL(DP),INTENT(out) ::    Af(3,input%nconf,S%nat3,basis%grid%nq)
    !
    REAL(DP) :: sigma_ry(input%nconf)
    !
    REAL(DP),ALLOCATABLE :: aux(:,:,:,:)
    INTEGER  :: iq, it, ix, nu
    !
    ALLOCATE(aux(3,input%nconf,S%nat3,basis%grid%nq))
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    ! Apply the first 1/sqrt(A_out)
    DO iq = 1,basis%grid%nq
      ! bang!
      DO nu = 1,S%nat3
        DO it = 1,input%nconf
        DO ix = 1,3
          aux(ix,it,nu,iq) = inv_sqrt_A_out(it,nu,iq)*f(ix,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
    ENDDO
      
    DO iq = 1,basis%grid%nq
      ! apply A_in
      Af(:,:,:,iq) = A_in_times_f_q(aux, iq, input%nconf, input%T, sigma_ry, S, basis, &
                                    fc2, fc3, input%isotopic_disorder)
      ! Apply 1+ and the second 1/sqrt(A)
      DO nu = 1,S%nat3
        DO it = 1,input%nconf
        DO ix = 1,3
          Af(ix,it,nu,iq) = f(ix,it,nu,iq) + inv_sqrt_A_out(it,nu,iq)*Af(ix,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
      !
    ENDDO
    !
    DEALLOCATE(aux)
    !
  END SUBROUTINE tilde_A_times_f
  ! \/o\________\\\_________________________________________/^>
  ! Compute 3*nat lines of the matrix A (those corresponding to the
  ! q-point iq0) and apply them to the vector f. In output it returns
  ! 3*nat elements of Af (for each configuration)
  FUNCTION A_in_times_f_q(f, iq0, nconf, T, sigma, S, basis, fc2, fc3, isotopic_disorder)
    USE q_grids,            ONLY : q_basis
    USE fc2_interpolate,    ONLY : forceconst2_grid, bose_phq, set_nu0, &
                                 freq_phq_safe, bose_phq, ip_cart2pat
    USE fc3_interpolate,    ONLY : forceconst3
    USE isotopes_linewidth, ONLY : sum_isotope_scattering_modes
    !USE constants, ONLY : RY_TO_CMM1
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_basis),INTENT(in)      :: basis
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    REAL(DP),INTENT(in) :: f(3,nconf,S%nat3,basis%grid%nq)
    INTEGER,INTENT(in)  :: iq0
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    LOGICAL,INTENT(in)  :: isotopic_disorder
    !
    ! FUNCTION RESULT:
    REAL(DP) :: A_in_times_f_q(3,nconf,S%nat3)
    REAL(DP) :: P3(S%nat3,S%nat3) !aux
    REAL(DP) :: P3_isot(S%nat3,S%nat3) !aux
    !
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:), V3Bsq(:,:,:)
    INTEGER :: iq, jq, mu, nu, it, nu0(4), ix
    !
    REAL(DP) :: freq(S%nat3,5), bose(S%nat3,5), xq(3,5), dq
    REAL(DP),PARAMETER :: epsq = 1.e-6_dp
    !
    ! xq(:,1) -> xq1
    ! xq(:,2) -> xq2
    ! xq(:,3) -> -xq1-xq2
    ! xq(:,4) -> -xq1+xq2
    ! xq(:,5) -> -xq2
    ! And accordingly for U(:,:,i).
    ! Note that U(5) = CONJG(U(2)) and freq(5)=freq(2)
    

    ALLOCATE(U(S%nat3, S%nat3,5))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(V3Bsq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    !
    A_in_times_f_q = 0._dp
    P3_isot = 0._dp
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
      xq(:,5) = -xq(:,2) ! => xq4 = -xq5-xq1
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
        ! xq(:,5) = 0._dp
        ! U(:,:,5) = U(:,:,2)
      ELSE
        !CALL freq_phq_safe(xq(:,5), S, fc2, freq(:,5), U(:,:,5))
        !IF(ANY(ABS(U(:,:,2)-CONJG(U(:,:,5)))>1.d-10)  ) STOP 999
        U(:,:,5) = CONJG(U(:,:,2))
        CALL fc3%interpolate(xq(:,5), xq(:,4), S%nat3, D3)
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,5), U(:,:,4))
        V3Bsq = REAL( CONJG(D3)*D3 , kind=DP)
      ENDIF
      !
      ! These loops are not in the optimal order, but it does not matter for performance
      DO it = 1,nconf
        bose(:,1) = basis%be(:,it,iq0)
        bose(:,2) = basis%be(:,it,iq)
        DO jq = 3,4
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
        !
        ! P3 is a 3*nat x 3*nat minor of the A matrix, the implicit indexes are 
        ! iq0 and iq, the matrix A has dimension (3*nat*nq x 3*nat*nq)
        ! Do not forget the minus sign!!
        P3 = - sum_scattering_modes( S%nat3, sigma(it), freq, bose, V3sq, V3Bsq, nu0 )
        !
        IF(isotopic_disorder)THEN
          P3_isot = sum_isotope_scattering_modes(S%nat3, S%nat, sigma(it), freq, &
                                              bose, S%ntyp, S%ityp, S%amass_variance, U)
          P3 = P3 + P3_isot
        ENDIF
        !
        ! 3*nat lines of the A matrix are applied now to f to produce 3*nat elements of Af
        DO mu = 1,S%nat3
        DO nu = 1,S%nat3
          DO ix = 1,3
            A_in_times_f_q(ix,it,nu) = A_in_times_f_q(ix,it,nu) + P3(nu,mu)*f(ix,it,mu,iq)*dq
          ENDDO
        ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDDO &
    QPOINT_INNER_LOOP
    !
    !A_in_times_f_q = p/basis%grid%nq
    !
    DEALLOCATE(U, V3sq, V3Bsq, D3)
    !
  END FUNCTION A_in_times_f_q
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
   
    !
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP          PRIVATE(i, j, k, bose_a, bose_b, bose_c, dom_a, dom_b, dom_c) &
!$OMP          PRIVATE(ctm_a, ctm_b, ctm_c, sum_a, sum_bc, norm_a, norm_bc) &
!$OMP          REDUCTION(+: p)
!$OMP DO  COLLAPSE(3)
    DO k = nu0(3),nat3
      DO j = nu0(2),nat3
        DO i = nu0(1),nat3
          !
          bose_a = bose(i,1) * bose(j,2) * (bose(k,3)+1)
          dom_a =  freq(i,1) + freq(j,2) - freq(k,3)
          ctm_a = bose_a *  f_gauss(dom_a, sigma)
          norm_a  = norm/(freq(i,1)*freq(j,2)*freq(k,3))
          sum_a  = norm_a  * ctm_a         * V3sq(i,j,k)
          p(i,j) = p(i,j) - sum_a
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO COLLAPSE(3)
    DO k = nu0(4),nat3
      DO j = nu0(2),nat3
        DO i = nu0(1),nat3
          !
          bose_b = bose(i,1) * (bose(j,2)+1) * bose(k,4)
          bose_c = (bose(i,1)+1) * bose(j,2) * bose(k,4)
          dom_b =  freq(i,1) - freq(j,2) + freq(k,4)
          dom_c = -freq(i,1) + freq(j,2) + freq(k,4)
          ctm_b = bose_b *  f_gauss(dom_b, sigma)
          ctm_c = bose_c *  f_gauss(dom_c, sigma)
          norm_bc = norm/(freq(i,1)*freq(j,2)*freq(k,4))
          sum_bc = norm_bc * (ctm_b+ctm_c) * V3Bsq(i,j,k)
          p(i,j) = p(i,j) + sum_bc
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL 
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  f(ix,it,nu,iq)*b(jx,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO      
    ENDDO
    ENDDO
!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
    CALL print_tk(tk, nconf, "simple")
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  0.5_dp*f(ix,it,nu,iq)*Af(jx,it,nu,iq) - f(ix,it,nu,iq)*b(jx,it,nu,iq)
        ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO      
    ENDDO
    ENDDO
!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = -2*Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
    CALL print_tk(tk, nconf, "variational")
    !
  END SUBROUTINE calc_tk_variational
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute thermal conductivity as -2\lambda \over { N T^2 } 1/2 ( f \dot g - f \dot b) )
  ! g = Af-b
  ! f.g = f.Af - b.f
  ! 1/2 (f.g -b.f ) = 1/2 f.Af - b.f
  FUNCTION calc_tk_gf(g, f, b, T, Omega, nconf, nat3, nq) RESULT(tk)
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) +  0.5_dp*( f(ix,it,nu,iq)*g(jx,it,nu,iq) &
                                                 -f(ix,it,nu,iq)*b(jx,it,nu,iq) )
        ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO      
    ENDDO
    ENDDO
!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = -2*Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
  END FUNCTION calc_tk_gf
  ! \/o\________\\\_________________________________________/^>
  ! Compute thermal conductivity as -2\lambda \over { N T^2 } 1/2 ( f \dot g - f \dot b) )
  ! g = Af-b
  ! f.g = f.Af - b.f
  ! 1/2 (f.g -b.f ) = 1/2 f.Af - b.f
  FUNCTION calc_tk_rxb(x, r, b, T, Omega, nconf, nat3, nq) RESULT(tk)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    USE timers
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: x(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: r(3, nconf, nat3, nq)
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix,jx)
    DO iq = 1,nq
    DO nu = 1,nat3
!$OMP DO COLLAPSE(3)    
      DO it = 1,nconf
        DO jx = 1,3
        DO ix = 1,3
          tk(ix,jx,it) = tk(ix,jx,it) -  0.5_dp*( x(ix,it,nu,iq)*r(jx,it,nu,iq) &
                                                 -x(ix,it,nu,iq)*b(jx,it,nu,iq) )
        ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO      
    ENDDO
    ENDDO
!$OMP END PARALLEL
    DO it = 1,nconf
      tk(:,:,it) = -2*Omega*tk(:,:,it)/nq/T(it)**2
    ENDDO
    !
  END FUNCTION calc_tk_rxb
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_tk(tk, nconf, name)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: tk(3,3,nconf)
    INTEGER,INTENT(in) :: nconf
    CHARACTER(len=*),INTENT(in) :: name
    !
    INTEGER :: it
    WRITE(*,*) name
    DO it = 1,nconf
      WRITE(*,'(i6,3(3e17.8,3x))') it, &
      tk(1,1,it)*RY_TO_WATTMM1KM1, &
      tk(2,2,it)*RY_TO_WATTMM1KM1, &
      tk(3,3,it)*RY_TO_WATTMM1KM1, &
      tk(1,2,it)*RY_TO_WATTMM1KM1, &
      tk(1,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,3,it)*RY_TO_WATTMM1KM1, &
      tk(2,1,it)*RY_TO_WATTMM1KM1, &
      tk(3,1,it)*RY_TO_WATTMM1KM1, &
      tk(3,2,it)*RY_TO_WATTMM1KM1

    ENDDO
    !
  END SUBROUTINE print_tk
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_gradmod2_tk(g, name, T, Omega, nconf, nat3, nq)
    USE more_constants,     ONLY : RY_TO_WATTMM1KM1
    IMPLICIT NONE
    REAL(DP),INTENT(in) :: g(3,nconf)
    CHARACTER(len=*),INTENT(in) :: name
    REAL(DP),INTENT(in) :: T(nconf)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    REAL(DP),INTENT(in) :: Omega ! cell volume (bohr^3)
    !
    INTEGER :: it
    WRITE(*,*) name
    DO it = 1,nconf
      WRITE(*,'(12x,i6,3(3e17.8,3x))') it, g(:,it)*RY_TO_WATTMM1KM1*Omega/nq/T(it)**2
    ENDDO
    !
  END SUBROUTINE print_gradmod2_tk
  !
!   ! \/o\________\\\_________________________________________/^>
!   FUNCTION b_over_A(b, A, nconf, nat3, nq) RESULT(f)
!     IMPLICIT NONE
!     !
!     REAL(DP):: f(3, nconf, nat3, nq)
!     !
!     REAL(DP),INTENT(in) :: A(nconf, nat3, nq)
!     REAL(DP),INTENT(in) :: b(3, nconf, nat3, nq)
!     INTEGER,INTENT(in)  :: nconf, nat3, nq
!     !
!     INTEGER  :: iq, it, ix, nu
!     !
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq,nu,it,ix) COLLAPSE(4)
!     DO iq = 1,nq
!     DO nu = 1,nat3
!       DO it = 1,nconf
!       DO ix = 1,3
!         IF(A(it,nu,iq)/=0._dp)THEN
!           f(ix,it,nu,iq) = b(ix,it,nu,iq)/A(it,nu,iq)
!         ELSE
!           f(ix,it,nu,iq) = 0._dp
!         ENDIF
!       ENDDO
!       ENDDO
!     ENDDO
!     ENDDO
! !$OMP END PARALLEL DO    
!     !
!   END FUNCTION b_over_A
  !
  !  
END MODULE variational_tk
