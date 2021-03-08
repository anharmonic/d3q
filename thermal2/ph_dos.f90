
MODULE ph_dos

  CONTAINS
  
!   ! \/o\________\\\_________________________________________/^>
  FUNCTION joint_dos_q(qgrid, sigma, T, S, fc, xq_i)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid !, q_basis, setup_grid, prepare_q_basis
!     USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    USE random_numbers,   ONLY : randy
    IMPLICIT NONE
    TYPE(q_grid),INTENT(in)  :: qgrid
    TYPE(ph_system_info)     :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP), INTENT(in) :: sigma, T, xq_i(3)
    REAL(DP) :: joint_dos_q
    !
!    TYPE(q_basis) :: qbasis, sbasis
    !
    REAL(DP) :: freqi(S%nat3), freqj(S%nat3), freqk(S%nat3)
    REAL(DP) :: bosej(S%nat3), bosek(S%nat3)
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    !
    !REAL(DP) :: xq0(3) = (/ 0._dp, 0._dp, 0._dp /)
    !REAL(DP) :: xq_random(3)
    !
    REAL(DP) :: jd_C(S%nat3), jd_X(S%nat3)
    REAL(DP) :: xq_j(3), xq_k(3)
    REAL(DP) :: jdq
    REAL(DP) :: dom_C(S%nat3), dom_X(S%nat3)
    REAL(DP) ::  ctm_C(S%nat3), ctm_X(S%nat3), bose_C, bose_X
    INTEGER :: jq, k,j,i
    !
    !sigma_ry = input%sigma(1)/RY_TO_CMM1
    
    CALL freq_phq(xq_i, S, fc, freqi, U)
      
    jd_C = 0._dp
    jd_X = 0._dp

    !xq_random  = (/ randy(), randy(), randy() /)
!     CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
!                 qgrid, scatter=.true.)
    
    DO jq = 1, qgrid%nq
      xq_j = qgrid%xq(:,jq)
      CALL freq_phq(xq_j, S, fc, freqj, U)
      bosej(:) = f_bose(freqj, T)
      !WRITE(20001, '(6f14.6)') freqj*RY_TO_CMM1
      
      xq_k = -(xq_i + xq_j)
      CALL freq_phq(xq_k, S, fc, freqk, U)
      bosek(:) = f_bose(freqk, T)
      !
      DO k = 1,S%nat3
        DO j = 1,S%nat3
          !
          bose_C = 2*(bosej(j) - bosek(k))
          dom_C(:) =freqi(:)+freqj(j)-freqk(k) ! cohalescence
          ctm_C = bose_C * f_gauss(dom_C, sigma) !delta 
          !
          bose_X = bosej(j) + bosek(k) + 1
          dom_X(:) =freqi(:)-freqj(j)-freqk(k) ! scattering/decay
          ctm_X = bose_X * f_gauss(dom_X, sigma) !delta
          !
          jd_C(:) = jd_C(:) + qgrid%w(jq)*ctm_C(:)
          jd_X(:) = jd_X(:) + qgrid%w(jq)*ctm_X(:)
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
    jdq = SUM(jd_C+jd_X)
    IF(qgrid%scattered) CALL mpi_bsum(jdq)
    joint_dos_q =jdq
    !print*, joint_dos_q
    !
  END FUNCTION joint_dos_q
!   ! \/o\________\\\_________________________________________/^>

END MODULE ph_dos
