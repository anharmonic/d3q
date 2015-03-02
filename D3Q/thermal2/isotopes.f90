!
! Written by Lorenzo Paulatto (2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE isotopes_linewidth

  USE kinds, ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid
  USE fc2_interpolate, ONLY : fftinterp_mat2, mat2_diag
  !
  !
  CONTAINS
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION isotopic_linewidth_q(xq0, nconf, T, sigma, S, grid, fc2) &
  RESULT(lw)
    !
    USE nanoclock
    !
    USE q_grid,                 ONLY : q_grid_type
    USE fc2_interpolate,             ONLY : freq_phq_safe, bose_phq
    USE nist_isotopes_db, ONLY : compute_gs


    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid_type),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    ! FUNCTION RESULT:
    REAL(DP) :: lw(S%nat3,nconf)
    !
    COMPLEX(DP) :: U(S%nat3, S%nat3,2)
    INTEGER  :: iq, jq, nu, it, ia, alpha
    !
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,2)
    REAL(DP) :: gs2(S%nat3) !, gm(S%nat3), auxm(S%ntyp), auxs(S%ntyp)
    !
    lw = 0._dp
!     ! Compute isotope-phonon cross section for the various atomic species
!     ! NOTE: move this code to input!
!     DO it = 1, S%ntyp
!       CALL compute_gs(auxm(it), auxs(it), S%atm(it), 0, 0)
!     ENDDO
    ! Rearrange the cross sections with the mode index, for simplicity
    nu = 0
    DO ia = 1, S%nat
      it = S%ityp(ia)
      DO alpha = 1,3
        nu = nu+1
        !gm(nu)  = auxm(it)
        gs2(nu) = S%amass_variance(it)
      ENDDO
    ENDDO
    !
!     print*, "gm and gs: "
!     print*, gm
!     print*, gs2
    
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      CALL freq_phq_safe(xq(:,2), S, fc2, freq(:,2), U(:,:,2))
      !
      DO it = 1,nconf
        !
        lw(:,it) = lw(:,it) + sum_isotope_modes( S%nat3, sigma(it), freq, gs2, U )
        !
      ENDDO
      !
    ENDDO
    !
    lw = lw/grid%nq
    !
  END FUNCTION isotopic_linewidth_q
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_isotope_modes(nat3, sigma, freq, gs2, zz)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3
    REAL(DP),INTENT(in) :: sigma
    REAL(DP),INTENT(in) :: gs2(nat3)
    REAL(DP),INTENT(in) :: freq(nat3,2)
!     REAL(DP),INTENT(in) :: bose(nat3,2)
    COMPLEX(DP),INTENT(in) :: zz(nat3,nat3,2)
    !
    REAL(DP) :: sum_isotope_modes(nat3)
    !
    REAL(DP) :: bose_f, freq_f, sum_zz2, dfreq
    COMPLEX(DP) :: sum_zz
    !
    !
    INTEGER :: i,j
    REAL(DP) :: lw(nat3)
    lw(:) = 0._dp
    !
    !
    DO j = 1,nat3
      DO i = 1,nat3
        !
        dfreq = f_gauss(freq(i,1)-freq(j,2), sigma)
        !
        !IF(dfreq>eps)THEN
        !bose_f = bose(i,1)*bose(j,2)+0.5_dp*(bose(i,1)+bose(j,2))
        freq_f = freq(i,1)*freq(j,2)
        !
        sum_zz = SUM( CONJG(zz(:,i,1)) * zz(:,j,2) )
        sum_zz2 = REAL(CONJG(sum_zz)*sum_zz, kind=DP)
        !
        lw(i) = lw(i) + gs2(i)*freq_f*dfreq*sum_zz2
        !ENDIF
        !
      ENDDO
    ENDDO
    !
    sum_isotope_modes = 0.5_dp*pi*lw
    !
  END FUNCTION sum_isotope_modes
  !
END MODULE isotopes_linewidth
!
