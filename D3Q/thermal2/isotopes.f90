!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE isotopes_linewidth

  USE kinds, ONLY : DP
  USE input_fc, ONLY : ph_system_info, forceconst2_grid
  USE fc2_interpolate, ONLY : mat2_diag
  !
  !
  CONTAINS
  !
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Returns the HALF width half maximum (note the 0.5 factor) of phonon modes
  ! due to isotope-isotope scattering
  FUNCTION isotopic_linewidth_q(xq0, nconf, T, sigma, S, grid, fc2, freq1, U1)
  !RESULT(lw)
    !
    USE nanoclock
    !
    USE q_grids,           ONLY : q_grid
    USE fc2_interpolate,   ONLY : freq_phq_safe, bose_phq
    USE nist_isotopes_db,  ONLY : compute_gs
    USE mpi_thermal,       ONLY : mpi_bsum

    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    !
    REAL(DP),OPTIONAL,INTENT(in) :: freq1(S%nat3)
    COMPLEX(DP),OPTIONAL,INTENT(in) :: U1(S%nat3,S%nat3)
    !
    ! FUNCTION RESULT:
    REAL(DP) :: isotopic_linewidth_q(S%nat3,nconf)
    REAL(DP) :: lw(S%nat3,nconf)
    !
    COMPLEX(DP) :: U(S%nat3, S%nat3,2)
    INTEGER  :: iq, jq, it
    !
    REAL(DP) :: freq(S%nat3,3), xq(3,2)
    !REAL(DP) :: gs2(S%nat3) !, gm(S%nat3), auxm(S%ntyp), auxs(S%ntyp)
    !
    lw = 0._dp
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
     IF(present(freq1) .and. present(U1)) THEN
      freq(:,1) = freq1
      U(:,:,1)    = U1
    ELSE
      CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    ENDIF
    !
    DO iq = 1, grid%nq
      !
      ! Compute eigenvalues and eigenmodes at q2 and q3
      xq(:,2) = grid%xq(:,iq)
      CALL freq_phq_safe(xq(:,2), S, fc2, freq(:,2), U(:,:,2))
      !
      DO it = 1,nconf
        !
        lw(:,it) = lw(:,it) &
                  +0.5_dp * grid%w(iq)                      &
                    * sum_isotope_linewidth_modes(          &
                        S%nat3, S%nat, sigma(it), freq,     &
                        S%ntyp, S%ityp, S%amass_variance, U )
        !
      ENDDO
      !
    ENDDO
    !
    IF(grid%scattered) CALL mpi_bsum(S%nat3,nconf,lw)
    isotopic_linewidth_q = lw
    !
  END FUNCTION isotopic_linewidth_q
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_isotope_linewidth_modes(nat3, nat, sigma, freq, ntyp, ityp, gs2, zz)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3, nat
    REAL(DP),INTENT(in) :: sigma
    INTEGER,INTENT(in)  :: ntyp
    INTEGER,INTENT(in)  :: ityp(nat)
    REAL(DP),INTENT(in) :: gs2(ntyp)
    REAL(DP),INTENT(in) :: freq(nat3,2)
    COMPLEX(DP),INTENT(in) :: zz(nat3,nat3,2)
    !
    REAL(DP) :: sum_isotope_linewidth_modes(nat3)
    !
    REAL(DP) :: freq_f, sum_zz2
    COMPLEX(DP) :: sum_zz
!    COMPLEX(DP) :: czz1(nat3,nat3), zz2(nat3,nat3)
    !
    !
    INTEGER :: i,j, ia, it, ix, nu
    REAL(DP) :: lw(nat3)
    lw(:) = 0._dp
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,ia,it,ix,nu,freq_f,sum_zz,sum_zz2) &
!$OMP REDUCTION(+:lw) COLLAPSE(2)
    DO i = 1,nat3
      DO j = 1,nat3
        !
        freq_f = freq(i,1)*freq(j,2) * f_gauss(freq(i,1)-freq(j,2), sigma)
        !
        !IF(freq_f > 1.d-8)THEN
        sum_zz2 = 0._dp
        nu = 0
        DO ia = 1,nat
          it = ityp(ia)
          sum_zz = 0._dp
          DO ix = 1,3
            nu = nu + 1
            !
            sum_zz =  sum_zz + CONJG(zz(nu,i,1)) * zz(nu,j,2)
            !
          ENDDO
          sum_zz2 = sum_zz2 + gs2(it)*REAL(CONJG(sum_zz)*sum_zz, kind=DP)
        ENDDO
        !ENDIF
        !
        lw(i) = lw(i) + freq_f*sum_zz2
      ENDDO
    ENDDO
!$OMP END PARALLEL DO 
    !
    sum_isotope_linewidth_modes = 0.5_dp*pi*lw
    !
  END FUNCTION sum_isotope_linewidth_modes
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_isotope_scattering_modes(nat3, nat, sigma, freq, bose, ntyp, ityp, gs2, zz)
    USE functions, ONLY : f_gauss
    USE constants, ONLY : pi
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nat3, nat
    REAL(DP),INTENT(in) :: sigma
    INTEGER,INTENT(in)  :: ntyp
    INTEGER,INTENT(in)  :: ityp(nat)
    REAL(DP),INTENT(in) :: gs2(ntyp)
    REAL(DP),INTENT(in) :: freq(nat3,2)
    REAL(DP),INTENT(in) :: bose(nat3,2)
    COMPLEX(DP),INTENT(in) :: zz(nat3,nat3,2)
    !
    REAL(DP) :: sum_isotope_scattering_modes(nat3,nat3)
    !
    REAL(DP) :: freq_f, bose_f, sum_zz2
    COMPLEX(DP) :: sum_zz
    !
    !
    INTEGER :: i,j, ia, it, ix, nu
    REAL(DP) :: P(nat3, nat3)
    P(:,:) = 0._dp
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,ia,it,ix,nu,bose_f,freq_f,sum_zz,sum_zz2) &
!$OMP REDUCTION(+:P) COLLAPSE(2)
    DO j = 1,nat3
      DO i = 1,nat3
        !
        freq_f = freq(i,1)*freq(j,2) * f_gauss(freq(i,1)-freq(j,2), sigma)
        bose_f = bose(i,1)*bose(j,2) + 0.5_dp*(bose(i,1)+bose(j,2))
        !
        !IF(freq_f > 1.d-8)THEN
        nu = 0
        sum_zz2 = 0._dp
        DO ia = 1,nat
          it = ityp(ia)
          !
          sum_zz = 0._dp
          DO ix = 1,3
            nu = nu + 1
            !
            sum_zz =  sum_zz + CONJG(zz(nu,i,1)) * zz(nu,j,2)
            !
          ENDDO
          sum_zz2 = sum_zz2 + gs2(it)*REAL(CONJG(sum_zz)*sum_zz, kind=DP)
        ENDDO
        !ENDIF
        !
        P(i,j) = P(i,j) + bose_f*freq_f*sum_zz2
      ENDDO
    ENDDO
!$OMP END PARALLEL DO 
    !
    sum_isotope_scattering_modes = 0.5_dp*pi*P
    !
  END FUNCTION sum_isotope_scattering_modes
  !
END MODULE isotopes_linewidth
!
