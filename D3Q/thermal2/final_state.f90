!
! Written by Lorenzo Paulatto (2014-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE final_state
  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : ph_system_info, forceconst2_grid
  USE fc2_interpolate, ONLY : ip_cart2pat
  USE fc3_interpolate, ONLY : forceconst3
  USE mpi_thermal,     ONLY : mpi_bsum, ionode
#include "mpi_thermal.h"
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  ! Full spectral function, computed as in eq. 1 of arXiv:1312.7467v1
  FUNCTION final_state_q(xq0, qpath, nconf, T, sigma, S, grid, fc2, fc3,&
                         ei, ne, ener, qresolved, qsummed, sigmaq, outdir, prefix)
    USE fc2_interpolate,     ONLY : bose_phq, freq_phq_safe
    USE linewidth,      ONLY : sum_selfnrg_modes
    USE q_grids,        ONLY : q_grid
    USE functions,      ONLY : refold_bz, refold_bz_mod, f_gauss
    USE constants,      ONLY : RY_TO_CMM1
    USE more_constants, ONLY : write_conf
    USE nanoclock,      ONLY : print_percent_wall
    USE timers
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(q_grid),INTENT(in) :: qpath
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    REAL(DP),INTENT(in) :: ei       ! energy to examine (cm^-1)
    !
    INTEGER,INTENT(in)  :: ne       ! number of final state energies
    REAL(DP),INTENT(in) :: ener(ne) ! the final state energies
    REAL(DP),INTENT(in) :: sigmaq ! gaussian width used to select the q points
    LOGICAL,INTENT(in) :: qresolved, qsummed
    CHARACTER(len=256),INTENT(in) :: outdir, prefix
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    ! FUNCTION RESULT:
    REAL(DP) :: final_state_q(ne,S%nat3,nconf)
    !
    INTEGER :: unit
    INTEGER, EXTERNAL :: find_free_unit
    ! To interpolate D2 and D3:
    INTEGER :: iq, jq, nu, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie
    REAL(DP) :: gamma, delta, omega, denom
    COMPLEX(DP),ALLOCATABLE :: fstate_q(:,:,:)
    !
    REAL(DP) :: sumaux(ne,S%nat3)
    INTEGER     :: iqpath
    REAL(DP),ALLOCATABLE    :: xqbar(:,:,:,:)
    REAL(DP),ALLOCATABLE    :: xqsum(:,:,:)
    REAL(DP) :: qbarweight
    REAL(DP) :: pl,dpl
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(fstate_q(ne,S%nat3,nconf))
    !
    fstate_q = (0._dp, 0._dp)
    !
    IF(qresolved .or. qsummed)THEN
      ioWRITE(*,'(2x,a,3f12.4,a,f12.6,a)') "Q-resolved final state", xq0, &
                                       " energy = ", ei*RY_TO_CMM1, "cm^-1"
    ENDIF
    !
    IF(qresolved)THEN
      timer_CALL t_qresolved_io%start()
      ALLOCATE(xqbar(ne,S%nat3,qpath%nq,nconf), stat=i)
      IF(i/=0)THEN
        ioWRITE(*,'(///,5x,a,f6.3,a,//,5x,a)') "Trying to allocate: ", &
          DBLE(ne)*DBLE(S%nat3)*DBLE(qpath%nq)*DBLE(nconf)*DBLE(DP)/2._dp**30, &
          " GB", "Use less configurations/q-points/energies "&
               //"to do q-resolved final state."
        CALL errore("final_state_q","out of memory",1)
      ENDIF
      xqbar = 0._dp
      timer_CALL t_qresolved_io%stop()
    ENDIF
    !
    IF(qsummed) THEN
      ALLOCATE(xqsum(S%nat3,qpath%nq,nconf))
      xqsum = 0._dp
    ENDIF
    !
    ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q1
    xq(:,1) = xq0
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    GRID_LOOP : &
    DO iq = 1, grid%nq
        timer_CALL t_spf%start()
      !
      CALL print_percent_wall(10._dp, 300._dp, iq, grid%nq, (iq==1))
      !
      xq(:,2) = grid%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
      DO jq = 2,3
        CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
      ENDDO
!/nope/!$OMP END PARALLEL DO
      !
      ! ------ start of CALL scatter_3q(S,fc2,fc3, xq(:,1),xq(:,2),xq(:,3), V3sq)
        timer_CALL t_fc3int%start()
      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
        timer_CALL t_fc3int%stop()
        timer_CALL t_fc3rot%start()
      CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
        timer_CALL t_fc3rot%stop()
        timer_CALL t_fc3m2%start()
      V3sq = REAL( CONJG(D3)*D3 , kind=DP)
        timer_CALL t_fc3m2%stop()
      !
      CONF_LOOP : &
      DO it = 1,nconf
        ! Compute eigenvalues, eigenmodes and bose-einstein occupation at q2 and q3
!/nope/!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(jq)
          timer_CALL t_bose%start()
        DO jq = 1,3
          CALL bose_phq(T(it),S%nat3, freq(:,jq), bose(:,jq))
        ENDDO
!/nope/!$OMP END PARALLEL DO
          timer_CALL t_bose%stop()
          timer_CALL t_sum%start()
        sumaux = grid%w(iq) * &
            sum_final_state_e( S, sigma(it), T(it),freq, bose, V3sq, ei, ne, ener )
        fstate_q(:,:,it) = fstate_q(:,:,it) + sumaux
          timer_CALL t_sum%stop()
        !
        IF(qresolved) THEN
          timer_CALL t_qresolved%start()
          DO iqpath = 1,qpath%nq
            qbarweight = & !qpath%w(iqpath)*&
              f_gauss(refold_bz_mod(qpath%xq(:,iqpath)-xq(:,2), S%bg), sigmaq)
            xqbar(:,:,iqpath,it) = xqbar(:,:,iqpath,it) - 0.5_dp * sumaux *  qbarweight
          ENDDO
          timer_CALL t_qresolved%stop()
        ENDIF
        !
        IF(qsummed)THEN
          timer_CALL t_qsummed%start()
          IF(qresolved)THEN
            DO iqpath = 1,qpath%nq
            DO nu = 1,S%nat3
            DO ie = 1,ne
              xqsum(nu,iqpath,it) = xqsum(nu,iqpath,it) + xqbar(ie,nu,iqpath,it)
            ENDDO
            ENDDO
            ENDDO
          ELSE
            DO iqpath = 1,qpath%nq
            qbarweight = & !qpath%w(iqpath)*&
              f_gauss(refold_bz_mod(qpath%xq(:,iqpath)-xq(:,2), S%bg), sigmaq)
              DO nu = 1,S%nat3
                xqsum(nu,iqpath,it) = xqsum(nu,iqpath,it) &
                              - 0.5_dp * SUM(sumaux(:,nu)) *  qbarweight
              ENDDO
            ENDDO
          ENDIF
          timer_CALL t_qsummed%stop()
        ENDIF
        !
      ENDDO &
      CONF_LOOP
      !
        timer_CALL t_spf%stop()
    ENDDO &
    GRID_LOOP
    !
    IF(qresolved)THEN
      timer_CALL t_qresolved_io%start()
      CALL mpi_bsum(ne,S%nat3,qpath%nq, nconf, xqbar)
      !
      ! To avoid messy graphs because many softwares do not understand 1.00-120 (i.e. 10^{-120})
      WHERE(ABS(xqbar)<1.d-99) xqbar = 0._dp
      IF(ionode)THEN
      unit = find_free_unit()
      DO it = 1,nconf
        OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                        "_qresolved_T"//TRIM(write_conf(it,nconf,T))//&
                        "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".out")
        WRITE(unit,'(a)') "# ener  path_l/q_weight  q_x q_y q_z total band1 band2 ..."
        DO iqpath = 1,qpath%nq
          DO ie = 1,ne
            WRITE(unit,'(100e15.5)') qpath%w(iqpath), &
                            ener(ie)*RY_TO_CMM1, qpath%xq(:,iqpath), &
                            SUM(xqbar(ie,:,iqpath,it)),  xqbar(ie,:,iqpath,it)
          ENDDO
          ioWRITE(unit,*) 
        ENDDO
        CLOSE(unit)
        !
      ENDDO
      ENDIF
      DEALLOCATE(xqbar)
      timer_CALL t_qresolved_io%stop()
    ENDIF
    !
    IF(qsummed)THEN
      timer_CALL t_qsummed_io%start()
      CALL mpi_bsum(S%nat3,qpath%nq, nconf, xqsum)
      !
      WHERE(ABS(xqsum)<1.d-99) xqsum = 0._dp
      !
      IF(ionode)THEN
      unit = find_free_unit()
      DO it = 1,nconf
        OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                        "_qsum_T"//TRIM(write_conf(it,nconf,T))//&
                        "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".out")
        WRITE(unit,'(a)') "# path_l/q_weight  q_x q_y q_z total band1 band2 ..."
        DO iqpath = 1,qpath%nq
          WRITE(unit,'(100e15.5)') &
                          qpath%w(iqpath), qpath%xq(:,iqpath), &
                          SUM(xqsum(:,iqpath,it)),  xqsum(:,iqpath,it)
        ENDDO
        CLOSE(unit)
      ENDDO
      ENDIF
      DEALLOCATE(xqsum)
      timer_CALL t_qsummed_io%stop()
    ENDIF
        
        
    IF(grid%scattered) CALL mpi_bsum(ne,S%nat3,nconf, fstate_q)
    final_state_q = -0.5_dp * fstate_q
    !
    DEALLOCATE(U, V3sq, D3, fstate_q)
    !
  END FUNCTION final_state_q
  !
  ! Sum the self energy at the provided ener(ne) input energies
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_final_state_e(S, sigma, T, freq, bose, V3sq, ei, ne, ener)
    USE functions,        ONLY : f_gauss, df_bose
    USE merge_degenerate, ONLY : merge_degen
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma, T   ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
    !
    REAL(DP),INTENT(in) :: ei ! the energy of the state under scrutiny
    !
    INTEGER,INTENT(in)  :: ne           ! number of energies on which to decompose the final state
    REAL(DP),INTENT(in) :: ener(ne)     ! the energies
    REAL(DP)            :: de
    !
    ! _P -> scattering, _M -> cohalescence
    REAL(DP) :: bose_P, bose_M      ! final/initial state populations 
    REAL(DP) :: freqtot, freqtotm1
    REAL(DP) :: omega_P,  omega_M   ! \delta\omega
    REAL(DP) :: omega_P2, omega_M2  ! \delta\omega
    REAL(DP) :: wfinal(ne)
    COMPLEX(DP) :: ctm_P, ctm_M, reg, num
    !
    INTEGER :: i,j,k, ie
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    REAL(DP) :: sum_final_state_e(ne,S%nat3)
    COMPLEX(DP),ALLOCATABLE :: fsdf(:,:) ! final state decay function
    !
    ALLOCATE(fsdf(ne,S%nat3))
    fsdf = (0._dp, 0._dp)
    !
    de = ABS(ener(2)-ener(ne))/100._dp ! FIXME
    !
    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(ie,i,j,k,bose_P,bose_M,omega_P,omega_M,omega_P2,omega_M2,&
!$OMP                     ctm_P,ctm_M,reg,freqtot,freqtotm1,wfinal) &
!$OMP             REDUCTION(+: fsdf) COLLAPSE(2)
    DO k = 1,S%nat3
      DO j = 1,S%nat3
        !
        FORALL(i=1:ne) wfinal(i) = f_gauss((ener(i)-freq(j,2)), de)
        !

        bose_P   = 1 + bose(j,2) + bose(k,3)
        omega_P  = freq(j,2)+freq(k,3)
        !
        bose_M   = bose(k,3) - bose(j,2)
        omega_M  = freq(j,2)-freq(k,3)
        !
        IF(sigma>0._dp)THEN
          omega_P2 = omega_P**2
          omega_M2 = omega_M**2
        ELSE IF(sigma<0._dp)THEN
          ctm_P = 2 * bose_P *omega_P/(omega_P**2+sigma**2)
          ctm_M = 2 * bose_M *omega_M/(omega_M**2+sigma**2)
        ELSE !IF (sigma==0._dp)THEN
          !
          IF(omega_P>0._dp)THEN
            ctm_P = 2 * bose_P /omega_P
          ELSE
            ctm_P = 0._dp
          ENDIF
          !
          IF(ABS(omega_M)>1.e-2_dp)THEN
            ctm_M = 2 * bose_M /omega_M
          ELSE
            IF(T>0._dp)THEN
              ctm_M = -2* df_bose(0.5_dp * omega_P, T)
            ELSE
              ctm_M = 0._dp
            ENDIF
          ENDIF
          !
        ENDIF
        !bose_P   = 1 + bose(j,2) + bose(k,3)
        !omega_P  = freq(j,2)+freq(k,3)
        !omega_P2 = omega_P**2
        !
        !bose_M   = bose(k,3) - bose(j,2)
        !omega_M  = freq(j,2)-freq(k,3)
        !omega_M2 = omega_M**2
        !
        DO i = 1,S%nat3
!           i = 6
          !
          ! This comes from the definition of u_qj, Ref. 1. in linewidth.f90
          freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
          !
          IF (freqtot/=0._dp) THEN
            freqtotm1 = 1 / freqtot
            !
            DO ie = 1, ne
            ! regularization:
              !
              ! regularization:
              IF(sigma>0._dp)THEN
                reg = CMPLX(ei, sigma, kind=DP)**2
                !reg = CMPLX(freq(i,1), sigma, kind=DP)**2
                ctm_P = 2 * bose_P *omega_P/(omega_P2-reg )
                ctm_M = 2 * bose_M *omega_M/(omega_M2-reg )
              ENDIF

              !ctm_P = 2 * bose_P *omega_P/(omega_P2-reg)
              !ctm_M = 2 * bose_M *omega_M/(omega_M2-reg)
              !
              fsdf(ie,i) = fsdf(ie,i) + wfinal(ie)*(ctm_P + ctm_M) * V3sq(i,j,k) * freqtotm1
            ENDDO
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    CALL merge_degen(ne, S%nat3, fsdf, freq(:,1))
    !
    sum_final_state_e = -DIMAG(fsdf)
    !sum_final_state_e = REAL(fsdf)
    DEALLOCATE(fsdf)
    !
  END FUNCTION sum_final_state_e
  
END MODULE final_state 
