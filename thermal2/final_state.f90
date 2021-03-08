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
  USE mpi_thermal,     ONLY : mpi_bsum, ionode, mpi_wbarrier
#include "mpi_thermal.h"
  INTEGER,PARAMETER :: TOT=1, C=2, X=3, LW=2, LS=3, NTERMS=3
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  FUNCTION final_state_q(xq0, qpath, nconf, T, sigma, S, grid, fc2, fc3,&
                         nui, ei, ne, ener, sigma_e, qresolved, qsummed, sigmaq, outdir, prefix)
    USE fc2_interpolate,      ONLY : bose_phq, freq_phq_safe, set_nu0
    USE linewidth,            ONLY : sum_selfnrg_modes
    USE q_grids,              ONLY : q_grid
    USE functions,            ONLY : refold_bz, refold_bz_mod, f_gauss, quicksort
    USE constants,            ONLY : RY_TO_CMM1
    USE more_constants,       ONLY : write_conf
    USE timers
    !
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq0(3)
    TYPE(q_grid),INTENT(in) :: qpath
    INTEGER,INTENT(in)  :: nconf
    REAL(DP),INTENT(in) :: T(nconf)     ! Kelvin
    !
    INTEGER,INTENT(in)  :: nui      ! mode of initial energy, 0=get the sum of all modes
    REAL(DP),INTENT(in) :: ei       ! energy to examine (cm^-1), if <0 and nui>0 get omega_nui
    !
    INTEGER,INTENT(in)  :: ne       ! number of final state energies
    REAL(DP),INTENT(in) :: ener(ne) ! the final state energies
    REAL(DP),INTENT(in) :: sigma_e  ! gaussian width to smooth plots in energy
    REAL(DP),INTENT(in) :: sigmaq ! gaussian width used to select the q points
    LOGICAL,INTENT(in)  :: qresolved, qsummed
    CHARACTER(len=256),INTENT(in) :: outdir, prefix
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: grid
    REAL(DP),INTENT(in) :: sigma(nconf) ! ry
    ! FUNCTION RESULT:
    REAL(DP) :: final_state_q(ne,S%nat3,NTERMS,nconf)
    !
    INTEGER :: unit
    INTEGER, EXTERNAL :: find_free_unit
    ! To interpolate D2 and D3:
    INTEGER :: iq, jq, nu, it, nu0(3)
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    REAL(DP),ALLOCATABLE    :: V3sq(:,:,:)
    REAL(DP) :: freq(S%nat3,3), bose(S%nat3,3), xq(3,3)
    !
    ! To compute the spectral function from the self energy:
    INTEGER  :: i, ie
    REAL(DP) :: gamma, delta, omega, denom
    REAL(DP),ALLOCATABLE :: fstate_q(:,:,:,:)
    !
    REAL(DP) :: sumaux(ne,S%nat3,NTERMS)
    INTEGER     :: iqpath, ibnd
    REAL(DP),ALLOCATABLE    :: xqbar(:,:,:,:)
    REAL(DP),ALLOCATABLE    :: xqsum(:,:,:), xqsumsum(:,:), aux_sumsum(:)
    LOGICAL :: found(4)
    REAL(DP) :: sum_aux, aux
    REAL(DP) :: qbarweight
    REAL(DP) :: pl,dpl
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    ALLOCATE(U(S%nat3, S%nat3,3))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
    ALLOCATE(fstate_q(ne,S%nat3,NTERMS,nconf))
    !
    fstate_q = 0._dp
    !
!     IF(qresolved .or. qsummed)THEN
!       ioWRITE(*,'(2x,a,3f12.4,a,f12.6,a)') "Q-resolved final state", xq0, &
!                                        " energy = ", ei*RY_TO_CMM1, "cm^-1"
!     ENDIF
    !
    IF(qresolved)THEN
      timer_CALL t_qresolved_io%start()
      ALLOCATE(xqbar(ne,S%nat3,qpath%nq,nconf), stat=i)
      IF(i/=0)THEN
        ioWRITE(*,'(///,5x,a,f6.3,a,//,5x,a)') "Trying to allocate: ", &
          DBLE(ne)*DBLE(S%nat3)*DBLE(qpath%nq)*DBLE(nconf)*DBLE(DP)/2._dp**30, &
          " GB", "Use less configurations/q-points/energies "&
               //"to do q-resolved final state!"
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
    nu0(1) = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc2, freq(:,1), U(:,:,1))
    !
    ! Next check moved to sum_final_state_e
!     IF(ei<0._dp)THEN
!       IF(nui>0 .and. nui<=S%nat3) ei=freq(nui,1)
!       IF(ei<0) CALL errore("final_state_q","you must set initial mode or initial energy (or both)",1)
!     ENDIF
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
        nu0(jq) = set_nu0(xq(:,jq), S%at)
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
            sum_final_state_e( S, sigma(it), T(it),freq, bose, V3sq, &
                               nui, ei, ne, ener, sigma_e, nu0 )
        !
        fstate_q(:,:,:,it) = fstate_q(:,:,:,it) + sumaux(:,:,:)
          timer_CALL t_sum%stop()
        !
        IF(qresolved) THEN
          timer_CALL t_qresolved%start()
          DO iqpath = 1,qpath%nq
            qbarweight = & !qpath%w(iqpath)*&
              f_gauss(refold_bz_mod(qpath%xq(:,iqpath)-xq(:,3), S%bg), sigmaq)
            xqbar(:,:,iqpath,it) = xqbar(:,:,iqpath,it) - 0.5_dp * sumaux(:,:,TOT) *  qbarweight
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
              f_gauss(refold_bz_mod(qpath%xq(:,iqpath)-xq(:,3), S%bg), sigmaq)
              DO nu = 1,S%nat3
                xqsum(nu,iqpath,it) = xqsum(nu,iqpath,it) &
                              - 0.5_dp * SUM(sumaux(:,nu,TOT)) *  qbarweight
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
      !WHERE(ABS(xqbar)<1.d-99) xqbar = 0._dp
      IF(ionode)THEN
      unit = find_free_unit()
      DO it = 1,nconf
        OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                        "_qresolved_T"//TRIM(write_conf(it,nconf,T))//&
                        "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".out")
        WRITE(unit,'(a)') "# ener  path_l/q_weight  q_x q_y q_z total band1 band2 ..."
        DO iqpath = 1,qpath%nq
          DO ie = 1,ne
            WRITE(unit,'(100ES27.15E3)') qpath%w(iqpath), &
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
    QSUMMED_IO : &
    IF(qsummed)THEN
      timer_CALL t_qsummed_io%start()
      CALL mpi_bsum(S%nat3,qpath%nq, nconf, xqsum)
      !
      !WHERE(ABS(xqsum)<1.d-99) xqsum = 0._dp
      ALLOCATE(xqsumsum(qpath%nq,nconf))
      DO it = 1,nconf
      DO iqpath = 1,qpath%nq
        xqsumsum(iqpath,it) = SUM(xqsum(:,iqpath,it))
      ENDDO
      ENDDO
      !
      IONODE_IO : &
      IF(ionode)THEN
        unit = find_free_unit()
        !
        XCRYSDEN_ELSE : &
        IF(qpath%type == 'bxsf')THEN
          !
          DO it = 1,nconf
            ! Find the value which defines an isovalue-surface that has 25%, 50% 
            ! 75% and 90% of contribution to linewidth inside
            !
            ALLOCATE(aux_sumsum(qpath%nq))
              aux_sumsum(:) = xqsumsum(:,it)
              sum_aux = SUM(aux_sumsum)
              CALL quicksort(aux_sumsum,1, qpath%nq)
              aux = 0._dp
              found(:) = .false.
              WRITE(*,"(a)") "Scanning for isovalues:"
              WRITE(*,"(a)") "% Gamma, isovalue, % points inside:"
              SEEK_90_PERCENT : &
              DO iqpath = qpath%nq,1,-1
                aux = aux + aux_sumsum(iqpath)
                IF(aux>.25_dp*sum_aux .and. .not. found(1)) THEN
                  found(1)=.true.
                  WRITE(*,'(4x,a,e12.3,f12.6)') "25%", aux_sumsum(iqpath), 100*(1-DBLE(iqpath)/qpath%nq)
                ENDIF
                IF(aux>.50_dp*sum_aux .and. .not. found(2)) THEN
                  found(2)=.true.
                  WRITE(*,'(4x,a,e12.3,f12.6)') "50%", aux_sumsum(iqpath), 100*(1-DBLE(iqpath)/qpath%nq)
                ENDIF
                IF(aux>.75_dp*sum_aux .and. .not. found(3)) THEN
                  found(3)=.true.
                  WRITE(*,'(4x,a,e12.3,f12.6)') "75%", aux_sumsum(iqpath), 100*(1-DBLE(iqpath)/qpath%nq)
                ENDIF
                IF(aux>.90_dp*sum_aux .and. .not. found(4)) THEN
                  found(4)=.true.
                  WRITE(*,'(4x,a,e12.3,f12.6)') "90%", aux_sumsum(iqpath), 100*(1-DBLE(iqpath)/qpath%nq)
                ENDIF
                IF(all(found)) EXIT SEEK_90_PERCENT
              ENDDO SEEK_90_PERCENT 
              WRITE(*,*)
            DEALLOCATE(aux_sumsum)
            !
            OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                            "_qsum_T"//TRIM(write_conf(it,nconf,T))//&
                            "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".bxsf")
            WRITE(unit,*) "BEGIN_INFO"
            WRITE(unit,'(x,a,1e12.3)') "Fermi Energy:", MAXVAL(xqsum(:,:,it))*.9_dp
            WRITE(unit,*) "END_INFO"
            WRITE(unit,*) "BEGIN_BLOCK_BANDGRID_3D"
            WRITE(unit,*) "final_state"
            WRITE(unit,*) "BEGIN_BANDGRID_3D_BANDS"
            WRITE(unit,*) S%nat3+2
            WRITE(unit,*) qpath%n
            WRITE(unit,*) qpath%xq0
            WRITE(unit,*) S%bg(:,1)
            WRITE(unit,*) S%bg(:,2)
            WRITE(unit,*) S%bg(:,3)
            DO ibnd = 1, S%nat3
              WRITE(unit,*) "BAND:", ibnd
              DO iqpath = 1,qpath%nq
                WRITE(unit,'(1ES27.15E3)') xqsum(ibnd,iqpath,it)
              ENDDO
            ENDDO
            WRITE(unit,*) "BAND:", S%nat3+1
            DO iqpath = 1,qpath%nq
              WRITE(unit,'(1ES27.15E3)') xqsumsum(iqpath,it)
            ENDDO
            WRITE(unit,*) "BAND:", S%nat3+2
            DO iqpath = 1,qpath%nq
              qbarweight = & !qpath%w(iqpath)*&
                f_gauss(refold_bz_mod(qpath%xq(:,iqpath)-xq0, S%bg), sigmaq)
              WRITE(unit,'(1ES27.15E3)') qbarweight 
            ENDDO
            WRITE(unit,*) "END_BANDGRID_3D"
            WRITE(unit,*) "END_BLOCK_BANDGRID_3D"
            !
            CLOSE(unit)
            !
          ENDDO
          !
        ELSE IF(qpath%type == 'xsf')THEN XCRYSDEN_ELSE
          !
          DO it = 1,nconf
            OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                            "_qsum_T"//TRIM(write_conf(it,nconf,T))//&
                            "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".xsf")
            WRITE(unit,'(a)') "CRYSTAL"
            WRITE(unit,'(a)') "PRIMVEC"
            WRITE(unit,'(3f12.6)') S%bg(:,1)
            WRITE(unit,'(3f12.6)') S%bg(:,2)
            WRITE(unit,'(3f12.6)') S%bg(:,3)
            WRITE(unit,'(a)') "PRIMCOORD"
            WRITE(unit,'(2i2)') 1, 1
            WRITE(unit,'(i4,3f12.6)') 1, xq0

            WRITE(unit,'(a)') "BEGIN_BLOCK_DATAGRID_3D"
            WRITE(unit,'(a)') "  final_state_conf_"//int_to_char(it)
            DO ibnd = 0, S%nat3
            !WRITE(unit,'(a)') "DATAGRID_3D_UNKNOWN"
            IF(ibnd==0)THEN
               WRITE(unit,'(a)') "  DATAGRID_3D_final_state_sum"
            ELSE
               WRITE(unit,'(a)') "  DATAGRID_3D_final_state_"//int_to_char(ibnd)
            ENDIF
            WRITE(unit,'(3i4)') qpath%n
            WRITE(unit,'(3f12.6)') 0.0, 0.0, 0.0
            WRITE(unit,'(3f12.6)') S%bg(:,1)
            WRITE(unit,'(3f12.6)') S%bg(:,2)
            WRITE(unit,'(3f12.6)') S%bg(:,3)
            DO iqpath = 1,qpath%nq
              IF(ibnd==0)THEN
                  WRITE(unit,'(1ES27.15E3)') SUM(xqsum(:,iqpath,it))
              ELSE
                  WRITE(unit,'(1ES27.15E3)') xqsum(ibnd,iqpath,it)
              ENDIF
              IF(mod(iqpath, qpath%n(1)*qpath%n(2))==0 .and. iqpath<qpath%nq) WRITE(unit, *) 
            ENDDO
            WRITE(unit,'(a)') "  END_DATAGRID_3D"
            ENDDO
            WRITE(unit,'(a)') "END_BLOCK_DATAGRID_3D"
            !
            CLOSE(unit)
            !
          ENDDO
          !
        ELSE XCRYSDEN_ELSE
          !
          DO it = 1,nconf
            OPEN(unit, file=TRIM(outdir)//"/"//TRIM(prefix)//&
                            "_qsum_T"//TRIM(write_conf(it,nconf,T))//&
                            "_s"//TRIM(write_conf(it,nconf,sigma*RY_TO_CMM1))//".out")
            WRITE(unit,'(a)') "# path_l/q_weight  q_x q_y q_z total band1 band2 ..."
            DO iqpath = 1,qpath%nq
              WRITE(unit,'(100ES27.15E3)') qpath%w(iqpath), qpath%xq(:,iqpath), &
                                      xqsumsum(iqpath,it),  xqsum(:,iqpath,it)
            ENDDO
            CLOSE(unit)
          ENDDO
        ENDIF &
        XCRYSDEN_ELSE
        !
      ENDIF &
      IONODE_IO
      !
      DEALLOCATE(xqsum, xqsumsum)
      CALL mpi_wbarrier()
      timer_CALL t_qsummed_io%stop()
      !
    ENDIF &
    QSUMMED_IO
    !
    IF(grid%scattered) CALL mpi_bsum(ne,S%nat3,NTERMS,nconf, fstate_q)
    final_state_q = -0.5_dp * fstate_q !-0.5_dp * fstate_q
    !
    DEALLOCATE(U, V3sq, D3, fstate_q)
    !
  END FUNCTION final_state_q
  !
  ! Sum the self energy at the provided ener(ne) input energies
  ! \/o\________\\\_________________________________________/^>
!   FUNCTION sum_final_iselfnrg_e(S, sigma, T, freq, bose, V3sq, nu_initial, e_initial, ne, ener, sigma_e, nu0)
! !  FUNCTION sum_final_state_e(S, sigma, T, freq, bose, V3sq, e_initial, ne, ener, sigma_e, nu0)
!     USE functions,        ONLY : f_gauss, df_bose
!     USE merge_degenerate, ONLY : merge_degen
!     IMPLICIT NONE
!     TYPE(ph_system_info),INTENT(in)   :: S
!     REAL(DP),INTENT(in) :: sigma, T   ! smearing (regularization) (Ry)
!     REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
!     REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
!     REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
!     !
!     REAL(DP),INTENT(in) :: e_initial ! the energy of the state under scrutiny
!     !
!     INTEGER,INTENT(in)  :: ne           ! number of energies on which to decompose the final state
!     REAL(DP),INTENT(in) :: ener(ne)     ! the energies
!     REAL(DP)            :: de
!     REAL(DP),INTENT(in) :: sigma_e      ! smearing for the energy axis
!     INTEGER,INTENT(in)  :: nu0(3)
!     !
!     ! _X -> scattering, _C -> cohalescence
!     REAL(DP) :: bose_X, bose_C      ! final/initial state populations 
!     REAL(DP) :: freqtot, freqtotm1
!     REAL(DP) :: omega_X,  omega_C   ! \delta\omega
!     REAL(DP) :: omega_X2, omega_C2  ! \delta\omega
!     REAL(DP) :: wfinal3(ne)
!     COMPLEX(DP) :: ctm_X, ctm_C, reg, num
!     !
!     INTEGER :: i,j,k, ie
!     !
!     ! Note: using the function result in an OMP reduction causes crash with ifort 14
!     REAL(DP) :: sum_final_iselfnrg_e(ne,S%nat3,NTERMS)
!     COMPLEX(DP),ALLOCATABLE :: fsdf(:,:) ! final state decay function
!     !
!     ALLOCATE(fsdf(ne,S%nat3))
!     fsdf = (0._dp, 0._dp)
!     !
!     de = ABS(ener(2)-ener(ne))/100._dp ! FIXME
!     !
!     !
! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP             PRIVATE(ie,i,j,k,bose_X,bose_C,omega_X,omega_C,omega_X2,omega_C2,&
! !$OMP                     ctm_X,ctm_C,reg,freqtot,freqtotm1,wfinal3) &
! !$OMP             REDUCTION(+: fsdf) COLLAPSE(2)
!     DO k = 1,S%nat3
!       DO j = 1,S%nat3
!         !
!         FORALL(i=1:ne) wfinal3(i) = f_gauss((ener(i)-freq(j,2)), de)
!         !FORALL(i=1:ne) wfinal3(i) = f_gauss((ener(i)-freq(j,3)), de)
!         !
! 
!         bose_X   = 1 + bose(j,2) + bose(k,3)
!         omega_X  = freq(j,2)+freq(k,3)
!         !
!         bose_C   = bose(k,3) - bose(j,2)
!         omega_C  = freq(j,2)-freq(k,3)
!         !
!         IF(sigma>0._dp)THEN
!           omega_X2 = omega_X**2
!           omega_C2 = omega_C**2
!         ELSE IF(sigma<0._dp)THEN
!           ctm_X = 2 * bose_X *omega_X/(omega_X**2+sigma**2)
!           ctm_C = 2 * bose_C *omega_C/(omega_C**2+sigma**2)
!         ELSE !IF (sigma==0._dp)THEN
!           !
!           IF(omega_X>0._dp)THEN
!             ctm_X = 2 * bose_X /omega_X
!           ELSE
!             ctm_X = 0._dp
!           ENDIF
!           !
!           IF(ABS(omega_C)>1.e-2_dp)THEN
!             ctm_C = 2 * bose_C /omega_C
!           ELSE
!             IF(T>0._dp)THEN
!               ctm_C = -2* df_bose(0.5_dp * omega_X, T)
!             ELSE
!               ctm_C = 0._dp
!             ENDIF
!           ENDIF
!           !
!         ENDIF
!         !bose_X   = 1 + bose(j,2) + bose(k,3)
!         !omega_X  = freq(j,2)+freq(k,3)
!         !omega_X2 = omega_X**2
!         !
!         !bose_C   = bose(k,3) - bose(j,2)
!         !omega_C  = freq(j,2)-freq(k,3)
!         !omega_C2 = omega_C**2
!         !
!         DO i = 1,S%nat3
! !           i = 6
!           !
!           ! This comes from the definition of u_qj, Ref. 1. in linewidth.f90
!           freqtot = 8*freq(i,1)*freq(j,2)*freq(k,3)
!           !
!           IF (freqtot/=0._dp) THEN
!             freqtotm1 = 1 / freqtot
!             !
!             DO ie = 1, ne
!             ! regularization:
!               !
!               ! regularization:
!               IF(sigma>0._dp)THEN
!                 reg = CMPLX(e_initial, sigma, kind=DP)**2
!                 !reg = CMPLX(freq(i,1), sigma, kind=DP)**2
!                 ctm_X = 2 * bose_X *omega_X/(omega_X2-reg )
!                 ctm_C = 2 * bose_C *omega_C/(omega_C2-reg )
!               ENDIF
! 
!               !ctm_X = 2 * bose_X *omega_X/(omega_X2-reg)
!               !ctm_C = 2 * bose_C *omega_C/(omega_C2-reg)
!               !
!               fsdf(ie,i) = fsdf(ie,i) + wfinal3(ie)*(ctm_X + ctm_C) * V3sq(i,j,k) * freqtotm1
!             ENDDO
!           ENDIF
!           !
!         ENDDO
!       ENDDO
!     ENDDO
! !$OMP END PARALLEL DO
!     !
! !    CALL merge_degen(ne, S%nat3, fsdf, freq(:,1))
!     !
!     sum_final_iselfnrg_e(:,:,LW) = -DIMAG(fsdf)
!     sum_final_iselfnrg_e(:,:,LS) = DBLE(fsdf)
!     !sum_final_iselfnrg_e = -DIMAG(fsdf)
!     !sum_final_iselfnrg_e = REAL(fsdf)
!     DEALLOCATE(fsdf)
!     !
!   END FUNCTION sum_final_iselfnrg_e

  !
  ! Compute the decomposition of gamma over the final states. 
  ! If nu_initial is provided, then it will compute d\gamma_{nu_initial} and will 
  ! set e_initial to omega_{nu_initial}, unless a different value of e_initial is required
  ! In this case the different columns of the return array will correspond to the different
  ! final states
  !
  ! If nu_initial is set to zero, it will return \sum_nu d\gamma_nu
  ! In this case, the different columns will correspond to \gamma_{1..nat3}
  !
  ! \/o\________\\\_________________________________________/^>
  FUNCTION sum_final_state_e(S, sigma, T, freq, bose, V3sq, nu_initial, e_initial, &
                             ne, ener, sigma_e, nu0)
    USE functions,        ONLY : f_gauss, df_bose
    USE merge_degenerate, ONLY : merge_degen
    USE constants,        ONLY : pi, RY_TO_CMM1
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    REAL(DP),INTENT(in) :: sigma, T   ! smearing (regularization) (Ry)
    REAL(DP),INTENT(in) :: freq(S%nat3,3)  ! phonon energies (Ry)
    REAL(DP),INTENT(in) :: bose(S%nat3,3)  ! bose/einstein distribution of freq
    REAL(DP),INTENT(in) :: V3sq(S%nat3,S%nat3,S%nat3) ! |D^3|**2 on the basis of phonons patterns
    !
    REAL(DP),INTENT(in) :: e_initial ! the initial energy
    !
    INTEGER,INTENT(in)  :: ne           ! number of energies on which to decompose the final state
    REAL(DP),INTENT(in) :: ener(ne)     ! the energies
    REAL(DP),INTENT(in) :: sigma_e      ! smearing for the energy axis
    INTEGER,INTENT(in)  :: nu0(3)
    INTEGER,INTENT(in),OPTIONAL  :: nu_initial
    !
    ! _X -> scattering, _C -> cohalescence
    REAL(DP) :: bose_X, bose_C      ! final/initial state populations 
    REAL(DP) :: freqtot, freqtotm1, freqtotm1_23
    REAL(DP) :: omega_X,  omega_C   ! \delta\omega
    REAL(DP) :: dom_C, dom_X   ! \delta\omega
    REAL(DP) :: ctm_C, ctm_X   !
    REAL(DP) :: pref_C, pref_X
    REAL(DP) :: wfinal3(ne), wfinal2(ne), wfinal(ne)
    REAL(DP) :: freqm1(S%nat3,3), freq_initial
    !
    INTEGER,TARGET :: i,j,k
    INTEGER :: ie, start_i, end_i
    INTEGER,POINTER :: i_final
    !
    ! Note: using the function result in an OMP reduction causes crash with ifort 14
    REAL(DP) :: sum_final_state_e(ne,S%nat3,3)
    REAL(DP),ALLOCATABLE :: fsdf(:,:,:) ! final state decay function
    !
    IF(nu_initial>0 .and. nu_initial<=S%nat3) THEN
      start_i = nu_initial
      end_i   = nu_initial
      i_final => j
      IF(e_initial<0._dp) THEN
        freq_initial = freq(nu_initial,1)
!        print*, freq_initial*RY_TO_CMM1
      ELSE
        freq_initial = e_initial
      ENDIF
    ELSE
      start_i = 1
      end_i   = S%nat3
      i_final => i
      freq_initial = e_initial
      IF(freq_initial<0._dp) CALL errore("sum_final_state_e","nu_initial was unset" &
                                  //" or invalid and e_initial was not specified", 1)
    ENDIF
    !
    !
    ALLOCATE(fsdf(ne,S%nat3,3))
    fsdf = (0._dp, 0._dp)
    !
    freqm1 = 0._dp
    DO i = 1,S%nat3
      IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
      IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
      IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
    ENDDO    !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ie,i,j,k,bose_X,bose_C,omega_X,omega_C,dom_X,dom_C,pref_X,pref_C,&
!$OMP         ctm_X,ctm_C,freqtot,freqtotm1,FREQTOTM1_23,wfinal3) &
!$OMP REDUCTION(+: fsdf) COLLAPSE(1)
    DO k = 1,S%nat3
      FORALL(ie=1:ne) wfinal3(ie) = f_gauss((ener(ie)-freq(k,3)), sigma_e)
      DO j = 1,S%nat3
        !
        !FORALL(ie=1:ne) wfinal2(ie) = f_gauss((ener(ie)-freq(j,2)), sigma_e)
        !
        ! I define a final energy for Coalescence process which is symmetric in the
        ! exchange of q2 with q3. 
        ! The scattering (X) part is already symmetric, I could just take wfinal3
        wfinal = wfinal3 !(wfinal2+wfinal3)/2
        !wfinal = (wfinal2+wfinal3)/2
        !
        bose_C = 2* (bose(j,2) - bose(k,3))
        bose_X = bose(j,2) + bose(k,3) + 1
        freqtotm1_23= freqm1(j,2) * freqm1(k,3)
        !
        DO i = start_i, end_i !1,S%nat3
          !
          freqtotm1 = freqm1(i,1) * freqtotm1_23
          !IF(freqtot/=0._dp)THEN
          !
          !dom_C =(freq(i,1)+freq(j,2)-freq(k,3))
          dom_C =(freq_initial+freq(j,2)-freq(k,3))
          ctm_C = bose_C * f_gauss(dom_C, sigma)
          !
          !dom_X =(freq(i,1)-freq(j,2)-freq(k,3))
          dom_X =(freq_initial-freq(j,2)-freq(k,3))
          ctm_X = bose_X * f_gauss(dom_X, sigma)
          !
          pref_C =  -pi*ctm_C * V3sq(i,j,k) * freqtotm1
          pref_X =  -pi*ctm_X * V3sq(i,j,k) * freqtotm1
          !DO ie = 1, ne
            !lw(i) = lw(i) - pi*freqtotm1 * (ctm_C + ctm_X) * V3sq(i,j,k)
            fsdf(:,i_final,C) = fsdf(:,i_final,C) +pref_C* wfinal(:)
            fsdf(:,i_final,X) = fsdf(:,i_final,X) +pref_X* wfinal(:)
           !fsdf(:,i_final,TOT) = fsdf(:,i_final,C)+fsdf(:,i_final,X)
          !ENDDO
          !
        ENDDO
        !
        !
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    
    fsdf(:,:,TOT) = fsdf(:,:,C)+fsdf(:,:,X)
    !
!    CALL merge_degen(ne, S%nat3, 3, fsdf, freq(:,1))
    !
    sum_final_state_e = fsdf
    DEALLOCATE(fsdf)
    !
  END FUNCTION sum_final_state_e
  
  
  
END MODULE final_state 
