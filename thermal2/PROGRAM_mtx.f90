!
! Written by Lorenzo Paulatto (2015-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE matrix_program
  !
#include "mpi_thermal.h"
  USE kinds,        ONLY : DP
  USE mpi_thermal,  ONLY : ionode
  USE posix_signal, ONLY : check_graceful_termination
  USE q_grids,      ONLY : q_grid
  USE timers
  !
  ! This type contains a list of active scattering channels 
  ! (i.e. that respect momentum and energy conservation)
  ! associated with a list of q-points.
  ! An active scattering channels is defined by two q-points
  ! in the list (the third is q3 = -q1-q2) and a triplet of bands
  ! nu1, nu2, nu3 that conserve energy at q1,q2,q3. 
  ! This information is stocked in a decently efficient way.
  !
  ! TYPE :: bands_triplet
  !   INTEGER,ALLOCATABLE :: idx1(:)
  !   INTEGER,ALLOCATABLE :: idx2(:)
  !   INTEGER,ALLOCATABLE :: idx3(:)
  ! END TYPE
  !
  ! TYPE, EXTENDS(q_grid) :: scattering_channels
  !   ! The number of active triplets at q-points q_i,q_j,-(q_i+q_j)
  !   INTEGER,ALLOCATABLE :: n_triplets(:,:)
  !   ! The list of active triplets at q-points i,j,(-i-j)
  !   TYPE(bands_triplet),ALLOCATABLE :: bands(:,:)!%idx(nu1,nu2,nu3)
  ! END TYPE

  TYPE, EXTENDS(q_grid) :: q_list
    INTEGER,ALLOCATABLE :: idx(:,:)
  END TYPE

  TYPE :: event_list
    INTEGER :: nev
    INTEGER,ALLOCATABLE :: idx1(:)
    INTEGER,ALLOCATABLE :: idx2(:)
    INTEGER,ALLOCATABLE :: bnd1(:)
    INTEGER,ALLOCATABLE :: bnd2(:)
    INTEGER,ALLOCATABLE :: bnd3(:)
    REAL(DP),ALLOCATABLE :: v3sq(:) 
  END TYPE

  CONTAINS
  !
  ! SUBROUTINE MTX_ELEM(input, grid, S, fc2, fc3)
  !   USE linewidth,          ONLY : linewidth_q
  !   USE constants,          ONLY : RY_TO_CMM1, K_BOLTZMANN_RY, tpi
  !   USE more_constants,     ONLY : RY_TO_WATTMM1KM1, write_conf, ryvel_si
  !   USE q_grids,            ONLY : q_grid, setup_grid
  !   USE fc3_interpolate,    ONLY : forceconst3
  !   USE isotopes_linewidth, ONLY : isotopic_linewidth_q
  !   USE casimir_linewidth,  ONLY : casimir_linewidth_vel, mfp_scatter_vel
  !   USE input_fc,           ONLY : ph_system_info
  !   USE code_input,         ONLY : code_input_type
  !   USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, bose_phq
  !   USE ph_velocity,        ONLY : velocity, velocity_operator
  !   !USE overlap,            ONLY : order_type
  !   USE timers
  !   IMPLICIT NONE
  !   !
  !   TYPE(code_input_type),INTENT(in)  :: input
  !   TYPE(forceconst2_grid),INTENT(in) :: fc2
  !   CLASS(forceconst3),INTENT(in)     :: fc3
  !   TYPE(ph_system_info),INTENT(in)   :: S
  !   TYPE(q_grid),INTENT(in)      :: grid
  !   !
  !   REAL(DP) :: sigma_ry(input%nconf), xq(3,3), freq(S%nat3,3), freqm1(S%nat3,3)
  !   COMPLEX(DP) :: U(S%nat3,S%nat3,3)
  !   REAL(DP) :: max_ctm, avg_ctm, avg_ctm_q
  !   !
  !   INTEGER  :: iq1, iq2, iqq, nu1, nu2, nu3, i, j, k, PASS
  !   !
  !   INTEGER :: n_active_bands(grid%nq, grid%nq)
  !   ! The list of active triplets at q-points i,j,(-i-j)
  !   TYPE(bands_triplet) :: active_bands(grid%nq, grid%nq)!%idx(nu1,nu2,nu3)
  !   !
  !   max_ctm=0._dp
  !   avg_ctm=0._dp
  !   avg_ctm_q=0._dp
  !   n_active_bands = 0

  !   DO PASS = 1,3
  !   QPOINT_LOOP : &
  !   DO iq2 = 1,grid%nq
  !     DO iq1 = 1,grid%nq
  !       !
  !       IF(PASS==2)THEN
  !         threshold = 0.0001_dp * max_ctm
  !       ENDIF

  !       ! on third pass, we know how many bands we have at each triplet
  !       IF(PASS==3) THEN
  !         ALLOCATE(active_bands(iq1,iq2)%idx1(n_active_bands(iq1,iq2)))
  !         ALLOCATE(active_bands(iq1,iq2)%idx2(n_active_bands(iq1,iq2)))
  !         ALLOCATE(active_bands(iq1,iq2)%idx3(n_active_bands(iq1,iq2)
  !       ENDIF
  !       !
  !       freqm1 = 0._dp
  !       DO i = 1,nat3
  !         IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
  !         IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
  !         IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
  !       ENDDO
  !       !
  !       ! check if event is Normal or Umklapp
  !       !refold_bz(xq(:,1), S%bg) +  refold_bz(xq(:,2), S%bg) + refold_bz(xq(:,3), S%bg)
  !       !
  !       DO k = 1,nat3
  !         DO j = 1,nat3
  !           DO i = 1,nat3
  !             !
  !             bose_a = bose(i,1) * bose(j,2) * (bose(k,3)+1)
  !             !
  !             dom_a =  freq(i,1) + freq(j,2) - freq(k,3)
  !             !
  !             ctm_a = bose_a *  f_gauss(dom_a, sigma)
  !             !
              
  !             IF(PASS==1)THEN
  !               max_ctm = MAX(max_ctm, ctm_a)
  !               avg_ctm_q = avg_ctm_q + ctm_a
  !             ELSE IF(max_ctm>threshold)
  !               IF (PASS==2) THEN
  !                 n_active_bands(iq1,iq2) = n_active_bands(iq1,iq2)+1
  !               ELSE
  !                 norm_a = tpi*freqm1(i,1)*freqm1(j,2)*freqm1(k,3)
  !                 sum_a = norm_a * ctm_a * V3sq(i,j,k)
  !                 !WRITE(*,'(3i6,1e24.15)') i,j,k, V3sq(i,j,k)
  !               ENDIF
  !             ENDIF
  !             !
  !           ENDDO
  !         ENDDO
  !       ENDDO
  !       !
  !       avg_ctm = avg_ctm + avg_ctm_q
  !       !
  !     ENDDO ! inner Q_POINT
  !   ENDDO QPOINT_LOOP
  ! ENDDO ! PASS

  ! !
  ! !
  ! END SUBROUTINE MTX_ELEM
  !

  SUBROUTINE MTX_ELEM(grid, events, S, fc2, fc3)
    USE q_grids,            ONLY : q_grid
    USE fc3_interpolate,    ONLY : forceconst3, ip_cart2pat
    USE input_fc,           ONLY : ph_system_info
    USE fc2_interpolate,    ONLY : forceconst2_grid, freq_phq_safe, set_nu0
    USE timers
    USE mpi_thermal
    IMPLICIT NONE
    !
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    CLASS(forceconst3),INTENT(in)     :: fc3
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_list)              :: grid
    TYPE(event_list)          :: events
    !
    REAL(DP) :: xq(3,3), freq(S%nat3,3), faux !freqm1(S%nat3,3),
    COMPLEX(DP) :: U(S%nat3,S%nat3,3), aux, D3(S%nat3,S%nat3,S%nat3)
    !
    INTEGER :: iq1, iq2, iev, jq, nu0(3), i
    INTEGER :: start_iev, end_iev
    !
    ! The list of active triplets at q-points i,j,(-i-j)
    !TYPE(bands_triplet) :: active_bands(grid%nq, grid%nq)!%idx(nu1,nu2,nu3)
    !
    ! for parallelism, we divide in blocks, in order to recycle repeated 
    ! identical triplets
    IF(num_procs>1) events%v3sq = 0._dp
    CALL mpi_block_divide(events%nev, start_iev, end_iev)
    iq1 = -1
    iq2 = -1
    DO iev = start_iev, end_iev
      !
      CALL print_percent_wall(10._dp, 300._dp, iev, end_iev-start_iev+1, (iev==start_iev))
      !
      ! Only recompute the D3 matrix if this point is different from the previous one
      IF ( iq1 /= events%idx1(iev) .or. &
           iq2 /= events%idx2(iev) ) THEN
        !
        iq1 = events%idx1(iev)
        iq2 = events%idx2(iev)
        !
        xq(:,1) = grid%xq(:,iq2)
        xq(:,2) = grid%xq(:,iq2)
        xq(:,3) = -(xq(:,2)+xq(:,1))
        !
        DO jq = 1,3
          nu0(jq) = set_nu0(xq(:,jq), S%at)
          CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
        ENDDO
        !
        ! freqm1 = 0._dp
        ! DO i = 1,S%nat3
        !   IF(i>=nu0(1)) freqm1(i,1) = 0.5_dp/freq(i,1)
        !   IF(i>=nu0(2)) freqm1(i,2) = 0.5_dp/freq(i,2)
        !   IF(i>=nu0(3)) freqm1(i,3) = 0.5_dp/freq(i,3)
        ! ENDDO
        !
        CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
        !V3sq = REAL( CONJG(D3)*D3 , kind=DP)
      ENDIF
      !
      IF( events%bnd1(iev)>=nu0(1) .and. &
          events%bnd2(iev)>=nu0(2) .and. &
          events%bnd3(iev)>=nu0(3) ) THEN

        aux = D3( events%bnd1(iev), &
                  events%bnd2(iev), &
                  events%bnd3(iev)  )

        faux = ( 8 * freq(events%bnd1(iev),1) &
                   * freq(events%bnd2(iev),2) &
                   * freq(events%bnd3(iev),3) )
        !
        events%v3sq(iev) = REAL( CONJG(aux)*aux, kind=DP) / faux
        !
      ENDIF
  ENDDO
  !
  CALL mpi_bsum(events%nev, events%v3sq)
  IF(ionode)THEN
    DO i = 1, events%nev
      WRITE(999,*) events%v3sq(iev) 
    ENDDO
  ENDIF
  !
  END SUBROUTINE MTX_ELEM  
  ! {{{{{{{{{{{{}}}}}}}}}}}}}}
  !
  SUBROUTINE read_qpoints(filename, qpts, S)
    USE kinds, ONLY : DP
    USE input_fc,         ONLY : ph_system_info
    USE mpi_thermal
    TYPE(q_list),INTENT(out)  :: qpts
    TYPE(ph_system_info)       :: S
    INTEGER :: ios, u, i, idx(3)
    REAL(DP) :: xq(3)

    CHARACTER(len=*) :: filename
    OPEN(newunit=u, file=filename, status='old')
    !CALL setup_path(xq, naux, qpts, S%at)

    IF (ionode) THEN
      i=0
      DO WHILE(ios==0)
        READ(u,*,iostat=ios) idx, xq
        IF(ios/=0) EXIT
        i = i+1
      ENDDO
      WRITE(*,*) "Read", i, "q-points"
    ENDIF
    qpts%nq = i

    CALL mpi_broadcast(qpts%nq)
    allocate(qpts%xq(3,qpts%nq))
    allocate(qpts%idx(3,qpts%nq))
    allocate(qpts%w(qpts%nq))
    qpts%nqtot = qpts%nq

    qpts%n  = 0
    qpts%scattered = .false.
    qpts%shifted =  .false.
    qpts%iq0 = 0
    qpts%xq0 = 0._dp    

    IF (ionode) THEN
      REWIND(u)
      !
      DO i = 1, qpts%nq
        READ(u,*,iostat=ios) qpts%idx(:,i), qpts%xq(:,i)
        qpts%w(i) = 1._dp
      ENDDO

      CLOSE(u)
    ENDIF
    !
    CALL mpi_broadcast(3, qpts%nq, qpts%idx)
    CALL mpi_broadcast(3, qpts%nq, qpts%xq)
    CALL mpi_broadcast(qpts%nq, qpts%w)
    !
  END SUBROUTINE
  !
  ! {{{{{{{{{{{{}}}}}}}}}}}}}}
  !
  FUNCTION qlist_get_idx(i, idx, n, n0)
    INTEGER, INTENT(in) :: n, i(3), idx(3,n)
    INTEGER, INTENT(inout) :: n0
    INTEGER :: j, qlist_get_idx
    DO j = n0, n
      IF(ALL(i == idx(:,j))) THEN
        qlist_get_idx = j
        n0 = j
        RETURN
      ENDIF
    ENDDO
    DO j = 1, n0-1
      IF(ALL(i == idx(:,j))) THEN
        qlist_get_idx = j
        n0 = j
        RETURN
      ENDIF
    ENDDO

    IF(j>n) CALL errore("idx","did not find index",1)
  END FUNCTION
  !
  ! {{{{{{{{{{{{}}}}}}}}}}}}}}
  !
  SUBROUTINE read_events_list(filename, qpts, S, events)
    USE kinds,      ONLY : DP
    USE input_fc,   ONLY : ph_system_info
    USE mpi_thermal
    TYPE(q_list),INTENT(in)         :: qpts
    TYPE(ph_system_info),INTENT(in) :: S
    TYPE(event_list)          :: events
    INTEGER :: ios, u, i, i1, i2, idx1(3), idx2(3), bnd(3), nev, n01, n02
    REAL(DP) :: xq(3)

    CHARACTER(len=*) :: filename
    !CALL setup_path(xq, naux, qpts, S%at)
    IF(ionode) THEN
      OPEN(newunit=u, file=filename, status='old')
      i=0
      DO WHILE(ios==0)
        READ(u,*,iostat=ios) idx1, idx2, bnd
        IF(ios/=0) EXIT
        i = i+1
      ENDDO
      WRITE(*,*) "Found", i, "events from file ", filename
      nev = i
    ENDIF

    CALL mpi_broadcast(nev)

    events%nev = nev
    ALLOCATE(events%idx1(nev))
    ALLOCATE(events%idx2(nev))
    ALLOCATE(events%bnd1(nev))
    ALLOCATE(events%bnd2(nev))
    ALLOCATE(events%bnd3(nev))
    ALLOCATE(events%v3sq(nev))

    IF(ionode) THEN
      REWIND(u)
      !
      n01 = 1
      n02 = 1
      DO i = 1, nev
        READ(u,*,iostat=ios) idx1, idx2, bnd
        events%idx1(i) = qlist_get_idx(idx1, qpts%idx, qpts%nq, n01)
        events%idx2(i) = qlist_get_idx(idx2, qpts%idx, qpts%nq, n02)
        events%bnd1(i) = bnd(1)
        events%bnd2(i) = bnd(2)
        events%bnd3(i) = bnd(3)
      ENDDO
    ENDIF

    !write(*,*) "broadcast", nev, my_id
    CALL mpi_broadcast(nev, events%idx1)
    CALL mpi_broadcast(nev, events%idx2)
    CALL mpi_broadcast(nev, events%bnd1)
    CALL mpi_broadcast(nev, events%bnd2)
    CALL mpi_broadcast(nev, events%bnd3)
    !write(*,*) "broadcast done", my_id

    CLOSE(u)
    !
  END SUBROUTINE

  END MODULE matrix_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM matrix
  USE kinds,            ONLY : DP
  USE matrix_program
  USE input_fc,           ONLY : same_system, read_fc2, aux_system, &
                                 forceconst2_grid, ph_system_info, div_mass_fc2
  USE fc3_interpolate,  ONLY : read_fc3, forceconst3
  !USE code_input,       ONLY : READ_INPUT, code_input_type, 
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE nanoclock,        ONLY : init_nanoclock
  USE more_constants,   ONLY : print_citations_linewidth
  USE asr2_module,      ONLY : impose_asr2
  !
  USE posix_signal,       ONLY : set_TERMINATE_GRACEFULLY
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid)     :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)       :: S, S3
  !TYPE(code_input_type)      :: mtxinput
  TYPE(q_list)              :: qpts
  !TYPE(q_grid)              :: qpts2
  TYPE(event_list)          :: events
  ! TYPE(event_list)          :: events_1, events_2, events_1u, events_2u

!   CALL mp_world_start(world_comm)
!   CALL environment_start('TK')
  CALL init_nanoclock()
  CALL start_mpi()
  CALL remove_stack_limit()
  CALL print_citations_linewidth()
  CALL set_TERMINATE_GRACEFULLY() !print_timers_and_die)

  ! READ_INPUT also reads force constants from disk, using subroutine READ_DATA
  !CALL READ_INPUT("TK", mtxinput, grid, S, fc2, fc3)
  CALL read_fc2("mat2R", S,  fc2)
  !
  CALL read_qpoints("qpoints_GaAs_30x15x15", qpts, S)
  !
  !CALL qpts%copy(qpts2)
  !CALL qpts2%scatter()

  CALL read_events_list("interface_30x15x15_GaAs_class1u",  qpts, S, events)
  ! CALL read_events_list("interface_30x15x15_GaAs_class1",  qpts, S, events_1)
  ! CALL read_events_list("interface_30x15x15_GaAs_class2",  qpts, S, events_2)
  ! CALL read_events_list("interface_30x15x15_GaAs_class1u", qpts, S, events_1u)
  ! CALL read_events_list("interface_30x15x15_GaAs_class2u", qpts, S, events_2u)
  !
  !
  CALL aux_system(S)
  CALL impose_asr2('simple', S%nat, fc2, S%zeu)
  CALL div_mass_fc2(S, fc2)
  fc3 => read_fc3("mat3R.asr.sparse", S3)
  CALL fc3%div_mass(S)  !
  CALL MTX_ELEM(qpts, events, S, fc2, fc3)
  !
  !CALL MTX_ELEM(mtxinput, grid, S, fc2, fc3)
  !
  IF(ionode) CALL print_citations_linewidth()
  CALL stop_mpi()

END PROGRAM matrix
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

