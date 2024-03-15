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
  ! Extend the list of q-point, inclue also the index (i,j,k) in a grid
  TYPE, EXTENDS(q_grid) :: q_list
    INTEGER,ALLOCATABLE :: idx(:,:)
  END TYPE

  ! 3-phonon scattering event
  TYPE :: event_list
    INTEGER :: nev
    INTEGER,ALLOCATABLE :: idx1(:) ! first phonon from a list
    INTEGER,ALLOCATABLE :: idx2(:) ! second phonon
    INTEGER,ALLOCATABLE :: bnd1(:) ! bands1
    INTEGER,ALLOCATABLE :: bnd2(:)
    INTEGER,ALLOCATABLE :: bnd3(:)
    REAL(DP),ALLOCATABLE :: v3sq(:) ! matrix element squared |V_3|^2
  END TYPE

  CONTAINS

  SUBROUTINE MTX_ELEM(grid, events, S, fc2, fc3, filename)
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
    CHARACTER(len=*),INTENT(in) :: filename
    TYPE(q_list)                :: grid
    TYPE(event_list)            :: events
    !
    REAL(DP) :: xq(3,3), freq(S%nat3,3), faux !freqm1(S%nat3,3),
    COMPLEX(DP) :: U(S%nat3,S%nat3,3), aux, D3(S%nat3,S%nat3,S%nat3)
    !
    INTEGER :: iq1, iq2, iev, jq, nu0(3), i, unit
    INTEGER :: start_iev, end_iev
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
        xq(:,1) = grid%xq(:,iq1)
        xq(:,2) = grid%xq(:,iq2)
        xq(:,3) = -(xq(:,2)+xq(:,1))
        !
        DO jq = 1,3
          nu0(jq) = set_nu0(xq(:,jq), S%at)
          CALL freq_phq_safe(xq(:,jq), S, fc2, freq(:,jq), U(:,:,jq))
        ENDDO
        !
       !
        CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
        CALL ip_cart2pat(D3, S%nat3, U(:,:,1), U(:,:,2), U(:,:,3))
      ENDIF
      !
      IF( events%bnd1(iev)>=nu0(1) .and. &
          events%bnd2(iev)>=nu0(2) .and. &
          events%bnd3(iev)>=nu0(3) ) THEN

        aux = D3( events%bnd1(iev), &
                  events%bnd2(iev), &
                  events%bnd3(iev)  )

        ! the factor 1/sqrt(mass) is put directly in the force constants
        ! here we put sqrt(1/2w), squared
        faux = ( 8 * freq(events%bnd1(iev),1) &
                   * freq(events%bnd2(iev),2) &
                   * freq(events%bnd3(iev),3) )
        !
        events%v3sq(iev) = REAL( CONJG(aux)*aux, kind=DP) / faux
        !
      ENDIF
  ENDDO
  !
  ! Collect in the stupidest way possible
  IF(ionode .and. num_procs>1) WRITE(*,*) "Collecting..."
  CALL mpi_bsum(events%nev, events%v3sq)
  !
  IF(ionode)THEN
    OPEN(newunit=unit, file=filename, status='new')
    WRITE(*,*) "Writing to disk..."
    DO iev = 1, events%nev
      WRITE(unit,*) events%v3sq(iev) 
    ENDDO
    CLOSE(unit)
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
      WRITE(*,*) "Found", i, "q-points"
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
      WRITE(*,*) "Read", i, "q-points"

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
  ! Find a q-point in a list, start looking from n0
  ! Optimization: because q-points are in incremental order,
  ! we start looking from the last found one, then wrap around
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
    ! wrap
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
      WRITE(*,*) "Found", i, "events from file ", TRIM(filename)
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

    CALL mpi_broadcast(nev, events%idx1)
    CALL mpi_broadcast(nev, events%idx2)
    CALL mpi_broadcast(nev, events%bnd1)
    CALL mpi_broadcast(nev, events%bnd2)
    CALL mpi_broadcast(nev, events%bnd3)

    CLOSE(u)
    !
  END SUBROUTINE

  END MODULE matrix_program
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM matrix
  USE kinds,            ONLY : DP
  USE input_fc,           ONLY : same_system, read_fc2, aux_system, &
                                 forceconst2_grid, ph_system_info, div_mass_fc2
  USE fc3_interpolate,  ONLY : read_fc3, forceconst3
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi
  USE nanoclock,        ONLY : init_nanoclock
  USE more_constants,   ONLY : print_citations_linewidth
  USE asr2_module,      ONLY : impose_asr2
  !
  USE matrix_program
  USE cmdline_param_module
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid)     :: fc2
  CLASS(forceconst3),POINTER :: fc3
  TYPE(ph_system_info)       :: S, S3
  TYPE(q_list)               :: qpts
  TYPE(event_list)           :: events
  CHARACTER(len=1024) :: q_filename, i_filename, o_filename, filemat3, filemat2
  CHARACTER(len=8) :: asr_type
  INTEGER :: event_type 

  q_filename = cmdline_param_char("q", "qpoints_GaAs_30x15x15")
  i_filename = cmdline_param_char("i", "interface_30x15x15_GaAs_class1")
  o_filename = cmdline_param_char("o", TRIM(i_filename)//"_mtx")
  IF(o_filename == i_filename .or. o_filename == q_filename) CALL errore("mtx", "bad filename choice, I wont overwrite input!",1)
  event_type = cmdline_param_int("e", 1)
  asr_type = cmdline_param_char("a", "simple")
  filemat2 = cmdline_param_char("2", "mat2R")
  filemat3 = cmdline_param_char("3", "mat3R.asr.sparse")

  CALL init_nanoclock()
  CALL start_mpi()
  CALL print_citations_linewidth()
  !CALL set_TERMINATE_GRACEFULLY() !print_timers_and_die)
  !
  CALL read_fc2(filemat2, S,  fc2)
  !
  CALL read_qpoints(q_filename, qpts, S)
  !
  !CALL qpts%copy(qpts2)
  !CALL qpts2%scatter()

  CALL read_events_list(i_filename,  qpts, S, events)
  !
  CALL aux_system(S)
  CALL impose_asr2(asr_type, S%nat, fc2, S%zeu)
  CALL div_mass_fc2(S, fc2)
  !
  fc3 => read_fc3(filemat3, S3)
  CALL fc3%div_mass(S) 
  !
  CALL MTX_ELEM(qpts, events, S, fc2, fc3, o_filename)
  !
  IF(ionode) CALL print_citations_linewidth()
  CALL stop_mpi()

END PROGRAM matrix
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

