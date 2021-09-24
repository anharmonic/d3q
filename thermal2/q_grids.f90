!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE q_grids

  USE kinds,     ONLY : DP
  USE mpi_thermal,      ONLY : ionode
#include "mpi_thermal.h"
  
  TYPE q_grid
    CHARACTER(len=9) :: basis = ''
    CHARACTER(len=9) :: type = ''
    INTEGER :: n(3) = -1
    INTEGER :: nq = 0
    INTEGER :: nqtot = 0
    LOGICAL :: scattered = .false.
    INTEGER :: iq0 = 0
    LOGICAL :: shifted = .false.
    REAL(DP) :: xq0(3) = 0._dp ! the shift applied to the grid (if shifted)
    REAL(DP),ALLOCATABLE :: xq(:,:) ! coordinates of the q-point
    REAL(DP),ALLOCATABLE :: w(:)    ! weight for integral of the BZ
    CONTAINS
      procedure :: scatter => q_grid_scatter
      procedure :: destroy => q_grid_destroy
  END TYPE q_grid
  !
  TYPE  q_basis
    INTEGER :: nbnd  ! number of phonon bands
    INTEGER :: nconf ! number of temperatures
    INTEGER :: nq    ! number of points 
    REAL(DP),ALLOCATABLE :: w(:,:)    ! phonon frequencies
    REAL(DP),ALLOCATABLE :: c(:,:,:)  ! phonon group velocity, \nabla_q w
    REAL(DP),ALLOCATABLE :: be(:,:,:) ! bose-einstein occupation (one per configuration)
    REAL(DP),ALLOCATABLE :: b(:,:,:,:)  ! b vector, as in PRB 88, 045430 (2013)
!     TYPE(q_grid),POINTER :: grid
!     CONTAINS
!       procedure :: B  => B_right_hand_side
  END TYPE q_basis
  !
  CONTAINS
!   ! \/o\________\\\_________________________________________/^>

  SUBROUTINE q_grid_destroy(grid)
    IMPLICIT NONE
    CLASS(q_grid),INTENT(inout) :: grid
    IF(allocated(grid%xq)) DEALLOCATE(grid%xq)
    IF(allocated(grid%w)) DEALLOCATE(grid%w)
    grid%n  = 0
    grid%scattered = .false.
    grid%shifted =  .false.
    grid%nq = 0
    grid%nqtot = 0
    grid%iq0 = 0
    grid%xq0 = 0._dp
  END SUBROUTINE

  SUBROUTINE q_grid_scatter(grid, quiet)
    USE mpi_thermal, ONLY : scatteri_vec, scatteri_mat
    USE functions,   ONLY : default_if_not_present
    IMPLICIT NONE
    CLASS(q_grid),INTENT(inout) :: grid
    LOGICAL,OPTIONAL,INTENT(in) :: quiet
    INTEGER :: nq
    IF(grid%scattered) RETURN
    nq = grid%nq
    CALL scatteri_mat(3,nq, grid%xq)
    nq = grid%nq
    CALL scatteri_vec(nq, grid%w, grid%iq0)
    grid%nq = nq
    grid%scattered = .true.
    IF (.not.default_if_not_present(.false.,quiet) .and. ionode) & 
          WRITE(stdout,"(2x,a)") "Grid scattered with MPI"
  END SUBROUTINE
!   !
!   ! Nasty subroutine that sets some global variables of QE.
!
!       NOTE: all the at, bg, nat, tau and so on must be set in 
!             the QE global modules for this to work
!
!   SUBROUTINE setup_symmetry(S)
!     USE input_fc,  ONLY : ph_system_info 
!     USE symm_base, ONLY : nrot, nsym, set_sym_bl, find_sym
!     IMPLICIT NONE
!     TYPE(ph_system_info),INTENT(in)   :: S ! = System
!     REAL(DP) :: m_loc(3,S%nat) ! fake starting magnetisation
!     m_loc = 0._dp
!     ! ######################### symmetry setup #########################
!     ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
!     CALL set_sym_bl ( )
!     WRITE(*, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
!     !
!     ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
!     CALL find_sym ( S%nat, S%tau, S%ityp, 6,6,6, .false., m_loc )
!     WRITE(*, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
!    END SUBROUTINE setup_symmetry

  SUBROUTINE setup_plane_grid(grid_type, bg, n1,n2, e1, e2, xq0, grid)
    USE input_fc,         ONLY : ph_system_info
    USE random_numbers,   ONLY : randy
    USE functions,        ONLY : default_if_not_present
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)     :: grid_type
    REAL(DP),INTENT(in)   :: bg(3,3)
    INTEGER,INTENT(in) :: n1,n2
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: e1(3), e2(3), xq0(3)

    INTEGER :: i,j,k, idx
    
    grid%type = "plane"
    grid%n(1) = 0
    grid%n(2) = 0
    grid%n(3) = 0
    grid%nq = n1*n2
    !
    IF(allocated(grid%xq)) &
        CALL errore("setup_plane_grid", "grid is already allocated", 1)
    ALLOCATE(grid%xq(3,grid%nq))
    ALLOCATE(grid%w(grid%nq))
    grid%w = 1._dp
    !
    idx = 0
    DO i = 0, n1-1
    DO j = 0, n2-1
      !
      idx = idx+1
      grid%xq(:,idx) = REAL(i,kind=DP)/REAL(n1,kind=DP)*e1 &
                      +REAL(j,kind=DP)/REAL(n2,kind=DP)*e2 &
                      +xq0
      !
    ENDDO
    ENDDO
    
    IF(idx/=grid%nq) CALL errore("setup_plane_grid","wrong nq",1)
    
  END SUBROUTINE

! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_grid(grid_type, bg, n1,n2,n3, grid, xq0, scatter, quiet)
    USE input_fc,         ONLY : ph_system_info
    USE random_numbers,   ONLY : randy
    USE functions,        ONLY : default_if_not_present
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)     :: grid_type
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    LOGICAL,OPTIONAL,INTENT(in) :: scatter
    LOGICAL,OPTIONAL,INTENT(in) :: quiet
    INTEGER :: ipol
    !
    LOGICAL  :: do_scatter
    do_scatter = .false.
    IF(present(scatter)) do_scatter = scatter
    IF(present(xq0)) THEN
      IF(SUM(ABS(xq0))>0._dp) THEN
        grid%xq0 = xq0
        grid%shifted = .true.
        IF(grid_type=="random" .or. grid_type=="randws") THEN
          ioWRITE(stdout, "(2x,a,3f12.6)") &
            "WARNING! Random grid: ignoring grid shift from input", grid%xq0
        ELSE
          IF (.not.default_if_not_present(.false.,quiet) .and. ionode) & 
          WRITE(stdout, "(2x,a,3f12.6)") "Applying grid shift", grid%xq0
        ENDIF
      ENDIF
    ENDIF
    !
    IF(grid_type=="simple" .or. grid_type=="grid")THEN
      IF(.not.do_scatter)THEN
        ! If the grid is not mpi-scattered I use the simple subroutine
        CALL setup_simple_grid(bg, n1,n2,n3, grid, xq0)
      ELSE
        ! otherwise, I use this subroutine instead, which directly scatters over MPI:
        CALL setup_scattered_grid(bg, n1,n2,n3, grid, xq0)
      ENDIF
    ELSE IF(grid_type=="bxsf" .or. grid_type=="xsf")THEN
      CALL setup_xcrysden_grid(grid_type, bg, n1,n2,n3, grid, xq0)
      IF(do_scatter) CALL grid%scatter()
    ELSE IF(grid_type=="random" .or. grid_type=="randws")THEN
      ! Do not add a random shift in the direction where 
      ! only on point is requested.
      ! The random shift is inside the first grid cell
      grid%xq0 = 0._dp
      IF(n1>1) grid%xq0(1) = randy()/DBLE(n1)
      IF(n2>1) grid%xq0(2) = randy()/DBLE(n2)
      IF(n3>1) grid%xq0(3) = randy()/DBLE(n3)
      CALL cryst_to_cart(1,grid%xq0,bg, +1)
      grid%shifted = .true.
      IF (.not.default_if_not_present(.false.,quiet) .and. ionode) & 
        WRITE(stdout, "(2x,a,3f12.6)") "Random grid shift", grid%xq0
      !
      IF(.not.do_scatter)THEN
        ! If the grid is not mpi-scattered I use the simple subroutine
        CALL setup_simple_grid(bg, n1,n2,n3, grid, grid%xq0)
      ELSE
        ! otherwise, I use this subroutine instead, which directly scatters over MPI:
        CALL setup_scattered_grid(bg, n1,n2,n3, grid, grid%xq0)
      ENDIF
    ELSE IF (grid_type=="ws" .or. grid_type=="bz")THEN
      ! This grid has to be generated in its entirety and then scatterd, 
      ! mind the memory bottleneck!
      !print*,1
      CALL setup_bz_grid(bg, n1,n2,n3, grid, grid%xq0)
      !print*,2
      IF(do_scatter) CALL grid%scatter()
      !print*,3
    ELSE
      CALL errore("setup_grid", "wrong grid type '"//TRIM(grid_type)//"'", 1)
    ENDIF
    IF (.not.default_if_not_present(.false.,quiet) .and. ionode) & 
       WRITE(stdout,'(2x,"Setup a ",a," grid of",i9," q-points")') grid_type, grid%nqtot

  END SUBROUTINE setup_grid

  ! \/o\________\\\_________________________________________/^>
  ! Setup a grid that respects the symmetry of the Brillouin zone
  SUBROUTINE setup_bz_grid(bg, n1,n2,n3, grid, xq0)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    TYPE(q_grid) :: sg
    INTEGER :: i,j,k, iq, nqtot
    REAL(DP) :: g(3), qg(3), qbz(3,12), qgmod2, qgmod2_min
    INTEGER :: nqbz
    INTEGER,PARAMETER :: far=2
    REAL(DP),PARAMETER :: eps = 1.d-4

    grid%type = 'bz'
    ! Do not put the optional xq0 shift here, I'll add it at the end to
    ! have a grid that is centered around xq0
    CALL setup_simple_grid(bg, n1,n2,n3, sg, xq0)

    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    !
    IF(allocated(grid%xq)) CALL errore("setup_simple_grid", "grid is already allocated", 1)
    !
    nqtot = 0
    grid%nq = 0

    ! first count the points
    DO iq = 1, sg%nq
      qgmod2_min = SUM(sg%xq(:,iq)**2)
      qbz = 0._dp
      nqbz = 0
      DO i = -far,far
      DO j = -far,far
      DO k = -far,far
        g = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
        qg = sg%xq(:,iq)+g
        qgmod2 = SUM(qg**2)
        IF(ABS(qgmod2-qgmod2_min)< eps)THEN
          ! same distance: add this vector to the list f degeneracies
          nqbz = nqbz+1
          qbz(:,nqbz) = qg
          qgmod2_min = qgmod2
        ELSEIF(qgmod2<qgmod2_min)THEN
          ! shorter distance: reinit the list of degeneracies
          qbz = 0._dp
          nqbz = 1 
          qbz(:,nqbz) = qg
          qgmod2_min = qgmod2
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      !
      nqtot = nqtot + nqbz
    ENDDO    

    ALLOCATE(grid%xq(3,nqtot))
    ALLOCATE(grid%w(nqtot))
    
    DO iq = 1, sg%nq
      !print*, 2000, iq
      qgmod2_min = SUM(sg%xq(:,iq)**2)
      qbz = 0._dp
      nqbz = 0
      DO i = -far,far
      DO j = -far,far
      DO k = -far,far
        g = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
        qg = sg%xq(:,iq)+g
        qgmod2 = SUM(qg**2)
        IF(ABS(qgmod2-qgmod2_min)< eps)THEN
          ! same distance: add this vector to the list f degeneracies
          nqbz = nqbz+1
          qbz(:,nqbz) = qg
          qgmod2_min = qgmod2
        ELSEIF(qgmod2<qgmod2_min)THEN
          ! shorter distance: reinit the list of degeneracies
          qbz = 0._dp
          nqbz = 1 
          qbz(:,nqbz) = qg
          qgmod2_min = qgmod2
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      !
      grid%xq(:,grid%nq+1:grid%nq+nqbz) = qbz(:,1:nqbz)
      grid%w(grid%nq+1:grid%nq+nqbz) = 1._dp/nqbz
      grid%nq = grid%nq + nqbz
      ! here we have a complete list of equivalent q-points
      ! we add them to the grid with their weight
      !print*, "found", nqbz, qbz(:,1:nqbz)
!       CALL expand_wlist(grid%nq, grid%w,  nqbz)
!       CALL expand_qlist(grid%nq, grid%xq, nqbz, qbz)
    ENDDO
    IF(grid%nq /= nqtot) CALL errore('setup_bz_grid', 'nq problem', 1)

    grid%basis = 'cartesian'
    IF(NINT(sum(grid%w)) /= sg%nq) CALL errore('setup_bz_grid','wrong weight',1)
    ! renormalize the weights
    grid%w = grid%w / DBLE(sg%nq)
    !
!    IF(present(xq0)) THEN
!      DO iq = 1,grid%nq
!        grid%xq(:,iq) = grid%xq(:,iq) + xq0
!      ENDDO
!    ENDIF
    !
    grid%nqtot = grid%nq
    !
  END SUBROUTINE setup_bz_grid
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_simple_grid(bg, n1,n2,n3, grid, xq0)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    INTEGER :: i,j,k, idx
    
    grid%type = 'simple'
    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    grid%nq = n1*n2*n3
    !
    IF(allocated(grid%xq)) CALL errore("setup_simple_grid", "grid is already allocated", 1)
    ALLOCATE(grid%xq(3,grid%nq))
    ALLOCATE(grid%w(grid%nq))
    grid%w = 1._dp/grid%nq
    !
    idx = 0
    DO i = 0, n1-1
    DO j = 0, n2-1
    DO k = 0, n3-1
      !
      idx = idx+1
      grid%xq(1,idx) = REAL(i,kind=DP)/REAL(n1,kind=DP)
      grid%xq(2,idx) = REAL(j,kind=DP)/REAL(n2,kind=DP)
      grid%xq(3,idx) = REAL(k,kind=DP)/REAL(n3,kind=DP)
      !
    ENDDO
    ENDDO
    ENDDO
    !
    CALL cryst_to_cart(grid%nq,grid%xq,bg, +1)
    grid%basis = 'cartesian'
    !
    IF(present(xq0)) THEN
      DO idx = 1,grid%nq
        grid%xq(:,idx) = grid%xq(:,idx) + xq0
      ENDDO
      grid%xq0 = xq0
    ELSE
      grid%xq0 = 0._dp
    ENDIF
    !
    grid%nqtot = grid%nq
    !
  END SUBROUTINE setup_simple_grid
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_xcrysden_grid(grid_type, bg, n1,n2,n3, grid, xq0)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: grid_type
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    INTEGER :: i,j,k, idx
    grid%n(1) = n1+1
    grid%n(2) = n2+1
    grid%n(3) = n3+1
    grid%nq = (n1+1)*(n2+1)*(n3+1)
    !
    IF(allocated(grid%xq)) CALL errore("setup_xcrysden_grid", "grid is already allocated", 1)
    ALLOCATE(grid%xq(3,grid%nq))
    ALLOCATE(grid%w(grid%nq))
    grid%w = 1._dp/(n1*n2*n3) ! Note: not 1/grid%nq
    grid%type = grid_type
    !
    idx = 0
    DO k = 0, n3
    DO j = 0, n2
    DO i = 0, n1
      !
      idx = idx+1
      grid%xq(1,idx) = REAL(i,kind=DP)/REAL(n1,kind=DP)
      grid%xq(2,idx) = REAL(j,kind=DP)/REAL(n2,kind=DP)
      grid%xq(3,idx) = REAL(k,kind=DP)/REAL(n3,kind=DP)
      !
    ENDDO
    ENDDO
    ENDDO
    !
    CALL cryst_to_cart(grid%nq,grid%xq,bg, +1)
    grid%basis = 'cartesian'
    !
    IF(present(xq0)) THEN
      DO idx = 1,grid%nq
        grid%xq(:,idx) = grid%xq(:,idx) + xq0
      ENDDO
      grid%xq0 = xq0
    ELSE
      grid%xq0 = 0._dp
    ENDIF
    !
    grid%nqtot = grid%nq
    !
  END SUBROUTINE setup_xcrysden_grid
  ! \/o\________\\\_________________________________________/^>
  ! Create a scattered grid without generating a full local grid first
  ! which becomes a memory bottleneck in some cases
  SUBROUTINE setup_scattered_grid(bg, n1,n2,n3, grid, xq0)
    USE input_fc,    ONLY : ph_system_info
    USE mpi_thermal
    IMPLICIT NONE
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    INTEGER :: nqme, nqtot, nqdiv, nqresidual, nqfirst
    INTEGER :: i,j,k, idx, my_idx
    !
    IF(.not.mpi_started) CALL errore("setup_scattered_grid", "mpi must be started first",1)
    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    grid%type = 'simple'
    ! find the number of points for this cpu, the first cpus can have 
    ! an extra point if num_procs is not a divisor of nqtot
    nqtot = n1*n2*n3
    grid%nqtot = nqtot
    !
    nqdiv = nqtot/num_procs
    !print*, nqtot, num_procs, nqdiv
    nqresidual = nqtot-nqdiv*num_procs
    nqme = nqdiv
    IF(my_id<nqresidual) nqme = nqme+1
    grid%nq=nqme
    
    nqfirst = 1
    DO idx = 0,my_id-1
      nqfirst = nqfirst+nqdiv
      IF(idx<nqresidual) nqfirst = nqfirst+1
    ENDDO
    grid%iq0 = nqfirst-1
    !print*, "me", my_id, nqfirst, nqme, nqresidual, nqdiv, num_procs
    !
    IF(allocated(grid%xq)) CALL errore("setup_simple_grid", "grid is already allocated", 1)
    ALLOCATE(grid%xq(3,grid%nq))
    ALLOCATE(grid%w(grid%nq))
    grid%w = 1._dp/grid%nqtot
    grid%scattered = .true.
    !
    idx = 0
    NQ_LOOP : &
    DO i = 0, n1-1
    DO j = 0, n2-1
    DO k = 0, n3-1
      !
      idx = idx+1
      my_idx = idx -nqfirst+1
      IF(my_idx>nqme) EXIT NQ_LOOP
      IF(my_idx>0)THEN
        !print*, "yayaya", my_id, idx, my_idx
        grid%xq(1,my_idx) = REAL(i,kind=DP)/REAL(n1,kind=DP)
        grid%xq(2,my_idx) = REAL(j,kind=DP)/REAL(n2,kind=DP)
        grid%xq(3,my_idx) = REAL(k,kind=DP)/REAL(n3,kind=DP)
      ENDIF
      !
    ENDDO
    ENDDO
    ENDDO &
    NQ_LOOP
    !
    CALL cryst_to_cart(grid%nq,grid%xq,bg, +1)
    grid%basis = 'cartesian'
    !
    IF(present(xq0)) THEN
      DO idx = 1,grid%nq
        grid%xq(:,idx) = grid%xq(:,idx) + xq0
        grid%xq0 = xq0
      ENDDO
    ENDIF
    !
  END SUBROUTINE setup_scattered_grid
  !
  ! \/o\________\\\_________________________________________/^>
  ! Create a line of nq q-point from the last point in the path end xqi 
  ! weights are improperly used as path length
  ! SPECIAL: 
  ! 1. If the new point xqi is equivalent to the last point in the line, the path
  ! length is not incremented (only works if nq_new == 1)
  ! 2. If nq_new is -1 the new point is added to the line and the path length 
  ! is reset to zero, also lw.x will print an empy line in the output
  SUBROUTINE setup_path(xqi, nq_new, path, at)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nq_new
    REAL(DP),INTENT(in) :: xqi(3)
    REAL(DP),INTENT(in),OPTIONAL :: at(3,3)
    TYPE(q_grid),INTENT(inout) :: path
    !
    INTEGER :: i, nq_old
    REAL(DP),PARAMETER:: eps = 1.d-8
    REAL(DP),ALLOCATABLE :: auxq(:,:), auxw(:)
    REAL(DP) :: dq(3), dql
    LOGICAL ::  equiv
    !
    IF(nq_new==0) THEN
      IF(path%nq==0) &
        CALL errore('setup_path', 'cannot bootstrap the path', 2)
      ioWRITE(*,'(2x,"Path skip point:    ",f12.6)') path%w(path%nq)
      !RETURN
    ENDIF
    !
    IF(nq_new>1 .and. path%nq==0) &
      CALL errore('setup_path', 'cannot bootstrap the path', 1)
    !
    IF(nq_new<-1) CALL errore("setup_path", "nq_new must be >= -1",1)
    !
    ! If I add a point to a grid, I get something strange
    IF( path%nq==0)THEN
      path%type = 'path'
    ELSE
      IF(path%type /= 'path') path%type = 'unknown'
    ENDIF
    !
    ADD_TO_PATH : &
    IF(path%nq >0)THEN
      ALLOCATE(auxq(3,path%nq),auxw(path%nq))
      auxq = path%xq
      auxw = path%w
      DEALLOCATE(path%xq, path%w)
      ! If xqi is the same as the last point in the list, 
      ! do not add this distance to the path length
      nq_old = path%nq
      IF(nq_new>0)THEN
        path%nq = path%nq + nq_new
      ELSE
        path%nq = path%nq + 1
      ENDIF
      !
      ALLOCATE(path%xq(3,path%nq),path%w(path%nq))
      path%xq(:,1:nq_old) = auxq(:,1:nq_old)
      path%w(1:nq_old) = auxw(1:nq_old)
      ! check for turning point
      !
      DEALLOCATE(auxq,auxw)
    ELSE ADD_TO_PATH
      nq_old=0
      path%nq = nq_new
      ALLOCATE(path%xq(3,path%nq),path%w(path%nq))
      path%w(1)=0._dp
    ENDIF &
    ADD_TO_PATH
    !
    !
    !
    equiv = .false.
    COMPUTE_NEW_PATH : &
    IF(nq_new>1)THEN
      dq = (xqi-path%xq(:,nq_old))/(nq_new)
      dql = DSQRT(SUM(dq**2))
      ! build the actual path
      DO i = nq_old+1, path%nq
        path%xq(:,i) = path%xq(:,nq_old) + dq * (i-nq_old)
        IF(i>1) path%w(i) = path%w(i-1) + DSQRT(SUM(dq**2))
      ENDDO
    ELSE &
    COMPUTE_NEW_PATH
      IF(nq_old>0)THEN
        dq = xqi-path%xq(:,nq_old)
        ! compute path length before taking it to crystal axes
        dql = DSQRT(SUM(dq**2))
        ! if cell vectors are available, we can check for periodic image equivalence
        IF(present(at)) THEN
          CALL cryst_to_cart(1,dq,at,-1)
          !ioWRITE(*,'(a,2(3f12.6,5x))') "equiv a", dqG, REAL(NINT(dq))
          dq = dq-NINT(dq)
        ENDIF
        equiv=SUM(dq**2)<1.d-8
        !
        !IF(equiv.and.nq_new==1)THEN
        IF(equiv)THEN
          ioWRITE(*,*) "equiv", xqi
          path%w(nq_old+1) = path%w(nq_old)
        ELSE
          path%w(nq_old+1) = path%w(nq_old)+dql
        ENDIF
        IF(nq_new==-1) path%w(nq_old+1) = 0._dp
        IF(nq_new==0)  path%w(nq_old+1) = path%w(nq_old)
      ENDIF
      !
      !
      path%xq(:,path%nq) = xqi
    ENDIF  &
    COMPUTE_NEW_PATH
    ! 
    IF(nq_old>1)THEN
    IF(SUM( ( (path%xq(:,nq_old-1) - path%xq(:,nq_old)) &
             -(path%xq(:,nq_old) - path%xq(:,nq_old+1)))**2 ) > 1.d-6 ) THEN
      ioWRITE(*,'(2x,"Path turning point: ",f12.6)') path%w(nq_old)
    ENDIF
    ENDIF
    !
  END SUBROUTINE setup_path
  ! Prepare a derived type that contains phonon frequencies, group velocities
  ! Bose-Einstein distribution for the input grid and configurations.
  ! Also allocate the space to store the perturbed phonon population "b"
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE prepare_q_basis(qgrid, qbasis, nconf, T, S, fc2)
    USE fc2_interpolate,    ONLY : freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity
    USE input_fc,           ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    TYPE(q_grid),INTENT(in)   :: qgrid
    TYPE(q_basis),INTENT(out) :: qbasis
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    INTEGER             :: nconf
    REAL(DP),INTENT(in) :: T(nconf)
    !
    INTEGER :: ix, nu, iq, it
    COMPLEX(DP),ALLOCATABLE :: U(:,:)
    !
    IF(qgrid%nq<=0) CALL errore("prepare_q_basis", "grid must be set up first", 1)
    !
    qbasis%nbnd = S%nat3
    qbasis%nconf = nconf
    qbasis%nq    = qgrid%nq
    !qbasis%grid => qgrid
    !
    ALLOCATE(qbasis%w(S%nat3, qbasis%nq))
    ALLOCATE(qbasis%c(3, S%nat3, qbasis%nq))
    ALLOCATE(qbasis%be(S%nat3, qbasis%nconf, qbasis%nq))
    ALLOCATE(qbasis%b(3,nconf,S%nat3,qbasis%nq))
    !
    ALLOCATE(U(S%nat3,S%nat3))
    
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix, U)
!$OMP DO
    DO iq = 1,qbasis%nq
      qbasis%c(:,:,iq) = velocity(S, fc2, qgrid%xq(:,iq))
      CALL  freq_phq_safe(qgrid%xq(:,iq), S, fc2, qbasis%w(:,iq), U)
      DO it = 1, nconf
        CALL  bose_phq(T(it), S%nat3, qbasis%w(:,iq), qbasis%be(:,it,iq))
      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO
    DO iq = 1,qbasis%nq
    DO nu = 1,S%nat3
      DO it = 1,nconf
      DO ix = 1,3
        qbasis%b(ix,it,nu,iq) = - qbasis%c(ix,nu,iq) * qbasis%w(nu,iq) &
                                 *qbasis%be(nu,it,iq)*(qbasis%be(nu,it,iq)+1)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
    DEALLOCATE(U)
    
  END SUBROUTINE prepare_q_basis
  ! \/o\________\\\_________________________________________/^>
  FUNCTION qbasis_dot(x, y, nconf, nat3, nq) RESULT(z)
    IMPLICIT NONE
    !
    REAL(DP) :: z(3, nconf)
    !
    REAL(DP),INTENT(in) :: x(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: y(3, nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, ix, nu
    !
    z = 0._dp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix)
!$OMP DO COLLAPSE(4)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
      DO ix = 1,3
        z(ix,it) = z(ix,it)+x(ix,it,nu,iq)*y(ix,it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
  END FUNCTION qbasis_dot
  ! \/o\________\\\_________________________________________/^>
  FUNCTION qbasis_ax(a, x, nconf, nat3, nq) RESULT(y)
    IMPLICIT NONE
    !
    REAL(DP) :: y(3, nconf, nat3, nq)
    !
    REAL(DP),INTENT(in) :: a(3, nconf)
    REAL(DP),INTENT(in) :: x(3, nconf, nat3, nq)
    INTEGER,INTENT(in)  :: nconf, nat3, nq
    !
    INTEGER  :: iq, it, ix, nu
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix)
!$OMP DO COLLAPSE(4)
    DO iq = 1,nq
    DO nu = 1,nat3
      DO it = 1,nconf
      DO ix = 1,3
        y(ix,it,nu,iq) = a(ix,it)*x(ix,it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
  END FUNCTION qbasis_ax
  ! \/o\________\\\_________________________________________/^>
  FUNCTION qbasis_a_over_b(a, b, nconf) RESULT(c)
    IMPLICIT NONE
    !
    REAL(DP) :: c(3, nconf)
    !
    REAL(DP),INTENT(in) :: a(3, nconf)
    REAL(DP),INTENT(in) :: b(3, nconf)
    INTEGER,INTENT(in)  :: nconf
    !
    INTEGER  :: it, ix
    !
    c = 0._dp
    DO it = 1,nconf
    DO ix = 1,3
      IF(b(ix,it)/=0._dp)THEN
        c(ix,it) = a(ix,it)/b(ix,it)
      ELSE
        IF(a(ix,it)/=0._dp) WRITE(*,'(3x,a,2i4)') "WARNING! Division by zero qbasis a/b", ix, it
      ENDIF
    ENDDO
    ENDDO
    !
  END FUNCTION qbasis_a_over_b
  !
  ! \/o\________\\\_________________________________________/^>
  REAL(DP) FUNCTION B_right_hand_side(self, ix, it, nu, iq) &
           RESULT (B)
    IMPLICIT NONE
    CLASS(q_basis),INTENT(in) :: self
    INTEGER :: ix, nu, iq, it
    !
    B = -self%c(ix,nu,iq) * self%w(nu,iq) * self%be(nu,it,iq) * (self%be(nu,it,iq)+1)
    !
  END FUNCTION B_right_hand_side

END MODULE q_grids





