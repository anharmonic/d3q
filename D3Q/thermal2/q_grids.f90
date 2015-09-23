!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE q_grids

  USE kinds,     ONLY : DP
  USE mpi_thermal,      ONLY : ionode
#include "mpi_thermal.h"
  
  TYPE q_grid
    CHARACTER(len=9) :: basis = ''
    INTEGER :: n(3) = -1
    INTEGER :: nq = 0
    LOGICAL :: scattered = .false.
    INTEGER :: iq0 = 0
    REAL(DP),ALLOCATABLE :: xq(:,:) ! coordinates of the q-point
    REAL(DP),ALLOCATABLE :: w(:)    ! weight for integral of the BZ
    CONTAINS
      procedure :: scatter => q_grid_scatter
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
! !   !
! !   ! A scalar in this special space (i.e. an array of size 3,nconf)
! !   TYPE  q_grid_real_scalar
! !     INTEGER :: nconf ! number of temperatures
! !     REAL(DP),ALLOCATABLE :: r(:,:)
! !     CONTAINS
! !       procedure :: assign_q_grid_real_scalar
! !   END TYPE q_grid_real_scalar
! !   INTERFACE q_grid_real_scalar
! !     MODULE PROCEDURE real_scalar_on_q_grid
! !   END INTERFACE
! !   !
! !   ! A vector in this special space (i.e. an array of size 3,nconf,3*nat,nqpts)
! !   TYPE  q_grid_real_vector
! !     INTEGER :: nconf ! number of temperatures
! !     REAL(DP),ALLOCATABLE :: r(:,:,:,:)
! !     CONTAINS
! !       procedure :: assign_q_grid_real_scalar
! !   END TYPE q_grid_real_vector
! !   INTERFACE q_grid_real_vector
! !     MODULE PROCEDURE real_vector_on_q_grid
! !   END INTERFACE
  
  
  CONTAINS
!   ! \/o\________\\\_________________________________________/^>

  SUBROUTINE q_grid_scatter(grid)
    USE mpi_thermal, ONLY : scatteri_vec, scatteri_mat
    IMPLICIT NONE
    CLASS(q_grid),INTENT(inout) :: grid
    INTEGER :: nq
    nq = grid%nq
    CALL scatteri_mat(3,nq, grid%xq)
    nq = grid%nq
    CALL scatteri_vec(nq, grid%w, grid%iq0)
    grid%nq = nq
    grid%scattered = .true.
    ioWRITE(stdout,"(2x,a)") "Grid scattered with MPI"
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
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE setup_grid(grid_type, bg, n1,n2,n3, grid, xq0)
    USE input_fc, ONLY : ph_system_info
    IMPLICIT NONE
    CHARACTER(len=6),INTENT(in)     :: grid_type
    REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    !
    IF(grid_type=="simple")THEN
      CALL setup_simple_grid(bg, n1,n2,n3, grid, xq0)
    ELSE IF (grid_type=="bz")THEN
      CALL setup_bz_grid(bg, n1,n2,n3, grid, xq0)
    ELSE
      CALL errore("setup_grid", "wrong grid type", 1)
    ENDIF
    ioWRITE(stdout,'(2x,"Setup a ",a," grid of",i6," q-points")') grid_type, grid%nq
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
    INTEGER :: i,j,k, iq
    REAL(DP) :: g(3), qg(3), qbz(3,12), qgmod2, qgmod2_min
    INTEGER :: nqbz
    INTEGER,PARAMETER :: far=2
    REAL(DP),PARAMETER :: eps = 1.d-4

    ! Do not put the optional xq0 shift here, I'll add it at the end to
    ! have a grid that is centered around xq0
    CALL setup_simple_grid(bg, n1,n2,n3, sg)

    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    !
    IF(allocated(grid%xq)) CALL errore("setup_simple_grid", "grid is already allocated", 1)
    !
    grid%nq = 0
    
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
      ! here we have a complete list of equivalent q-points
      ! we add them to the grid with their weight
      !print*, "found", nqbz, qbz(:,1:nqbz)
      CALL expand_wlist(grid%nq, grid%w,  nqbz)
      CALL expand_qlist(grid%nq, grid%xq, nqbz, qbz)
      grid%nq = grid%nq + nqbz
    ENDDO

    grid%basis = 'cartesian'
    IF(NINT(sum(grid%w)) /= sg%nq) CALL errore('setup_bz_grid,''wrong weight',1)
    grid%w = grid%w / DBLE(sg%nq)
    !
    IF(present(xq0)) THEN
      DO iq = 1,grid%nq
        grid%xq(:,iq) = grid%xq(:,iq) + xq0
      ENDDO
    ENDIF
    !
  END SUBROUTINE setup_bz_grid
  !
  SUBROUTINE expand_qlist(nq, xq, nq_add, xq_add)
    INTEGER,INTENT(inout) :: nq 
    INTEGER,INTENT(in) :: nq_add
    REAL(DP),ALLOCATABLE,INTENT(inout) :: xq(:,:)
    REAL(DP),INTENT(in) :: xq_add(3,nq_add)
    !
    REAL(DP),ALLOCATABLE :: xq_aux(:,:)
    IF(allocated(xq)) THEN
      ALLOCATE(xq_aux(3,nq))
      xq_aux = xq
      DEALLOCATE(xq)
      ALLOCATE(xq(3,nq+nq_add))
      xq(:, 1:nq)           = xq_aux
      xq(:, nq+1:nq+nq_add) = xq_add(:, 1:nq_add)
      !nq = nq+nq_add
      DEALLOCATE(xq_aux)
    ELSE
      IF(nq/=0) CALL errore("expand_qlist", "i do not understand", 1)
      ALLOCATE(xq(3,nq_add))
      xq = xq_add(:, 1:nq_add)
      !nq = nq_add
    ENDIF
    RETURN
  END SUBROUTINE
  SUBROUTINE expand_wlist(nq, w, nq_add)
    INTEGER,INTENT(in) :: nq 
    INTEGER,INTENT(in) :: nq_add
    REAL(DP),ALLOCATABLE,INTENT(inout) :: w(:)
    !
    REAL(DP),ALLOCATABLE :: w_aux(:)
    IF(allocated(w)) THEN
      !print*,"got", size(w), nq, nq_add 
      IF(size(w)/=nq) CALL errore('expand_wlist', 'wrong size', 1)
      ALLOCATE(w_aux(nq))
      w_aux(1:nq) = w(1:nq)
      DEALLOCATE(w)
      ALLOCATE(w(nq+nq_add))
      w(1:nq)           = w_aux
      w(nq+1:nq+nq_add) = 1._dp/DBLE(nq_add)
      DEALLOCATE(w_aux)
    ELSE
      !print*,"got", size(w), nq, nq_add 
      IF(nq/=0) CALL errore("expand_qlist", "i do not understand", 1)
      ALLOCATE(w(nq_add))
      w = 1._dp/DBLE(nq_add)
    ENDIF
    RETURN
  END SUBROUTINE


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
    ENDIF
    !
  END SUBROUTINE setup_simple_grid
  ! \/o\________\\\_________________________________________/^>
  ! Create a line of nq q-point from xqi to xqf, the list is appended
  ! to the grid 
  SUBROUTINE setup_path(xqi, xqf, nq, path)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nq
    REAL(DP),INTENT(in) :: xqi(3), xqf(3)
    TYPE(q_grid),INTENT(inout) :: path
    !
    INTEGER :: i, n0
    REAL(DP),PARAMETER:: eps = 1.d-8
    REAL(DP),ALLOCATABLE :: auxq(:,:)
    REAL(DP) :: dq(3)
    !
    IF(nq==0) RETURN
    !
    IF(path%nq >0)THEN
      ALLOCATE(auxq(3,path%nq))
      auxq = path%xq
      DEALLOCATE(path%xq)
      ! If xqi is the same as the last point in the list, do not
      ! add it again. FIXME: check for q = q'+G, use optional input S info
      IF(SUM(ABS(auxq(:,path%nq)-xqi))<eps) THEN
        n0 = path%nq 
        path%nq = path%nq + nq -1
      ELSE
        n0 = path%nq +1
        path%nq = path%nq + nq
      ENDIF
      !
      ALLOCATE(path%xq(3,path%nq))
      path%xq(:,1:size(auxq,2)) = auxq(:,1:size(auxq,2))
      DEALLOCATE(auxq)
    ELSE
      n0 = 1
      path%nq = nq
      ALLOCATE(path%xq(3,path%nq))
    ENDIF
    !
    IF(nq>1)THEN
      dq = (xqf-xqi)/(nq-1)
      DO i = n0, path%nq
        path%xq(:,i) = xqi + dq * (i-n0)
      ENDDO
    ELSE
        ! With nq=1 xqf is appended to the list, xqi is ignored
        ! note that repeated points are dropped
        path%xq(:,path%nq) = xqf
    ENDIF
    !
    ! paths are not suitable for integrals, set al weights to zero just in case
    IF(allocated(path%w)) DEALLOCATE(path%w)
    ALLOCATE(path%w(path%nq))
    path%w = 0._dp 
    !
  END SUBROUTINE setup_path
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE prepare_q_basis(qgrid, qbasis, nconf, T, S, fc2)
    USE fc2_interpolate,    ONLY : freq_phq_safe, bose_phq
    USE ph_velocity,        ONLY : velocity_proj
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
    !
    ALLOCATE(U(S%nat3,S%nat3))
    
    DO iq = 1,qbasis%nq
      qbasis%c(:,:,iq) = velocity_proj(S, fc2, qgrid%xq(:,iq))
      CALL  freq_phq_safe(qgrid%xq(:,iq), S, fc2, qbasis%w(:,iq), U)
      DO it = 1, nconf
        CALL  bose_phq(T(it), S%nat3, qbasis%w(:,iq), qbasis%be(:,it,iq))
      ENDDO
    ENDDO

    ALLOCATE(qbasis%b(3,nconf,S%nat3,qbasis%nq))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iq,nu,it,ix)
!$OMP DO COLLAPSE(4)
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
  SUBROUTINE qbasis_x_times_y(z, x, y, nconf, nat3, nq)
    IMPLICIT NONE
    !
    REAL(DP),INTENT(out):: z(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: x(3, nconf, nat3, nq)
    REAL(DP),INTENT(in) :: y(3, nconf, nat3, nq)
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
        z(ix,it,nu,iq) = x(ix,it,nu,iq)*y(ix,it,nu,iq)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !
  END SUBROUTINE qbasis_x_times_y
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

  !
END MODULE q_grids





