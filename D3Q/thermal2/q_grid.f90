!
! Copyright Lorenzo Paulatto, Giorgia Fugallo 2013 - released under the CeCILL licence v 2.1
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE q_grid

  USE kinds,     ONLY : DP
  
  TYPE q_grid_type
    INTEGER :: n(3), nq
    REAL(DP),ALLOCATABLE :: xq(:,:)
    !REAL(DP),ALLOCATABLE :: w(:)
  END TYPE
  !
!   TYPE q_3grid_type
!     INTEGER :: n(3), nq
!     TYPE(q_grid_type) :: sub(3)
!     !REAL(DP),ALLOCATABLE :: w(:)
!   END TYPE

  CONTAINS
!   ! \/o\________\\\_________________________________________/^>
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
!   END SUBROUTINE setup_symmetry
!   ! \/o\________\\\_________________________________________/^>
  
  SUBROUTINE setup_simple_grid(S, n1,n2,n3, grid)
    USE input_fc, ONLY : ph_system_info 
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    INTEGER,INTENT(in) :: n1,n2,n3
    TYPE(q_grid_type),INTENT(inout) :: grid
!     REAL(DP),OPTIONAl,INTENT(in) :: xq0
    !
    INTEGER :: i,j,k, idx
    grid%n(1) = n1
    grid%n(2) = n2
    grid%n(3) = n3
    grid%nq = n1*n2*n3
    !
    ALLOCATE(grid%xq(3,grid%nq))
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
    CALL cryst_to_cart(grid%nq,grid%xq,S%bg, +1)
    !
    !IF(present(xq0)) grid%xq = grid%xq + xq0
    !
  END SUBROUTINE setup_simple_grid
  !
END MODULE q_grid





