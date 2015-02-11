!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! This is a didactical example on how to use object oriented fortran 2003 features.
! It is quite pointless.
MODULE vector_class
  USE kinds, ONLY : DP
  !
  PRIVATE
  PUBLIC :: vector
  ! self-growing vector
  INTERFACE vector
    MODULE PROCEDURE new_vector
  END INTERFACE vector
  !
  TYPE vector
    ! q points
    INTEGER :: size_ = 0   ! number of elements in data array
    INTEGER,PRIVATE :: space_ = 0  ! currently allocated space
    REAL(DP),POINTER :: x(:) => null()
    CONTAINS
      procedure,nopass :: new => new_vector         ! create a new vector, optional initial size
      procedure        :: init => initialize_vector ! initialize this vector, optional initial size
      procedure :: free    ! clean up data
      procedure :: append  ! add one element at the end
      procedure :: insert  ! insert(idx,item): insert item before idx
      procedure :: pop     ! remove last element and return it
      procedure :: reverse ! reverse(i,j) exchange value of element i with element j
      procedure :: get => vector_get       ! return all data
      procedure :: size => vector_size     ! return number of elements
      procedure :: space__ => vector_space ! return allocated size (should not be used explicitly)
      procedure :: assign_vector           ! assign_vector(to,from)  copy "from" over "to"
      procedure :: write => write_vector   ! writes content of vector with specified format
      procedure :: extend          ! increase available space
      !procedur :: compact ! reduce allocation to minimal possible size (TODO)
      !
      generic :: assignment(=) => assign_vector !
      final   :: destroy                        ! calls self%free
      generic :: write(formatted) => write      ! calls self%write (unsupported by intel compiler!)
    END TYPE vector
    !
    ! When running out of space inrcease the allocation size by this factor:
    INTEGER,PARAMETER,PRIVATE :: grow_factor = 2
    ! Initial allocated size:
    INTEGER,PARAMETER,PRIVATE :: initial_size = 128

  CONTAINS
  ! <<^self^\\=========================================//-//-//========//O\\//
  FUNCTION vector_get(self) RESULT(X)
    CLASS(vector),INTENT(in) :: self
    REAL(DP) :: X(self%size_)
    X = self%x(1:self%size_)
    RETURN
  END FUNCTION vector_get
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE destroy(self)
    IMPLICIT NONE
    TYPE(vector),INTENT(inout) :: self
    CALL self%free()
    RETURN
  END SUBROUTINE destroy
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE free(self)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    IF(associated(self%x)) DEALLOCATE(self%x)
    self%x => null()
    self%space_ = 0
    self%size_ = 0
    RETURN
  END SUBROUTINE free
  ! \/o\________\\\_________________________________________/^>
  ! aka %init
  SUBROUTINE initialize_vector(self, space_)
    IMPLICIT NONE
    CLASS(vector),INTENT(out) :: self
    INTEGER,INTENT(in),OPTIONAL :: space_
    self%size_  = 0
    self%space_ = 0
    self%x     => null()
    IF(present(space_)) THEN
      CALL self%extend(space_)
    ELSE
      CALL self%extend(initial_size)
    ENDIF
    !print*, "vector initialized to size", self%space_
    RETURN
  END SUBROUTINE initialize_vector
  ! \/o\________\\\_________________________________________/^>
  ! aka %new and default constructor
  TYPE(vector) FUNCTION new_vector(space_) RESULT(self)
    IMPLICIT NONE
    ! Note: the "optional" status is passed down to %init
    INTEGER,INTENT(in),OPTIONAL :: space_
    CALL self%init(space_)
    RETURN
  END FUNCTION new_vector
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_vector(self, unit, iotype, v_list, iostat, iomsg)
    CLASS(vector),INTENT(in) :: self
    INTEGER, INTENT(IN) :: unit
    CHARACTER (LEN=*), INTENT(IN) :: iotype
    INTEGER, INTENT(in) :: v_list(:)
    INTEGER, INTENT(out) :: iostat
    CHARACTER (LEN=*), INTENT(inout) :: iomsg
    
    iostat = 0
    CALL vector_error("write_vector: this would work if Intel supported it")
    !WRITE(unit, iotype, iostat=iostat, iomsg=iomsg) self%x
!     iomsg="write_vector unsupported by intel compiler"
!     iostat=1
    
  END SUBROUTINE
  ! \/o\________\\\_________________________________________/^>
  ! aka the operator "="
  SUBROUTINE assign_vector(to, from)
    IMPLICIT NONE
    CLASS(vector),INTENT(out) :: to
    CLASS(vector),INTENT(in) :: from
    !
    !IF(to%space_ < from%space_) 
    !print*, "copying from to", from%space_, to%space_
    ! NOTE: we can use from%size_ or from%space_, I think the latter is neater
    CALL to%extend(from%size_)
    to%size_ = from%size_
    to%x(1:to%size_) = from%x(1:from%size_)
    RETURN
  END SUBROUTINE
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE append(self, item)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    REAL(DP),INTENT(in) :: item
    !
    self%size_ = self%size_+1
    IF(self%size_>self%space_) CALL self%extend(grow_factor*self%size_)
    self%x(self%size_) = item
    RETURN
      
  END SUBROUTINE append
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE reverse(self, i,j)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    INTEGER,INTENT(in) :: i,j
    REAL(DP) :: aux
    !
    IF(i<1.or.j<1.or.i>self%size_.or.j>self%size_) &
      CALL vector_error("reverse: bad index")
    aux = self%x(i)
    self%x(i) = self%x(j)
    self%x(j) = aux
    RETURN
      
  END SUBROUTINE reverse
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE extend(self, newsize)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    INTEGER,INTENT(in) :: newsize
    REAL(DP),POINTER   :: aux(:)
    !
    IF(newsize<self%space_) RETURN
    !
    IF(self%space_>0)THEN
      aux => self%x
      ALLOCATE(self%x(newsize))
      IF(self%space_>0) self%x(1:self%space_) = aux (1:self%space_)
      DEALLOCATE(aux)
    ELSE
      ALLOCATE(self%x(newsize))
    ENDIF
    self%space_ = newsize
    !print*, "vector resized to size", self%space_
    RETURN
  END SUBROUTINE extend
  ! \/o\________\\\_________________________________________/^>
  INTEGER FUNCTION vector_space(self)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    vector_space = self%space_
    RETURN
  END FUNCTION vector_space
  ! \/o\________\\\_________________________________________/^>
  INTEGER FUNCTION vector_size(self)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    vector_size = self%size_
    RETURN
  END FUNCTION vector_size
  ! \/o\________\\\_________________________________________/^>
  REAL(DP) FUNCTION pop(self, idx)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    INTEGER,OPTIONAL,INTENT(in) :: idx
    INTEGER :: item,i
    !
    IF(.not.present(idx))THEN
      pop = self%x(self%size_)
    ELSE
      IF(idx>self%size_ .or. idx<1) CALL vector_error("pop:bad index")
      pop = self%x(idx)
      IF(.not.(idx==self%size_))THEN
        self%x(idx:self%size_-1) = self%x(idx+1:self%size_)
      ENDIF
    ENDIF
    self%size_ = self%size_ -1
    RETURN
      
  END FUNCTION pop
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE insert(self, idx, item)
    IMPLICIT NONE
    CLASS(vector),INTENT(inout) :: self
    INTEGER,OPTIONAL,INTENT(in) :: idx
    REAL(DP),INTENT(in)         :: item
    !
    IF(idx<1) CALL vector_error("insert: bad index")
    self%size_ = MAX(self%size_+1,idx)
    IF(self%size_>self%space_) CALL self%extend(grow_factor*self%size_)
    !
    IF(.not.(idx==self%size_))THEN
      self%x(idx+1:self%size_) = self%x(idx:self%size_-1)
    ENDIF
    self%x(idx) = item
    RETURN
      
  END SUBROUTINE insert
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE vector_error(msg)
    CHARACTER(len=*) :: msg
    WRITE(0,"(3x,a)") "\/o\________\\\_________________________________________/^>"
    WRITE(0,"(3x,a)") "\/  selfector library encountered an error:"
    WRITE(0,"(3x,2a)") "\/  ", msg
    WRITE(0,"(3x,a)") "\/o\________\\\_________________________________________/^>"
    STOP
  END SUBROUTINE vector_error
  
END MODULE vector_class

! \\\///oOo\\\===========\.\.\==================================================//self^>>
! MODULE vector_extension
!   USE vector_class
!   USE kinds, ONLY : DP
! 
!   TYPE,EXTENDS(vector) :: matrix
!     INTEGER :: n
!     REAL(DP),POINTER :: x3(:,:) => null()
!     CONTAINS
!     procedure :: extend        ! increase available space
!     procedure :: append        ! increase available space
! 
!   END TYPE
!   
!   CONTAINS
!   ! \/o\________\\\_________________________________________/^>
!   SUBROUTINE extend(self, newsize)
!     IMPLICIT NONE
!     CLASS(matrix),INTENT(inout) :: self
!     INTEGER,INTENT(in) :: newsize
!     REAL(DP),POINTER   :: aux(:)
!     !
!     CALL self%vector%extend(self%n*newsize)
!     self%x3(1:self%n, 1:newsize) => self%x(1:self%n*newsize)
!     RETURN
!   END SUBROUTINE extend
!   ! \/o\________\\\_________________________________________/^>
!   SUBROUTINE append(self, item)
!     IMPLICIT NONE
!     CLASS(vector),INTENT(inout) :: self
!     REAL(DP),INTENT(in) :: item
!     !
!     self%size_ = self%size_+1
!     IF(self%size_>self%space_) CALL self%extend(grow_factor*self%size_)
!     self%x(self%size_) = item
!     RETURN
!       
!   END SUBROUTINE append
!   
! END MODULE vector_extension



