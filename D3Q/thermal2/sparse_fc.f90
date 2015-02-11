!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! Module that uses object oriented features of Fortran 2003 to deal with both
! regular grid and sparse representations of Force constants in a transparent way.
! May not compile with some less maintained compilers, should work at least with 
! gfortran, xlf, ifort and pgi compiler. Does not compile with g95.
!
! It is not written defensively, e.g. writing and interpolating FCs do not check 
! if the data has been read. May be improved in the future.
!
MODULE sparse_fc
  USE kinds,            ONLY : DP
  USE parameters,       ONLY : ntypx
  USE io_global,        ONLY : stdout
  USE interp_fc,        ONLY : ph_system_info
  ! Base implementation: all methods stop with error
  TYPE,ABSTRACT :: forceconst3
    ! q points
    INTEGER :: n_R = 0
    INTEGER :: i_00 = -1 ! set to the index where both R2 and R3 are zero
    INTEGER :: nq(3) ! initial q-point grid size (unused)
    INTEGER,ALLOCATABLE  :: yR2(:,:), yR3(:,:) ! crystalline coords  3*n_R
    REAL(DP),ALLOCATABLE :: xR2(:,:), xR3(:,:) ! cartesian coords    3*n_R
    CONTAINS
      procedure(fftinterp_mat3_error),deferred :: interpolate
      procedure :: destroy      => destroy_fc3_
      procedure(read_fc3_error),deferred :: read
      procedure(write_fc3_error),deferred :: write
      procedure(div_mass_fc3_error),deferred :: div_mass
  END TYPE forceconst3
  ! \/o\________\\\_________________________________________/^>
  ! Interfaces to the deferred subroutines:
  ABSTRACT INTERFACE
  SUBROUTINE fftinterp_mat3_error(fc, xq2,xq3, nat3, D)
    USE kinds,     ONLY : DP
    IMPORT forceconst3 
    CLASS(forceconst3),INTENT(in) :: fc
    INTEGER,INTENT(in)   :: nat3
    REAL(DP),INTENT(in) :: xq2(3), xq3(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
  END SUBROUTINE fftinterp_mat3_error
  END INTERFACE
  !
  ABSTRACT INTERFACE
  SUBROUTINE read_fc3_error(fc, filename, S)
    USE interp_fc,        ONLY : ph_system_info
    IMPORT forceconst3 
    CHARACTER(len=*),INTENT(in)        :: filename
    TYPE(ph_system_info),INTENT(inout) :: S ! = System
    CLASS(forceconst3),INTENT(inout)   :: fc
  END SUBROUTINE read_fc3_error
  END INTERFACE
  !
  ABSTRACT INTERFACE
  SUBROUTINE write_fc3_error(fc, filename, S)
    USE interp_fc,        ONLY : ph_system_info
    IMPORT forceconst3 
    CHARACTER(len=*),INTENT(in)     :: filename
    TYPE(ph_system_info),INTENT(in) :: S ! = System
    CLASS(forceconst3),INTENT(in)   :: fc
  END SUBROUTINE write_fc3_error
  END INTERFACE
  !
  ABSTRACT INTERFACE
  SUBROUTINE div_mass_fc3_error(fc, S)
    USE interp_fc,        ONLY : ph_system_info
    IMPORT forceconst3 
    CLASS(forceconst3),INTENT(inout) :: fc
    TYPE(ph_system_info),INTENT(in)  :: S ! = System
  END SUBROUTINE div_mass_fc3_error
  END INTERFACE
  ! \/o\________\\\_________________________________________/^>
  !
  ! First implementation: not sparse
  TYPE,EXTENDS(forceconst3) :: grid
    REAL(DP),ALLOCATABLE :: FC(:,:,:,:) ! 3*nat,3*nat,3*nat, n_R
    !
    CONTAINS
      procedure :: interpolate  => fftinterp_mat3_grid
      procedure :: destroy      => destroy_fc3_grid
      procedure :: read         => read_fc3_grid
      procedure :: write        => write_fc3_grid
      procedure :: div_mass     => div_mass_fc3_grid
  END TYPE grid
  INTERFACE grid
    MODULE PROCEDURE create_fc3_grid
  END INTERFACE
  ! \/o\________\\\_________________________________________/^>
  ! Second implementation: sparse
  TYPE forceconst3_sparse_helper
    REAL(DP),ALLOCATABLE :: fc(:)    ! the matrix elements
    INTEGER,ALLOCATABLE  :: idx(:,:) ! the index 3*nat,3*nat,3*nat
  END TYPE forceconst3_sparse_helper
  !
  TYPE,EXTENDS(forceconst3) :: sparse
    INTEGER,ALLOCATABLE  :: n_terms(:) ! number of non-zero elements for each R2,R3
    TYPE(forceconst3_sparse_helper),ALLOCATABLE :: dat(:)
    !
    CONTAINS
      procedure :: interpolate  => fftinterp_mat3_sparse
      procedure :: destroy      => destroy_fc3_sparse
      procedure :: read         => read_fc3_sparse
      procedure :: write        => write_fc3_sparse
      procedure :: div_mass     => div_mass_fc3_sparse
  END TYPE sparse
  INTERFACE sparse
    MODULE PROCEDURE create_fc3_sparse
  END INTERFACE
  ! \/o\________\\\_________________________________________/^>
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  !
  ! Returns a pointer to FC matrices, call as:
  !    TYPE(ph_system_info),INTENT(in)   :: S
  !    CLASS(forceconst3),POINTER :: fc
  !     ...
  !    fc => read_fc3("file_mat3", S)
  ! This will read both sparse and grid files, there is room for improvement but not much.
  !
  FUNCTION read_fc3(filename, S) RESULT(fc)
    USE input_fc, ONLY : read_system
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(out)   :: S ! = System
    CLASS(forceconst3),POINTER :: fc    
    !
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: unit, ios
    CHARACTER(len=32) :: buf
    !
    ! Scan the type of file:
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore("read_fc3","opening '"//TRIM(filename)//"'", 1)
    READ(unit, '(a32)') buf
    CLOSE(unit)
    !
    IF(buf=="sparse representation") THEN
      ALLOCATE(sparse:: fc)
    ELSE
      ALLOCATE(grid::   fc)
    ENDIF
    !
    ! Use class specific read method:
    CALL fc%read(filename, S)
    !
    RETURN
    !
  END FUNCTION read_fc3  
  !
  ! Placeholder creators:
  TYPE(grid) FUNCTION create_fc3_grid()
    RETURN
  END FUNCTION
  TYPE(sparse) FUNCTION create_fc3_sparse()
    RETURN
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fc3_grid_to_sparse(nat,fc, sfc, thr)
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(grid),INTENT(in)    :: fc
    TYPE(sparse),INTENT(out) :: sfc
    REAL(DP),OPTIONAL,INTENT(in) :: thr
    !
    REAL(DP),PARAMETER :: eps12 = 1.d-12
    INTEGER :: iR, n_found
    INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
    REAL(DP) :: eps
    !
    eps = eps12
    IF(present(thr)) eps = thr
    !
    sfc%n_R  = fc%n_R
    sfc%i_00 = fc%i_00
    sfc%nq   = fc%nq
    !
    ALLOCATE(sfc%yR2(3,sfc%n_R), sfc%yR3(3,sfc%n_R))
    ALLOCATE(sfc%xR2(3,sfc%n_R), sfc%xR3(3,sfc%n_R))
    sfc%yR2 = fc%yR2;    sfc%yR3 = fc%yR3
    sfc%xR2 = fc%xR2;    sfc%xR3 = fc%xR3
    ALLOCATE(sfc%dat(sfc%n_R))
    ALLOCATE(sfc%n_terms(sfc%n_R))
    
    DO iR = 1, sfc%n_R
      sfc%n_terms(iR) = COUNT(  ABS(fc%FC(:,:,:,iR))>eps )
      ALLOCATE(sfc%dat(iR)%idx(3,sfc%n_terms(iR)))
      ALLOCATE(sfc%dat(iR)%fc(sfc%n_terms(iR)))
    ENDDO
    
    DO iR = 1, sfc%n_R
      !
      n_found = 0
      !
      DO na3=1,nat 
      DO na2=1,nat 
      DO na1=1,nat
        DO j3=1,3     
        jn3 = j3 + (na3-1)*3
        DO j2=1,3     
        jn2 = j2 + (na2-1)*3
        DO j1=1,3
        jn1 = j1 + (na1-1)*3
          !
          IF(ABS(fc%FC(jn1,jn2,jn3,iR))>eps)THEN
            n_found = n_found + 1
            IF(n_found > sfc%n_terms(iR)) CALL errore("fc3_grid_to_sparse", "too many terms", iR)
            sfc%dat(iR)%idx(1,n_found) = jn1
            sfc%dat(iR)%idx(2,n_found) = jn2
            sfc%dat(iR)%idx(3,n_found) = jn3
            sfc%dat(iR)%FC(n_found) = fc%FC(jn1,jn2,jn3,iR)
          ENDIF
          !
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      !
      IF(n_found /= sfc%n_terms(iR)) CALL errore("fc3_grid_to_sparse", "wrong number of terms", iR)
      !
    ENDDO
    !
  END SUBROUTINE fc3_grid_to_sparse
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat3_grid(fc, xq2,xq3, nat3, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)   :: nat3
    CLASS(grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq2(3), xq3(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,arg,phase) REDUCTION(+: D)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:,:,:) = D(:,:,:) + phase * fc%fc(:,:,:, i)
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE fftinterp_mat3_grid
  !
  SUBROUTINE fftinterp_mat3_sparse(fc, xq2,xq3, nat3, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)   :: nat3
    CLASS(sparse),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq2(3), xq3(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i,j
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,arg,phase) REDUCTION(+: D)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      DO j = 1, fc%n_terms(i)
        D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
          = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
            + phase * fc%dat(i)%fc(j)
      ENDDO
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE fftinterp_mat3_sparse
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE read_fc3_sparse(fc, filename, S)
    USE input_fc, ONLY : read_system
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    CLASS(sparse),INTENT(inout) :: fc
    !
    CHARACTER(len=13),PARAMETER :: sub = "read_fc3_sparse"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(len=32) :: buf
    !
    INTEGER :: i, j, i_, j_
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    READ(unit, '(a32)') buf
    IF(buf/="sparse representation") CALL errore(sub, "cannot read this format", 1)
    !
    WRITE(stdout,*) "** Reading sparse FC3 file"
    CALL read_system(unit, S)
    !
    READ(unit, *) fc%nq
    READ(unit, *) fc%n_R
    WRITE(stdout,*) "   Original FC3 grid:", fc%nq
    WRITE(stdout,*) "   Number of R:      ", fc%n_R
    !
    ALLOCATE(fc%n_terms(fc%n_R))
    ALLOCATE(fc%yR2(3,fc%n_R), fc%yR3(3,fc%n_R))
    DO i = 1,fc%n_R
      READ(unit, *) fc%n_terms(i),  fc%yR2(:,i), fc%yR3(:,i)
    ENDDO
    !
    ALLOCATE(fc%dat(fc%n_R))
    !
    DO i = 1, fc%n_R
      !
      ALLOCATE(fc%dat(i)%idx(3,fc%n_terms(i)))
      ALLOCATE(fc%dat(i)%fc(fc%n_terms(i)))
      !
      DO j = 1, fc%n_terms(i)
        READ(unit,*) fc%dat(i)%idx(:,j), fc%dat(i)%fc(j)
      ENDDO
    ENDDO
    !
    CLOSE(unit)
    ! Prepare the lists of R in cartesian coords
    ALLOCATE(fc%xR2(3,fc%n_R), fc%xR3(3,fc%n_R))
    fc%xR2 = REAL(fc%yR2, kind=DP)
    CALL cryst_to_cart(fc%n_R,fc%xR2,S%at,1)
    fc%xR3 = REAL(fc%yR3, kind=DP)
    CALL cryst_to_cart(fc%n_R,fc%xR3,S%at,1)
    !
  END SUBROUTINE read_fc3_sparse
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_fc3_sparse(fc, filename, S)
    USE input_fc, ONLY : write_system
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    CLASS(sparse),INTENT(in) :: fc
    !
    CHARACTER(len=14),PARAMETER :: sub = "write_fc3_sparse"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    !
    INTEGER :: i, j, n_digits(3)
    CHARACTER(len=64) :: cformat
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    WRITE(unit, '(a)') "sparse representation"
    !
    CALL write_system(unit, S)
    WRITE(unit, '(3i9)') fc%nq
    WRITE(unit, '(3i9)') fc%n_R
    !
    n_digits(1) = CEILING(LOG10(DBLE(MAXVAL(fc%n_terms)))) ! no need to leave space in front
    n_digits(2) = CEILING(LOG10(DBLE(MAXVAL(fc%yR2))))+2 ! +2 to account minus signs
    n_digits(3) = CEILING(LOG10(DBLE(MAXVAL(fc%yR3))))+2
    cformat = '(i'//int_to_char(n_digits(1))// &
              ',3i'//int_to_char(n_digits(2))// &
              ',3i'//int_to_char(n_digits(3))//')'
    DO i = 1,fc%n_R
      WRITE(unit, cformat) fc%n_terms(i),  fc%yR2(:,i), fc%yR3(:,i)
    ENDDO
    !
    ! Writing with the minimum number of white spaces saves loads of disk space (unless you zip)
    ! keeping everything aligned for readability
    n_digits(2) = CEILING(LOG10(DBLE(3*S%nat)))+1
    n_digits(1) = n_digits(2)-1
    !
    cformat = '(i'//int_to_char(n_digits(1))// &
              ',2i'//int_to_char(n_digits(2))//',1pe25.15)'
    !
    DO i = 1, fc%n_R
      DO j = 1, fc%n_terms(i)
        !
        WRITE(unit,cformat) fc%dat(i)%idx(:,j), fc%dat(i)%fc(j)
        !
      ENDDO
    ENDDO
    !
    CLOSE(unit)
    !
  END SUBROUTINE write_fc3_sparse
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE destroy_fc3_sparse(fc)
    USE input_fc, ONLY : write_system
    IMPLICIT NONE
    CLASS(sparse),INTENT(inout) :: fc
    INTEGER :: i
    !
    DEALLOCATE(fc%yR2)
    DEALLOCATE(fc%yR3)
    DEALLOCATE(fc%xR2)
    DEALLOCATE(fc%xR3)
    DO i = 1, fc%n_R
      !
      DEALLOCATE(fc%dat(i)%idx)
      DEALLOCATE(fc%dat(i)%fc)
      !
    ENDDO
    DEALLOCATE(fc%dat)
    !
  END SUBROUTINE destroy_fc3_sparse
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE div_mass_fc3_sparse (fc,S)
    USE kinds, only : DP
    IMPLICIT NONE
    CLASS(sparse),INTENT(inout) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j
    !
    IF(.not.ALLOCATED(S%sqrtmm1)) &
      call errore('div_mass_fc3_sparse', 'missing sqrtmm1, call aux_system first', 1)
    
    DO i = 1, fc%n_R
      DO j = 1, fc%n_terms(i)
        fc%dat(i)%fc(j) = fc%dat(i)%fc(j) &
                         *S%sqrtmm1( fc%dat(i)%idx(1,j) ) &
                         *S%sqrtmm1( fc%dat(i)%idx(2,j) ) &
                         *S%sqrtmm1( fc%dat(i)%idx(3,j) )
      ENDDO
    ENDDO
    !
  END SUBROUTINE div_mass_fc3_sparse
  ! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE read_fc3_grid(fc, filename, S)
    USE input_fc, ONLY : read_system
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    CLASS(grid),INTENT(inout) :: fc
    !
    CHARACTER(len=13),PARAMETER :: sub = "read_fc3_grid"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
    INTEGER :: na1_, na2_, na3_, j1_, j2_, j3_
    INTEGER :: n_R, i
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    WRITE(stdout,*) "** Reading full grid FC3 file"
    !
    CALL read_system(unit, S)
    !
    READ(unit, *) fc%nq
    WRITE(stdout,*) "   Original FC3 grid:", fc%nq
    WRITE(stdout,*) "   Number of R:      ", fc%n_R
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
    DO na3=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3
      DO j3=1,3     
      jn3 = j3 + (na3-1)*3
          !
          READ(unit,*) j1_, j2_, j3_, na1_, na2_, na3_
          IF ( ANY((/na1,na2,na3,j1,j2,j3/) /= (/na1_,na2_,na3_,j1_,j2_,j3_/)) ) THEN
            print*, (/na1,na2,na3,j1,j2,j3/)
            print*, (/na1_,na2_,na3_,j1_,j2_,j3_/)
            CALL errore(sub,'not matching na1,na2,na3,j1,j2,j3 in file "'//TRIM(filename)//"'",1)
          ENDIF
          !
          READ(unit,*) n_R
          IF ( fc%n_R == 0) THEN 
            IF(allocated(fc%yR2) .or. allocated(fc%xR2) .or. &
               allocated(fc%yR3) .or. allocated(fc%xR3) ) &
              CALL errore(sub, 'some element are already allocated', 1)
            !
            ALLOCATE(fc%yR2(3,n_R), fc%yR3(3,n_R))
            ALLOCATE(fc%xR2(3,n_R), fc%xR3(3,n_R))
            ALLOCATE(fc%FC(3*S%nat,3*S%nat,3*S%nat,n_R))
            fc%n_R = n_R
          ELSE
            IF(n_R/=fc%n_R) CALL errore(sub, "cannot read this fc format",1)
          ENDIF
          !
          DO i = 1, fc%n_R
            READ(unit,*) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i)
            ! also, find the index of R=0
            IF( ALL(fc%yR2(:,i)==0) .and. ALL(fc%yR3(:,i)==0) ) fc%i_00 = i
          ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    ! Prepare the lists of R in cartesian coords
    fc%xR2 = REAL(fc%yR2, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR2,S%at,1)
    fc%xR3 = REAL(fc%yR3, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR3,S%at,1)
    !
  END SUBROUTINE read_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_fc3_grid(fc, filename, S)
    USE input_fc, ONLY : write_system
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    CLASS(grid),INTENT(in) :: fc
    !
    CHARACTER(len=14),PARAMETER :: sub = "write_fc3_grid"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    ! format:
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    INTEGER :: n_digits(3)
    CHARACTER(len=64) :: cformat, cformat2
    !
    INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
    INTEGER :: i
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    CALL write_system(unit, S)
    !
    WRITE(unit, '(3i9)') fc%nq
    !
    n_digits(1) = CEILING(LOG10(DBLE(S%nat)))+1
    cformat='(i1,2i2,x,3i'//int_to_char(n_digits(1))//')'
    !
    n_digits(2) = CEILING(LOG10(DBLE((MAXVAL(fc%yR2)))))+2 ! +2 to account minus signs
    n_digits(1) = n_digits(2)-1 ! no space at lines beginning
    n_digits(3) = CEILING(LOG10(DBLE((MAXVAL(fc%yR3)))))+2
    !
    cformat2 = '( i'//int_to_char(n_digits(1))// &
               ',2i'//int_to_char(n_digits(2))// &
               ',3i'//int_to_char(n_digits(3))//',1pe25.15)'
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
    DO na3=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3
      DO j3=1,3     
      jn3 = j3 + (na3-1)*3
          !
          WRITE(unit,cformat) j1, j2, j3, na1, na2, na3
          !
          WRITE(unit,'(i9)') fc%n_R
          !
          DO i = 1, fc%n_R
            WRITE(unit,cformat2) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i)
          ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    !
  END SUBROUTINE write_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE destroy_fc3_grid(fc)
    USE input_fc, ONLY : write_system
    IMPLICIT NONE
    CLASS(grid),INTENT(inout) :: fc
    !
    DEALLOCATE(fc%yR2)
    DEALLOCATE(fc%yR3)
    DEALLOCATE(fc%xR2)
    DEALLOCATE(fc%xR3)
    DEALLOCATE(fc%FC)
    !
  END SUBROUTINE destroy_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE div_mass_fc3_grid (fc,S)
    USE kinds, only : DP
    IMPLICIT NONE
    CLASS(grid),INTENT(inout) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j, k, i_R
    !
    IF(.not.ALLOCATED(S%sqrtmm1)) &
      call errore('div_mass_fc3_grid', 'missing sqrtmm1, call aux_system first', 1)
    
    DO i_R = 1, fc%n_R
      DO k = 1, S%nat3
      DO j = 1, S%nat3
      DO i = 1, S%nat3
        fc%FC(i, j, k, i_R) = fc%FC(i, j, k, i_R) &
                    * S%sqrtmm1(i)*S%sqrtmm1(j)*S%sqrtmm1(k)
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE div_mass_fc3_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE destroy_fc3_(fc)
    USE input_fc, ONLY : write_system
    IMPLICIT NONE
    CLASS(forceconst3),INTENT(inout) :: fc
    !
    DEALLOCATE(fc%yR2)
    DEALLOCATE(fc%yR3)
    DEALLOCATE(fc%xR2)
    DEALLOCATE(fc%xR3)
    !
  END SUBROUTINE destroy_fc3_
  !
!   ! \/o\________\\\_________________________________________/^>
!   SUBROUTINE fftinterp_mat3_error(fc, xq2,xq3, nat3, D)
!     USE constants, ONLY : tpi
!     IMPLICIT NONE
!     INTEGER,INTENT(in)   :: nat3
!     CLASS(forceconst3),INTENT(in) :: fc
!     REAL(DP),INTENT(in) :: xq2(3), xq3(3)
!     COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
!     CALL errore("fftinterp_mat3_error", "not implemented for this class", 1)
!     D = 0._dp
!   END SUBROUTINE fftinterp_mat3_error
!   ! \/o\________\\\_________________________________________/^>
!   SUBROUTINE read_fc3_error(fc, filename, S)
!     IMPLICIT NONE
!     CHARACTER(len=*),INTENT(in)        :: filename
!     TYPE(ph_system_info),INTENT(inout) :: S ! = System
!     CLASS(forceconst3),INTENT(inout)   :: fc
!     CALL errore("read_fc3_error", "not implemented for this class", 1)
!   END SUBROUTINE read_fc3_error
!   !
!   SUBROUTINE write_fc3_error(fc, filename, S)
!     IMPLICIT NONE
!     CHARACTER(len=*),INTENT(in)     :: filename
!     TYPE(ph_system_info),INTENT(in) :: S ! = System
!     CLASS(forceconst3),INTENT(in)   :: fc
!     CALL errore("write_fc3_error", "not implemented for this class", 1)
!   END SUBROUTINE write_fc3_error
!   !
!   SUBROUTINE div_mass_fc3_error(fc, S)
!     IMPLICIT NONE
!     CLASS(forceconst3),INTENT(inout) :: fc
!     TYPE(ph_system_info),INTENT(in)  :: S ! = System
!     CALL errore("div_mass_fc3_error", "not implemented for this class", 1)
!   END SUBROUTINE div_mass_fc3_error
  
END MODULE
