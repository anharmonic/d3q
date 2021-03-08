!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!

! <<^V^\\=========================================//-//-//========//O\\//
MODULE input_fc
  !
  USE kinds,            ONLY : DP
!   USE parameters,       ONLY : ntypx
  USE mpi_thermal,      ONLY : ionode
  USE ph_system
#include "mpi_thermal.h"
  !
    ! \/o\________\\\_________________________________________/^>
  TYPE forceconst2_grid
    ! q points
    INTEGER :: n_R = 0, i_0 = -1
    INTEGER,ALLOCATABLE  :: yR(:,:) ! crystalline coords  3*n_R
    REAL(DP),ALLOCATABLE :: xR(:,:) ! cartesian coords    3*n_R
    REAL(DP),ALLOCATABLE :: FC(:,:,:) ! 3*nat,3*nat, n_R
    INTEGER :: nq(3) ! initial grid size, only kept for reconstruction purposes
  END TYPE forceconst2_grid
  !
  CONTAINS
  ! <<^V^\\=========================================//-//-//========//O\\//
  
  ! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE read_fc2(filename, S, fc)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)          :: filename
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    !
    CHARACTER(len=8),PARAMETER :: sub = "read_fc2"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER :: na1,  na2,  j1,  j2, jn1, jn2, jn3
    INTEGER :: na1_, na2_, j1_, j2_
    INTEGER :: n_R, i
    !
    ioWRITE(stdout,*) "** Reading FC2 file: ", TRIM(filename)
    fc%i_0 = -1
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    CALL read_system(unit, S)
    !
    READ(unit, *) jn1, jn2, jn3
    ioWRITE(stdout,*) "Original FC2 grid:", jn1, jn2, jn3
    fc%nq(1) = jn1; fc%nq(2) = jn2; fc%nq(3) = jn3
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3            
          !
          READ(unit,*) j1_, j2_, na1_, na2_
          IF ( ANY((/na1,na2,j1,j2/) /= (/na1_,na2_,j1_,j2_/)) ) &
            CALL errore(sub,'not matching na1,na2,j1,j2',1)
          !
          READ(unit,*) n_R
          IF ( fc%n_R == 0) CALL allocate_fc2_grid(n_R, S%nat, fc)
          !
          DO i = 1, fc%n_R
            READ(unit,*) fc%yR(:,i), fc%FC(jn1,jn2,i)
            ! also, find the index of R=0
            IF( ALL(fc%yR(:,i)==0) ) fc%i_0 = i
          ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    IF(fc%i_0 < 0) CALL errore("read_fc2","could not find i_0",1)
    ! Compute list of R in cartesian coords
    fc%xR = REAL(fc%yR, kind=DP)
    CALL cryst_to_cart(n_R,fc%xR,S%at,1)
    !
  END SUBROUTINE read_fc2
  ! <<^V^\\=========================================//-//-//========//O\\//
  SUBROUTINE write_fc2(filename, S, fc)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)       :: filename
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    CHARACTER(len=8),PARAMETER :: sub = "write_fc2"
    !
    INTEGER :: unit, ios
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER :: na1,  na2,  j1,  j2, jn1, jn2, jn3
    INTEGER :: n_R, i
    !
    unit = find_free_unit()
    OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
    IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
    !
    CALL write_system(unit, S)
    !
    WRITE(unit, '(3i9)') fc%nq
    !
    DO na1=1,S%nat
    DO na2=1,S%nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3            
          !
          WRITE(unit,'(2i9,3x,2i9)') j1, j2, na1, na2
          WRITE(unit,'(i9)') fc%n_R
          !
          DO i = 1, fc%n_R
            WRITE(unit,'(3i4,1pe25.15)') fc%yR(:,i), fc%FC(jn1,jn2,i)
!             WRITE(unit,'(3i4,2x,1pe18.11)') fc%yR(:,i), fc%FC(jn1,jn2,i)
          ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    !
    CLOSE(unit)
    !
  END SUBROUTINE write_fc2
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE allocate_fc2_grid(n_R, nat, fc)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n_R, nat
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    CHARACTER(len=16),PARAMETER :: sub = "allocate_fc2_grid"
    !
    IF(allocated(fc%yR) .or. allocated(fc%xR) .or. allocated(fc%FC)) &
      CALL errore(sub, 'some element is already allocated', 1)
    !
    ALLOCATE(fc%yR(3,n_R))
    ALLOCATE(fc%xR(3,n_R))
    ALLOCATE(fc%FC(3*nat,3*nat,n_R))
    fc%n_R = n_R
    !
  END SUBROUTINE allocate_fc2_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE deallocate_fc2_grid(fc)
    IMPLICIT NONE
    TYPE(forceconst2_grid),INTENT(inout) :: fc
    CHARACTER(len=18),PARAMETER :: sub = "deallocate_fc2_grid"
    !
    IF(.NOT.(allocated(fc%yR) .and. allocated(fc%xR) .and. allocated(fc%FC))) &
      CALL errore(sub, 'some element is already deallocated', 1)
    !
    DEALLOCATE(fc%yR)
    DEALLOCATE(fc%xR)
    DEALLOCATE(fc%FC)
    fc%n_R = 0
    fc%i_0 = -1
    !
  END SUBROUTINE deallocate_fc2_grid
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE div_mass_fc2 (S,fc)
    USE kinds, only : DP
    IMPLICIT NONE
    TYPE(forceconst2_grid) :: fc
    TYPE(ph_system_info)   :: S
    !
    INTEGER :: i, j, i_R
    !
    IF(.not.ALLOCATED(S%sqrtmm1)) &
      call errore('div_mass_fc2', 'missing sqrtmm1, call aux_system first', 1)
    
    DO i_R = 1, fc%n_R
      DO j = 1, S%nat3
      DO i = 1, S%nat3
        fc%FC(i, j, i_R) = fc%FC(i, j, i_R) * S%sqrtmm1(i)*S%sqrtmm1(j)
      ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE div_mass_fc2
  FUNCTION multiply_mass_dyn (S,dyn) RESULT(dyn_mass)
    USE kinds, only : DP
    !USE input_fc,         ONLY : ph_system_info
    IMPLICIT NONE
    TYPE(ph_system_info)   :: S
    COMPLEX(DP),INTENT(in) :: dyn(S%nat3, S%nat3)
    COMPLEX(DP) :: dyn_mass(S%nat3, S%nat3)
    !
    INTEGER :: i, j
    !
    DO j = 1, S%nat3
    DO i = 1, S%nat3
      dyn_mass(i,j) = dyn(i,j)/(S%sqrtmm1(i)*S%sqrtmm1(j))
    ENDDO
    ENDDO
    !
  END FUNCTION multiply_mass_dyn
  FUNCTION div_mass_dyn (S,dyn) RESULT(dyn_mass)
    USE kinds, only : DP
    !USE input_fc,         ONLY : ph_system_info
    IMPLICIT NONE
    TYPE(ph_system_info)   :: S
    COMPLEX(DP),INTENT(in) :: dyn(S%nat3, S%nat3)
    COMPLEX(DP) :: dyn_mass(S%nat3, S%nat3)
    !
    INTEGER :: i, j
    !
    DO j = 1, S%nat3
    DO i = 1, S%nat3
      dyn_mass(i,j) = dyn(i,j)*(S%sqrtmm1(i)*S%sqrtmm1(j))
    ENDDO
    ENDDO
    !
  END FUNCTION div_mass_dyn
  ! write dynamical matrix on file in the ph.x format
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_dyn (filename, xq, d2, S)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: filename
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    COMPLEX(DP),INTENT(in) :: d2(S%nat3, S%nat3)!  the dynamical matrix
    REAL(DP),INTENT(in) :: xq (3)! the q vector
    !
    COMPLEX(DP) :: phi(3,3,S%nat,S%nat)!  the dynamical matrix
    INTEGER :: iudyn  ! unit number, number of atom in the unit cell
    INTEGER :: i, j, na, nb, icar, jcar
    INTEGER,EXTERNAL :: find_free_unit
    CHARACTER(len=9) :: cdate, ctime
 
    !
    DO i = 1, S%nat3
      na = (i - 1) / 3 + 1
      icar = i - 3 * (na - 1)
      DO j = 1, S%nat3
          nb = (j - 1) / 3 + 1
          jcar = j - 3 * (nb - 1)
          phi(icar, jcar, na, nb) = d2 (i, j)
      ENDDO
    ENDDO
    
    iudyn = find_free_unit()
    OPEN(unit=iudyn, file=filename, status='UNKNOWN', action="WRITE")
    WRITE(iudyn, *) "Dynamical matrix file from thermal2"
    WRITE(iudyn, *)
    CALL write_system(iudyn, S, matdyn=.true.)
    CALL write_dyn_on_file_ph(xq, phi, S%nat, iudyn)
    CALL date_and_tim( cdate, ctime )
    WRITE(iudyn, *)
    WRITE(iudyn, *) "File generated by thermal2 codes Lorenzo Paulatto ", cdate, " ", ctime
    WRITE(iudyn, *)
    CLOSE(iudyn)
    
  END SUBROUTINE write_dyn
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_dyn_on_file_ph(xq, phi, nat, iudyn)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! input variables
    INTEGER,INTENT(in) :: iudyn, nat ! unit number, number of atom in the unit cell
    COMPLEX(DP),INTENT(in) :: phi (3, 3, nat, nat)!  the dynamical matrix
    REAL(DP),INTENT(in) :: xq (3)! the q vector
    !
    INTEGER :: na, nb, icar, jcar
    WRITE (iudyn, "(/,5x,'Dynamical  Matrix in cartesian axes', &
        &       //,5x,'q = ( ',3f14.9,' ) ',/)") (xq (icar), icar = 1, 3)
        
    DO na = 1, nat
    DO nb = 1, nat
      WRITE (iudyn, '(2i5)') na, nb
      DO icar = 1, 3
!             WRITE (iudyn, '(3e24.12)') (phi(icar,jcar,na,nb), jcar=1,3)
          write (iudyn, '(3(2f14.8,2x))') (phi(icar,jcar,na,nb), jcar=1,3)
      ENDDO
    ENDDO
    ENDDO
  END SUBROUTINE write_dyn_on_file_ph
  ! \/o\________\\\_________________________________________/^>
END MODULE input_fc








