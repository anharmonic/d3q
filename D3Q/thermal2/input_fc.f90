!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE input_fc
  !
  USE kinds,            ONLY : DP
  USE parameters,       ONLY : ntypx
  USE mpi_thermal,      ONLY : ionode
#include "para_io.h"
  !
  ! \/o\________\\\_________________________________________/^>
  TYPE ph_system_info
    ! atoms
    INTEGER              :: ntyp
    REAL(DP)             :: amass(ntypx)
    REAL(DP)             :: amass_variance(ntypx)
    CHARACTER(len=3  )   :: atm(ntypx)
    ! atoms basis
    INTEGER              :: nat
    REAL(DP),ALLOCATABLE :: tau(:,:), zeu(:,:,:)
    INTEGER, ALLOCATABLE :: ityp(:)
    ! unit cell, and reciprocal
    INTEGER              :: ibrav
    CHARACTER(len=9)     :: symm_type
    REAL(DP)             :: celldm(6), at(3,3), bg(3,3)
    REAL(DP)             :: omega
    ! phonon switches (mostly unused here)
    REAL(DP)             :: epsil(3,3)
    LOGICAL              :: lrigid
    ! ''''''''''''''''''''''''''''''''''''''''
    ! auxiliary quantities:
    REAL(DP),ALLOCATABLE :: sqrtmm1(:) ! 1/sqrt(amass)
    INTEGER :: nat3, nat32, nat33
  END TYPE
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
  
  FUNCTION same_system(S,Z) RESULT (same)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in) :: S,Z
    LOGICAL :: same
    REAL(DP),PARAMETER :: eps = 1.d-6
    LOGICAL :: verbose
    !
    verbose = ionode.and..FALSE. !just shut up
    
    !
    ! NOT checking : atm, amass, symm_type
    !
    same = .true.
    same = same .and. (S%ntyp == Z%ntyp)
    IF(.not.same.and.verbose) WRITE(stdout,*) "ntyp", S%ntyp, Z%ntyp
    same = same .and. (S%nat == Z%nat)
    IF(.not.same.and.verbose) WRITE(stdout,*) "nat", S%nat, Z%nat
    same = same .and. (S%ibrav == Z%ibrav)
    IF(.not.same.and.verbose) WRITE(stdout,*) "ibrav", S%ibrav, Z%ibrav
    same = same .and. ALL( S%ityp(1:S%ntyp) == Z%ityp(1:Z%ntyp))
    IF(.not.same.and.verbose) WRITE(stdout,*) "ityp", S%ityp, Z%ityp
    
    IF(allocated(S%tau).and.allocated(Z%tau)) THEN
      same = same .and. ALL( ABS(S%tau -Z%tau) < eps)
      IF(.not.same.and.verbose) WRITE(stdout,*) "tau", S%tau, Z%tau
    ENDIF
    IF(allocated(S%zeu).and.allocated(Z%zeu)) THEN
      same = same .and. ALL( ABS(S%zeu -Z%zeu) < eps)
      IF(.not.same.and.verbose) WRITE(stdout,*) "zeu", S%zeu, Z%zeu
    ENDIF

    same = same .and. ALL( ABS(S%celldm -Z%celldm) < eps)
    IF(.not.same.and.verbose) WRITE(stdout,*) "celldm", S%celldm, Z%celldm
    same = same .and. ALL( ABS(S%at -Z%at) < eps)
    IF(.not.same.and.verbose) WRITE(stdout,*) "at", S%at, Z%at
    same = same .and. ALL( ABS(S%bg -Z%bg) < eps)
    IF(.not.same.and.verbose) WRITE(stdout,*) "bg", S%bg, Z%bg
    same = same .and. ( ABS(S%omega -Z%omega) < eps)
    IF(.not.same.and.verbose) WRITE(stdout,*) "omega", S%omega, Z%omega

!     same = same .and. (S%lrigid .or. Z%lrigid)
!     same = same .and. ALL( ABS(S%epsil -Z%epsil) < eps)
!     IF(.not.same) ioWRITE(stdout,*) "epsil", S%epsil, Z%epsil
    
  END FUNCTION same_system
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE read_system(unit, S)
    USE more_constants, ONLY : MASS_DALTON_TO_RY
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(out)   :: S ! = System
    INTEGER,INTENT(in) :: unit
    !
    CHARACTER(len=11),PARAMETER :: sub = "read_system"
    !
    INTEGER :: ios, dummy
    !
    INTEGER :: nt, na
    !
    READ(unit,*,iostat=ios) S%ntyp, S%nat, S%ibrav, S%celldm(1:6)
    IF(ios/=0) CALL errore(sub,"reading S%ntyp, S%nat, S%ibrav, S%celldm(1:6)", 1)
    !
    IF(S%ibrav==0)THEN
      READ(unit,*,iostat=ios) S%at(1:3,1:3)
      IF(ios/=0) CALL errore(sub,"reading S%at(1:3,1:3)", 1)
    ENDIF
    !
    ! generate at, bg, volume
    IF (S%ibrav /= 0) THEN
      CALL latgen(S%ibrav, S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega)
      S%at = S%at / S%celldm(1)  !  bring at from bohr to units of alat 
    ENDIF
    CALL volume(S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega)
    CALL recips(S%at(:,1), S%at(:,2), S%at(:,3), S%bg(:,1), S%bg(:,2), S%bg(:,3))
    !
    DO nt = 1, S%ntyp
      READ(unit,*,iostat=ios) dummy, S%atm(nt), S%amass(nt)
      IF(ios/=0) CALL errore(sub,"reading nt, S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
    IF(any(S%amass(1:S%ntyp)<500._dp)) THEN
        ioWRITE(stdout,*) "WARNING: Masses seem to be in Dalton units: rescaling"
        ioWRITE(stdout,*) "old:", S%amass(1:S%ntyp)
        S%amass(1:S%ntyp) = S%amass(1:S%ntyp)* MASS_DALTON_TO_RY
        ioWRITE(stdout,*) "new:", S%amass(1:S%ntyp)
    ENDIF
    S%amass_variance = 0._dp
    !
    ALLOCATE(S%ityp(S%nat), S%tau(3,S%nat))
    DO na = 1, S%nat
      READ(unit,*,iostat=ios) dummy, S%ityp(na), S%tau(:,na)
      IF(ios/=0) CALL errore(sub,"reading na, S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
  END SUBROUTINE read_system
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_system(unit, S)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    INTEGER,INTENT(in) :: unit
    !
    CHARACTER(len=11),PARAMETER :: sub = "write_system"
    !
    INTEGER :: ios
    !
    INTEGER :: nt, na
    !
    WRITE(unit,"(3i9,6f25.16)",iostat=ios) S%ntyp, S%nat, S%ibrav, S%celldm(1:6)
    IF(ios/=0) CALL errore(sub,"writing S%ntyp, S%nat, S%ibrav, S%celldm(1:6)", 1)
    !
    IF(S%ibrav==0)THEN
      WRITE(unit,"(3f25.16)",iostat=ios) S%at(1:3,1:3)
      IF(ios/=0) CALL errore(sub,"writing S%at(1:3,1:3)", 1)
    ENDIF
    !
    ! generate at, bg, volume
    DO nt = 1, S%ntyp
      WRITE(unit,'(i9,2x,a5,f25.16)',iostat=ios) nt, "'"//S%atm(nt)//"'", S%amass(nt)
      IF(ios/=0) CALL errore(sub,"writing nt, S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
    IF (.not.allocated(S%ityp) .or. .not. allocated(S%tau))&
      CALL errore(sub, 'missing ityp, tau', 1)
    DO na = 1, S%nat
      WRITE(unit,'(2i9,3f25.16)',iostat=ios) na, S%ityp(na), S%tau(:,na)
      IF(ios/=0) CALL errore(sub,"writing na, S%atm(nt), S%amass(nt)", nt)
    ENDDO
    !
  END SUBROUTINE write_system
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE destroy_system(unit, S)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    INTEGER,INTENT(in) :: unit
    !
    S%ntyp = 0
    S%nat  = 0
    S%ibrav = 0
    S%celldm(1:6) = 0._dp
    S%atm(:) = ""
    S%amass(:) = 0._dp
    IF(allocated(S%ityp)) DEALLOCATE(S%ityp)
    IF(allocated(S%tau))  DEALLOCATE(S%tau)
    !
  END SUBROUTINE destroy_system
  ! \/o\________\\\_________________________________________/^>
  ! compute redundant but useful quantities
  SUBROUTINE aux_system(S)
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(inout)   :: S ! = System
    INTEGER :: i, na
    S%nat3 = 3*S%nat
    S%nat32 = S%nat3**2
    S%nat33 = S%nat3**3
    
    IF(allocated(S%sqrtmm1)) CALL errore("aux_system","should not be called twice",1)
    ALLOCATE(S%sqrtmm1(S%nat3))
    DO i = 1,S%nat3
      na = (i-1)/3 +1
      S%sqrtmm1(i) = 1._dp/SQRT(S%amass(S%ityp(na)))
    ENDDO
    
  END SUBROUTINE aux_system
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
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE print_citations_linewidth
    WRITE(stdout,*)
    WRITE(stdout,'(5x,a)') &
        " ",&
        "For any third order calculation please cite:",&
        "  Lorenzo Paulatto, Francesco Mauri, and Michele Lazzeri",&
        "  Phys. Rev. B 87, 214303 – Published 7 June 2013"
    WRITE(stdout,*)
    WRITE(stdout,'(5x,a)') &
        " ",&
        "For thermal transport calculations please cite:",&
        "  Giorgia Fugallo, Michele Lazzeri, Lorenzo Paulatto, and Francesco Mauri",&
        "  Phys. Rev. B 88, 045430 – Published 17 July 2013"
    WRITE(stdout,*)

  END SUBROUTINE print_citations_linewidth
  ! \/o\________\\\_________________________________________/^>
END MODULE input_fc








