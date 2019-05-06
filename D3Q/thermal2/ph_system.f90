!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE ph_system
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
#include "mpi_thermal.h"  
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
    REAL(DP)             :: omega, alat, tpiba
    ! phonon switches (mostly unused here)
    REAL(DP)             :: epsil(3,3)
    LOGICAL              :: lrigid
    ! ''''''''''''''''''''''''''''''''''''''''
    ! auxiliary quantities:
    REAL(DP),ALLOCATABLE :: sqrtmm1(:) ! 1/sqrt(amass)
    INTEGER :: nat3, nat32, nat33
  END TYPE
  !
  CONTAINS
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
    
    
    IF(allocated(S%ityp).and.allocated(Z%ityp)) THEN
      same = same .and. ALL( S%ityp(1:S%ntyp) == Z%ityp(1:Z%ntyp))
      IF(.not.same.and.verbose) WRITE(stdout,*) "ityp", S%ityp, Z%ityp
    ELSE
      WRITE(stdout,*) "One system is missing ions types!"
    ENDIF
    
    IF(allocated(S%tau).and.allocated(Z%tau)) THEN
      same = same .and. ALL( ABS(S%tau -Z%tau) < eps)
      IF(.not.same.and.verbose) WRITE(stdout,*) "tau", S%tau, Z%tau
    ELSE
      WRITE(stdout,*) "One system is missing ions base!"
    ENDIF
    !
    ! Not checking: lrigid and co., because they are not in mat3R files
!     IF(S%lrigid.or.Z%lrigid)THEN
!       IF(.not.S%lrigid) THEN
!         !WRITE(stdout,*) "Only Z has rigid"
!       ELSE IF(.not.Z%lrigid) THEN
!         !WRITE(stdout,*) "Only S has rigid"
!       ELSE
!         IF(allocated(S%zeu).and.allocated(Z%zeu)) THEN
!           same = same .and. ALL( ABS(S%zeu -Z%zeu) < eps)
!           IF(.not.same.and.verbose) WRITE(stdout,*) "zeu", S%zeu, Z%zeu
!         ENDIF
!       ENDIF
!     ENDIF

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
    USE constants,      ONLY : tpi
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(out)   :: S ! = System
    INTEGER,INTENT(in) :: unit
    !
    CHARACTER(len=11),PARAMETER :: sub = "read_system"
    !
    INTEGER :: ios, dummy
    !
    INTEGER :: nt, na
    CHARACTER(len=256) cdummy
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
!    S%at = S%at*S%celldm(1)
!    S%bg = tpi*S%bg/S%celldm(1)
!    S%celldm(1) = 1._dp
    !
    S%alat  = S%celldm(1)
    S%tpiba = tpi/S%alat
    ioWRITE(*,'(5x,a,f12.6)') "alat (bohr)", S%alat
    ioWRITE(*,'(5x,a,f12.6)') "tpiba (1/bohr)", S%tpiba
    ioWRITE(*,'(5x,a,3f12.6,a)') "at_1 = (", S%at(:,1)," )"
    ioWRITE(*,'(5x,a,3f12.6,a)') "at_2 = (", S%at(:,2)," )"
    ioWRITE(*,'(5x,a,3f12.6,a)') "at_3 = (", S%at(:,3)," )"
    ioWRITE(*,'(5x,a,3f12.6,a)') "bg_1 = (", S%bg(:,1)," )"
    ioWRITE(*,'(5x,a,3f12.6,a)') "bg_2 = (", S%bg(:,2)," )"
    ioWRITE(*,'(5x,a,3f12.6,a)') "bg_3 = (", S%bg(:,3)," )"
    ioWRITE(*,'(5x,a,f12.6)') "volume (bohr^3)", S%omega
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
    READ(unit,*,iostat=ios) S%lrigid
    !print*, "lrigid", S%lrigid
    IF(ios/=0) CALL errore(sub,"reading rigid", 1)
    IF(S%lrigid)THEN
!       READ(unit,*,iostat=ios) cdummy
      READ(unit,*,iostat=ios) S%epsil
      IF(ios/=0) CALL errore(sub,"reading epsilon (2)", 1)
!    print*, "epsil", S%epsil
!       READ(unit,*) cdummy
!       READ(unit,*) cdummy
      ALLOCATE(S%zeu(3,3,S%nat))
      DO na = 1, S%nat
        READ(unit,*,iostat=ios) S%zeu(:,:,na)
!      print*, "zeu", na, S%zeu(:,:,na)
        IF(ios/=0) CALL errore(sub,"reading zeu (2)", na)
!         READ(unit,*) cdummy
      ENDDO
     
!     ELSE
!       ALLOCATE(S%zeu(0,0,0))
    ENDIF
    !
  END SUBROUTINE read_system
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE write_system(unit, S, matdyn)
    USE functions, ONLY : default_if_not_present
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S ! = System
    INTEGER,INTENT(in) :: unit
    LOGICAL,OPTIONAL,INTENT(in) :: matdyn 
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
      IF(present(matdyn))THEN
        IF(matdyn)  WRITE(unit, *) "Basis vectors"
      ENDIF
      WRITE(unit,"(3f25.16)",iostat=ios) S%at(1:3,1:3)
      IF(ios/=0) CALL errore(sub,"writing S%at(1:3,1:3)", 1)
    ENDIF
    !
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
   
!    WRITE(*,*) "present matdyn", present(matdyn), matdyn, default_if_not_present(matdyn,.false.) 
    IF(.not. default_if_not_present(.false.,matdyn)) THEN
    WRITE(unit,'(5x,l)',iostat=ios) S%lrigid
    IF(ios/=0) CALL errore(sub,"writing rigid", 1)
    IF(S%lrigid)THEN
      WRITE(unit,'(3(3f25.16,/))',iostat=ios) S%epsil
      IF(ios/=0) CALL errore(sub,"writing epsilon (2)", 1)
      DO na = 1, S%nat
        WRITE(unit,'(3(3f25.16,/))',iostat=ios) S%zeu(:,:,na)
      IF(ios/=0) CALL errore(sub,"writing zeu (2)", na)
      ENDDO
    ENDIF
    ENDIF
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
      S%sqrtmm1(i) = 1._dp/DSQRT(S%amass(S%ityp(na)))
    ENDDO

  END SUBROUTINE aux_system
  
END MODULE 
