!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE import3py_module
  USE kinds, ONLY : DP
#include "mpi_thermal.h"
  IMPLICIT NONE
  !
  CONTAINS
  !
  ! Read FORCE_CONSTANTS_3RD files created by ShengBTE "thirdorder" 
  ! (and possibly Phono3py) as documented at 
  ! <https://bitbucket.org/sousaw/shengbte/src/master/README.md>
  SUBROUTINE read_3py(filename, fc, S)
    USE fc3_interpolate,  ONLY : grid
    USE input_fc,         ONLY : ph_system_info
    USE constants,        ONLY : ANGSTROM_AU, RYTOEV
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in) :: filename
    TYPE(grid), INTENT(inout) :: fc
    TYPE(ph_system_info),INTENT(in) :: S
    CHARACTER(len=8) :: sub = "read_3py"
    INTEGER :: u, iR, iR2, n_R
    REAL(DP) :: R2(3), R3(3), F
    INTEGER :: i1,i2,i3, j1,j2,j3, na1,na2,na3, jn1,jn2,jn3
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    REAL(DP),PARAMETER :: F_FACTOR= 1._dp/RYTOEV/ANGSTROM_AU**3
    !
    !n_R = NINT(SQRT(DBLE(n_R)))
    !IF(n_R**2/=n_R) CALL errore(sub, "problem with R and n_R", 1)
    
    u = find_free_unit()
    OPEN(unit=u, file=filename, status='OLD', form='FORMATTED')
      READ(u,*) n_R
      ALLOCATE(fc%yR2(3,n_R), fc%yR3(3,n_R))
      ALLOCATE(fc%xR2(3,n_R), fc%xR3(3,n_R))
      ALLOCATE(fc%FC(3*S%nat,3*S%nat,3*S%nat,n_R))
      fc%xR2 = 0._dp
      fc%xR3 = 0._dp
      fc%yR2 = 0
      fc%yR3 = 0
      fc%FC  = 0._dp
      fc%n_R = n_R
      WRITE(stdout,*) "reading "//TRIM(int_to_char(n_R))//" blocks"
      DO iR = 1, n_R
        READ(u,*)
        READ(u,*) iR2
        IF(iR/=iR2) CALL errore(sub,"i does not match i2", 1)
        READ(u,*) R2
        READ(u,*) R3
        R2 = R2*ANGSTROM_AU/S%celldm(1)
        R3 = R3*ANGSTROM_AU/S%celldm(1)
        !WRITE(stdout,*) iR
        !WRITE(stdout,'(2(3f12.6,5x))') R2,R3
        fc%xR2(:,iR) = R2
        fc%xR3(:,iR) = R3
        CALL cryst_to_cart(1,R2,S%bg,-1)
        CALL cryst_to_cart(1,R3,S%bg,-1)
        !WRITE(stdout,'(2(3f12.6,5x))') R2,R3
        fc%yR2(:,iR) = NINT(R2)
        fc%yR3(:,iR) = NINT(R3)
        
        READ(u,*) na1, na2, na3
        DO j1=1,3
        jn1 = j1 + (na1-1)*3
        DO j2=1,3     
        jn2 = j2 + (na2-1)*3
        DO j3=1,3     
        jn3 = j3 + (na3-1)*3
          READ(u,*) i1,i2,i3, F
          IF(i1/=j1 .or. i2/=j2 .or. i3/=j3) CALL errore(sub, "unexpected i/=j", 1)
          fc%FC(jn1,jn2,jn3,iR) = F*F_FACTOR
        ENDDO
        ENDDO
        ENDDO
      ENDDO
    CLOSE(u)
  END SUBROUTINE
  ! \/o\________\\\________________\\/\_________________________/^>
  FUNCTION guess_nq(nR, yR2, yR3) RESULT(nq)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nR
    INTEGER,INTENT(in) :: yR2(3,nR), yR3(3,nR)
    INTEGER :: nq(3)
    nq(1) = MAXVAL( (/ MAXVAL(ABS(yR2(1,:))),MAXVAL(ABS(yR3(1,:))), &
                       MAXVAL(ABS(yR3(1,:)+yR3(1,:))) /) )
    nq(2) = MAXVAL( (/ MAXVAL(ABS(yR2(2,:))),MAXVAL(ABS(yR3(2,:))), &
                       MAXVAL(ABS(yR3(2,:)+yR3(2,:))) /) )
    nq(3) = MAXVAL( (/ MAXVAL(ABS(yR2(3,:))),MAXVAL(ABS(yR3(3,:))), &
                       MAXVAL(ABS(yR3(3,:)+yR3(3,:))) /) )
    RETURN
  END FUNCTION
  !
END MODULE


PROGRAM import3py
  USE fc3_interpolate, ONLY : grid
  USE input_fc,        ONLY : read_system, ph_system_info
  USE ph_system,       ONLY : aux_system
  USE import3py_module
  USE f3_bwfft
  IMPLICIT NONE
  INTEGER :: narg, i
  CHARACTER(len=512) :: argv
  CHARACTER(len=256) :: filename, fileout
  TYPE(grid) :: fc, fcb
  TYPE(ph_system_info) :: S
  !
  INTEGER :: nq(3), nq_trip
  TYPE(d3_list),ALLOCATABLE :: d3grid(:)
  !
  filename = "FORCE_CONSTANTS_3RD"
  fileout  = "mat3R.3py"

  OPEN(unit=999,file="mat2R",action='READ',status='OLD')
  CALL read_system(999, S)
  CALL aux_system(S)
  CLOSE(999)
  !
  !CALL scan_3py(filename, S)
  CALL read_3py(filename, fc, S)
  CALL fc%write(fileout, S)
  !
  nq = guess_nq(fc%n_R, fc%yR2, fc%yR3)
  WRITE(*,*) "Guessing grid size", nq
  !
  narg = command_argument_count()
  IF(narg>=3)THEN
    DO i = 1,3
      CALL get_command_argument(i,argv)
      READ(argv,*) nq(i)
    ENDDO
    WRITE(*,*) "Grid size from input:", nq
  ELSE
    WRITE(*,*) "To change the grid size by hand use syntax:"
    WRITE(*,*) "  d3_import3py.x nqx nqy nqz"
  ENDIF
  WRITE(*,*)
  
!  print*, "nq:", nq
!   IF(ANY(fc%yR(1,:)>nq(1)) .or. ANY(fc%yR(1,:)<-nq(1))) PRINT*, "problem nq1"
!   IF(ANY(fc%yR(1,:)>nq(1)) .or. ANY(fc%yR(1,:)<-nq(1))) PRINT*, "problem nq1"
!   IF(ANY(fc%yR(1,:)>nq(1)) .or. ANY(fc%yR(1,:)<-nq(1))) PRINT*, "problem nq1"
  
!  nq = (/2, 2, 2/)
  nq_trip = (nq(1)*nq(2)*nq(3))**2
  ALLOCATE(d3grid(nq_trip))
  CALL regen_fwfft_d3(nq, nq_trip, S, d3grid, fc)
  
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fcb)
  CALL fcb%write("mat3R.3py_b", S)
  !
  !
END PROGRAM
!
