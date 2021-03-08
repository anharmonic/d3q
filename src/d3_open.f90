!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto - Nov 2011
MODULE d3_open
  PUBLIC :: diropn_d3, close_d3, reopn_d3
  PUBLIC :: listu_d3 ! debug
  !
  PRIVATE
  !
  TYPE d3_opened_file
    INTEGER :: u, recl
    CHARACTER(len=256) :: ext  ! should be enough...
    CHARACTER(len=256) :: tmp
  END TYPE d3_opened_file
  INTEGER,PARAMETER    :: d3_max_opened_files = 200
  INTEGER              :: d3_opened_files     = 0
  TYPE(d3_opened_file) :: d3_files(d3_max_opened_files)

CONTAINS
!
SUBROUTINE diropn_d3(unit, extension, recl, exst, tmp_dir)
  USE io_files,   ONLY : diropn
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(in) :: extension ! input: name of the file to open
  CHARACTER(len=*),INTENT(in) :: tmp_dir   ! if present it is used as tmp_dir
  INTEGER,INTENT(in)  :: unit, recl ! unit of the file to open, length of the records
  LOGICAL,INTENT(out) :: exst
  !
  CHARACTER(len=9),PARAMETER :: sub = 'diropn_d3'
  INTEGER :: idx
  !
  CALL diropn(unit, extension, recl, exst, tmp_dir)
  idx = scanu_d3(unit)
  IF(idx>0) CALL errore(sub, 'Unit was not closed properly, or still opened', 2)

  d3_opened_files = d3_opened_files+1
  IF(d3_opened_files>d3_max_opened_files) THEN
    d3_opened_files = d3_max_opened_files
    CALL listu_d3(-1)
    CALL errore(sub, 'Too many opened files, cannot open: "'//TRIM(extension)//'"',1)
  ENDIF
  !
  d3_files(d3_opened_files)%u    = unit
  d3_files(d3_opened_files)%recl = recl
  d3_files(d3_opened_files)%ext  = extension
  d3_files(d3_opened_files)%tmp  = tmp_dir
  !
  RETURN
  !
END SUBROUTINE diropn_d3

SUBROUTINE close_d3(unit, keep, safe)
  IMPLICIT NONE
  INTEGER,INTENT(in)  :: unit
  LOGICAL,OPTIONAl,INTENT(in) :: keep, safe
  INTEGER :: idx
  LOGICAL :: opnd, safe_
  CHARACTER(len=8),PARAMETER :: sub = 'close_d3'
  CHARACTER(len=6) :: stat
  !
  safe_=.false.
  IF(present(safe)) safe_=safe
  !
  stat="KEEP"
  IF(present(keep)) THEN
    IF(.not.keep) stat='DELETE'
  ENDIF
  !
  INQUIRE(UNIT = unit, OPENED = opnd)
  !
  IF(.not.opnd) THEN
    IF(safe_) RETURN
    CALL errore(sub, 'Unit was not opened', unit)
  ENDIF
  !
  idx = popu_d3(unit)
  IF(idx<0) THEN
    CALL errore(sub, 'Unit was not registered',unit)
  ENDIF
  CLOSE(unit, status=stat)
  !
  RETURN
  !
END SUBROUTINE close_d3

SUBROUTINE reopn_d3(unit) !, extension, recl, exst, tmp_dir)
  USE io_files,   ONLY : diropn
  IMPLICIT NONE
  INTEGER,INTENT(in)  :: unit
  INTEGER :: idx
  LOGICAL :: exst
  CHARACTER(len=8),PARAMETER :: sub = 'reopn_d3'
  !
  idx = scanu_d3(unit)
  IF(idx<0) CALL errore(sub, 'Cannot reopen a file that is not opened', unit)
  !
  CLOSE(unit)
  CALL diropn(d3_files(idx)%u, d3_files(idx)%ext, d3_files(idx)%recl, exst, d3_files(idx)%tmp)
  IF(.not.exst) CALL errore(sub,'File must exist "'//TRIM(d3_files(idx)%ext)//'"', unit)
  !
  RETURN
  !
END SUBROUTINE reopn_d3
!
SUBROUTINE listu_d3(unit)
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER,INTENT(in)  :: unit
  INTEGER :: i
  !
  IF(unit<0) THEN
    DO i = 1,d3_opened_files
      CALL print_d3_file(stdout, d3_files(i))
    ENDDO
  ELSE
    CALL print_d3_file(stdout, d3_files(scanu_d3(unit)))
  ENDIF
  !
  RETURN
  !
END SUBROUTINE listu_d3
FUNCTION scanu_d3(unit)
  IMPLICIT NONE
  INTEGER :: scanu_d3
  INTEGER,INTENT(in)  :: unit
  !
  INTEGER :: i
  !
  scanu_d3 = -1
  !
  DO i = 1,d3_opened_files
    IF( d3_files(i)%u == unit) THEN
      scanu_d3 = i
      RETURN
    ENDIF
  ENDDO
  !
  RETURN
  !
END FUNCTION scanu_d3
FUNCTION popu_d3(unit)
  IMPLICIT NONE
  INTEGER :: popu_d3
  INTEGER,INTENT(in)  :: unit
  !
  INTEGER :: i,LAST
  !
  i = scanu_d3(unit)
  popu_d3 = i
  IF(i<0) RETURN
  !
  CALL destroy_d3_file(d3_files(i))
  LAST = d3_opened_files
  d3_opened_files = d3_opened_files-1
  !
  IF(d3_opened_files==0 .or. i==LAST) RETURN
  !
  CALL copy_d3_file(d3_files(LAST), d3_files(i))
  CALL destroy_d3_file(d3_files(LAST))
  !
  RETURN
  !
END FUNCTION popu_d3
SUBROUTINE destroy_d3_file(f)
  IMPLICIT NONE
  TYPE(d3_opened_file) :: f
  f%u    =-1
  f%recl =-1
  f%ext  =''
  f%tmp  =''
  RETURN
END SUBROUTINE destroy_d3_file
SUBROUTINE copy_d3_file(i, o)
  IMPLICIT NONE
  TYPE(d3_opened_file) :: i,o
  o%u    = i%u
  o%recl = i%recl
  o%ext  = i%ext
  o%tmp  = i%tmp
  RETURN
END SUBROUTINE copy_d3_file
SUBROUTINE print_d3_file(u, f)
  IMPLICIT NONE
  INTEGER :: u
  TYPE(d3_opened_file) :: f
  WRITE(u,'(5x,"u:",1i4,3x,"recl:",1i10,3x,"ext: ",a30,3x,"tmp: ",a)') f%u, f%recl, TRIM(f%ext(1:30)), TRIM(f%tmp)
  RETURN
END SUBROUTINE print_d3_file

END MODULE d3_open
