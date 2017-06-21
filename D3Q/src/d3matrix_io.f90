!
! Copyright (C) 2001-2011 Quantum-ESPRSSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the ionode_id directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3matrix_io
!-----------------------------------------------------------------------
  CHARACTER(len=5),PARAMETER :: format_version = "1.1.0"
!
CONTAINS
!
!
!-----------------------------------------------------------------------
FUNCTION d3matrix_filename(xq1, xq2, xq3, at, basename)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE dfile_autoname, ONLY : dfile_generate_name
  IMPLICIT NONE
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  REAL(DP),INTENT(in) :: at(3,3)
  CHARACTER(len=*),INTENT(in) :: basename
  CHARACTER(len=512) :: d3matrix_filename

  d3matrix_filename =  TRIM(basename) &
            // TRIM(dfile_generate_name(xq1(:), at, "_Q1")) &
            // TRIM(dfile_generate_name(xq2(:), at, "_Q2")) &
            // TRIM(dfile_generate_name(xq3(:), at, "_Q3"))
  RETURN
  !-----------------------------------------------------------------------
END FUNCTION d3matrix_filename
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_d3dyn_xml(basename, xq1,xq2,xq3, d3, ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  USE d3com,      ONLY : code_version => version
  USE io_global,  ONLY : stdout
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! input variables
  !
  CHARACTER(len=*),INTENT(in) :: basename
  REAL(DP),INTENT(in)         :: xq1(3), xq2(3), xq3(3)
  COMPLEX(DP),INTENT(in)      :: d3 (3, 3, 3, nat, nat, nat)
  !
  INTEGER,INTENT(in)          :: nat, ntyp, ibrav
  !
  REAL(DP),INTENT(in)         :: celldm(6), tau(3,nat), amass(ntypx), at(3,3)
  INTEGER,INTENT(in)          :: ityp(nat)
  CHARACTER(len=3),INTENT(in) :: atm(ntypx)

  !
  ! local variables
  !
  INTEGER :: i,j,k, u
  CHARACTER(len=iotk_attlenx) :: attr
  CHARACTER(len=iotk_namlenx) :: filename
  CHARACTER(LEN=9)  :: cdate, ctime

  filename = d3matrix_filename(xq1, xq2, xq3, at, basename)
  WRITE(stdout,"(5x,' -->',a)") TRIM(filename)
  !
  CALL iotk_free_unit(u)
  CALL date_and_tim( cdate, ctime )
  !
    CALL iotk_write_attr(attr, "code_version",   code_version, FIRST=.TRUE.)
    CALL iotk_write_attr(attr, "format_version", format_version)
    CALL iotk_write_attr(attr, "date",    cdate)
    CALL iotk_write_attr(attr, "time",    ctime)
  CALL iotk_open_write(u, filename, root="d3", attr=attr, skip_head=.true.)
    !___________________________________________________________
    !
      CALL iotk_write_attr(attr, "bravais_index", ibrav, FIRST=.TRUE.)
    CALL iotk_write_begin(u, "lattice", attr=attr)
    IF(ibrav==0)THEN
      CALL iotk_write_comment(u, "In units of alat which is in bohr units")
        CALL iotk_write_attr(attr, "alat", celldm(1), FIRST=.TRUE.)
      CALL iotk_write_dat(u, "unit_cell", at, columns=3, attr=attr)
    ELSE
      CALL iotk_write_comment(u, "see Doc/INPUT_PW.txt for description of ibrav and lattice parameters")
      CALL iotk_write_dat(u, "bravais_parameters", celldm, columns=3)
    ENDIF
    !
  CALL iotk_write_end(u, "lattice")
  !___________________________________________________________
  !
      CALL iotk_write_attr(attr, "number_of_species", ntyp,  newline=.true., FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "number_of_atoms",   nat,   newline=.true.)
  CALL iotk_write_begin(u, "atoms", attr=attr)
    !
      CALL iotk_write_comment(u, "positions are in alat units, cartesian coordinates")
    CALL iotk_write_dat(u, "atomic_positions", tau,   columns=3)
    CALL iotk_write_dat(u, "atomic_types",     ityp,  columns=1)
    CALL iotk_write_dat(u, "species_names",    atm(1:ntyp),   columns=ntyp)
      CALL iotk_write_comment(u, "masses are in Dalton atomic mass units (i.e. mass C^12=12)")
    CALL iotk_write_dat(u, "species_masses",   amass(1:ntyp), columns=ntyp)
    !
  CALL iotk_write_end(u, "atoms")
  !___________________________________________________________
  !
  CALL iotk_write_begin(u, "perturbation")
    CALL iotk_write_comment(u, "in cartesian units of 2pi/alat")
    CALL iotk_write_dat(u, "q1", xq1, columns=3)
    CALL iotk_write_dat(u, "q2", xq2, columns=3)
    CALL iotk_write_dat(u, "q3", xq3, columns=3)
  CALL iotk_write_end(u, "perturbation")
  !___________________________________________________________
  !
  CALL iotk_write_begin(u, 'd3matrix')
  CALL iotk_write_comment(u, "Each block contains the 3x3x3 tensor of cartesian displacements for atoms # i,j and k")
  CALL iotk_write_comment(u, "i is associated with q1, j with q2 and k with q3")
  CALL iotk_write_comment(u, "Matrix is NOT divided by the square root of masses")
  DO k = 1,nat
  DO j = 1,nat
  DO i = 1,nat
      CALL iotk_write_attr(attr, "atm_i", i, FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "atm_j", j)
      CALL iotk_write_attr(attr, "atm_k", k)
    CALL iotk_write_dat(u, "matrix", d3(:,:,:,i,j,k), attr=attr, columns=3)
  ENDDO
  ENDDO
  ENDDO
  CALL iotk_write_end(u, 'd3matrix')
  !
  CALL iotk_close_write(u)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE write_d3dyn_xml
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_d3dyn_xml(basename, xq1,xq2,xq3, d3, ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass,found,seek,file_format_version)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE parameters, ONLY : ntypx
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! input variables
  !
  CHARACTER(len=*),INTENT(in)           :: basename
  ! 1) if seek is present and .true. only the first part of the name (as specified in d3 input in fild3dyn)
  ! 2) if seek is not present or .false. the full name of the file to read
  REAL(DP),OPTIONAL,INTENT(inout)                   :: xq1(3), xq2(3), xq3(3)
  ! q vectors can be used in 2 ways:
  ! 1) if present and seek=.true. they are used to generate the complete filename
  ! 2) if present but seek=.false. or not present, they will be read from the file
  COMPLEX(DP),ALLOCATABLE,OPTIONAL,INTENT(inout)      :: d3 (:,:,:, :,:,:) ! (3,3,3, nat,nat,nat) the D3 matrix 
  !
  INTEGER,OPTIONAL,INTENT(out)          :: nat, ntyp       ! number of atoms and atomic species
  INTEGER,OPTIONAL,INTENT(out)          :: ibrav           ! lattice type
  !
  REAL(DP),OPTIONAL,INTENT(out)             :: celldm(6)   ! cell parameters depending on ibrav
  REAL(DP),OPTIONAL,INTENT(out)             :: at(3,3)     ! unit cell vectors NOTE: only filled if ibrav=0
  REAL(DP),ALLOCATABLE,OPTIONAL,INTENT(out) :: tau(:,:)    ! (3,nat) atomic positions, cartesians and alat
  REAL(DP),OPTIONAL,INTENT(out)             :: amass(ntypx)! (ntyp) mass of ions
  INTEGER,ALLOCATABLE,OPTIONAL,INTENT(out)  :: ityp(:)     ! (nat)  index of atomic types 
  CHARACTER(len=3),OPTIONAL,INTENT(out)     :: atm(ntypx)  ! (ntyp) atomic labels (es. Si)
  ! The version of the written file
  CHARACTER(len=5),OPTIONAL,INTENT(out) :: file_format_version
  ! Switches:
  LOGICAL,OPTIONAL,INTENT(out)              :: found
  ! if present and file is not found set it to false and return, instead of issuing an error
  LOGICAL,OPTIONAL,INTENT(in)               :: seek
  ! if present, xq1, xq2 and xq3 MUST be present, use them to generate the file name and find the file (if possible)

  !
  ! local variables
  !
  INTEGER :: i,j,k, l,m,n, u, ierr
  INTEGER :: i_nat, i_ntyp
  CHARACTER(len=iotk_attlenx) :: attr
  CHARACTER(len=iotk_namlenx) :: filename
  REAL(DP) :: xp1(3), xp2(3), xp3(3) ! dummy variables
  CHARACTER(len=14) :: sub='read_d3dyn_xml'
  LOGICAL :: do_seek
  LOGICAL,ALLOCATABLE :: d3ck(:,:,:)
  !
  IF(present(found)) found=.true.
  !
  do_seek = .false.
  IF(present(seek)) do_seek = seek
  !
  IF (do_seek) THEN
    IF(.not.(present(xq1).and.present(xq2).and.present(xq3).and.present(at))) &
      CALL errore(sub, '"seek" requires xq1, xq2, xq3 and at to be present!', 7)
    filename = d3matrix_filename(xq1, xq2, xq3, at, basename)
  ELSE
    filename = TRIM(basename)
  ENDIF
  !
  CALL iotk_free_unit(u)
  CALL iotk_open_read(u, filename, ierr=ierr, attr=attr)
  IF(present(file_format_version)) &
      CALL iotk_scan_attr(attr, "format_version", file_format_version)
  !
  IF (ierr/=0) THEN
    IF(present(found)) THEN
      found=.false.
      RETURN
    ELSE
      CALL errore(sub, 'D3 matrix file not found "'//TRIM(filename)//'"',1)
    ENDIF
  ENDIF
  !___________________________________________________________
  !
  IF(present(ibrav)) THEN
    IF(.not.present(celldm)) &
      CALL errore(sub, '"ibrav", "celldm" and "at" go together', 1)
      !
      CALL iotk_scan_begin(u, "lattice", attr=attr)
        CALL iotk_scan_attr(attr, "bravais_index", ibrav)
      IF(ibrav==0)THEN
        celldm = 0._dp
        IF(.not. present(at)) &
          CALL errore(sub, 'ibrav=0 and at was not present',8)
        CALL iotk_scan_dat(u, "unit_cell", at, attr=attr)
          CALL iotk_scan_attr(attr, "alat", celldm(1))
      ELSE
        CALL iotk_scan_dat(u, "bravais_parameters", celldm)
!         CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
      ENDIF
    !
    CALL iotk_scan_end(u, "lattice")
  ENDIF
  !___________________________________________________________
  !
  CALL iotk_scan_begin(u, "atoms", attr=attr)
    CALL iotk_scan_attr(attr, "number_of_species", i_ntyp)
    CALL iotk_scan_attr(attr, "number_of_atoms",   i_nat)
    IF(present(ntyp)) ntyp = i_ntyp
    IF(present(nat))  nat  = i_nat
    !
    IF(present(tau)) THEN
      IF(allocated(tau)) DEALLOCATE(tau)
      ALLOCATE(tau(3,i_nat))
      CALL iotk_scan_dat(u, "atomic_positions", tau)
    ENDIF
    IF(present(ityp)) THEN
      IF(allocated(ityp)) DEALLOCATE(ityp)
      ALLOCATE(ityp(i_nat))
      CALL iotk_scan_dat(u, "atomic_types",     ityp)
    ENDIF
    ! the next two have fixed size ntypx
    IF(present(atm)) &
      CALL iotk_scan_dat(u, "species_names",    atm(1:i_ntyp))
    IF(present(amass)) &
      CALL iotk_scan_dat(u, "species_masses",   amass(1:i_ntyp))
    !
  CALL iotk_scan_end(u, "atoms")
  !___________________________________________________________
  !
  CALL iotk_scan_begin(u, "perturbation")
    CALL iotk_scan_dat(u, "q1", xp1)
    CALL iotk_scan_dat(u, "q2", xp2)
    CALL iotk_scan_dat(u, "q3", xp3)
  CALL iotk_scan_end(u, "perturbation")
  IF(do_seek)THEN
    IF( SUM(ABS(xp1-xq1))+SUM(ABS(xp2-xq2))+SUM(ABS(xp3-xq3)) > 1.e-6_dp) THEN
      IF(present(found)) THEN
        found=.false.
        RETURN
      ELSE
        CALL errore(sub, 'File does not contain the correct q1, q2 and/or q3', 3)
      ENDIF
    ENDIF
  ELSE
    IF(present(xq1)) xq1 = xp1
    IF(present(xq2)) xq2 = xp2
    IF(present(xq3)) xq3 = xp3
  ENDIF

  !___________________________________________________________
  !
  IF(present(d3)) THEN
    IF(ALLOCATED(d3)) DEALLOCATE(d3)
    ALLOCATE(d3(3,3,3, i_nat,i_nat,i_nat))
    ALLOCATE(d3ck(i_nat,i_nat,i_nat))
    !
    d3ck = .true.
    CALL iotk_scan_begin(u, 'd3matrix')
    DO k = 1,i_nat
    DO j = 1,i_nat
    DO i = 1,i_nat
      d3ck(i,j,k) = .false.
      CALL iotk_scan_dat(u, "matrix", d3(:,:,:,i,j,k), attr=attr)
        CALL iotk_scan_attr(attr, "atm_i", l)
        CALL iotk_scan_attr(attr, "atm_j", m)
        CALL iotk_scan_attr(attr, "atm_k", n)
      IF(i/=l .or. m/=j .or. n/=k) &
        CALL errore(sub, 'problem with d3 matrix indexes', 4)
    ENDDO
    ENDDO
    ENDDO
    CALL iotk_scan_end(u, 'd3matrix')
    IF(ANY(d3ck)) CALL errore(sub,'could not read all the d3 matrix',40)
  ENDIF
  !
  CALL iotk_close_read(u)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE read_d3dyn_xml
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3matrix_io
!-----------------------------------------------------------------------
