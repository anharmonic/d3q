!
! Copyright (C) 2001-2011 Quantum-ESPRSSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the ionode_id directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3matrix_io2
!-----------------------------------------------------------------------
  CHARACTER(len=5),PARAMETER :: format_version = "1.1.0"
!
CONTAINS
!
!
!-----------------------------------------------------------------------
FUNCTION d3matrix_filename2(xq1, xq2, xq3, at, basename)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE dfile_autoname, ONLY : dfile_generate_name
  IMPLICIT NONE
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  REAL(DP),INTENT(in) :: at(3,3)
  CHARACTER(len=*),INTENT(in) :: basename
  CHARACTER(len=512) :: d3matrix_filename2

  d3matrix_filename2 =  TRIM(basename) &
            // TRIM(dfile_generate_name(xq1(:), at, "_Q1")) &
            // TRIM(dfile_generate_name(xq2(:), at, "_Q2")) &
            // TRIM(dfile_generate_name(xq3(:), at, "_Q3")) !//".2"
  RETURN
  !-----------------------------------------------------------------------
END FUNCTION d3matrix_filename2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_d3dyn_xml2(basename, xq1,xq2,xq3, d3, ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  USE d3com,      ONLY : code_version => version
  USE io_global,  ONLY : stdout
  USE xmltools
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
  CHARACTER(len=512) :: filename, atms
  CHARACTER(LEN=9)  :: cdate, ctime

  filename = d3matrix_filename2(xq1, xq2, xq3, at, basename)
  WRITE(stdout,"(5x,' -->',a)") TRIM(filename)
  !
  !CALL iotk_free_unit(u)
  u = xml_open_file(filename)
  IF(u==-1) CALL errore("write_d3dyn_xml", "cannot open "//TRIM(filename),1)
  CALL date_and_tim( cdate, ctime )
  !
    CALL add_attr("code_version",   code_version)
    CALL add_attr("format_version", format_version)
    CALL add_attr("date",    cdate)
    CALL add_attr("time",    ctime)
  CALL xmlw_opentag("d3")
!  CALL iotk_open_write(u, filename, root="d3", skip_head=.true.)
    !___________________________________________________________
    !
      CALL add_attr("bravais_index", ibrav)
    CALL xmlw_opentag("lattice")
    IF(ibrav==0)THEN
      !CALL iotk_write_comment(u, "In units of alat which is in bohr units")
        CALL add_attr("alat", celldm(1))
      CALL xmlw_writetag("unit_cell", at)
    ELSE
!      CALL iotk_write_comment(u, "see Doc/INPUT_PW.txt for description of ibrav and lattice parameters")
      CALL xmlw_writetag("bravais_parameters", celldm)
    ENDIF
    !
  CALL xmlw_closetag() !"lattice"
  !___________________________________________________________
  !
      CALL add_attr("number_of_species", ntyp)
      CALL add_attr("number_of_atoms",   nat)
  CALL xmlw_opentag("atoms")
    !
!      CALL iotk_write_comment(u, "positions are in alat units, cartesian coordinates")
    CALL xmlw_writetag("atomic_positions", tau)
    CALL xmlw_writetag("atomic_types",     ityp)
    WRITE(atms,*) atm(1:ntyp)
    CALL xmlw_writetag("species_names",    atms)
!      CALL iotk_write_comment(u, "masses are in Dalton atomic mass units (i.e. mass C^12=12)")
    CALL xmlw_writetag("species_masses",   amass(1:ntyp))
    !
  CALL xmlw_closetag() !"atoms"
  !___________________________________________________________
  !
  CALL xmlw_opentag("perturbation")
!    CALL iotk_write_comment(u, "in cartesian units of 2pi/alat")
    CALL xmlw_writetag("q1", xq1)
    CALL xmlw_writetag("q2", xq2)
    CALL xmlw_writetag("q3", xq3)
  CALL xmlw_closetag()! "perturbation"
  !___________________________________________________________
  !
  CALL xmlw_opentag('d3matrix')
!  CALL iotk_write_comment(u, "Each block contains the 3x3x3 tensor of cartesian displacements for atoms # i,j and k")
!  CALL iotk_write_comment(u, "i is associated with q1, j with q2 and k with q3")
!  CALL iotk_write_comment(u, "Matrix is NOT divided by the square root of masses")
  DO k = 1,nat
  DO j = 1,nat
  DO i = 1,nat
      CALL add_attr("atm_i", i)
      CALL add_attr("atm_j", j)
      CALL add_attr("atm_k", k)
    !CALL xmlw_writetag("matrix", RESHAPE(d3(:,:,:,i,j,k),[27]))
    CALL xmlw_writetag("matrix", d3(:,:,:,i,j,k))
  ENDDO
  ENDDO
  ENDDO
  CALL xmlw_closetag() !'d3matrix'
  CALL xmlw_closetag() !'d3'
  !
  CALL xml_closefile ( )
  !CALL iotk_close_write(u)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE write_d3dyn_xml2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_d3dyn_xml2(basename, xq1,xq2,xq3, d3, ntyp, nat, ibrav, celldm, at, &
                          ityp, tau, atm, amass,found,seek,file_format_version)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE parameters, ONLY : ntypx
  USE xmltools
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
  CHARACTER(len=512) :: filename
  REAL(DP) :: xp1(3), xp2(3), xp3(3) ! dummy variables
  CHARACTER(len=14) :: sub='read_d3dyn_xml2'
  LOGICAL :: do_seek
  LOGICAL,ALLOCATABLE :: d3ck(:,:,:)
  !COMPLEX(DP) :: d3block(27)
  CHARACTER(len=512) :: atms
  !
  IF(present(found)) found=.true.
  !
  do_seek = .false.
  IF(present(seek)) do_seek = seek
  !
  IF (do_seek) THEN
    IF(.not.(present(xq1).and.present(xq2).and.present(xq3).and.present(at))) &
      CALL errore(sub, '"seek" requires xq1, xq2, xq3 and at to be present!', 7)
    filename = d3matrix_filename2(xq1, xq2, xq3, at, basename)
  ELSE
    filename = TRIM(basename)
  ENDIF
  !
  !CALL iotk_free_unit(u)
  !CALL iotk_open_read(u, filename, ierr=ierr)
  u =  xml_open_file ( filename )
      IF ( u == -1 ) CALL errore('read_d3dyn', 'cannot open file '//TRIM(filename),1)
  call xmlr_opentag("d3", ierr)
  
  IF(present(file_format_version) .and. ierr==0) &
      CALL get_attr("format_version", file_format_version)
  !
  IF (ierr/=0) THEN
    CALL xml_closefile()
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
      CALL xmlr_opentag("lattice")
        CALL get_attr("bravais_index", ibrav)
      IF(ibrav==0)THEN
        celldm = 0._dp
        IF(.not. present(at)) &
          CALL errore(sub, 'ibrav=0 and at was not present',8)
        CALL xmlr_readtag("unit_cell", at)
          CALL get_attr("alat", celldm(1))
      ELSE
        CALL xmlr_readtag("bravais_parameters", celldm)
!         CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
      ENDIF
    !
    CALL xmlr_closetag() !"lattice")
  ENDIF
  !___________________________________________________________
  !
  CALL xmlr_opentag("atoms")
    CALL get_attr("number_of_species", i_ntyp)
    CALL get_attr("number_of_atoms",   i_nat)
    IF(present(ntyp)) ntyp = i_ntyp
    IF(present(nat))  nat  = i_nat
    !
    IF(present(tau)) THEN
      IF(allocated(tau)) DEALLOCATE(tau)
      ALLOCATE(tau(3,i_nat))
      CALL xmlr_readtag("atomic_positions", tau)
    ENDIF
    IF(present(ityp)) THEN
      IF(allocated(ityp)) DEALLOCATE(ityp)
      ALLOCATE(ityp(i_nat))
      CALL xmlr_readtag("atomic_types",     ityp)
    ENDIF
    ! the next two have fixed size ntypx
    IF(present(atm)) THEN
      CALL xmlr_readtag("species_names",    atms)
      READ(atms,*) atm(1:ntyp)
    ENDIF
    IF(present(amass)) &
      CALL xmlr_readtag("species_masses",   amass(1:i_ntyp))
    !
  CALL xmlr_closetag() !"atoms")
  !___________________________________________________________
  !
  CALL xmlr_opentag("perturbation")
    CALL xmlr_readtag("q1", xp1)
    CALL xmlr_readtag("q2", xp2)
    CALL xmlr_readtag("q3", xp3)
  CALL xmlr_closetag() !"perturbation")
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
    CALL xmlr_opentag('d3matrix')
    DO k = 1,i_nat
    DO j = 1,i_nat
    DO i = 1,i_nat
      d3ck(i,j,k) = .false.
      !CALL xmlr_readtag("matrix", d3block)
      CALL xmlr_readtag("matrix", d3(:,:,:,i,j,k))
        CALL get_attr("atm_i", l)
        CALL get_attr("atm_j", m)
        CALL get_attr("atm_k", n)
      IF(i/=l .or. m/=j .or. n/=k) &
        CALL errore(sub, 'problem with d3 matrix indexes', 4)
      !d3(:,:,:,i,j,k)=RESHAPE(d3block, [3,3,3])
    ENDDO
    ENDDO
    ENDDO
    CALL xmlr_closetag() !'d3matrix'
    IF(ANY(d3ck)) CALL errore(sub,'could not read all the d3 matrix',40)
  ENDIF
  !
  CALL xmlr_closetag() !'d3'
  !
  call xml_closefile( )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE read_d3dyn_xml2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3matrix_io2
!-----------------------------------------------------------------------
