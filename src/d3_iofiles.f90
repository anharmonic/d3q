!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE d3_iofiles
  !
  USE phcom, ONLY : lrdrho
  ! To be made private sooner of later: unit numbers.
  PUBLIC :: iu_dwfc, iu_psi_dH_psi, iu_dpsi_dH_psi, iudpdvp
  PUBLIC :: iu_drho_q, iu_drho_cc_q ! FIXME: these should be private (fix close_files)
  PUBLIC :: fild1rho, fild2rho, fild3rho, fildrho_dir
  ! Subroutines:
  PUBLIC :: openfild3, closefild3, openfile_drho
  PUBLIC :: d3_add_rho_core
  PUBLIC :: addcore_d3
  PUBLIC :: read_drho
  PUBLIC :: fildrho_q, fildrho_q_names
  PUBLIC :: setup_d3_iofiles
  ! ugly:
  PUBLIC :: tmp_dir_d3
  PUBLIC :: d3_basename
  ! LET'S STRESS THIS (use read_drho instead):
  PRIVATE :: davcio_drho_d3

  PRIVATE
  ! prefix of D3 calculation, it depends on the q point for grid calculations:
  CHARACTER(len=512) :: tmp_dir_d3 = ''
  ! units connected to the variation of the ground-state density,
  INTEGER :: iu_drho_q(1:3)  = (/ -1, -1, -1 /)

  ! fildrho names as read from input (i.e. without fildrho_autoname mangling)
  CHARACTER(len=256) :: fild1rho, fild2rho, fild3rho
  ! directory where fildrho files are stored e.g. the same as in ph input,
  ! it is the same as outdir (a.k.a. tmp_dir) by default
  CHARACTER(len=256) :: fildrho_dir = ''
  ! fildrho names after mangling (in a grid calculation they will change at each step)
  TYPE fildrho_q_names
    CHARACTER(len=256) :: name
  END TYPE fildrho_q_names
  TYPE(fildrho_q_names) :: fildrho_q(3)
  ! Same as iu_drho_q, but with added core charge variation
  ! it will be 2000+iu_drho_q, if there is core charge, just
  ! iu_drho_q otherwise:
  INTEGER :: iu_drho_cc_q(1:3) = (/ -1, -1, -1 /)
  LOGICAL :: cc_added_to_drho(1:3) = .false.
  LOGICAL :: drho_changed_q(1:3) = .false.
  !
  INTEGER :: iu_dwfc(-3:3,-3:3)
  !
  ! The wavefunction derivatives projected on conduction bands
  ! Pc_(q_prj) |d_(q_prt) psi_(k+q_wfc)> are saved in the following units:
  !   iu_dwfc(iq_wfc, iq_prt) = 5000+MODULO(10*iq_prt,100)+MODULO(iq_wfc,10)
  ! Where the values of iq_X are -3,-2,-1,0,1,2,3 corresponding to
  ! -q3,-q2,-q1,Gamma,q1,q2,q3 respectively.
  !
  ! Consequently, units will be in the form 50YX where X describes
  ! iq_wfc and Y iq_prt in the following way:
  ! 0 --> Gamma
  ! 1,2,3 --> q1,q2,q3
  ! 9,8,7 --> -q1,-q2,-q3
  ! q_prj is defined as q_prj = q_prt+q_wfc hence it is not necessary to select the unit.
  !
  INTEGER :: iu_psi_dH_psi(-3:3,-3:3) ! = iu_dwfc + 1000
  ! The perturbation of the eigenvalues epsilon_ij = <psi_(k+q_prj)| d_(q_prt) H |psi_(k+q_wfc)> will
  ! be stored following the exact same scheme but +1000, i.e.:
  ! iu_psi_dH_psi(iq_wfc, iq_prt) = 6000+MODULO(10*iq_prt,100)+MODULO(iq_wfc,10)
  !                        ^^^^^^
  ! When 2 or more q vectors are equal, or opposite, or equal to Gamma we do not need all
  ! these files, in these cases the precedence order is q1,q2,q3,Gamma,-q1,-q2,-q3 (determined
  ! in module set_kplus3q) some examples:
  !  1) if q2=-q3 then iu_psi_dH_psi(-3,-2) = iu_dH(2,3) = 5032
  !  2) if q1=-q1=Gamma then iu_psi_dH_psi(-1,1) = iu_dH(0,0) = 5011
  ! This should be totally transparent for i/o as direct access is always used (i.e. no need to rewind).
  !
  INTEGER :: iu_dpsi_dH_psi(-3:3,-3:3) ! = iu_dwfc + 2000
  ! This will contain the units for terms like <dpsi|dH|psi>, more precisely:
  ! iu_dpsi_dH_psi(iq_wfc,iq_prt) 
  !  --> <d^(q_wfc) \psi_k+q_wfc| d^(q_prt) H | \psi_k-q_shift >
  ! Where q_shift is -(q_wfc+q_prt), i.e. the third one
  INTEGER :: iudpdvp(3) = (/ -1, -1, -1 /)
  ! Similar to the previous, only used for metals when one of the q is zero;
  ! it contains the diagonal (i,i) terms of the form:
  !   <d^-q \psi_k,i| d^q H | \psi_k,i >
  ! In particular, for iudpdvp(1) -> (q1,q2,q3) = (0,q,-q)
  ! iudpdvp(2) -> (q,0,-q), iudpdvp(2) -> (q,-q,0), 
  !
  INTEGER,PUBLIC :: &
      iuef,           &! unit with fermi energy shift
      lrpdqvp,        &! length of <psi| dV  |psi>
      lrdpdvp          ! length of <dpsi |  dV |psi> records

CONTAINS
!
! This subroutine open the files with variation of the charge density with
! respect to a perturbation at generic q. If necessary it also shifts drho from
! its q to an equivalent q+G adding a phase, the shifted drho is then saved to 
! another file which is actually used.
!-----------------------------------------------------------------------
SUBROUTINE setup_d3_iofiles(xq1, xq2, xq3)
  !-----------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : at
  USE dfile_autoname,   ONLY : dfile_name
  USE io_files,         ONLY : prefix, check_tempdir
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE mp,               ONLY : mp_bcast, mp_barrier
  USE d3matrix_io2,     ONLY : d3matrix_filename2
  USE d3_control,       ONLY : d3dir
  USE clib_wrappers,    ONLY : f_mkdir
  USE io_files,         ONLY : tmp_dir, wfc_dir
  ! to copy rho:
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : rho
  USE io_rho_xml,       ONLY : write_scf
  USE mp_world,         ONLY : world_comm
  !USE xml_io_base,     ONLY : write_rho

  !
  IMPLICIT NONE
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  INTEGER :: iq !, ios=-1
  CHARACTER(len=16),PARAMETER :: sub = "setup_d3_iofiles"
  !CHARACTER(len=512) :: dirname
  LOGICAL :: exists, parallel_fs
  !
  ! Generate a prefix for the D3 calculation which depends on the q point, to prevent overwrites:
  IF(ionode)THEN
    tmp_dir_d3 = TRIM(d3matrix_filename2(xq1, xq2, xq3, at, TRIM(d3dir)//"/D3"))//"/"
    WRITE(stdout, '(5x,a,/,7x,a)') "Temporary directory set to:", TRIM(tmp_dir_d3)
    tmp_dir = tmp_dir_d3
    wfc_dir = tmp_dir_d3
  ENDIF
  CALL mp_bcast(tmp_dir,    ionode_id, world_comm)
  CALL mp_bcast(wfc_dir,    ionode_id, world_comm)
  CALL mp_bcast(tmp_dir_d3, ionode_id, world_comm)
  !
  CALL check_tempdir(tmp_dir_d3, exists, parallel_fs)
!  CALL write_rho( rho, nspin )
  CALL write_scf(rho, nspin)

  !dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
  !CALL check_tempdir(dirname, exists, parallel_fs)
  !CALL write_rho( dirname, rho%of_r, nspin )

  WRITE(stdout,'(5x,a)') "SCF rho copied to D3 tmp dir"
  !
  ! Set up fildrho names must be done AFTER read_file
  !
  IF(ionode)THEN
    WRITE(stdout, '(5x,a)') "Scanning for fildrho files..."
    fildrho_q(1)%name = dfile_name(xq1, at, fild1rho,TRIM(fildrho_dir)//TRIM(prefix),.false.,0)
    fildrho_q(2)%name = dfile_name(xq2, at, fild2rho,TRIM(fildrho_dir)//TRIM(prefix),.false.,0)
    fildrho_q(3)%name = dfile_name(xq3, at, fild3rho,TRIM(fildrho_dir)//TRIM(prefix),.false.,0)
    DO iq = 1,3
      WRITE(stdout,'(9x,a,i1,3a)') "--> drho file for q", iq, " found as '",TRIM(fildrho_q(iq)%name),"'"
    ENDDO
  ENDIF
  !
  CALL mp_barrier(world_comm)
  !
  DO iq=1,3
    CALL mp_bcast(fildrho_q(iq)%name, ionode_id, world_comm)
  ENDDO
  !
  !
  !-----------------------------------------------------------------------
END SUBROUTINE setup_d3_iofiles
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE openfile_drho()
  !-----------------------------------------------------------------------
  USE constants,  ONLY : eps8
  USE io_global,  ONLY : ionode, stdout
  USE d3_open,    ONLY : diropn_d3
  USE fft_base,   ONLY : dfftp
!   USE mp_global,  ONLY : me_pool, root_pool
  USE kplus3q,    ONLY : kplusq, q_names
  USE lsda_mod,   ONLY : nspin
  USE uspp,       ONLY : nlcc_any
  USE mp,         ONLY : mp_barrier
  USE mp_world,   ONLY : world_comm
  !
  IMPLICIT NONE
  CHARACTER (len=256) :: filint
  INTEGER :: iq
  LOGICAL :: exst, opened(3)
  !
  WRITE(stdout,'(/,5x,"Opening files of charge density derivative.")')
  !
  ! FIXME:
  !tmp_dir = tmp_dir_ph
  !
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin
  !
  OPENED = .false.
  !
  !   is opened only by the first task of each pool
  !
  ! 3 more files are necessary when there is core-charge, this is to prevent
  ! drho from being overwritten (the code would work anyway, but any subsequent
  ! execution would produce bad results)
  !
  iu_drho_q(1:3)  = (/ 301, 302, 303 /)
  !
  IF(nlcc_any) THEN
     iu_drho_cc_q = 2000+iu_drho_q
  ELSE
     iu_drho_cc_q = iu_drho_q
  ENDIF
  !
  ONLY_ROOT_PROCESSOR : &
  IF ( ionode ) THEN !me_pool == root_pool ) THEN
    !
    ! If we have drho for an equivalent q+G point we shift it, this also opens the new file
    ! containg the shifted drho and modifies fildrho_q(iq)%name.
    ! NOTE: the file with the original unshifted drho MUST BE immediately closed to avoid
    !       double-opening when q_X = q_Y+G for X,Y = 1,2,3; X /= Y
    ! NOTE2: fildrho_q(iq)%name must be modified to avoid double-opening of drho+drho_core file
    !
    DO iq = 1,3
      IF(kplusq(iq)%ldrho_is_mine .and. &
        SUM(ABS(kplusq(iq)%xq_drho-kplusq(iq)%xq)) > eps8 ) THEN
        !
!         filint = TRIM(fildrho_q(iq)%name)//'.u'
        CALL drho_change_q(iq, fildrho_q(iq)%name, kplusq(iq)%xq_drho, kplusq(iq)%xq)
        ! Reopen the newly created file (it was initially opened only by ionode)
        CALL diropn_d3(iu_drho_q(iq), TRIM(fildrho_q(iq)%name), lrdrho, exst, tmp_dir_d3)
        OPENED(iq) = .true.
        !
      ENDIF
    ENDDO
    !
    ! Now we open the files that do not need to be shifted, and the files for drho+drho_core
    ! for equal q point we set them to use the same file
    !
    DO iq = 1,3
      !
      filint = TRIM(fildrho_q(iq)%name)
      !
      ONLY_OPEN_IF_HAS_DRHO : &
      IF(kplusq(iq)%ldrho_is_mine) THEN
        IF(.not. opened(iq) ) THEN
            WRITE(stdout,'(9x,3a)') q_names(iq),&
                      "--> opening: ", TRIM(filint)
            CALL diropn_d3(iu_drho_q(iq), filint, lrdrho, exst, fildrho_dir)
            !
            IF (.not.exst) &
              CALL errore ('openfile_drho', 'file ' // TRIM(filint) //' not found', 1)
            !
        ELSE
            WRITE(stdout,'(9x,3a)') q_names(iq),&
                      "--> just opened: ", TRIM(filint)
        ENDIF
        !
        ! Also open the file for d(rho+core), if necessary
        IF (nlcc_any) &
          CALL diropn_d3(iu_drho_cc_q(iq), TRIM(filint)//'+nlcc', lrdrho, exst, tmp_dir_d3)
        !
      ELSE ONLY_OPEN_IF_HAS_DRHO
        WRITE(stdout,'(9x,3a)') q_names(iq),"--> not opening: ", TRIM(filint)
        ! set the unit of this skipped file to use the file which it copies
        iu_drho_q(iq)    = iu_drho_q( ABS(kplusq(iq)%copy_of) )
        iu_drho_cc_q(iq) = iu_drho_cc_q( ABS(kplusq(iq)%copy_of) )
      ENDIF &
      ONLY_OPEN_IF_HAS_DRHO
      !
    ENDDO
    !
    !
  END IF &
  ONLY_ROOT_PROCESSOR
  !
  CALL mp_barrier(world_comm)
  !
  ! FIXME:
  !tmp_dir = tmp_dir_save
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE openfile_drho
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE drho_change_q(iq, filint, xq_old, xq_new)
  !-----------------------------------------------------------------------
  !  Reads the variation of the charge from iudrho_x, adds the variation
  !  of the core_charge, than saves it to iudrho_cc_x
  !
  USE ions_base,  ONLY : nat!, ityp, ntyp => nsp, tau
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE fft_base,   ONLY : dfftp
  USE phcom,      ONLY : lrdrho
  USE uspp_param, ONLY : upf
  USE mp,         ONLY : mp_barrier
  USE io_global,  ONLY : stdout, ionode
  USE drho_add_phase_module, ONLY : drho_add_phase
  USE d3_open,    ONLY : diropn_d3, close_d3

  IMPLICIT NONE

  INTEGER,INTENT(in)   :: iq
  CHARACTER(len=*),INTENT(inout) :: filint
  REAL(DP),INTENT(in)  :: xq_new(3), xq_old(3) ! q point
  !
  INTEGER :: ipert, iu_tmp
  COMPLEX (DP), ALLOCATABLE :: drho_full(:)
  LOGICAL :: exst
  CHARACTER(len=16) :: postfix
  INTEGER,EXTERNAL :: find_free_unit
  !
  CALL start_clock('drho_change_q')
  !
  IF( drho_changed_q(iq)) CALL errore('drho_change_q', 'This rho is already translated', iq)
  drho_changed_q(iq) = .true.
  !
  IF(iq<=0) CALL errore('drho_change_q', 'This subroutine can only work for iq>0',5)
  !
  iu_tmp = find_free_unit()
  ALLOCATE(drho_full(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  !
  IF(.not. ionode) RETURN
  !
  WRITE(stdout, '(7x,a,3f8.4,a,3f10.5,a,i1)') &
                "Translating drho from xq=(",xq_old,") to xq=(",xq_new,") for point q",iq
  !
  ! open original file
  CALL diropn_d3(iu_tmp, TRIM(filint), lrdrho, exst, fildrho_dir)
  !
  ! change the file name and open file for traslated drho
  WRITE(postfix, '(a,i1)') "_newq", iq
  filint = TRIM(filint)//TRIM(postfix)
  CALL diropn_d3(iu_drho_q(iq), TRIM(filint), lrdrho, exst, tmp_dir_d3)
  !
  ! This operation is done entirely by the ionode, so we use a large array to contain
  ! the entire drho and use davcio (instead of davcio_drho_d3) to read/write it
  DO ipert = 1, 3 * nat
     ! read from filint (iu_tmp)
!      CALL davcio_drho_d3(drho, lrdrho, iu_tmp, ipert, -1, pool_only=.true.)
     CALL davcio(drho_full, lrdrho, iu_tmp, ipert, -1)
     ! shift
     CALL drho_add_phase(drho_full, xq_old, xq_new)
     ! write to filint+postfix (iu_drho_q)
!      CALL davcio_drho_d3(drho, lrdrho, iu_drho_q(iq), ipert, +1)
     CALL davcio(drho_full, lrdrho, iu_drho_q(iq), ipert, +1)
     !
  ENDDO
  !
  WRITE(stdout, '(9x,3a)')  "--> new drho stored to '",TRIM(filint),"'"
  CALL close_d3(iu_tmp)
  CALL close_d3(iu_drho_q(iq))
  !
!   CALL mp_barrier()
  !
  DEALLOCATE(drho_full)
  !
  CALL stop_clock('drho_change_q')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE drho_change_q
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE closefild3(cleanup)
  !-----------------------------------------------------------------------
  USE units_lr,  ONLY : iuwfc
  USE d3_open,   ONLY : close_d3
  !
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: cleanup
  !
  INTEGER :: iq1, iq2!, io_level_save
  LOGICAL :: opnd, keep=.true.
  !
  keep = .not. cleanup
  !
!  ! delete wavefunctions and rho from D3 directory
!  IF(cleanup)THEN
!     io_level_save = io_level
!     io_level = -1
!     CALL close_files()
!     io_level = io_level_save
!  ENDIF
  !
  CALL close_d3(iuwfc, keep=keep)
  !
  INQUIRE(UNIT = iuef, OPENED = opnd)
  IF( opnd ) CLOSE(iuef)

  DO iq1 = 1,3
    ! never delete the original drho file from ph:
    CALL close_d3(iu_drho_q(iq1), &
                  keep=(keep.or..not.drho_changed_q(iq1)), &
                  safe=.true.) 
    ! the others are safe to delete
    CALL close_d3(iudpdvp(iq1),      keep=keep, safe=.true.)
    CALL close_d3(iu_drho_cc_q(iq1), keep=keep, safe=.true.)
    !
  ENDDO
  drho_changed_q = .false.
  !
  DO iq1 = -3, 3
    !
    DO iq2 = -3,3
      !
      CALL close_d3(iu_dwfc(iq1, iq2),        keep=keep, safe=.true.)
      CALL close_d3(iu_psi_dH_psi(iq1, iq2),  keep=keep, safe=.true.)
      CALL close_d3(iu_dpsi_dH_psi(iq1, iq2), keep=keep, safe=.true.)
      !
!       INQUIRE (UNIT = iu_dwfc(iq1, iq2), OPENED = opnd)
!       IF( opnd ) CLOSE(iu_dwfc(iq1, iq2))
!       !
!       INQUIRE (UNIT = iu_psi_dH_psi(iq1, iq2), OPENED = opnd)
!       IF( opnd ) CLOSE(iu_psi_dH_psi(iq1, iq2))
!       !
!       INQUIRE (UNIT = iu_dpsi_dH_psi(iq1, iq2), OPENED = opnd)
!       IF( opnd ) CLOSE(iu_dpsi_dH_psi(iq1, iq2))
      !
    ENDDO
    !
  ENDDO
  !
  cc_added_to_drho = .false.
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE openfild3
  !-----------------------------------------------------------------------
  !
  !     This subroutine opens all the files necessary for the
  !     third derivative calculation.
  !
  USE kinds,           ONLY : DP
  USE pwcom,           ONLY : npwx, nbnd, degauss
  USE units_lr,        ONLY : iuwfc, lrwfc, lrdwf
  USE io_files,        ONLY : prefix, seqopn
  USE io_global,       ONLY : ionode
  USE kplus3q,         ONLY : kplusq, q_sum_rule, q_names
  USE d3_open,         ONLY : diropn_d3
  !
  IMPLICIT NONE
  !
  INTEGER :: iq1, iq2, iq_wfc, iq_prt
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint!, tmp_dir_save2
  ! the name of the file
  LOGICAL :: exst, lmetal ! logical variable to check file exists
  !
  !twfcollect=.FALSE.
  lmetal = (degauss /= 0._dp)

  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfild3', 'wrong prefix', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20
  lrwfc = 2 * nbnd * npwx
  CALL diropn_d3 (iuwfc, 'wfc', lrwfc, exst, tmp_dir_d3)
!   IF (.NOT.exst) THEN
!      CALL errore ('openfild3', 'file ' // TRIM(prefix) //'.wfc not found', 1)
!   END IF
  !
  !    The file with deltaV_{bare} * psi
  !
  !
  !    The file with the solution delta psi
  !
  lrdwf = 2 * nbnd * npwx
  lrpdqvp = 2 * nbnd * nbnd
  lrdpdvp = 2 * nbnd * nbnd
  !
  ! We need a certain number of files for the variation of
  ! the wavefunctions: each wavefunction can belong to a certain k+q_rgt grid,
  ! be derived w.r.t a perturbatiom defined by q_prt and projected on the
  ! conduction manifold of a k+q_lft grid. Furthermore, 
  !
  iu_dwfc = -1
  iu_psi_dH_psi   = -1
  iu_dpsi_dH_psi = -1
  !
  DO iq1 = -3,3 ! perturbation index
    iq_prt = kplusq(iq1)%copy_of
    !
    DO iq2 = -3,3 ! k+q_iq2
      !
      iq_wfc = kplusq(iq2)%copy_of
      !
      ! Only consider wavefunctions derivatives that have a meaningful periodicity
!       iq_shift = q_sum_rule(iq2,iq1,exst)
      !                             ^^^^
!       GOOD_COMBINATION : &
!       IF(exst) THEN
        !
        ! If this file has NOT yet been opened for some equivalent q: open it now
        IF( iu_dwfc(iq_wfc, iq_prt) < 0 ) THEN
          !
          iu_dwfc(iq_wfc, iq_prt) = 5000+MODULO(10*iq_prt,100)+MODULO(iq_wfc,10)
          iu_dwfc(iq2, iq1) = iu_dwfc(iq_wfc, iq_prt)
          !
          iu_psi_dH_psi(iq_wfc, iq_prt)  = iu_dwfc(iq_wfc, iq_prt) + 1000
          iu_psi_dH_psi(iq2, iq1)        = iu_dwfc(iq_wfc, iq_prt) + 1000
          !
          iu_dpsi_dH_psi(iq_wfc, iq_prt) = iu_dwfc(iq_wfc, iq_prt) + 2000
          iu_dpsi_dH_psi(iq2, iq1)       = iu_dwfc(iq_wfc, iq_prt) + 2000
          !
          ! open the |psi_(q_prj)><psi_(q_prj)|d_(q_prt) psi_(q_wfc)> file:
!           WRITE(filint, '("d_",a,"__psi_",a,"__k_pool")')  TRIM(ADJUSTL(q_names(iq_prt))), TRIM(ADJUSTL(q_names(iq_wfc)))
          WRITE(filint, '("d",a,"p",a,".")')  TRIM(ADJUSTL(q_names(iq_prt))), TRIM(ADJUSTL(q_names(iq_wfc)))
          CALL diropn_d3 (iu_dwfc(iq2, iq1), TRIM(filint), lrdwf, exst, tmp_dir_d3)
          !
          ! and the file for <psi_(k+q_prj)| d_(q_prt) H |psi_(k+q_wfc)>
!           WRITE(filint, '("psi__dH_",a,"__psi_",a,"__pool")')  &
          WRITE(filint, '("pdv",a,"p",a,".")')  &
            TRIM(ADJUSTL(q_names(iq_prt))), TRIM(ADJUSTL(q_names(iq_wfc)))
          CALL diropn_d3 (iu_psi_dH_psi(iq2, iq1), TRIM(filint), lrpdqvp, exst, tmp_dir_d3)
          !
          ! and, if we have a metal, the file for  <d^(q_wfc) \psi_k+q_wfc| d^(q_prt) H | \psi_k-q_shift >
          IF (lmetal) THEN
!             WRITE(filint, '("dpsi_",a,"__dH_",a,"__psi__pool")')  &
            WRITE(filint, '("dp",a,"dv",a,"p.")')  &
              TRIM(ADJUSTL(q_names(iq_wfc))), TRIM(ADJUSTL(q_names(iq_prt)))
            CALL diropn_d3 (iu_dpsi_dH_psi(iq2, iq1), TRIM(filint), lrdpdvp, exst, tmp_dir_d3)
          ENDIF
          !
        ELSE
          ! otherwise, skip it
          iu_dwfc(iq2, iq1)        = iu_dwfc(iq_wfc, iq_prt)
          iu_psi_dH_psi(iq2, iq1)  = iu_dwfc(iq_wfc, iq_prt) + 1000
          iu_dpsi_dH_psi(iq2, iq1) = iu_dwfc(iq_wfc, iq_prt) + 2000
        ENDIF
!       ENDIF &
!       GOOD_COMBINATION
    ENDDO
  ENDDO
  !
  !
!   iudpdvp = -1
  iudpdvp(1:3) = (/ 8001, 8002, 8003 /)

  METAL : &
  IF (lmetal) THEN
    !
    !    The file with   <dqpsi| dqV |psi> (only in the metallic case)
    !
    IF (kplusq(1)%lgamma .and. &
        kplusq(2)%lgamma .and. &
        kplusq(3)%lgamma ) &
    THEN
      iudpdvp(:) =  iu_dpsi_dH_psi(0,0)! they are equal in this case
    ELSE
      ! Only one of these files will be used, the others will remain empty:
      CALL diropn_d3(iudpdvp(1), 'dpvp1.' , lrdpdvp, exst, tmp_dir_d3)
      CALL diropn_d3(iudpdvp(2), 'dpvp2.' , lrdpdvp, exst, tmp_dir_d3)
      CALL diropn_d3(iudpdvp(3), 'dpvp3.' , lrdpdvp, exst, tmp_dir_d3)
    ENDIF
    !
    ! The file containing the variation of the FermiEnergy ef_sh
    !
    ! opened only by the first task of the first pool
    !
    iuef = 41
    IF(ionode) CALL seqopn(iuef, 'efs', 'unformatted', exst)
    !
  ENDIF &
  METAL

  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE openfild3
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_add_rho_core (scalef)
  !-----------------------------------------------------------------------
  !
  !   Used when non_linear_core_correction are present to change the files
  !   containing the variation of the charge
  !   iflag = +1 :
  !       adds the variation of the core charge to the variation of the
  !       valence charge ( both for xq.eq.0 and xq.ne.0 )
  !
  !   iflag = -1 :
  !       subtracts the variation of the core charge to the variation of
  !       the total charge --used to set drho and d0rho as they were
  !       before the first call of drho_cc--
  !
  USE kinds,            ONLY : DP
  USE uspp,             ONLY : nlcc_any
  USE kplus3q,          ONLY : kplusq, q_names
  USE d3_basis,         ONLY : patq
  USE d3com,            ONLY : d3c
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE ions_base,        ONLY : nat


  IMPLICIT NONE
  REAL(DP),INTENT(in) :: scalef
  INTEGER :: nu, iq, ierr
  COMPLEX (DP), ALLOCATABLE :: drhov(:)

!   print*, "calling the drho_drc with", iud0rho, iud0rhoc, nlcc_any, lgamma
  IF (.NOT.nlcc_any) RETURN
  !
  CALL start_clock('d3_add_rho_core')
  !
  ALLOCATE(drhov(dfftp%nnr), stat=ierr)
  IF(ierr/=0) CALL errore("d3_add_rho_core", "cannot allocate tmp space",1)
  DO iq = 1,3
    IF(kplusq(iq)%ldrho_is_mine) THEN
      WRITE(stdout, "(7x,a,a)") "Adding variation of core charge for ", q_names(iq)
      cc_added_to_drho(iq) = .true.
      !
      DO nu = 1, 3 * nat
        !
        IF(.not.ALLOCATED(patq(iq)%u)) CALL errore("d3_add_rho_core", "error 1", 1)
        IF(.not.ASSOCIATED(d3c(iq)%drc)) CALL errore("d3_add_rho_core", "error 2", 2)
        CALL davcio_drho_d3(drhov, lrdrho, iu_drho_q(iq), nu, -1)
        CALL addcore_d3(kplusq(iq)%xq, patq(iq)%u, nu, d3c(iq)%drc, drhov, +1._dp)
        CALL davcio_drho_d3(drhov, lrdrho, iu_drho_cc_q(iq), nu, +1)
        !
      ENDDO      
    ENDIF
  ENDDO
  DEALLOCATE(drhov)
  !
  CALL stop_clock('d3_add_rho_core')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_add_rho_core
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE read_drho(drho, iq, ipert, with_core, pool_only)
  !-----------------------------------------------------------------------
  USE kinds,    ONLY : DP
  USE phcom,    ONLY : lrdrho
  USE uspp,     ONLY : nlcc_any
  USE kplus3q,  ONLY : kplusq
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE
  INTEGER,INTENT(in)          :: iq, ipert
  COMPLEX(DP),INTENT(inout)   :: drho(dfftp%nnr)
  LOGICAL,INTENT(in)          :: with_core
  LOGICAL,OPTIONAL,INTENT(in) :: pool_only
  !
  INTEGER :: iudrho, jq = -99
  LOGICAL :: do_complex_conjg = .false.
  CHARACTER(len=9),PARAMETER :: sub="read_drho"
!     LOGICAL  :: ldrho    ! this point has it's own drho file from phonon, if ldrho is false:
!     INTEGER  :: drho_from ! get drho from this other q-point
!     LOGICAL  :: ldrho_cc    ! this point has it's own drho file from phonon, if ldrho is false:
!     INTEGER  :: drho_cc_from ! drho can obtained from this other point by complex-conjugate
  !
  CALL start_clock('read_drho')
  !
  IF ( kplusq(iq)%ldrho ) THEN
    jq = kplusq(iq)%drho_from
    do_complex_conjg = .false.
  ELSE IF(kplusq(iq)%ldrho_cc) THEN
    jq = kplusq(iq)%drho_cc_from
    do_complex_conjg = .true.
  ELSE
    CALL errore(sub, 'I have no drho for this q-point', 10+iq)
  ENDIF
  !
  IF(with_core) THEN
    IF(.not. cc_added_to_drho(jq) .and. nlcc_any) THEN
      CALL errore(sub, 'drho with core correction not available at this q, '//&
                       'call "d3_add_rho_core" before!', 20+iq )
    ENDIF
    iudrho = iu_drho_cc_q(jq)
  ELSE
    iudrho = iu_drho_q(jq)
  ENDIF
  !
!   write(*,'(7x,a,i3,a,i3,a,l2)') "\\\ input iq:",iq," used jq:",jq," with cc:", do_complex_conjg
  !
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  CALL davcio_drho_d3(drho, lrdrho, iudrho, ipert, -1, pool_only)
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !
  IF(do_complex_conjg) THEN
!    print*, "   read_drho: reading",iq,"from",iudrho,jq,"with cc"
    drho = CONJG(drho)
!  ELSE
!    print*, "   read_drho: reading",iq,"from",iudrho,jq
  ENDIF
  !
  CALL stop_clock('read_drho')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE read_drho
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE davcio_drho_d3( drho, lrec, iunit, nrec, isw, pool_only )
  !----------------------------------------------------------------------------
  !
  ! reads/writes variation of the charge with respect to a perturbation
  ! on a file.
  ! isw = +1 : gathers data from the nodes and writes on a single file
  ! isw = -1 : reads data from a single file and distributes them
  !
  USE kinds,        ONLY : DP
  USE fft_base,  ONLY : dfftp
!   USE phcom
  USE io_global,    ONLY : ionode_id, ionode
  USE mp_pools,  ONLY : inter_pool_comm, intra_pool_comm, me_pool, root_pool
  !USE mp_global,    ONLY : intra_pool_comm, inter_pool_comm, me_pool, root_pool 
  USE mp,           ONLY : mp_bcast, mp_barrier
  USE scatter_mod,  ONLY : gather_grid, scatter_grid
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in)          :: iunit, lrec, nrec, isw
  COMPLEX(DP),INTENT(inout)   :: drho (dfftp%nnr)
  LOGICAL,OPTIONAL,INTENT(in) :: pool_only
#ifdef __MPI
  !
  ! local variables
  !
  INTEGER :: itmp, proc
  COMPLEX(DP), ALLOCATABLE :: ddrho (:)
  LOGICAL :: intra_pool_only
  !
  intra_pool_only = .false.
  IF(present(pool_only)) intra_pool_only = pool_only
  !
  ALLOCATE (ddrho( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  !
  IF (isw == 1) THEN
     !
     ! First task of the pool gathers and writes in the file
     !
     !CALL cgather_sym (drho, ddrho)
     CALL gather_grid (dfftp, drho, ddrho)
     !
     IF ( ionode ) CALL davcio (ddrho, lrec, iunit, nrec, + 1)
  ELSEIF (isw < 0) THEN
     !
     ! First task of the pool reads ddrho, and broadcasts to all the
     ! processors of the pool
     !
     !print*, "davcioing", intra_pool_only
     IF(.not. intra_pool_only) THEN
        IF ( ionode ) CALL davcio (ddrho, lrec, iunit, nrec, -1)
        CALL mp_bcast( ddrho, ionode_id, inter_pool_comm )
     ELSE
        IF ( me_pool == 0 ) CALL davcio (ddrho, lrec, iunit, nrec, -1)
     ENDIF

     CALL scatter_grid ( dfftp, ddrho(:), drho(:) )

!     CALL mp_bcast( ddrho, root_pool, intra_pool_comm )
!     !
!     ! Distributes ddrho between between the tasks of the pool
!     itmp = 1
!     DO proc = 1, me_pool
!        !itmp = itmp + dfftp%nnp * dfftp%npp (proc)
!        itmp = itmp + dfftp%nnp * dfftp%my_nr3p
!     ENDDO
!     drho (:) = (0.d0, 0.d0)
!     !CALL zcopy (dfftp%nnp * dfftp%npp (me_pool+1), ddrho (itmp), 1, drho, 1)
!     CALL zcopy (dfftp%nnp * dfftp%my_nr3p, ddrho (itmp), 1, drho, 1)
  ENDIF

  DEALLOCATE(ddrho)
#else
  CALL davcio(drho, lrec, iunit, nrec, isw)
#endif
  !
  RETURN
  !
  !----------------------------------------------------------------------------
END SUBROUTINE davcio_drho_d3
!----------------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE addcore_d3(xq, u, mode, drc, drhoc, factor)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the change of the core charge
  !    when the atoms moves along the given mode
  !    drhoc is in real space in input and in output
  !
  USE constants, only : tpi
  USE kinds, only : DP
  use uspp_param, only: upf
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau
  use cell_base, only: tpiba
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  use gvect, only: ngm, mill, eigts1, eigts2, eigts3, g !, nl
  use uspp, only: nlcc_any
  implicit none

  !
  !   The dummy variables
  INTEGER, INTENT(IN) :: mode
  REAL(DP), INTENT(IN) :: xq(3)
  COMPLEX(DP),INTENT(IN) :: u(3*nat, 3*nat), &  ! the transformation modes patterns
                            drc(ngm, ntyp)      ! contain the rhocore (without structure factor)
  
  COMPLEX(DP), intent(INOUT) :: drhoc (dfftp%nnr)
  COMPLEX(DP),ALLOCATABLE :: coreg (:)
  REAL(DP),INTENT(IN) :: factor
  ! input: rho / output: rho+core
  !
  !   Local variables
  !
  INTEGER :: nt, ig, mu, na
  COMPLEX(DP) :: fact, gu, gu0, u1, u2, u3, gtau
  COMPLEX(DP) :: eigqts(nat)
  real(DP)    :: arg
  !
  !
  if (.not.nlcc_any) return
  !
  ALLOCATE(coreg(dfftp%nnr))
  !
  ! compute the derivative of the core charge  along the given mode
  !
  DO na = 1, nat
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     eigqts(na) = CMPLX( COS( arg ), -SIN( arg ) ,kind=DP)
  END DO  
  !
  coreg(:) = (0.d0, 0.d0)
  do na = 1, nat
     nt = ityp (na)
     if (upf(nt)%nlcc) then
        fact = tpiba * (0.d0, -1.d0) * eigqts (na)
        mu = 3 * (na - 1)
        if ( abs(u(mu + 1, mode) ) + &
             abs(u(mu + 2, mode) ) + &
             abs(u(mu + 3, mode) ) > 1.0d-12) then
           u1 = u(mu + 1, mode)
           u2 = u(mu + 2, mode)
           u3 = u(mu + 3, mode)
           gu0 = xq(1)*u1 + xq(2)*u2 + xq(3)*u3
           do ig = 1, ngm
              gtau = eigts1 (mill (1,ig), na) &
                   * eigts2 (mill (2,ig), na) &
                   * eigts3 (mill (3,ig), na)
              gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
              !coreg (nl (ig) ) = coreg(nl(ig)) &
              coreg (dfftp%nl (ig) ) = coreg(dfftp%nl(ig)) &
                               + factor * drc (ig, nt) * gu * fact * gtau
           enddo
        endif
     endif
  enddo
  !
  !   transform to real space
  !
  CALL invfft ('Rho', coreg, dfftp)
  
  drhoc = drhoc+coreg
  DEALLOCATE(coreg)
  !
  return

!-----------------------------------------------------------------------
end subroutine addcore_d3
!-----------------------------------------------------------------------
!
FUNCTION d3_basename(filename) RESULT(basename)
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(in) :: filename
  CHARACTER(len=LEN(filename)) basename
  !
  INTEGER :: i, l
  l = LEN(filename)
  basename = ""
  IF(l<= 1) RETURN
  DO i = l,2,-1
    !print*, i, filename(i:i)
    IF(filename(i:i) == "/") THEN
      IF(filename(i-1:i-1) /= "/" ) THEN
        basename = filename(1:i-1)
        !print*, "basename", TRIM(basename)
        RETURN
      ENDIF
    ENDIF
  ENDDO
END FUNCTION
!
!-----------------------------------------------------------------------
END MODULE d3_iofiles
!-----------------------------------------------------------------------

