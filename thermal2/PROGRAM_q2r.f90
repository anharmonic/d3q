!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
! on top of the PHonon/PH/q2r.f90 from the Quantum-ESPRESSO package
!  GPLv2 licence and following, <http://www.gnu.org/copyleft/gpl.txt>
!
!----------------------------------------------------------------------------
PROGRAM q2r
  !----------------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE dynamicalq, ONLY : phiq, tau, ityp, zeu
  USE fft_scalar, ONLY : cfft3d
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail, &
                         write_dyn_mat_header, write_ifc
  USE environment, ONLY : environment_start, environment_end
  USE quter_module,       ONLY : quter
  USE input_fc,    ONLY : ph_system_info, forceconst2_grid, write_fc2
  USE rigid_d3, ONLY : rgd_blk_d3
  USE cmdline_param_module
  !
  IMPLICIT NONE
  !
  INTEGER,       PARAMETER :: ntypx = 10
  REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
  INTEGER                  :: nr1, nr2, nr3, nr(3)
  !     dimensions of the FFT grid formed by the q-point grid
  !
  CHARACTER(len=20)  :: crystal
  CHARACTER(len=256) :: fildyn, filin, filj, filf, flfrc, inputf
  CHARACTER(len=6)   :: atm(ntypx)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1, lnogridinfo, xmldyn
  CHARACTER (LEN=10) :: zasr
  INTEGER :: m1, m2, m3, m(3), l1, l2, l3, i, j, j1, j2, na1, na2, ipol, nn
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, ifile, nqs, nq_log
  INTEGER :: na, nt
  !
  INTEGER :: gid, ibrav, ierr, nspin_mag, ios, nfar, input_unit
  !
  INTEGER, ALLOCATABLE ::  nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phid(:,:,:,:,:), matq(:,:,:,:,:)
  INTEGER :: nqtot
  REAL(DP), ALLOCATABLE :: m_loc(:,:), gridq(:,:)
  !
  REAL(DP) :: celldm(6), at(3,3), bg(3,3)
  REAL(DP) :: q(3,48), omega, xq, amass(ntypx), resi
  REAL(DP) :: epsil(3,3), esum
  !
  logical           :: la2F
  LOGICAL, EXTERNAL :: has_xml
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc
  !
  NAMELIST / input / fildyn, flfrc, zasr, la2F, nfar
  !
  CALL mp_startup()
  CALL environment_start('D3_Q2R')
  !
  !If there are no command-line parameters, read from stdin
  IF(cmdline_check_any())THEN
    fildyn  = cmdline_param_char("d", "dyn")
    flfrc   = cmdline_param_char("o", "mat2R")
    zasr    = cmdline_param_char("z", "crystal")
    nfar    = cmdline_param_int ("f", 2)
    inputf  = cmdline_param_char("i", "///")
    !
    IF (cmdline_param_logical('h')) THEN
        WRITE(*,'(2x,a)') "Syntax: d3_q2r.x [-d FILDYN_prefix] [-o FILE_fc] [-z ZASR] [-f NFAR] "
        WRITE(*,'(2x,a)') ""
        WRITE(*,'(2x,a)') "  FILDYN_prefix (default: dyn): Input prefix of dynamical matrix files "
        WRITE(*,'(2x,a)') "  FILE_fc (mat2R): Output force constants"
        WRITE(*,'(2x,a)') "  ZASR (crystal): Sum rule for effective charges"
        WRITE(*,'(2x,a)') "  NFAR (2): 0->produce periodic force constants (like q2r)"
        WRITE(*,'(2x,a)') "            2->produce centered force constants with minimal image distance for interpolation"
        WRITE(*,'(2x,a)') " Note: can also read the same input variables from a &input namelist with switch -i INPUTFILE"
        STOP 1
    ENDIF
    !
!  cmdline = cmdline_residual()
    CALL cmdline_check_exausted()
  ELSE
    fildyn = ' '
    flfrc = 'mat2R'
    zasr = 'no'
    nfar = 2
    
    WRITE(*,*) "No command line option, reading from stdin"
    WRITE(*,*) "WARNING! THIS IS DEPRECATED! Future versions will only accepted command line parameter"
    WRITE(*,*) "use d3_q2r.x -h to get the full help"

    !IF (ionode) CALL input_from_file ( )
    inputf='-'
  ENDIF
     !
     !
  la2F=.false.
     !
     !
  IF (ionode .and. inputf /= '///')  THEN
   IF(TRIM(inputf)=='-')THEN
     input_unit=5
   ELSE
     OPEN(newunit=input_unit, file=inputf, status="OLD", action="READ",iostat=ios)
   ENDIF
   IF(ios==0) READ ( input_unit, input, IOSTAT =ios )
  ENDIF
 
  !CALL mp_bcast(ios, ionode_id, world_comm)
  !CALL errore('q2r','error reading input namelist', abs(ios))

  CALL mp_bcast(fildyn, ionode_id, world_comm)
  CALL mp_bcast(flfrc, ionode_id, world_comm)
  CALL mp_bcast(zasr, ionode_id, world_comm)
  CALL mp_bcast(la2f, ionode_id, world_comm)
  CALL mp_bcast(nfar, ionode_id, world_comm)
     !
     ! check input
     !
  !IF (flfrc == ' ')  CALL errore ('q2r',' bad flfrc',1)
     !
  xmldyn=has_xml(fildyn)

  IF (ionode) THEN
     OPEN (unit=1, file=TRIM(fildyn)//'0', status='old', form='formatted', &
          iostat=ierr)
     lnogridinfo = ( ierr /= 0 )
     IF (lnogridinfo) THEN
        WRITE (stdout,*)
        WRITE (stdout,*) ' file ',TRIM(fildyn)//'0', ' not found'
        WRITE (stdout,*) ' reading grid info from input'
        READ (5, *) nr1, nr2, nr3
        READ (5, *) nfile
     ELSE
        WRITE (stdout,'(/,4x," reading grid info from file ",a)') &
                                                          TRIM(fildyn)//'0'
        READ (1, *) nr1, nr2, nr3
        READ (1, *) nfile
        CLOSE (unit=1, status='keep')
     END IF
  ENDIF
  CALL mp_bcast(nr1, ionode_id, world_comm)
  CALL mp_bcast(nr2, ionode_id, world_comm)
  CALL mp_bcast(nr3, ionode_id, world_comm)
  CALL mp_bcast(nfile, ionode_id, world_comm)
  CALL mp_bcast(lnogridinfo, ionode_id, world_comm)
     !
     IF (nr1 < 1 .OR. nr1 > 1024) CALL errore ('q2r',' nr1 wrong or missing',1)
     IF (nr2 < 1 .OR. nr2 > 1024) CALL errore ('q2r',' nr2 wrong or missing',1)
     IF (nr3 < 1 .OR. nr2 > 1024) CALL errore ('q2r',' nr3 wrong or missing',1)
     IF (nfile < 1 .OR. nfile > 1024) &
        CALL errore ('q2r','too few or too many file',MAX(1,nfile))
     !
     ! copy nrX -> nr(X)
     !
     nr(1) = nr1
     nr(2) = nr2
     nr(3) = nr3
     !
     ! D matrix (analytical part)
     !
     ntyp = ntypx ! avoids spurious out-of-bound errors
     !
     ALLOCATE ( nc(nr1,nr2,nr3) )
     !
     nc = 0
     !
     ! Force constants in reciprocal space read from file
     !
     nqtot = 0
     !
     NFILE_LOOP : &
     DO ifile=1,nfile
        IF (lnogridinfo) THEN
           IF (ionode) READ(5,'(a)') filin
           call mp_bcast(filin, ionode_id, world_comm)
        ELSE
           filin = TRIM(fildyn) // TRIM( int_to_char( ifile ) )
        END IF
        WRITE (stdout,*) ' reading force constants from file ',TRIM(filin)

        IF (xmldyn) THEN
           CALL read_dyn_mat_param(filin,ntyp,nat)
           IF (ifile==1) THEN
              ALLOCATE (m_loc(3,nat))
              ALLOCATE (tau(3,nat))
              ALLOCATE (ityp(nat))
              ALLOCATE (zeu(3,3,nat))
           ENDIF
           IF (ifile==1) THEN
              CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass, tau, ityp, &
                 m_loc, nqs, lrigid, epsil, zeu )
           ELSE
              CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
           ENDIF
           ALLOCATE (phiq(3,3,nat,nat,nqs) )
           DO iq=1,nqs
              CALL read_dyn_mat(nat,iq,q(:,iq),phiq(:,:,:,:,iq))
           ENDDO
           CALL read_dyn_mat_tail(nat)
        ELSE
           IF (ionode) &
           OPEN (unit=1, file=filin,status='old',form='formatted',iostat=ierr)
           CALL mp_bcast(ierr, ionode_id, world_comm)
           IF (ierr /= 0) CALL errore('q2r','file '//TRIM(filin)//' missing!',1)
           CALL read_dyn_from_file (nqs, q, epsil, lrigid,  &
                ntyp, nat, ibrav, celldm, at, atm, amass)
           IF (ionode) CLOSE(unit=1)
        ENDIF
        IF (ifile == 1) THEN
           ! it must be allocated here because nat is read from file
           !ALLOCATE (phid(nr1*nr2*nr3,3,3,nat,nat) )
           ALLOCATE(matq(3,3,nat,nat,nr1*nr2*nr3))
           ALLOCATE(gridq(3,nr1*nr2*nr3))
           !
           lrigid1=lrigid

           CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
           at = at / celldm(1)  !  bring at in units of alat

           CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
           CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
           IF (lrigid .AND. (zasr.NE.'no')) THEN
              CALL set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
           END IF
        END IF
        IF (lrigid.AND..NOT.lrigid1) CALL errore('q2r', &
           & 'file with dyn.mat. at q=0 should be first of the list',ifile)
        !
        WRITE (stdout,*) ' nqs= ',nqs
        DO nq = 1,nqs
!           write(998,'(i6,6f12.6)') nq, q(:,nq)
!           write(998,'(6f12.6)') phiq(:,:,:,:,nq) 
            IF (lrigid) THEN
              !WRITE(stdout,*) "quite rigid"
              ! Remove non-analytic part before doing the Fourier transform
              !CALL rgd_blk_d3 (nr1,nr2,nr3,nat,phiq(:,:,:,:,nq),q(:,nq), &
              !              tau,epsil,zeu,bg,omega,-1.d0)
              CALL rgd_blk_d3 (2,2,2,nat,phiq(:,:,:,:,nq),q(:,nq), &
                    tau,epsil,zeu,bg,omega,celldm(1), .false.,-1.d0) ! 2D added celldm and flag
            END IF
!            write(998,'(6f12.6)') phiq(:,:,:,:,nq) 
            nqtot = nqtot+1
            matq(:,:,:,:,nqtot) = phiq(:,:,:,:,nq)
            gridq(:,nqtot)      = q(:,nq)
            WRITE(stdout,"(5x,3f12.6)") gridq(:,nqtot)
        END DO
        IF (xmldyn) DEALLOCATE(phiq)
     END DO &
     NFILE_LOOP
     !
     ! Check grid dimension
     !
     IF (nqtot == nr1*nr2*nr3) THEN
        WRITE (stdout,'(/5x,a,i4)') ' q-space grid ok, #points = ',nqtot
     ELSE
        CALL errore('init',' missing q-point(s)!',1)
     END IF
     !
     ! dyn.mat. FFT (use serial version)
     !
!      DO j1=1,3
!         DO j2=1,3doing
!            DO na1=1,nat
!               DO na2=1,nat
!                  CALL cfft3d ( phid (:,j1,j2,na1,na2), &
!                       nr1,nr2,nr3, nr1,nr2,nr3, 1 )
!                  phid(:,j1,j2,na1,na2) = &
!                       phid(:,j1,j2,na1,na2) / DBLE(nr1*nr2*nr3)
!               END DO
!            END DO
!         END DO
!      END DO
     !
!     ! auxiliary quantities:
!     REAL(DP),ALLOCATABLE :: sqrtmm1(:) ! 1/sqrt(amass)
!     INTEGER :: nat3, nat32, nat33     
     ALLOCATE(S%tau(3,nat), S%ityp(nat), S%zeu(3,3,nat))
     S%ntyp  = ntyp
     S%amass(1:ntyp) = amass(1:ntyp)
     S%atm(1:ntyp)   = atm(1:ntyp)
     S%nat   = nat
     S%tau   = tau
     
     ! simple sum rule on effective charges
!     do i=1,3
!        do j=1,3
!           esum=0.0d0
!           do na=1,nat
!              esum = esum + zeu(i,j,na)
!           end do
!           do na=1,nat
!              zeu(i,j,na) = zeu(i,j,na) - esum/nat
!           end do
!        end do
!     end do
     !
     S%zeu   = zeu
     S%ityp  = ityp
     S%ibrav = ibrav
     S%symm_type = "unknown"
     S%celldm  = celldm
     S%at      = at
     S%bg      = bg
     S%omega   = omega
     S%epsil   = epsil
     S%lrigid  = lrigid
!     write(998,*) nr1, nr2, nr3, nat
!     write(998,*) tau
!     write(998,*) at
!     write(998,*) bg
!     write(998,'(6f12.6)') matq
!     write(998,'(333f12.6)') gridq
     CALL quter(nr1, nr2, nr3, nat,tau,at,bg, matq, gridq, fc, nfar)
     CALL write_fc2(flfrc, S, fc)
!      STOP 0
!      !
!      ! Real space force constants written to file (analytical part)
!      !
!      IF (xmldyn) THEN
!         IF (lrigid) THEN
!            CALL write_dyn_mat_header( flfrc, ntyp, nat, ibrav, nspin_mag,  &
!                 celldm, at, bg, omega, atm, amass, tau, ityp,   &
!                 m_loc, nqs, epsil, zeu)
!         ELSE
!            CALL write_dyn_mat_header( flfrc, ntyp, nat, ibrav, nspin_mag,  &
!                 celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
!         ENDIF
!         CALL write_ifc(nr1,nr2,nr3,nat,phid)
!      ELSE IF (ionode) THEN
!      OPEN(unit=2,file=flfrc,status='unknown',form='formatted')
!      WRITE(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
!      if (ibrav==0) then
!         write (2,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
!      end if
!      DO nt = 1,ntyp
!         WRITE(2,*) nt," '",atm(nt),"' ",amass(nt)
!      END DO
!      DO na=1,nat
!         WRITE(2,'(2i5,3f18.10)') na,ityp(na),(tau(j,na),j=1,3)
!      END DO
!      WRITE (2,*) lrigid
!      IF (lrigid) THEN
!         WRITE(2,'(3f15.7)') ((epsil(i,j),j=1,3),i=1,3)
!         DO na=1,nat
!            WRITE(2,'(i5)') na
!            WRITE(2,'(3f15.7)') ((zeu(i,j,na),j=1,3),i=1,3)
!         END DO
!      END IF
!      WRITE (2,'(4i4)') nr1, nr2, nr3
!      DO j1=1,3
!         DO j2=1,3
!            DO na1=1,nat
!               DO na2=1,nat
!                  WRITE (2,'(4i4)') j1,j2,na1,na2
!                  nn=0
!                  DO m3=1,nr3
!                     DO m2=1,nr2
!                        DO m1=1,nr1
!                           nn=nn+1
!                           WRITE (2,'(3i4,2x,1pe18.11)')   &
!                                m1,m2,m3, DBLE(phid(nn,j1,j2,na1,na2))
!                        END DO
!                     END DO
!                  END DO
!               END DO
!            END DO
!         END DO
!      END DO
!      CLOSE(2)
!      ENDIF
!      resi = SUM ( ABS (AIMAG ( phid ) ) )
!      IF (resi > eps12) THEN
!         WRITE (stdout,"(/5x,' fft-check warning: sum of imaginary terms = ',es12.6)") resi
!      ELSE
!         WRITE (stdout,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
!      END IF
!      !
!      DEALLOCATE(phid, zeu, nc)
!      IF (.NOT.xmldyn) DEALLOCATE(phiq)
!      !
!      IF(la2F) CALL gammaq2r ( nfile, nat, nr1, nr2, nr3, at )
     !
     DEALLOCATE (tau, ityp)
     !
  !
  CALL environment_end('D3_Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2r
!
! !----------------------------------------------------------------------------
! SUBROUTINE gammaq2r( nqtot, nat, nr1, nr2, nr3, at )
!   !----------------------------------------------------------------------------
!   !
!   USE kinds, ONLY : DP
!   USE fft_scalar, ONLY : cfft3d
!   USE io_global, ONLY : ionode, ionode_id, stdout
!   USE mp_images, ONLY : intra_image_comm
!   USE mp,        ONLY : mp_bcast
!   USE mp_world,  ONLY : world_comm
!   !
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: nqtot, nat, nr1, nr2, nr3
!   REAL(DP), INTENT(IN) :: at(3,3)
!   !
!   INTEGER, ALLOCATABLE :: nc(:,:,:)
!   COMPLEX(DP), ALLOCATABLE :: gaminp(:,:,:,:,:), gamout(:,:,:,:,:)
!   !
!   REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
!   INTEGER  :: nsig = 10, isig, filea2F, nstar, count_q, nq, nq_log, iq, &
!        icar, ipol, m1,m2,m3, m(3), nr(3), j1,j2, na1, na2, nn
!   LOGICAL :: lq
!   REAL(DP) :: deg, ef, dosscf
!   REAL(DP) :: q(3,48), xq, resi
!   character(len=256) :: name
!   CHARACTER(LEN=256) :: elph_dir
!   CHARACTER(LEN=6) :: int_to_char
!   LOGICAL :: exst
!   INTEGER :: ios
! 
!   !
!   ALLOCATE (gaminp(3,3,nat,nat,48), gamout(nr1*nr2*nr3,3,3,nat,nat) )
!   ALLOCATE ( nc (nr1,nr2,nr3) )
!   write (stdout,*)
!   write (stdout,*) '  Preparing gamma for a2F '
!   write (stdout,*)
!   !
!   nr(1) = nr1
!   nr(2) = nr2
!   nr(3) = nr3
!   elph_dir='elph_dir/'
! !  IF (ionode) INQUIRE(FILE=TRIM(elph_dir), EXIST=exst)
! !  CALL mp_bcast(exst, ionode_id, intra_image_comm)
! !  IF (.NOT. exst) CALL errore('gammaq2r','elph_dir directory not exists',1)
!   !
!   DO isig=1, nsig
!      filea2F = 50 + isig
!      nc = 0
!      DO count_q=1,nqtot
!         name= TRIM(elph_dir) // 'a2Fq2r.' // TRIM(int_to_char(filea2F)) &
!                           // '.' // TRIM(int_to_char(count_q))
!         IF (ionode) open(unit=filea2F, file=name, STATUS = 'old', &
!                                  FORM = 'formatted', IOSTAT=ios)
!         CALL mp_bcast(ios, ionode_id, intra_image_comm)
!         IF (ios /= 0) CALL errore('gammaq2r','problem opening file' &
!                                    //TRIM(name), 1)
!      !
!      ! to pass to matdyn, for each isig, we read: degauss, Fermi energy and DOS
!      !
!         !
!         IF (ionode) THEN
!            READ(filea2F,*) deg, ef, dosscf
!            READ(filea2F,*) nstar
!         ENDIF
!         CALL mp_bcast(deg, ionode_id, world_comm)
!         CALL mp_bcast(ef, ionode_id, world_comm)
!         CALL mp_bcast(dosscf, ionode_id, world_comm)
!         CALL mp_bcast(nstar, ionode_id, world_comm)
!         !
!         CALL read_gamma ( nstar, nat, filea2F, q, gaminp )
!         !
!         do nq = 1,nstar
!            lq = .true.
!            do ipol=1,3
!               xq = 0.0d0
!               do icar=1,3
!                  xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
!               end do
!               lq = lq .AND. (ABS(NINT(xq) - xq) < eps)
!               iq = NINT(xq)
!               !
!               m(ipol)= mod(iq,nr(ipol)) + 1
!               if (m(ipol) < 1) m(ipol) = m(ipol) + nr(ipol)
!            end do !ipol
!            IF (.NOT.lq) CALL errore('gammaq2r','q not allowed',1)
!            !
!            if(nc(m(1),m(2),m(3)) == 0) then
!               nc(m(1),m(2),m(3)) = 1
!               CALL TRASL( gamout, gaminp, nq, nr1, nr2, nr3, nat, m(1), m(2), m(3) )
!            else
!               call errore('gammaq2r',' nc already filled: wrong q grid or wrong nr',1)
!            end if
!         enddo ! stars for given q-point
!      ENDDO ! q-points
!      !
!      nq_log = SUM (nc)
!      if (nq_log == nr1*nr2*nr3) then
!         write (stdout,*)
!         write (stdout,'(" Broadening = ",F10.3)') deg
!         write (stdout,'(5x,a,i4)') ' q-space grid ok, #points = ',nq_log
!      else
!         call errore('gammaq2r',' missing q-point(s)!',1)
!      end if
!      do j1=1,3
!         do j2=1,3
!            do na1=1,nat
!               do na2=1,nat
!                  call cfft3d ( gamout(:,j1,j2,na1,na2), &
!                       nr1,nr2,nr3, nr1,nr2,nr3, 1 )
!               end do
!            end do
!         end do
!      end do
!      gamout = gamout / DBLE (nr1*nr2*nr3)
!      !
!      IF (ionode) close(filea2F)
!      !
!      filea2F = 60 + isig
!      name = TRIM(elph_dir) // 'a2Fmatdyn.'// TRIM(int_to_char(filea2F)) 
!      IF (ionode) THEN
!      open(unit=filea2F, file=name, STATUS = 'unknown')
!      !
!      WRITE(filea2F,*) deg, ef, dosscf
!      write(filea2F,'(3i4)') nr1, nr2, nr3
! 
!      do j1=1,3
!         do j2=1,3
!            do na1=1,nat
!               do na2=1,nat
!                  write(filea2F,'(4i4)') j1,j2,na1,na2
!                  nn=0
!                  DO m3=1,nr3
!                     DO m2=1,nr2
!                        DO m1=1,nr1
!                           nn=nn+1
!                           write(filea2F,'(3i4,2x,1pe18.11)')   &
!                                m1,m2,m3, DBLE(gamout(nn,j1,j2,na1,na2))
!                        END DO
!                     END DO
!                  END DO
!               end do  ! na2
!            end do  ! na1
!         end do   !  j2
!      end do   ! j1
!      close(filea2F)
!      ENDIF  ! ionode
! 
!      resi = SUM ( ABS ( AIMAG( gamout ) ) )
! 
!      IF (resi > eps12) THEN
!         WRITE (stdout,"(/5x,' fft-check warning: sum of imaginary terms = ',es12.6)") resi
!      ELSE
!         WRITE (stdout,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
!      END IF
! 
!   ENDDO
!   !
!   DEALLOCATE (gaminp, gamout )
!   !
! END SUBROUTINE gammaq2r
! !
! !-----------------------------------------------------------------------
! subroutine read_gamma (nqs, nat, ifn, xq, gaminp)
!   !-----------------------------------------------------------------------
!   !
!   USE kinds, ONLY : DP
!   USE io_global, ONLY : ionode, ionode_id, stdout
!   USE mp,        ONLY : mp_bcast
!   USE mp_world,  ONLY : world_comm
!   implicit none
!   !
!   ! I/O variables
!   integer, intent(in) :: nqs, nat, ifn
!   real(DP), intent(out) :: xq(3,48)
!   complex(DP), intent(out) :: gaminp(3,3,nat,nat,48)
!   !
!   logical :: lrigid
!   integer :: i, j, na, nb, nt, iq
!   real(DP) :: phir(3),phii(3)
!   CHARACTER(LEN=75) :: line
!   !
!   !
!   Do iq=1,nqs
!      IF (ionode) THEN
!         READ(ifn,*)
!         READ(ifn,*)
!         READ(ifn,*)
!         READ(ifn,'(11X,3F14.9)')  (xq(i,iq),i=1,3)
!      !     write(*,*) 'xq    ',iq,(xq(i,iq),i=1,3)
!         READ(ifn,*)
!      END IF
!      CALL mp_bcast(xq(:,iq), ionode_id, world_comm)
!      do na=1,nat
!         do nb=1,nat
!            IF (ionode) read(ifn,*) i,j
!            CALL mp_bcast(i, ionode_id, world_comm)
!            CALL mp_bcast(j, ionode_id, world_comm)
!            if (i.ne.na) call errore('read_gamma','wrong na read',na)
!            if (j.ne.nb) call errore('read_gamma','wrong nb read',nb)
!            do i=1,3
!               IF (ionode) read (ifn,*) (phir(j),phii(j),j=1,3)
!               CALL mp_bcast(phir, ionode_id, world_comm)
!               CALL mp_bcast(phii, ionode_id, world_comm)
!               do j = 1,3
!                  gaminp(i,j,na,nb,iq) = CMPLX(phir(j),phii(j),kind=DP)
!               end do
!               !           write(*,*) 'gaminp  ',(gaminp(i,j,na,nb,iq),j=1,3)
!            end do
!         end do
!      end do
!      !
!   ENDDO
!   RETURN
!   !
! end subroutine read_gamma
! !
!----------------------------------------------------------------------------
SUBROUTINE trasl( phid, phiq, nq, nr1, nr2, nr3, nat, m1, m2, m3 )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, intent(in) ::  nr1, nr2, nr3, m1, m2, m3, nat, nq
  COMPLEX(DP), intent(in) :: phiq(3,3,nat,nat,48)
  COMPLEX(DP), intent(out) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  INTEGER :: j1,j2,  na1, na2
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              phid(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5d0 * (      phiq(j1,j2,na1,na2,nq) +  &
                          CONJG(phiq(j2,j1,na2,na1,nq)))
           END DO
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE trasl
!----------------------------------------------------------------------
subroutine set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
  !-----------------------------------------------------------------------
  !
  ! Impose ASR - refined version by Nicolas Mounet
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  implicit none
  character(len=10) :: zasr
  integer ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
  integer n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
  !
  real(DP) sum, zeu(3,3,nat)
  real(DP) tau(3,nat), zeu_new(3,3,nat)
  !
  real(DP) zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer zeu_less(6*3),nzeu_less,izeu_less
  ! indices of vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
  ! temporary vectors and parameters

  ! Initialization.
  ! n is the number of sum rules to be considered (if zasr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2'
  ! and (Oz) if axis='3')
  !
  if((zasr.ne.'simple').and.(zasr.ne.'crystal').and.(zasr.ne.'one-dim') &
                       .and.(zasr.ne.'zero-dim')) then
      call errore('q2r','invalid Acoustic Sum Rulei for Z*:' // zasr, 1)
  endif
  if(zasr.eq.'crystal') n=3
  if(zasr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('q2r','too many directions of &
             &   periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'zasr: rotational axis may be wrong'
     endif
     write(stdout,'("zasr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(zasr.eq.'zero-dim') n=6

  ! Acoustic Sum Rule on effective charges
  !
  if(zasr.eq.'simple') then
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
               sum = sum + zeu(i,j,na)
            end do
            do na=1,nat
               zeu(i,j,na) = zeu(i,j,na) - sum/nat
            end do
         end do
      end do
   else
      ! generating the vectors of the orthogonal of the subspace to project
      ! the effective charges matrix on
      !
      zeu_u(:,:,:,:)=0.0d0
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu_new(i,j,na)=zeu(i,j,na)
            enddo
         enddo
      enddo
      !
      p=0
      do i=1,3
         do j=1,3
            ! These are the 3*3 vectors associated with the
            ! translational acoustic sum rules
            p=p+1
            zeu_u(p,i,j,:)=1.0d0
            !
         enddo
      enddo
      !
      if (n.eq.4) then
         do i=1,3
            ! These are the 3 vectors associated with the
            ! single rotational sum rule (1D system)
            p=p+1
            do na=1,nat
               zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
               zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
            enddo
            !
         enddo
      endif
      !
      if (n.eq.6) then
         do i=1,3
            do j=1,3
               ! These are the 3*3 vectors associated with the
               ! three rotational sum rules (0D system - typ. molecule)
               p=p+1
               do na=1,nat
                  zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                  zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
               enddo
               !
            enddo
         enddo
      endif
      !
      ! Gram-Schmidt orthonormalization of the set of vectors created.
      !
      nzeu_less=0
      do k=1,p
         zeu_w(:,:,:)=zeu_u(k,:,:,:)
         zeu_x(:,:,:)=zeu_u(k,:,:,:)
         do q=1,k-1
            r=1
            do izeu_less=1,nzeu_less
               if (zeu_less(izeu_less).eq.q) r=0
            enddo
            if (r.ne.0) then
               call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
               zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
            endif
         enddo
         call sp_zeu(zeu_w,zeu_w,nat,norm2)
         if (norm2.gt.1.0d-16) then
            zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
         else
            nzeu_less=nzeu_less+1
            zeu_less(nzeu_less)=k
         endif
      enddo
      !
      ! Projection of the effective charge "vector" on the orthogonal of the
      ! subspace of the vectors verifying the sum rules
      !
      zeu_w(:,:,:)=0.0d0
      do k=1,p
         r=1
         do izeu_less=1,nzeu_less
            if (zeu_less(izeu_less).eq.k) r=0
         enddo
         if (r.ne.0) then
            zeu_x(:,:,:)=zeu_u(k,:,:,:)
            call sp_zeu(zeu_x,zeu_new,nat,scal)
            zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
         endif
      enddo
      !
      ! Final substraction of the former projection to the initial zeu, to get
      ! the new "projected" zeu
      !
      zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
      call sp_zeu(zeu_w,zeu_w,nat,norm2)
      write(stdout,'("Norm of the difference between old and new effective ", &
           &  "charges: " , F25.20)') DSQRT(norm2)
      !
      ! Check projection
      !
      !write(6,'("Check projection of zeu")')
      !do k=1,p
      !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
      !  call sp_zeu(zeu_x,zeu_new,nat,scal)
      !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
      !enddo
      !
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu(i,j,na)=zeu_new(i,j,na)
            enddo
         enddo
      enddo
   endif
   !
   !
   return
 end subroutine set_zasr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY : DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
