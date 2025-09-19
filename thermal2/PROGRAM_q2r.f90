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
  CHARACTER(len=256) :: fildyn, filin, filj, filf, flfrc, non_periodic
  CHARACTER(len=6)   :: atm(ntypx)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1, lnogridinfo, xmldyn, nopbc(3)
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
  LOGICAL,EXTERNAL    :: matches
  !
  NAMELIST / input / fildyn, flfrc, zasr, la2F, nfar
  !
  CALL mp_startup()
  CALL environment_start('D3_Q2R')
  !
  !If there are no command-line parameters, read from stdin
  !IF(cmdline_check_any())THEN
    fildyn  = cmdline_param_char("d", "dyn")
    flfrc   = cmdline_param_char("o", "mat2R")
    zasr    = cmdline_param_char("z", "crystal")
    nfar    = cmdline_param_int ("f", 2)
    non_periodic  = cmdline_param_char("n", "x")
    !
    IF (cmdline_param_logical('h')) THEN
        WRITE(*,'(2x,a)') "Syntax: d3_q2r.x [-d FILDYN_prefix] [-o FILE_fc] [-z ZASR] [-f NFAR] [-n NOPBC]"
        WRITE(*,'(2x,a)') ""
        WRITE(*,'(2x,a)') "  FILDYN_prefix (default: dyn): Input prefix of dynamical matrix files "
        WRITE(*,'(2x,a)') "  FILE_fc (mat2R): Output force constants"
        WRITE(*,'(2x,a)') "  ZASR (crystal): Sum rule for effective charges"
        WRITE(*,'(2x,a)') "  NFAR (2): 0->produce periodic force constants (like q2r)"
        WRITE(*,'(2x,a)') "            2->produce centered force constants with minimal image distance for interpolation"
        WRITE(*,'(2x,a)') "  NOPBC: specify along which cell axis the system is non periodic. Can be 1, 2 or 3 or"
        WRITE(*,'(2x,a)') "         any combination."
        WRITE(*,'(2x,a)') "     I.e.   3 -> slab with vacuum along the third axis (typically z). "
        WRITE(*,'(2x,a)') "          123 -> isolated molecule."
        WRITE(*,'(2x,a)') ""
        STOP 1
    ENDIF
    !
!  cmdline = cmdline_residual()
    CALL cmdline_check_exausted()

    WRITE(*,*) "Use d3_q2r.x -h to get the full help"

     !
     !
  la2F=.false.
     !


  CALL mp_bcast(fildyn, ionode_id, world_comm)
  CALL mp_bcast(flfrc, ionode_id, world_comm)
  CALL mp_bcast(zasr, ionode_id, world_comm)
  CALL mp_bcast(la2f, ionode_id, world_comm)
  CALL mp_bcast(nfar, ionode_id, world_comm)
  CALL mp_bcast(non_periodic, ionode_id, world_comm)

  nopbc(:) = .false. ! periodic boundary conditions unless stated otherwise
  IF(matches("1",trim(non_periodic))) nopbc(1) = .true.
  IF(matches("2",trim(non_periodic))) nopbc(2) = .true.
  IF(matches("3",trim(non_periodic))) nopbc(3) = .true.
  

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
     IF(nopbc(1).and.(nr(1)>1) ) CALL errore('q2r','required non-periodic along direction 1, but more than one q-point',1)
     IF(nopbc(2).and.(nr(2)>1) ) CALL errore('q2r','required non-periodic along direction 2, but more than one q-point',2)
     IF(nopbc(3).and.(nr(3)>1) ) CALL errore('q2r','required non-periodic along direction 3, but more than one q-point',3)     !
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
            IF (lrigid) THEN
              CALL rgd_blk_d3 (nopbc,nat,phiq(:,:,:,:,nq),q(:,nq), &
                    tau,epsil,zeu,bg,omega,celldm(1), -1.d0) ! 2D added celldm and flag
            END IF
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
     ALLOCATE(S%tau(3,nat), S%ityp(nat), S%zeu(3,3,nat))
     S%ntyp  = ntyp
     S%amass(1:ntyp) = amass(1:ntyp)
     S%atm(1:ntyp)   = atm(1:ntyp)
     S%nat   = nat
     S%tau   = tau
     
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
     S%nopbc     = nopbc

     
     CALL quter(nr1, nr2, nr3, nat,tau,at,bg, matq, gridq, fc, nfar)
     CALL write_fc2(flfrc, S, fc)

     IF(nr(1)==1 .and. .not. nopbc(1)) WRITE(*,*) "WARNING! please use '-n 1' if system is isolated along direction 1"
     IF(nr(2)==1 .and. .not. nopbc(2)) WRITE(*,*) "WARNING! please use '-n 2' if system is isolated along direction 2"
     IF(nr(3)==1 .and. .not. nopbc(3)) WRITE(*,*) "WARNING! please use '-n 3' if system is isolated along direction 3"     
     !
     DEALLOCATE (tau, ityp)
     !
  !
  CALL environment_end('D3_Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2r
!
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
