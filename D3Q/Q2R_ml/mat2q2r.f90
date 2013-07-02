!GF ricordati di riallocare le variabili
!
MODULE constants
  !----------------------------------------------------------------------------
  !
  !
  ! ... The constants needed everywhere
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... Mathematical constants
   INTEGER, PARAMETER :: DP = selected_real_kind(14,200) 
  ! 
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP 
  REAL(DP), PARAMETER :: tpi    = 2.0_DP * pi
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  REAL(DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP 
  REAL(DP), PARAMETER :: sqrttpi = 2.506628275
  REAL(DP), PARAMETER :: sqrtpm1= 1.0_DP / sqrtpi
  REAL(DP), PARAMETER :: sqrt2  = 1.41421356237309504880_DP
  !
  ! ... Physical constants, SI (NIST CODATA 2006), Web Version 5.1
  !     http://physics.nist.gov/constants
  REAL(DP), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_DP   ! J s
  REAL(DP), PARAMETER :: H_BAR_ERG        = 1.05457172647E-27_DP ! erg s
  REAL(DP), PARAMETER :: H_BAR_SI         = 1.05457172647E-34_DP ! J s
  REAL(DP), PARAMETER :: H_BAR_RY         = 4.837764376E-17_DP     ! Ry s
  REAL(DP), PARAMETER :: K_BOLTZMANN_SI   = 1.3806504E-23_DP    ! J K^-1 
  REAL(DP), PARAMETER :: ELECTRON_SI      = 1.602176487E-19_DP  ! C
  REAL(DP), PARAMETER :: ELECTRONVOLT_SI  = 1.602176487E-19_DP  ! J  
  REAL(DP), PARAMETER :: ELECTRONMASS_SI  = 9.10938215E-31_DP   ! Kg
  REAL(DP), PARAMETER :: HARTREE_SI       = 4.35974394E-18_DP   ! J
  REAL(DP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_DP   ! J
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_DP ! m
  REAL(DP), PARAMETER :: AMU_SI           = 1.660538782E-27_DP  ! Kg
  REAL(DP), PARAMETER :: C_SI             = 2.99792458E+8_DP    ! m sec^-1
  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
  !
  REAL(DP), PARAMETER :: K_BOLTZMANN_AU   = K_BOLTZMANN_SI / HARTREE_SI
  REAL(DP), PARAMETER :: K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
  !
  ! ... Unit conversion factors: energy and masses
  !
  REAL(DP), PARAMETER :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
  REAL(DP), PARAMETER :: RYTOEV           = AUTOEV / 2.0_DP
  REAL(DP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(DP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_DP
  !
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(DP), PARAMETER :: AU_PS            = AU_SEC * 1.0E+12_DP
  !
  ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
  !
  REAL(DP), PARAMETER :: AU_GPA           = HARTREE_SI / BOHR_RADIUS_SI ** 3 &
                                            / 1.0E+9_DP 
  REAL(DP), PARAMETER :: RY_KBAR          = 10.0_DP * AU_GPA / 2.0_DP
  !
  ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm 
  ! ...                                  = 3.3356409519*10^-30 C*m 
  ! ...                                  = 0.208194346 e*A
  ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
  !
  REAL(DP), PARAMETER :: DEBYE_SI         = 3.3356409519_DP * 1.0E-30_DP ! C*m 
  REAL(DP), PARAMETER :: AU_DEBYE         = ELECTRON_SI * BOHR_RADIUS_SI / &
                                            DEBYE_SI
  !
  REAL(DP), PARAMETER :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
  REAL(DP), PARAMETER :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
  REAL(DP), PARAMETER :: ry_to_cm1 =   13.6058d0 * 8065.5d0


  !
  ! .. Unit conversion factors: Energy to wavelength
  !
  REAL(DP), PARAMETER :: EVTONM = 1E+9_DP * H_PLANCK_SI * C_SI / &
                                  &ELECTRONVOLT_SI
  REAL(DP), PARAMETER :: RYTONM = 1E+9_DP * H_PLANCK_SI * C_SI / RYDBERG_SI
  !
  !  Speed of light in atomic units
  !
  REAL(DP), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  !
  ! ... zero up to a given accuracy
  !
  REAL(DP), PARAMETER :: eps4  = 1.0E-4_DP
  REAL(DP), PARAMETER :: eps6  = 1.0E-6_DP
  REAL(DP), PARAMETER :: eps8  = 1.0E-8_DP
  REAL(DP), PARAMETER :: eps12 = 1.0E-12_DP
  REAL(DP), PARAMETER :: eps14 = 1.0E-14_DP
  REAL(DP), PARAMETER :: eps16 = 1.0E-16_DP
  REAL(DP), PARAMETER :: eps24 = 1.0E-24_DP
  REAL(DP), PARAMETER :: eps32 = 1.0E-32_DP
  !
  REAL(DP), PARAMETER :: gsmall = 1.0E-12_DP
  !
  REAL(DP), PARAMETER :: e2 = 2.0_DP      ! the square of the electron charge
  REAL(DP), PARAMETER :: degspin = 2.0_DP ! the number of spins per level
  !
  !!!!!! COMPATIBIILITY
  !
  REAL(DP), PARAMETER :: amconv = AMU_RY
  REAL(DP), PARAMETER :: uakbar = RY_KBAR
  REAL(DP), PARAMETER :: bohr_radius_cm = bohr_radius_si * 100.0_DP
  REAL(DP), PARAMETER :: BOHR_RADIUS_ANGS = bohr_radius_cm * 1.0E8_DP
  REAL(DP), PARAMETER :: ANGSTROM_AU = 1.0_DP/BOHR_RADIUS_ANGS
  REAL(DP), PARAMETER :: DIP_DEBYE = AU_DEBYE
  REAL(DP), PARAMETER :: AU_TERAHERTZ  = AU_PS
  REAL(DP), PARAMETER :: AU_TO_OHMCMM1 = 46000.0_DP ! (ohm cm)^-1
  REAL(DP), PARAMETER :: RY_TO_THZ = 1.0_DP / AU_TERAHERTZ / FPI
  REAL(DP), PARAMETER :: RY_TO_CMM1 = 1.E+10_DP * RY_TO_THZ / C_SI
  !

END MODULE constants

PROGRAM q2r
  USE constants
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: natx=4, nqtx=6000, nqsx=24, nRoutx=10000, nRbigxx=1000000
  INTEGER, PARAMETER::  nrwsx=200

  INTEGER :: nqt, nqs, nr(3), nr_(3)

  real(kind=8), parameter :: eps=1.0d-6
  integer :: nrws
  integer :: nrtot, ii, iR, iiq


  character(len=20) :: crystal
  character(len=80) :: title
  CHARACTER(len=80) :: file_out, file_in0
  CHARACTER(len=80), ALLOCATABLE :: file_in(:)
  CHARACTER(len=80) :: filename
  CHARACTER(len=1) :: chlab1
  CHARACTER(len=4) :: chlab4

  INTEGER :: iu0

  !
  logical :: lq, lrigid, lrigid_save, idiot, asr, lread_tau
  integer :: iq, m1, m2, m3, l1, l2, l3, i, j, j1, j2, na1, na2, ipol, ia, na1_, na2_
  integer :: nat, ntyp, ibrav, icar, nfile, ifile
  integer :: na, nt, ic, jc, lr,il, nRbig
  integer :: nRout, iout, m(3), n_notall
  !

  real(DP) :: celldm(6), at(3,3), bg(3,3)
  real(DP) :: omega, xq, resi,sum, totalweight,wg
  real(DP) :: epsil(3,3), rrbig(3)

   REAL(DP) :: dist(3), wrk(3)
   REAL(DP) ::  atws(3,3), rws(0:3,nrwsx)
    real(DP) :: arg
    REAL(DP), EXTERNAL :: wsweight
  complex(DP) :: aus


  CHARACTER(len=3)  :: atm(natx)
  INTEGER  :: ityp(natx), icorr(natx)
  REAL(DP) :: qqs(3,nqsx), qqt(3,nqtx)
  REAL(DP) :: tau(3,natx), amass(natx), zeu(3,3,natx)
  REAL(DP) :: tau_(3,natx)
  REAL(DP) :: Rbig_(3,nRbigxx), Rbig(3), Rout(3,nRoutx)
  INTEGER, ALLOCATABLE :: nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: matq_in(:,:,:,:,:), matq_s(:,:,:,:), matq(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: matR(:,:,:,:,:)

  namelist / input / nr, file_out, idiot, file_in0, asr, lread_tau, tau_, icorr

  ALLOCATE ( matq_in(3,3,natx,natx,nqsx)   )
  ALLOCATE ( matq_s (3,3,natx,natx)        )
  ALLOCATE ( matq   (3,3,natx,natx,nqtx)   )
  ALLOCATE ( matR   (3,3,natx,natx,nRoutx) )

  nr(:) = 0
  idiot = .FALSE.
  file_out = ' '
  file_in0 = ' '
  asr = .FALSE.
  tau_(:,:) = 0.D0
  lread_tau = .FALSE.
  DO ia = 1, natx
     icorr(ia) = ia
  END DO

  READ (5,input)

  IF (file_in0 .NE. ' ') THEN
     iu0 = 23
     WRITE( chlab1(1:1), '(i1)' ) 0
     filename = TRIM(file_in0) // TRIM(chlab1)
     OPEN(UNIT=iu0, FILE = filename, STATUS = 'unknown')
     READ(iu0,*) nr_(1), nr_(2), nr_(3)
     READ(iu0,*) nfile
     IF (nfile.GT.9999) call errore('reading dyn0','too many files',1)
     ALLOCATE (file_in(nfile))
     DO ifile = 1, nfile
        CALL labels ( ifile, 9999, chlab4, 2 )
        file_in(ifile) = TRIM(file_in0) // TRIM(chlab4)
     ENDDO
     CLOSE(iu0)
     IF ( nr(1)==0 .AND. nr(2)==0 .AND. nr(3)==0 ) nr(:) = nr_(:)
  ELSE
     READ (5,*) nfile
     ALLOCATE (file_in(nfile))
     DO ifile = 1, nfile
        READ(5,'(a)') file_in(ifile)
     ENDDO
  ENDIF

  ALLOCATE ( nc( nr(1), nr(2), nr(3) ) )

  nrtot = nr(1)*nr(2)*nr(3)

  DO l1 = 1, nr(1)
  DO l2 = 1, nr(2)
  DO l3 = 1, nr(3)
     nc(l1,l2,l3) = 0
  END DO
  END DO
  END DO

  nqt = 0
  n_notall = 0
  DO ifile = 1, nfile
     write (6,*) ' reading dyn.mat. from file ',file_in(ifile)
     open(unit=1,file=file_in(ifile),status='old',form='formatted')
     call read_file(nqs,nqsx,qqs,matq_in,natx,epsil,zeu,lrigid,  &
                    ntyp,nat,ibrav,celldm,atm,amass,ityp,tau,at)
     IF (ifile.EQ.1) THEN
        lrigid_save=lrigid
        call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
        if (ibrav.ne.0) at = at / celldm(1)  !  bring at in units of alat 
        call volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
        call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
     END IF
     close(unit=1)

     do iq = 1,nqs
        write(6,'(a,3f12.8)') ' q= ',(qqs(i,iq),i=1,3)

        lq = .true.
        do ipol=1,3
           xq = 0.0
           do icar=1,3
              xq = xq + at(icar,ipol) * qqs(icar,iq) * nr(ipol)
           end do
           lq = lq .and. (abs(nint(xq) - xq) .lt. eps)
           iiq = nint(xq)
           !
           m(ipol)= mod(iiq,nr(ipol)) + 1
           if (m(ipol) .lt. 1) m(ipol) = m(ipol) + nr(ipol)
        end do
        if (lq) then
           if(nc(m(1),m(2),m(3)).eq.0) then
              nc(m(1),m(2),m(3))=1
              CALL trasl ( nat, natx, matq_in(1,1,1,1,iq), matq_s )
              nqt = nqt + 1
              IF (nqt.GT.nqtx) CALL errore('main','nqtx exceeded',1)
              matq(:,:,:,:,nqt) = matq_s(:,:,:,:)
              qqt(:,nqt) = qqs(:,iq)
           else
              write (*,'(3i4)') (m(i),i=1,3)
              call errore('init',' nc already filled: wrong q grid or wrong nr',1)
           end if
        else
           n_notall = n_notall + 1
        endif
     end do
  end do
  if (n_notall.gt.0) then
     write(*,*) 'points excluded from the list:', n_notall
     call errore('init','q not allowed',-1)
  end if


  !  ! Check grid dimension
  !
  if (nqt .eq. nr(1)*nr(2)*nr(3)) then
     write (6,'(/5x,a,i4)') ' q-space grid ok, #points = ',nqt
  else
     call errore('init',' missing q-point(s)!',1)
  end if

  IF (lread_tau) THEN
     IF (tau_(1,1).EQ.0.D0 .AND. tau_(2,1).EQ.0.D0 .AND. tau_(3,1).EQ.0.D0 ) THEN
        DO ia = 2, nat
           IF (tau_(1,ia).EQ.0.D0 .AND. tau_(2,ia).EQ.0.D0 .AND. tau_(3,ia).EQ.0.D0 ) &
              call errore('','',1)
        END DO
     END IF
  ELSE
     tau_(:,:) = tau(:,:)
     DO ia = 1, natx
        icorr(ia) = ia
     END DO
  ENDIF


  atws(:,1) = at(:,1)*DFLOAT(nr(1))
  atws(:,2) = at(:,2)*DFLOAT(nr(2))
  atws(:,3) = at(:,3)*DFLOAT(nr(3))
  ! initialize WS r-vectors
  CALL wsinit(rws,nrwsx,nrws,atws)

  nRbig=0
  do l1=-2*nr(1),2*nr(1)
     do l2=-2*nr(2),2*nr(2)
        do l3=-2*nr(3),2*nr(3)
          nRbig=nRbig+1
           if (nRbig.gt.nRbigxx) call errore('main','nRbigxx exceeded',1)
           Rbig_(:, nRbig) = at(:,1)*l1+at(:,2)*l2 +at(:,3)*l3
         end do
     end do
  end do

   

  open(unit=12,file='matF')
  !
  ! dyn.mat. FFT
  !
  nRout=0
  do na1=1,nat
     do na2=1,nat
        na1_ = icorr(na1)
        na2_ = icorr(na2)
        do j1=1,3
           do j2=1,3
              write(12,'(4i4)') j1, j2, na1_, na2_
              totalweight = 0.d0
              do iR= 1,nRbig
                 aus=cmplx(0.d0,0.d0)
                 dist(:) = Rbig_(:,iR) + tau_(:,na1_) - tau_(:,na2_)
                 Rbig(:) = dist(:) - tau(:,na1) + tau(:,na2)
                 wg = wsweight(dist,rws,nrws)
                 if(wg.ne.0) then 

                    do iiq=1,nqt
                       !   
                       arg=tpi*(qqt(1,iiq)*Rbig(1)+qqt(2,iiq)*Rbig(2)+ qqt(3,iiq)*Rbig(3)) 
                       aus = aus + cmplx(cos(arg),sin(arg),kind=DP)*matq(j1,j2,na1,na2,iiq)
                       !write(*,*) arg,aus,q(1,iiq),iiq    
                    end do
                    wrk(:) = Rbig_(:,iR)
                    call cryst_to_cart(1,wrk,bg,-1)
                    CALL find_in_list (nRoutx,Rout,wrk,iout,nRout)
                    matR(j1,j2,na1_,na2_,iout) = aus * wg  / dfloat(nrtot)
                 end if
                 totalweight=totalweight+wg
               end do
              if(abs(totalweight-nrtot) .gt.eps) call errore('main','wrong totalweight',1)  
           enddo
        enddo
     enddo
  end do

  !
  ! Applies sum rules
  !
  CALL add_asr (nrout,natx,nat,matR,Rout,asr)

  !
  ! Real space force constants written to file (analytical part)
  !
  resi = 0.0d0
  open(unit=2,file=file_out,status='unknown',form='formatted')
  write(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
  if (ibrav.eq.0) then
     write(2,'(3e24.12)') ((at(ic,jc),ic=1,3),jc=1,3)
  endif
  do nt = 1,ntyp
     write(2,*) nt," '",atm(nt),"' ",amass(nt)
  end do
  do na=1,nat
     write(2,'(2i5,3f15.7)') na,ityp(na),(tau_(j,na),j=1,3)
  end do
!
  write (2,'(4i4)') nr(1), nr(2), nr(3) 
  
  do na1=1,nat
     do na2=1,nat 
        do j1=1,3
           do j2=1,3
              !              do ir=1,nrtot
              do ir = 1, nRout
                 resi = resi + dabs(dimag(matR(j1,j2,na1,na2,ir)  ))
              end do
              write (2,'(4i4)') j1,j2,na1,na2
              write(2,'(4i4)') nRout
              do il = 1, nRout
                 write(2,'(3i4,2x,1pe18.11)') (nint(Rout(ic,il)),ic=1,3), real(matR(j1,j2,na1,na2,il))
              enddo
           end do
        end do
     end do
  end do
  close(2)
  write (6,"(/5x,' fft-check: imaginary sum = ',e15.7)") resi
  !


  stop
 
end program q2r
!
!-----------------------------------------------------------------------
subroutine read_file(nqs,nqx,xq,phi,nax,epsil,zeu,lrigid,             &
     &           ntyp,nat,ibrav,celldm,atm,amass,ityp,tau,at)
  !-----------------------------------------------------------------------
  !
  implicit none
  !
  ! I/O variables
  logical :: lrigid
  integer :: nqs, nqx, nax, ntyp, nat, ibrav, ityp(nax)
  real(kind=8) :: epsil(3,3),zeu(3,3,nax), at(3,3)
  real(kind=8) :: xq(3,nqx), celldm(6), amass(nax), tau(3,nax)
  complex(kind=8) :: phi(3,3,nax,nax,nqx)
  character(len=3) atm(nax)
  ! local variables
  integer :: ntyp1,nat1,ibrav1,ityp1
  integer :: i, j, na, nb, nt, ic, jc, ios
  real(kind=8) :: tau1(3), amass1, celldm1(6),q2
  real(kind=8) :: phir(3), phii(3), at1(3,3)
  complex(kind=8) dcmplx
  character(len=75) :: line
  character(len=3)  :: atm1
  logical :: first
  data first/.true./
  save first
  !
  read(1,*) 
  read(1,*) 
  if (first) then
     !
     ! read cell information from file
     !
     read(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     if (nat.gt.nax) call errore('read_f','nax too small',nat)
     if (ntyp.gt.nat) call errore('read_f','ntyp.gt.nat!!',ntyp)
     if ( ibrav.eq.0 ) then
        read(1,*) 
        read(1,*,iostat=ios) ((at(ic,jc),ic=1,3),jc=1,3)
        IF(ios/=0) call errore("read_file", "cannot read basis vectors", 1)
!        at(:,:) = at(:,:) * celldm(1)
     endif
     do nt = 1,ntyp
        read(1,*) i,atm(nt),amass(nt)
        if (i.ne.nt) call errore('read_f','wrong data read',nt)
     end do
     do na=1,nat
        read(1,*) i,ityp(na),(tau(j,na),j=1,3)
        if (i.ne.na) call errore('read_f','wrong data read',na)
     end do
     !
     first=.false.
     lrigid=.false.
     !
  else
     !
     ! check cell information with previous one
     !
     read(1,*) ntyp1,nat1,ibrav1,(celldm1(i),i=1,6)
     if (ntyp1.ne.ntyp) call errore('read_f','wrong ntyp',1)
     if (nat1.ne.nat) call errore('read_f','wrong nat',1)
     if (ibrav1.ne.ibrav) call errore('read_f','wrong ibrav',1)
     do i=1,6
        if(celldm1(i).ne.celldm(i)) call errore('read_f','wrong celldm',i)
     end do

     if ( ibrav.eq.0 ) then
        read(1,*) 
        read(1,*) ((at1(ic,jc),ic=1,3),jc=1,3)
!        at1(:,:) = at1(:,:) * celldm(1)
        do ic = 1, 3
        do jc = 1, 3
           if ( at(ic,jc) .ne. at1(ic,jc) ) then
              write(*,'(3f12.6)') at(:,:)
              write(*,*) ' '
              write(*,'(3f12.6)') at1(:,:)
              call errore('read_f','wrong at',ic)
           endif
        enddo
        enddo
     endif

     do nt = 1,ntyp
        read(1,*) i,atm1,amass1
        if (i.ne.nt) call errore('read_f','wrong data read',nt)
        if (atm1.ne.atm(nt)) call errore('read_f','wrong atm',nt)
        if (amass1.ne.amass(nt)) call errore('read_f','wrong amass',nt)
     end do
     do na=1,nat
        read(1,*) i,ityp1,(tau1(j),j=1,3)
        if (i.ne.na) call errore('read_f','wrong data read',na)
        if (ityp1.ne.ityp(na)) call errore('read_f','wrong ityp',na)
        if (tau1(1).ne.tau(1,na)) call errore('read_f','wrong tau1',na)
        if (tau1(2).ne.tau(2,na)) call errore('read_f','wrong tau2',na)
        if (tau1(3).ne.tau(3,na)) call errore('read_f','wrong tau3',na)
     end do
  end if
  !
  !
  nqs = 0
100 continue
  read(1,*)
  read(1,'(a)') line

  if (line(6:14).ne.'Dynamical') then
     if (nqs.eq.0) call errore('read',' stop with nqs=0 !!',1)
     q2 = xq(1,nqs)**2 + xq(2,nqs)**2 + xq(3,nqs)**2
     if (q2.ne.0.d0) return
     do while (line(6:15).ne.'Dielectric') 
        read(1,'(a)',err=200, end=200) line
     end do
     lrigid=.true.
     read(1,*) ((epsil(i,j),j=1,3),i=1,3)
     read(1,*)
     read(1,*)
     read(1,*)
     write (*,*) 'macroscopic fields =',lrigid
     write (*,'(3f10.5)') ((epsil(i,j),j=1,3),i=1,3)
     do na=1,nat
        read(1,*)
        read(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
        write (*,*) ' na= ', na
        write (*,'(3f10.5)') ((zeu(i,j,na),j=1,3),i=1,3)
     end do
     return
200  write (*,*) ' Dielectric Tensor not found'
     lrigid=.false.     
     return
  end if
  !
  nqs = nqs + 1
  if (nqs.gt.nqx) call errore('read_f','nqx exceeded',1)
  read(1,*) 
  read(1,'(a)') line
  read(line(11:75),*) (xq(i,nqs),i=1,3)
  read(1,*) 
  !
  do na=1,nat
     do nb=1,nat
        read(1,*) i,j
        if (i.ne.na) call errore('read_f','wrong na read',na)
        if (j.ne.nb) call errore('read_f','wrong nb read',nb)
        do i=1,3
           read (1,*) (phir(j),phii(j),j=1,3)
           do j = 1,3
              phi(i,j,na,nb,nqs) = dcmplx(phir(j),phii(j))
           end do
        end do
     end do
  end do
!        write(*,*) 'POPOPOPOP',nqs,phi(1,1,2,2,nqs) 


 !
  go to 100
  !
end subroutine read_file
!
!---------------------------------------------------------------------
subroutine trasl ( nat, natx, matin, matou)
  !---------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nat, natx
  COMPLEX(kind=8) :: matin(3,3,natx,natx), matou(3,3,natx,natx)

  INTEGER :: j1, j2, na1, na2
  DO j1 = 1,3
  DO j2 = 1,3
     DO na1 = 1, nat
     DO na2 = 1, nat
        matou (j1,j2,na1,na2) = &
          0.5d0 * (    matin(j1,j2,na1,na2) +  &
                 conjg(matin(j2,j1,na2,na1)))
     END DO
     END DO
  END DO
  END DO

  RETURN
END subroutine trasl


!-------------------------------------------------------------------
subroutine idiotic_cft3 (fin, nr1, nr2, nr3, nrx1, nrx2, nrx3)
  !-------------------------------------------------------------------
  !
  implicit none
  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3
  complex(kind=8) :: fin (nrx1,nrx2,nrx3)

  integer :: nrxl
  parameter (nrxl = 10)
  real(kind=8) :: twopi
  parameter (twopi= 2.0d0 * 3.14159265358979d0)
  real(kind=8) :: arg
  complex(kind=8) fou(nrxl,nrxl,nrxl), fct, tmp
  integer :: i1, i2, i3, j1, j2, j3

  if ( nr1.gt.nrxl .or. nr2.gt.nrxl .or. nr3.gt.nrxl )  &
        call errore ('idiotic_cft3','wrong parameters',1)

  do i1 = 0, nr1-1
  do i2 = 0, nr2-1
  do i3 = 0, nr3-1
     tmp = cmplx(0.d0,0.d0)
     do j1 = 0, nr1-1
     do j2 = 0, nr2-1
     do j3 = 0, nr3-1
        arg= twopi*(dfloat(i1*j1)/dfloat(nr1) +  &
                    dfloat(i2*j2)/dfloat(nr2) +  &
                    dfloat(i3*j3)/dfloat(nr3) )
        fct = cmplx(cos(arg),sin(arg))
        tmp = tmp + fct*fin(j1+1,j2+1,j3+1)
     enddo
     enddo
     enddo
     fou(i1+1,i2+1,i3+1) = tmp
  enddo
  enddo
  enddo

  do i1 = 1, nr1
  do i2 = 1, nr2
  do i3 = 1, nr3
     fin(i1,i2,i3) = fou(i1,i2,i3)
  enddo
  enddo
  enddo

  return
end subroutine idiotic_cft3
!--------------------------------------------------------------
SUBROUTINE  find_in_list (nrpx,rlist,rinit,iout,nrp)
  !--------------------------------------------------------------
  ! Given the list of point rlist find to which one corresponds rinit
  ! if rinit does not belong to the initial list it is added at the end of te list
  IMPLICIT NONE
  INTEGER :: nrpx, nrp, iout
  REAL(8) :: rlist(3,nrpx), rinit(3)

  REAL(8) :: ddist
  INTEGER :: ip
  REAL(8), PARAMETER :: eps=1.0d-6
 
  DO ip = 1, nrp
     ddist = dsqrt((rinit(1)-rlist(1,ip))**2.0d0 +(rinit(2)-rlist(2,ip))**2.0d0 +(rinit(3)-rlist(3,ip))**2.0d0 ) 
     IF (ddist.le.eps) THEN
        iout = ip
        RETURN
     END IF
  END DO
  nrp  = nrp + 1
  IF (nrp.GT.nrpx) CALL errore('find in list','nRoutx exceeded',1)
  iout = nrp
  rlist(:,iout) = rinit(:)
 
  RETURN
END SUBROUTINE find_in_list

!
!-----------------------------------------------------------------------
subroutine labels ( ilab, imax, chlab, isw )
  implicit none
  integer :: ilab, imax, isw
  character(len=4) :: chlab

  if ( ilab.lt.0 .or. ilab.gt.9999 ) call errore('lalbels',' not implemented',1)
  chlab = '    '  
  SELECT CASE (isw)
  CASE (1)
     if ( imax.lt.10 ) then
        write( chlab(1:1), '(i1)' ) ilab
     else if ( imax.lt.100 ) then
        if ( ilab.lt.10 ) then
           chlab = '0'
           write( chlab(2:2), '(i1)' ) ilab
        else
           write( chlab(1:2), '(i2)' ) ilab
        endif
     else if ( imax.lt.1000 ) then
        if ( ilab.lt.10 ) then
           chlab = '00'
           write( chlab(3:3), '(i1)' ) ilab
        else if ( ilab.lt.100 ) then
           chlab = '0'
           write( chlab(2:3), '(i2)' ) ilab
        else
           write( chlab(1:3), '(i3)' ) ilab
        endif
     else if ( imax.lt.10000 ) then
        if ( ilab.lt.10 ) then
           chlab = '000'
           write( chlab(4:4), '(i1)' ) ilab
        else if ( ilab.lt.100 ) then
           chlab = '00'
           write( chlab(3:4), '(i2)' ) ilab
        else if ( ilab.lt.1000 ) then
           chlab = '0'
           write( chlab(2:4), '(i3)' ) ilab
        else
           write( chlab(1:4), '(i4)' ) ilab
        endif
     else
        call errore('labels','wrong imax',1)
     endif
  CASE (2)
     write_chlab: SELECT CASE (ilab)
     CASE (1:9)
        WRITE( chlab(1:1), '(i1)' ) ilab
     CASE (10:99)
        WRITE( chlab(1:2), '(i2)' ) ilab
     CASE (100:999)
        WRITE( chlab(1:3), '(i3)' ) ilab
     CASE (1000:9999)
        WRITE( chlab(1:4), '(i4)' ) ilab

     END SELECT write_chlab
  CASE DEFAULT
     CALL errore('labels','wrong isw',1)
  END SELECT

  return
end subroutine labels
!--------------------------------------------------------------
SUBROUTINE add_asr (nRbig,natx,nat,mat,Rbig,asr)
  !--------------------------------------------------------------
  LOGICAL :: asr
  INTEGER :: nRbig, natx, nat
  REAL(kind=8) :: Rbig(3,nRbig)
  COMPLEX(kind=8) :: mat(3,3,natx,natx,nRbig)
  INTEGER :: iR0, ir, ic, jc, ia, ja
  REAL(kind=8) :: sum1

  IF (.not.asr) RETURN
  iR0 = -1
  DO ir = 1, nRbig
     IF ( NINT(Rbig(1,ir))==0 .and. NINT(Rbig(2,ir))==0 .and. NINT(Rbig(3,ir))==0 ) iR0 = ir
  END DO
  IF (iR0.eq.-1) call errore('add asr',' something wrong 123',1)

  DO ja = 1, nat
  DO jc = 1, 3
  DO ic = 1, 3
     sum1 = 0.0d0
     DO ia = 1, nat
        DO iR = 1, nRbig
           sum1 = sum1 + real(mat(ic,jc,ia,ja,iR))
        END DO
     END DO
     mat(ic,jc,ja,ja,iR0) = mat(ic,jc,ja,ja,iR0) - sum1
  END DO
  END DO
  END DO

  RETURN
END SUBROUTINE add_asr
!-----------------------------------------------------------------------
!

!!--------------------------------------------------------------
!SUBROUTINE  
!  !--------------------------------------------------------------
!  INTEGER ::
!  REAL(kid=8) ::
! 
!  RETURN
!END SUBROUTINE
!!-----------------------------------------------------------------------
!!
