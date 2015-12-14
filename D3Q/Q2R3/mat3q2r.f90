MODULE common_var
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: DP = kind(0.0d0)
  REAL(DP), PARAMETER :: pi    = 3.14159265358979323846_DP 
  REAL(DP), PARAMETER :: tpi   = 2.0_DP * pi

  INTEGER  :: ntyp, nat, ibrav, nat3
  INTEGER, ALLOCATABLE :: ityp(:)
  REAL(DP) :: celldm(6), at(3,3), bg(3,3), omega
  REAL(DP), ALLOCATABLE :: amass(:), tau(:,:)
  CHARACTER(len=3), ALLOCATABLE :: atm(:)

  INTEGER  :: nqptsx      ! maximum number of matrices (used only for the allocations)
  INTEGER  :: nqpts
  REAL(DP), ALLOCATABLE :: xq1(:,:), xq2(:,:), xq3(:,:)
  COMPLEX(DP), ALLOCATABLE :: matq(:,:,:,:)

END MODULE common_var
!
!-------------------------------------------------------------------------------------------
PROGRAM q2r
  USE common_var
  IMPLICIT NONE
  !
  integer, parameter ::  nrx1=4, nrx2=4, nrx3=4
  integer :: nqx=0,nRlistx=0,nperx=0
  real(DP), parameter :: eps=1.0d-5
  integer ::  jj,kk,ii
  integer:: kc
  INTEGER, PARAMETER::  nrwsx=200
  logical :: lq1,lq2, lacc
  integer :: l1, l2, l3, i, j, j1, j2, j3, na1, na2, ipol
  integer :: iq, icar, nqs
  integer :: na, nt, ic, jc,nRbigx,nRbig,iq1,iq2,nRsc
  integer :: nRscx,na3,jR,iR,iRs,jRs,iRlist,kcar,jcar
  integer :: m_1(3),m_2(3)
  integer ::  nc1(nrx1,nrx2,nrx3),  nc2(nrx1,nrx2,nrx3)
  real(DP) :: atsc(3,3),d1(3),d2(3),d3(3)
  real(DP) :: resi, norm,maxdiff
  real(DP) :: xq_1,xq_2,dd1,dd2,dd3,per,per0,diff
  real(DP) :: maxdiffq,diffq
  REAL(DP), ALLOCATABLE :: Rbig(:,:),dist(:),wrk(:),q1n(:,:),q2n(:,:),q3n(:,:)
  REAL(DP), ALLOCATABLE :: q1o(:,:),q2o(:,:),q3o(:,:),Rsc(:,:),wrk1(:),wrk2(:),q1n2(:,:),q2n2(:,:),q3n2(:,:)
   real(DP) :: arg
   real(dp),allocatable :: Rlist(:,:,:)
   REAL(DP), EXTERNAL :: wsweight
 integer,allocatable :: itable(:)
 integer nper
 integer,allocatable :: ijReq(:,:), irl(:)
 integer nRlist,iper,nqacc,qacc
  complex(DP), ALLOCATABLE :: phiqn(:,:,:,:),phiq2(:,:,:,:),phiqo(:,:,:,:),phiqn2(:,:,:,:)
  complex(DP) :: aus(3,3,3),phase,aus1(3,3,3)
  complex(DP),allocatable :: phiq1(:,:,:,:,:,:,:)
  complex(DP),allocatable :: mat(:,:,:,:,:,:,:)!(3,3,3,nat,nat,nat,nRlistx)


  INTEGER :: ifile

  INTEGER :: nr1, nr2, nr3, nr(3), nrtot
  INTEGER, PARAMETER :: nfilex = 100
  INTEGER :: nfile
  CHARACTER(len=256) :: filein(nfilex), fileout
  
  NAMELIST / input / nr1, nr2, nr3, nat, fileout, nqptsx


  !
  ! Set default values for the namelist
  !
  nr1=0
  nr2=0
  nr3=0
  nat=2
  !nqptsx = 400000   ! maximum number of q-points allowed (used for allocations)
  !nRbigx= 400000
  nRscx=  5**3            ! <-- check :lp:
  nperx= nRscx*6
  fileout = 'mat3R'
  !
  ! Reads input
  !
  READ (5,input)
  nqptsx=100*(nr1*nr2*nr3)**2
  WRITE(*,'(5x,"max number of q-points (est.)",i10)') nqptsx
  nqx = nqptsx !(nr1*nr2*nr3)**2            ! <-- check :lp:
  nRbigx = (nr1*nr2*nr3)**2          ! <-- check :lp:
  nRlistx = 32*(nr1*nr2*nr3+2)**2        ! <-- check :lp:
  
  READ(5,*) nfile
  IF ( nfile.GT.nfilex .OR. nfile.LE.0) CALL errore('reading input','wrong nfile',1)
  DO ifile = 1, nfile
    READ(5,'(a)')  filein(ifile)
  ENDDO

  IF (nr1 < 1) CALL errore ('q2r',' nr1 wrong or missing',1)
  IF (nr2 < 1) CALL errore ('q2r',' nr2 wrong or missing',2)
  IF (nr3 < 1) CALL errore ('q2r',' nr3 wrong or missing',3)
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3

  IF(nr1>nrx1) CALL errore('q2rd3', 'nr1x too small', 1)
  IF(nr2>nrx2) CALL errore('q2rd3', 'nr2x too small', 2)
  IF(nr3>nrx3) CALL errore('q2rd3', 'nr3x too small', 3)

  nrtot = nr(1)*nr(2)*nr(3)
  
  ! :lp: previously static:
  ALLOCATE(mat(3,3,3,nat,nat,nat,nRlistx))
  ALLOCATE(itable(nRlistx))
  ALLOCATE(Rlist(3,2,nRlistx))
  ALLOCATE( ijReq(2,nperx), irl(nperx))
  ALLOCATE(phiq1(3,3,3,nat,nat,nat,nqx))

  !
  ! Reads matrices from file and constructs
  ! the matrices equivalent for permutation of the indexes
  !
  CALL read_file ( nfile, filein )



 ALLOCATE(q1n(3,nqx),q2n(3,nqx),q3n(3,nqx))
 ALLOCATE(q1n2(3,nqx),q2n2(3,nqx),q3n2(3,nqx))
 ALLOCATE(q1o(3,nqx),q2o(3,nqx),q3o(3,nqx))
! ALLOCATE(q1c(3,nrtot),q2c(3,nrtot))
 ALLOCATE (dist(3),wrk(3),wrk1(3),wrk2(3))
 ALLOCATE (phiqo(3*nat,3*nat,3*nat,nqx))


    phiqo(1:3*nat,1:3*nat,1:3*nat,1:nqx) = matq(1:3*nat,1:3*nat,1:3*nat,1:nqx)
    q1o(:,1:nqx) = xq1(:,1:nqx)
    q2o(:,1:nqx) = xq2(:,1:nqx)
    q3o(:,1:nqx) = xq3(:,1:nqx)


 !if (nat.ne.nat) call errore('reading ',' nat .ne. nat',1)
 
 ALLOCATE(phiqn2(3*nat,3*nat,3*nat,nqx))
 !
 nqs=nqpts 
 write(*,*) ' calling qaccepted'
 call qaccepted(nqacc,nqx,nat,nqs,at,bg,q1o,q2o,q3o,phiqo,q1n2,q2n2,q3n2,phiqn2)
 write(*,*)'nqacc', nqacc
 !
 ! :lp: cleanup
 deallocate(q1o,q2o,q3o)
 deallocate(phiqo)
 ALLOCATE(phiqn(3*nat,3*nat,3*nat,nqx))

  !
  do l1=1,nr(1)
    do l2=1,nr(2)
      do l3=1,nr(3)
!        nc(l1,l2,l3)=0
        nc1(l1,l2,l3)=0
        nc2(l1,l2,l3)=0
      end do
    end do
  end do


! ALLOCATE ( Rbig(3,nRbigx),rout(3,nRbigx),Rsc(3,nRscx))
 ALLOCATE ( Rbig(3,nRbigx),Rsc(3,nRscx))

  qacc = 0
  DO j = 1, nqacc

    lq1 = .TRUE.
    lq2 = .true.
    DO ipol = 1, 3
      xq_1 = 0.0
      xq_2 = 0.0
      DO icar=1,3
        xq_1 = xq_1 + at(icar,ipol) * q1n2(icar,j) * nr(ipol)
        xq_2 = xq_2 + at(icar,ipol) * q2n2(icar,j) * nr(ipol)
      END DO
      lq1 = lq1 .and. (abs(nint(xq_1) - xq_1) .lt. eps)
      lq2 = lq2 .and. (abs(nint(xq_2) - xq_2) .lt. eps)
      iq1 = nint(xq_1)
      iq2 = nint(xq_2)
      m_1(ipol) = mod(iq1,nr(ipol)) + 1
      m_2(ipol) = mod(iq2,nr(ipol)) + 1
      IF (m_1(ipol) .LT. 1) m_1(ipol) = m_1(ipol) + nr(ipol)
      IF (m_2(ipol) .LT. 1) m_2(ipol) = m_2(ipol) + nr(ipol)
    END DO
    lacc = lq1 .and. lq2
    IF (.not.lacc )then !  call errore('init','q1 not allowed',1)
!      write(*,*) 'q1/q2 of 4x4x4'
    ELSE
      qacc = qacc + 1
      q1n(:,qacc) = q1n2(:,j)
      q2n(:,qacc) = q2n2(:,j)
      q3n(:,qacc) = q3n2(:,j)
      phiqn(:,:,:,qacc) = phiqn2(:,:,:,j)
    END IF
    if (lacc) then
      if(nc1(m_1(1),m_1(2),m_1(3)).lt.nrtot) then
        nc1(m_1(1),m_1(2),m_1(3))=nc1(m_1(1),m_1(2),m_1(3))+1
      else
        write(*,*) q2n2(:,j)
        write (*,*) (m_2(i),i=1,3),'nc1 already filled'
        !call errore('init',' nc1 already filled: wrong q grid or wrong nr',1)
      end if
     if (nc2(m_2(1),m_2(2),m_2(3)).lt.nrtot) then
       nc2(m_2(1),m_2(2),m_2(3))=nc2(m_2(1),m_2(2),m_2(3))+1
     else
       write(*,*) q2n2(:,j)
       write (*,*) (m_2(i),i=1,3),'nc2 already filled'
       !call errore('init',' nc2 already filled: wrong q grid or wrong nr',1)
     end if
    end if
  END DO
 ! :lp: cleanup
  deallocate(phiqn2)
  
  nqacc=qacc
  write(*,*)'---------------------------------------'
  write(*,*)'nqacc post grid', nqacc

  DO l1 = 1, nr(1)
  DO l2 = 1, nr(2)
  DO l3 = 1, nr(3)
    IF ( nc1(l1,l2,l3) .NE. nrtot ) THEN
     print*, nc1(l1,l2,l3), nrtot, l1,l2,l3, nr
     CALL errore('init', 'nc1 different from nrtot', 1)
    ENDIF
    IF ( nc2(l1,l2,l3) .NE. nrtot ) THEN
      print*, nc2(l1,l2,l3), nrtot, l1,l2,l3, nr
      CALL errore('init', 'nc2 different from nrtot', 1)
    ENDIF
  END DO
  END DO
  END DO

write(*,*) 'lab01'


    atsc(:,1) = at(:,1)*DFLOAT(nr(1))
    atsc(:,2) = at(:,2)*DFLOAT(nr(2))
    atsc(:,3) = at(:,3)*DFLOAT(nr(3))

  nRsc=0
  do l3 = -2, 2
    do l2 = -2, 2
      do l1 = -2, 2
        nRsc = nRsc+1
        Rsc(:, nRsc) = atsc(:,1)*l1 + atsc(:,2)*l2 + atsc(:,3)*l3
        if (nRsc.gt.nRscx) call errore('main','nRsc exceeded',1)
       end do
    end do
  end do
 
  nRbig=0
  do l3 = 0, nr(3)-1
    do l2 = 0, nr(2)-1
      do l1 = 0, nr(1)-1
        nRbig = nRbig + 1
        Rbig(:, nRbig) = at(:,1)*l1 + at(:,2)*l2 + at(:,3)*l3
        if (nRbig.gt.nRbigx) call errore('main','nRbig exceeded',1)
       end do
    end do
  end do
write(*,*) 'lab02'
 
open(unit=12,file=trim(fileout))
  nRlist = 0
  norm = 1.d0 / dfloat(nRbig*nRbig)
  mat = cmplx(0.d0,0.d0)
  do na3=1,nat
    do na2=1,nat
      do na1=1,nat
        
        itable = 0
        do jR=1,nRbig
          do iR=1,nRbig
            
            aus(:,:,:)=cmplx(0.d0,0.d0)
            do iq = 1, nqacc
              !write(*,'(3(2f10.4,3x))') q2n(:,iq), q3n(:,iq)
              arg=tpi*(q2n(1,iq)*Rbig(1,iR)+q2n(2,iq)*Rbig(2,iR)+ q2n(3,iq)*Rbig(3,iR) +&
                  q3n(1,iq)*Rbig(1,jR)+q3n(2,iq)*Rbig(2,jR)+ q3n(3,iq)*Rbig(3,jR) )
              phase = cmplx(cos(arg),sin(arg),kind=DP) ! cambiato il segno
              do ic = 1, 3
                ii = ic + 3*(na1-1)
                do jc = 1, 3
                  jj = jc + 3*(na2-1)
                  do kc = 1, 3
                    kk = kc + 3*(na3-1)
                    
                    aus(ic,jc,kc) = aus(ic,jc,kc) + phase * norm *phiqn(ii,jj,kk,iq) 
                  enddo
                enddo
              enddo
            end do
            !if(na2/=na3)then
              write(998,*)
              write(998,'(99i6)') na1,na2,na3, iR,jR
              write(998,'(3(2f10.4,3x))') aus
              !stop 10
            !endif
            !stop 10
            
            per0 = 1.d30
            nper = 0
            do jRs=1,nRsc
              do iRs=1,nRsc
                
                d1(:) = tau(:,na1)-tau(:,na2)-Rbig(:,iR)-Rsc(:,iRs)
                d2(:) = tau(:,na2)+Rbig(:,iR)+Rsc(:,iRs) -tau(:,na3)-Rbig(:,jR)-Rsc(:,jRs)
                d3(:) = tau(:,na3)+Rbig(:,jR)+Rsc(:,jRs) - tau(:,na1)
                
                dd1=dsqrt(d1(1)**2+d1(2)**2+d1(3)**2)
                dd2=dsqrt(d2(1)**2+d2(2)**2+d2(3)**2)
                dd3=dsqrt(d3(1)**2+d3(2)**2+d3(3)**2)
                
                per=dd1+dd2+dd3
                if (per.lt.per0-eps) then
                  nper = 1
                  ijReq(1,nper) = iRs
                  ijReq(2,nper) = jRs
                  per0 = per
                elseif(abs(per-per0).le.eps) then
                  nper = nper + 1
                  if (nper.gt.nperx) call errore('main','nperx exceeded',1)
                  ijReq(1,nper) = iRs
                  ijReq(2,nper) = jRs
                endif
               ! write(10,*) per,per0
                
              end do
            end do
            if (nper.le.0) call errore('main','unexpected',1)
!write(*,*) 'nat:', na1, na2, na3, ' nper:', nper
!write(*,'(2(i5,3f12.6,2x))') iR, Rbig(:,iR), jR, Rbig(:,jR)
!do iper = 1, nper
!write(*,'(2(i5,3f12.6,2x))') ijReq(1,iper), Rsc(:,ijReq(1,iper)), ijReq(2,iper), Rsc(:,ijReq(2,iper))
!enddo
!write(*,*) ' '
            call add_to_Rlist (nRbig,nRsc,Rbig,Rsc, nRlist, nRlistx,Rlist, iR, jR, nper, ijReq, irl)
!write(*,*) 'nRlist:', nRlist
!do iper = 1, nRlist
!write(*,'(3f12.6,2x,3f12.6)') Rlist(:,1,iper), Rlist(:,2,iper)
!enddo
!write(*,*) ' '
      

              ! write(10,*) 'nper matrici'
            do iper = 1, nper
              if (irl(iper).gt.nRlist .or. irl(iper).le.0 ) call errore('main','something wrong',1)
              if (itable(irl(iper)).ne.0) call errore('main','something wrong',1)
              do j3 =1,3
                do j2 = 1,3 
                  do j1 =1,3
                    mat(j1,j2,j3,na1,na2,na3,irl(iper)) = aus(j1,j2,j3)/ dfloat(nper)
                    itable(irl(iper)) = 1
                    

                  end do
                end do
              end do
             
            enddo

          end do
        end do

      end do
    end do
  end do

   ! set to exactly zero the very small terms
   where(abs(mat)<1.d-30) mat= 0.d0
  
  write(12,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
  if (ibrav.eq.0) then
    write(12,'(3e24.12)') ((at(ic,jc),ic=1,3),jc=1,3)
  endif
  do nt = 1,ntyp
    write(12,*) nt," '",atm(nt),"' ",amass(nt)
  end do
  do na=1,nat
    write(12,'(2i5,3f15.7)') na,ityp(na),(tau(j,na),j=1,3)
  end do
!
  write (12,'(4i4)') nr(1), nr(2), nr(3) 

  resi=0.0d0
  maxdiff=0.0d0

  do na1 =1,nat
    do na2 = 1,nat 
      do na3 =1,nat
        do j1 =1,3
          do j2 = 1,3 
            do j3 =1,3
              write(12,'(6i4)') j1,j2,j3,na1,na2,na3
              write(12,'(i8)') nRlist
              do iR=1,nRlist

                resi = resi + dabs(aimag(mat(j1,j2,j3,na1,na2,na3,iR)))
                diff=ABS(AIMAG(mat(j1,j2,j3,na1,na2,na3,iR)))
                IF(diff.GT.maxdiff) THEN
                  maxdiff=diff
                end if

                wrk1(:) = Rlist(:,1,iR)
                wrk2(:) = Rlist(:,2,iR)
                call cryst_to_cart(1,wrk1,bg,-1)
                call cryst_to_cart(1,wrk2,bg,-1)

!                 d1(:) = tau(:,na1) - tau(:,na2)-Rlist(:,1,iR)
!                 d2(:) = tau(:,na2)+Rlist(:,1,iR) - tau(:,na3)-Rlist(:,2,iR)
!                 d3(:) = tau(:,na3)+Rlist(:,2,iR) - tau(:,na1)
!                 dd1=dsqrt(d1(1)**2+d1(2)**2+d1(3)**2)
!                 dd2=dsqrt(d2(1)**2+d2(2)**2+d2(3)**2)
!                 dd3=dsqrt(d3(1)**2+d3(2)**2+d3(3)**2)
!                 per=dd1+dd2+dd3

                write(12,'(6i4,2x,1pe25.15,1pe12.3)')  &
                  (nint(wrk1(ic)),ic=1,3), &
                  (nint(wrk2(ic)),ic=1,3), &
                  real(mat(j1,j2,j3,na1,na2,na3,iR)),&
                  imag(mat(j1,j2,j3,na1,na2,na3,iR))
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  
  close(12)
  
  write (6,"(/5x,' check: imaginary sum = ',e15.7)") resi   
  write (6,"(/5x,' check: largest imaginary part = ',e15.7)") maxdiff




  IF (.true.) THEN

  write(*,*) "lab03"


! test to check the qmatrix.



  open(unit=4,file='dinQ')
!   ii=0
  do na1=1,nat
    do na2=1,nat
      do na3=1,nat 
        
        
        !
        do iq= 1,nqacc
          aus1=cmplx(0.d0,0.d0)
          do iRlist=1,nRlist
            
            
            
            
            arg=tpi*(q2n(1,iq)* Rlist(1,1,iRlist)+q2n(2,iq)* Rlist(2,1,iRlist)  + q2n(3,iq)* Rlist(3,1,iRlist) +&
                q3n(1,iq)* Rlist(1,2,iRlist)  +q3n(2,iq)*  Rlist(2,2,iRlist)  + q3n(3,iq)*  Rlist(3,2,iRlist)  )
            aus1 = aus1 + cmplx(cos(arg),-sin(arg),kind=DP)*REAL(mat(:,:,:,na1,na2,na3,iRlist),kind=DP)
            
            
          end do
          
          !             ii=ii+1
          phiq1(:,:,:,na1,na2,na3,iq) = aus1
          
 

!          write(4,'(1pe18.11)') real( phiq1(:,:,:,na1,na2,na3,iq))
!          write(4,*)  phiq1(:,:,:,na1,na2,na3,iq)

          !                      write(*,*)phiq(j1,j2,na1,na2,iiq)            
          !            phiq2(j1,j2,na1,na2)= phiq(j1,j2,na1,na2,iiq)
          
          !             CALL dyndiag(nat,ntyp,amass,ityp,phiq2,w2(1,iiq),z)
        end do
      enddo
    enddo
  enddo
if(.true.)then
do iq=1,nqacc
  do i = 1, 3 * nat
!    if (wrmode (i) ) then
      !
      write (4, '(/,12x,"modo:",i5,/)') i
      na3 = (i - 1) / 3 + 1
      kcar = i - 3 * (na3 - 1)
      do na1 = 1, nat
        do na2 = 1, nat
          write (4, '(2i3)') na1, na2
          do icar = 1, 3
            write (4, '(6e24.12)') (phiq1 (kcar, icar, jcar, na3, na1, na2,iq) &
               , jcar = 1, 3)
          enddo
        enddo
      enddo
      !
!    endif
  enddo
end do
end if




! test back ft

ALLOCATE (phiq2(3*nat,3*nat,3*nat,nqx))
maxdiffq=0.0d0
do iq=1,nqacc
  do na3=1,nat
    do kcar=1,3
      kk = kcar + 3*(na3-1)
      do na1=1,nat
        do na2=1,nat
          do icar=1,3
            ii = icar + 3*(na1-1)
            do jcar=1,3
              jj = jcar + 3*(na2-1)
              phiq2 (kk,ii,jj,iq) = phiq1(kcar,icar,jcar,na3,na1,na2,iq)
              diffq= abs(phiq2(kk,ii,jj,iq) - phiqn(kk,ii,jj,iq))
              if (maxdiffq.lt.diffq) maxdiffq=diffq

              !                write(4,'(1pe18.11)') real( phiq2(kk,ii,jj,iq))
              
            end do
          end do
        end do
      end do
    end do
  end do
end do
deallocate(phiq2)
!  
  write (6,"(/5x,' check: largest back phiq and original phiq  difference = ',e15.7)") maxdiffq   


END IF

CALL deallocations

 
end program q2r
!
!-----------------------------------------------------------------------
SUBROUTINE read_file ( nfile, filename )
  !-----------------------------------------------------------------------
  !
  USE common_var
  IMPLICIT NONE

  INTEGER :: nfile
  CHARACTER(len=256) :: filename(nfile)

  INTEGER :: nat_, ntyp_, ibrav_
  REAL(DP) :: celldm_(6)
  INTEGER :: nqp_, nqploc, ntmp
  INTEGER :: iun, ii, ifile, it, ia, ja, ka, ic, jc, kc, jj, kk
  INTEGER :: ic_,ia_,jc_,ja_,kc_,ka_
  CHARACTER(len=75) :: line
  CHARACTER(len=6) :: ch6
  REAL(DP) :: phir, phii

  iun = 1
  nqpts = 0
  DO ifile = 1, nfile
    WRITE(*,'(5x,a,i5,3x,a)') 'opening file #', ifile, trim(filename(ifile))
    OPEN ( UNIT=iun, FILE=filename(ifile), STATUS='old', FORM='formatted' )
    
    READ(iun,*) 
    READ(iun,*) 
    READ(iun,*) ntyp_, nat_, ibrav_, ( celldm_(ii), ii=1,6 )
    IF ( ifile .EQ. 1 ) THEN
      ntyp  = ntyp_
      nat  = nat_
      ibrav = ibrav_
      celldm(:) = celldm_(:)
      nat3 = 3 * nat
      ALLOCATE ( matq(nat3, nat3, nat3, nqptsx) )
      ALLOCATE ( xq1(3,nqptsx), xq2(3,nqptsx), xq3(3,nqptsx) )
      ALLOCATE ( amass(ntyp), ityp(nat), tau(3,nat)  )
      ALLOCATE ( atm(ntyp) )
    ELSE
      IF ( nat  .NE.  nat_  ) CALL errore('read_f', ' nat_', ifile)
      IF ( ntyp  .NE. ntyp_  ) CALL errore('read_f', ' nat_', ifile)
      IF ( ibrav .NE. ibrav_ ) CALL errore('read_f', ' nat_', ifile)
      DO ii = 1, 6
      IF ( celldm(ii) .NE. celldm(ii) ) CALL errore('read_f', ' celldm_', ifile)
      ENDDO
   ENDIF
    IF ( ibrav==0 .and. ifile==1) then
     READ(iun,*) ((at(ic,jc),ic=1,3),jc=1,3)
!     at(:,:) = at(:,:) * celldm(1)
    END IF
    DO it = 1, ntyp
      READ(iun,*) ii, atm(it), amass(it)
      IF ( ii .NE. it ) CALL errore('read_f','wrong data read',it)
    END DO
    DO ia = 1, nat
      READ(iun,*) ii, ityp(ia), ( tau(ic,ia), ic=1,3 )
      IF ( ii .NE. ia ) CALL errore('read_f','wrong data read',ia)
    END DO

    READ (iun,*)
    READ (iun,'(a)') line
    IF (line(6:21).ne.'Third derivative') then
      WRITE(*,*) line(6:21), 'wrong input thrid derivative file'
      CALL errore('','',1)
    END IF

    CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
    IF (ibrav .NE. 0) at = at / celldm(1)
    CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
    CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))

    nqploc = 0
    nqp_ = 0
110  continue
    nqploc = nqploc + 1    ! numero di matrice letta nel file
    nqpts  = nqpts + 1    ! numero totale raggiunto

    IF (nqpts.GT.nqptsx) call errore('read_f','nqptsx exceeded',1)
    READ (iun,*,END=120)
    READ (iun,*,END=120) ch6, (xq1(ic,nqpts), ic = 1, 3)
    READ (iun,*,END=120) ch6, (xq2(ic,nqpts), ic = 1, 3)
    READ (iun,*,END=120) ch6, (xq3(ic,nqpts), ic = 1, 3)

    DO ka = 1, nat
    DO ja = 1, nat
    DO ia = 1, nat
      DO kc = 1, 3
      DO jc = 1, 3
      DO ic = 1, 3
        READ(iun,*)ic_,ia_,jc_,ja_,kc_,ka_,phir,phii
        IF ( ic.NE.ic_ .or. ia.NE.ia_ .or.&
           jc.NE.jc_ .or. ja.NE.ja_ .or.&
           kc.NE.kc_ .or. ka.NE.ka_ ) THEN
           WRITE(*,*) ic_,ia_,jc_,ja_,kc_,ka_
           WRITE(*,*) ic,ia,jc,ja,kc,ka
           CALL errore('read mat','wrong indexing',1)
        END IF
        ii = ic + 3*(ia-1)
        jj = jc + 3*(ja-1)
        kk = kc + 3*(ka-1)
        matq(ii,jj,kk,nqpts) = dcmplx(phir,phii)  
      END DO
      END DO
      END DO
    END DO
    END DO
    END DO

    ntmp = nqpts
    CALL permutations (ntmp,nqptsx,nat,nqpts,at,bg, xq1, xq2, xq3, matq)
    nqp_ = nqp_ + nqpts - ntmp + 1
!    write(*,'(3(a,i5,2x))') 'matrix #:', nqploc, 'equivalents #:', nqpts - ntmp + 1, 'total read #', nqpts

    GOTO 110
120 continue
    nqploc = nqploc - 1
    nqpts = nqpts - 1

    WRITE(*,'(5x,3(a,i10,2x),/)') 'q-points read:', nqploc, 'after permutats:', nqp_, 'total #', nqpts

    CLOSE (iun)
  END DO

END SUBROUTINE read_file


!-----------------------------------------------------------------------------------------------------

subroutine  permutations (nnc,nqx,nat,nqs,at,bg, q1n,q2n,q3n,phiq)
  !
  ! Generates the matrices equivalent for permutation of the three q vectors
  ! xq1, xq2, xq3, mat are the full list of all the matrices and vectors
  ! in input this list contains nqs elements
  ! The routine generates the permutations of one matrix (with index nnc)
  ! appends the matrices to the end of the list 
  ! updates the value of nqs accordingly
  !
use common_var, only : DP
  implicit none
  ! I/O variables
  integer :: nnc, nqx, nat,  nqs
  !parameter (nqsoutx=100000)
  real(DP) :: q1n(3,nqx),q2n(3,nqx),q3n(3,nqx), at(3,3),bg(3,3)
  complex(DP) :: phiq(3*nat,3*nat,3*nat,nqx)

  ! local variables
  real(DP) :: q1(3),q2(3),q3(3), q1a(3),q2a(3),q3a(3)
  complex(DP) :: phi(3*nat,3*nat,3*nat)

  integer :: i,j,k
 logical :: equiv
 
  !equiv=.false.
  if (nnc.le.0) call errore('permutation','nnc negative',1)
  if (nnc.gt.nqs) call errore('permutation','nnc .lt. nqs',1)
  q1(:)  = q1n(:,nnc)
  q2(:)  = q2n(:,nnc)
  q3(:) = q3n(:,nnc)
  phi(:,:,:) = phiq(:,:,:,nnc)
  q1a(:)  = q1n(:,nnc)
  q2a(:)  = q2n(:,nnc)
  q3a(:) = q3n(:,nnc)

  if (equiv(q1a,q2a,at,bg) .and. equiv(q2a,q3a,at,bg).and. equiv(q1a,q3a,at,bg) ) then
    ! write (*,*) "q1 =q2 =q3"
    ! write(*,*)q1(:),q2(:),q3(:)

    ! return
  elseif ( equiv(q1a,q2a,at,bg) .or. equiv(q2a,q3a,at,bg) .or. equiv(q3a,q1a,at,bg)) then
  !  write (*,*) "qm =qn"
  !  write(*,*)q1(:),q2(:),q3(:)
    if (nqs+2.gt.nqx) call errore('permutation','nqx exceeded +2',nqs)

    ! 2 3 1 permutation
    nqs=nqs+1
    q1n(:,nqs)=q2(:)
    q2n(:,nqs)=q3(:)
    q3n(:,nqs)=q1(:)  
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(j,k,i,nqs )= phi(i,j,k)
        end do
      end do
    end do
    ! 3 1 2 permutation
    nqs=nqs+1
    q1n(:,nqs)=q3(:)
    q2n(:,nqs)=q1(:)
    q3n(:,nqs)=q2(:)
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(k,i,j,nqs)= phi(i,j,k)
        end do
      end do
    end do
  else
  !  write (*,*) "qm DIVERSO qn"
  !   write(*,*)q1(:),q2(:),q3(:)
    if (nqs+5.gt.nqx) call errore('permutation','nqx exceeded +5',nqs)
    ! 2 3 1 permutation
    nqs=nqs+1
    q1n(:,nqs)=q2(:)
    q2n(:,nqs)=q3(:)
    q3n(:,nqs)=q1(:)
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(j,k,i,nqs) = phi(i,j,k)
        end do
      end do
    end do
    ! 3 1 2 permutation
    nqs=nqs+1
    q1n(:,nqs)=q3(:)
    q2n(:,nqs)=q1(:)
    q3n(:,nqs)=q2(:)
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(k,i,j,nqs )= phi(i,j,k)
        end do
      end do
    end do
    ! 3 2 1 permutation
    nqs=nqs+1
    q1n(:,nqs)=q3(:)
    q2n(:,nqs)=q2(:)
    q3n(:,nqs)=q1(:)
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(k,j,i,nqs )= phi(i,j,k)
        end do
      end do
    end do
    ! 2 1 3 permutation
    nqs=nqs+1
    q1n(:,nqs)=q2(:)
    q2n(:,nqs)=q1(:)
    q3n(:,nqs)=q3(:)
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(j,i,k,nqs) = phi(i,j,k)
        end do
      end do
    end do
    ! 1 3 2 permutation
    nqs=nqs+1
    q1n(:,nqs)=q1(:)
    q2n(:,nqs)=q3(:)
    q3n(:,nqs)=q2(:)    
    do i=1,3*nat
      do j=1,3*nat
        do k=1,3*nat
          phiq(i,k,j,nqs )= phi(i,j,k)
        end do
      end do
    end do
  end if

  return
end subroutine permutations
!------------------------------------------------------------------------

function equiv (qi,qj,at,bg)
!-----------------------------------------------------------------------
!
use common_var, only : DP
implicit none
 real(kind=8) :: qi(3),qj(3),at(3,3),diff(3),bg(3,3)
 logical :: equiv
  equiv=.false.
  diff = qi - qj
  call cryst_to_cart(1,diff,at,-1)

  if ( abs(diff(1)-nint(diff(1))) + abs(diff(2)-nint(diff(2))) + abs(diff(3)-nint(diff(3))) .lt. 1.0d-6  ) equiv=.true.

  return
end function equiv
!----------------------------------------------------------------------------------------------
 subroutine add_to_Rlist(nRbig,nRsc,Rbig,Rsc, nRlist, nRlistx,Rlist, iR, jR, nper, ijReq, irl)

use common_var, only : DP
implicit none
 integer nRlist, nRlistx, nRbig, nRsc,ip,il 
 integer iR, jR, nper, ijReq(2,nper), irl(nper)
 real*8 Rbig(3,nRbig), Rlist(3,2,nRlistx), Rsc(3,nRsc)
 real*8 R1(3),R2(3)

! check if the two R vectors already belong to the list

do ip = 1, nper
  R1(:) = Rbig(:,iR) + Rsc(:,ijReq(1,ip))
  R2(:) = Rbig(:,jR) + Rsc(:,ijReq(2,ip))
  do il = 1, nRlist
    if ( abs(R1(1)-Rlist(1,1,il)) + abs(R1(2)-Rlist(2,1,il)) + abs(R1(3)-Rlist(3,1,il)) +&
       abs(R2(1)-Rlist(1,2,il)) + abs(R2(2)-Rlist(2,2,il)) + abs(R2(3)-Rlist(3,2,il)) .lt.1.0d-6) then
      irl(ip) = il
      goto 110
    endif
  enddo
  nRlist = nRlist + 1
  if (nRlist.gt.nRlistx) call errore('addtoRlist','nRlistx exceeded',nRlist)
  Rlist(:,1,nRlist) = R1(:)
  Rlist(:,2,nRlist) = R2(:)
  irl(ip) = nRlist
110 continue
enddo

end subroutine add_to_Rlist
!----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------

subroutine  qaccepted(nqacc,nqx,nat,nqs,at,bg,q1o,q2o,q3o,phiqo,q1n,q2n,q3n,phiqn)
  !
  ! Check whether the list contains equivalent matrices and deletes them
  ! from the list
  !
use common_var, only : DP

implicit none
  ! I/O variables
  integer :: nqx, nat, nqs,nqacc,nequiv
  !parameter (nqsoutx=100000)
  real(DP) :: q1n(3,nqx),q2n(3,nqx),q3n(3,nqx), at(3,3),bg(3,3)
  real(DP) :: q1o(3,nqx),q2o(3,nqx),q3o(3,nqx)
 complex(DP) :: phiqo(3*nat,3*nat,3*nat,nqx),phiqn(3*nat,3*nat,3*nat,nqx)

  ! local variables
  real(DP) :: q1b(3),q2b(3),q3b(3), q1a(3),q2a(3),q3a(3)
  complex(DP) :: phi(3*nat,3*nat,3*nat)

  integer :: i,j,j1,j2,j3
 logical :: equiv
! ALLOCATE(phiqn(3*nat,3*nat,3*nat,nqx))
  !equiv=.false.
 !  if (nnc.le.0) call errore('permutation','nnc negative',1)
 !  if (nnc.gt.nqs) call errore('permutation','nnc .lt. nqs',1)
 

 !  q1(:)  = q1n(:,nnc)
 !  q2(:)  = q2n(:,nnc)
 !  q3(:) = q3n(:,nnc)

nqacc=0
do i=1,nqs
  do j3=1,3*nat
    do j2=1,3*nat 
      do j1=1,3*nat
        phi(j1,j2,j3) = phiqo(j1,j2,j3,i)
      end do
    end do
  end do
  q1a(:)  = q1o(:,i)
  q2a(:)  = q2o(:,i)
  q3a(:)  = q3o(:,i)
  if (i>1) then
    nequiv=0
    look_for_equiv : &
    do j=1,i-1
      q1b(:)  = q1o(:,j)
      q2b(:)  = q2o(:,j)
      q3b(:)  = q3o(:,j)
      if(equiv(q1b,q1a,at,bg) .and. equiv(q2b,q2a,at,bg).and. equiv(q3b,q3a,at,bg) ) then
        nequiv=nequiv+1
        exit look_for_equiv
      end if
    end do &
    look_for_equiv
    
    if (nequiv.eq.0) then
      if(nqacc>nqx) call errore('qaccepted', 'nqacc too large', nqacc) 
      nqacc=nqacc+1
      q1n(:,nqacc)=q1a(:)
      q2n(:,nqacc)=q2a(:)
      q3n(:,nqacc)=q3a(:) 
      phiqn(:,:,:,nqacc)=  phi(:,:,:)
    end if
  else
    nqacc=1
    q1n(:,nqacc)=q1a(:)
    q2n(:,nqacc)=q2a(:)
    q3n(:,nqacc)=q3a(:) 
    phiqn(:,:,:,nqacc)=  phi(:,:,:)
  end if
  end do

 ! nqs = nqacc
 ! do i = 1, nqs
 !   phiqo(:,:,:,i) = phiqn(:,:,:,i)
 !   q1o(:,i) = q1n(:,1)
 !   q2o(:,i) = q2n(:,1)
 !   q3o(:,i) = q3n(:,1)
 ! enddo

  return
end subroutine qaccepted
!-----------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE deallocations
  !------------------------------------------------------------------------------
  USE common_var
  IMPLICIT NONE

  DEALLOCATE ( matq )
  DEALLOCATE ( xq1, xq2, xq3 )
  DEALLOCATE ( amass, ityp, tau )
  DEALLOCATE ( atm )

  RETURN
END SUBROUTINE deallocations
!------------------------------------------------------------------------------


!!------------------------------------------------------------------------------
!SUBROUTINE
!  !------------------------------------------------------------------------------
!  USE common_var, ONLY : dp
!  IMPLICIT NONE
!  INTEGER ::
!  REAL(DP) ::
!  RETURN
!END SUBROUTINE
!!------------------------------------------------------------------------------
