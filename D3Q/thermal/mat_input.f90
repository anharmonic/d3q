MODULE mat_input

 PRIVATE
 PUBLIC :: read_input
 PUBLIC :: read_mat2R, read_mat3R
 PUBLIC :: div_mass_mat2, div_mass_mat3
 !
 PUBLIC :: change_m2r_index
 PUBLIC :: asr2_s
 !
 PUBLIC :: change_m3r_index_asr3

 CONTAINS
!----------------------------------------------------------------------------------
subroutine read_input
  !----------------------------------------------------------------------------------
  USE mpi_base
  use constants, only : DP
  use common_variables, ONLY : &
                  filemat2, filemat3, &
                  ispeed2,ispeed3, asr2, asr3, nq1, nq2, deltaq1, &
                  deltaq2,sigma,ntemp,temp,alpha,&
                  ltest,asr3iters,x_iso,const_kuk,const_k,&
                  treshk,nxcryst,xcryst,Lcas,lbordo, &
                  ! following vars are not in input namelist:
                  ltobedone, dimrun, Fcas, iun1, &
                  nq1tot, nq2tot, tinit_0, T0, Lcasm1
  
  implicit none
  
  integer :: ix_,it,ix,irun,iun
  Real(DP),EXTERNAL :: nanosec

  namelist /input/ filemat2, filemat3, ispeed2,ispeed3, asr2, asr3, nq1, nq2, deltaq1, &
                   deltaq2,sigma,ntemp,temp,alpha,&
                   ltest,asr3iters,x_iso,const_kuk,const_k,&
                   treshk,nxcryst,xcryst,Lcas,lbordo
  open (unit=55,file='inp1.r2q')

  T0= nanosec(0.d0)
  tinit_0 = nanosec(T0)

  ! Default values
  filemat2 = ' '  ! file with the 2nd order materices
  filemat3 = ' '  ! file with the 3rd order matrices
  asr2=.true.      ! if true, imposes sum-rules on 2nd order matrices
  asr3=.true.      ! if true, imposes sum-rules on 3rd order matrices
  asr3iters = 1    ! a single iteration for asr3 (more are needed for very dense inner grids)
  ispeed2 = 0
  ispeed3 = 1
  nq1(:) = 1           ! grid dimension outer q-point loop
  nq2(:) = 30          ! grid dimension inner q-point loop
  deltaq1(:) = 0.0d0   ! grid shift     outer q-point loop
  deltaq2(:) = 0.0d0   ! grid shift     inner q-point loop
  ntemp =1
  temp(1) =9000.0d0
  sigma(1) = 2.5
  nxcryst=3
  xcryst(1)=1
  xcryst(2)=2
  xcryst(3)=3
  alpha =  0.5d0
  x_iso = 0.0d0  ! percentage of isotopes
  ltest = .false.
  const_kuk= 1.0d0 ! variable multiplying umklapp processes
  const_k = 1.0d0 ! variable multiplying klapp processes
  lbordo=.false.
  Lcas = 1889725.9885789 ! 0.01 cm 
  treshk=1.d-12

  ! reading namelist
  read (55,input)
  if (filemat2.eq.' ') call errore('reading input','wrong filename ', 1)
  if (filemat3.eq.' ') call errore('reading input','wrong filename ', 1)


  IF ( nq1(1).le.0 .or. nq1(2).le.0 .or. nq1(3).le.0 .or. &
       nq2(1).le.0 .or. nq2(2).le.0 .or. nq2(3).le.0 ) then
       write(*,*) ' nq1:', nq1
       write(*,*) ' nq2:', nq2
       call errore('main','wrong grid dimension',1)
  ENDIF  !
  nq1tot = nq1(1)*nq1(2)*nq1(3)
  nq2tot = nq2(1)*nq2(2)*nq2(3)

  dimrun = 3*ntemp
  allocate (ltobedone(dimrun))
  ltobedone(:) = .FALSE.
  DO it = 1, ntemp
     DO ix_ = 1, nxcryst
        ix = xcryst(ix_)
        irun = ix + 3*(it-1)
        ltobedone(irun) = .TRUE.
     END DO
  END DO


  iun1=29
!--------------------------------------------------------------------------------------
  ! border effects prefactor
  !  default Lcas= 0.01cm = 1889725.9885789 bohr

!  Lcasm1 = 1.0d0/ 56691779.657367d0 !0.3 cm
!   Lcasm1= 1.0d0/9448629.9428945   ! 0.05cm
   Lcasm1= 1.0d0/Lcas   !  default e' 0.01cm
  Fcas = 0.5d0
!-----------------------------------------------------------------------------------

iun=18


  return
end subroutine read_input

!-----------------------------------------------------------------------
subroutine read_mat2R
  !-----------------------------------------------------------------------
  use constants, only: DP,RY_TO_CMM1,K_BOLTZMANN_RY
  use common_variables, only : &
   ! system:
   ntypx, filemat2, &
   ntyp, ityp, at, tau, amass, bg, nat, ibrav, celldm, omega, atm, &
   ! thermal related stuff:
   const_condFI, sigmam1, sigma, nq1tot, temp,  ntemp, &
   const_iso, x_iso, const_cond, &
   ! aux variables:
   mat2_in, nRbig2_in, iRbig2_in, &
   fracmass2, nat3, nat33, nat32, sqrtm1, tempm1

  use mpi_base
  implicit none

  logical :: first
  integer :: iun, ic, jc, nt, ii, ia, ia_
  integer :: nr1, nr2, nr3
  integer :: na1, na2, j1, j2, na1_, na2_, j1_, j2_, jn1, jn2, iR, nta
  real(DP) :: raux

  real(DP) :: amass_b,deltam

  iun = 2
  open(unit=iun,file=filemat2)
  read(iun,'(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  if (ntyp.gt.ntypx) call errore('reading mat2','ntypx exceeded',1)
  if (ibrav.eq.0) then
    read(iun,*) ((at(ic,jc),ic=1,3),jc=1,3)
  endif

  allocate (ityp(nat))
  allocate (tau(3,nat))
  allocate (sqrtm1(3*nat)) 

  do nt = 1,ntyp
    read(iun,*) ii, atm(nt), amass(nt)
    if (ii.ne.nt) call errore('read_f','wrong data read matR',nt)
  end do
  do ia=1,nat
     read(iun,*) ia_, ityp(ia),(tau(ic,ia),ic=1,3)
     if (ia_.ne.ia) call errore('read_f','something wrong reading tau',ia)
  end do
  read (iun,*) nr1, nr2, nr3 
  call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  if (ibrav.ne.0) at = at / celldm(1)  !  bring at in units of alat 
  call volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
  call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))

  first = .true.
  do na1=1,nat
  do na2=1,nat 
     do j1=1,3
     do j2=1,3     
        read(iun,*) j1_, j2_, na1_, na2_
        if (j1.ne.j1_ .or. j2.ne.j2_ .or. na1.ne.na1_ .or. na2.ne.na2_ ) then
           write(*,*) j1, j2, na1, na2
           write(*,*) j1_, j2_, na1_, na2_
           call errore('reading matR','something wrong',1)
        endif
        jn1 = j1_ + (na1_-1)*3
        jn2 = j2_ + (na2_-1)*3            
        read(iun,*) nRbig2_in
        if (first) then
           allocate( iRbig2_in(3,          nRbig2_in)          )
           allocate(  mat2_in(3*nat,3*nat, nRbig2_in) )
           first = .false.
        endif
        do iR = 1, nRbig2_in
           read(iun,*) (iRbig2_in(ic,iR),ic=1,3), raux
           mat2_in(jn1,jn2,iR) = raux
        end do
     end do
     end do
  end do
  end do

  nat3  = 3*nat
  nat33 = 3*3*3*nat*nat*nat
  nat32 = 3*3*nat*nat
  do ia = 1, nat
      nta = ityp(ia)
      deltam= amass(nta)/12.0d0
      amass_b = amass(nta) + x_iso * deltam
      fracmass2  = (deltam /amass_b)**2.0d0

      do ic = 1, 3
         jn1 = ic + (ia-1)*3
         sqrtm1(jn1) = 1.d0/sqrt(amass(nta))
      end do
  end do

  const_iso = fracmass2 * x_iso *(1 - x_iso)

  do ii=1,ntemp
   sigma(ii) = sigma(ii)/RY_TO_CMM1 
 end do

  if(nrank.eq.1)  write(*,*) 'nq1tot',nq1tot,omega,celldm(1)
  if(nrank.eq.1)   write(*,*) 'sigma', ntemp, sigma(:)

  do ii = 1, ntemp
     tempm1(ii) = 1.d0/temp(ii)/K_BOLTZMANN_RY
     const_cond(ii)=1.0d0/(nq1tot*omega* temp(ii)*temp(ii)*K_BOLTZMANN_RY)
     sigmam1(ii) = 1.d0/sigma(ii)
  enddo

   const_condFI = 1.0d0/(nq1tot*omega)


  close(iun)
  return
end subroutine read_mat2R
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine read_mat3R
  !-----------------------------------------------------------------------------------------
  use constants, only : DP, eps8
  use common_variables, only : &
   ! system:
   ntypx, filemat3, &
   ntyp, ityp, at, tau, nat, ibrav, celldm, &
   ! 3rd order stuff
   nRbig3_in, iRbig3_in, mat3_in, nat3
   !
  use mpi_base
  implicit none

  integer :: ntyp_, nat_, ibrav_, ityp_(nat)
  real(DP):: celldm_(6), at_(3,3), tau_(3,nat), amass_(ntypx)
  character(len=3) :: atm_(ntypx)
  integer :: nr1_,nr2_,nr3_

  integer :: iun, ii, ic, jc, it, it_, ia, ia_
  integer :: na1, na2, na3, j1, j2, j3, na1_, na2_, na3_, j1_, j2_, j3_
  integer :: jn1, jn2, jn3, ir
  real(DP):: raux
  logical :: first

  iun = 3
  open(unit=iun,file=filemat3)
  read(iun,*) ntyp_,nat_,ibrav_,celldm_
  if (ntyp.ne.ntyp_ .or. nat.ne.nat_ .or. ibrav.ne.ibrav_ ) then
   if(nrank.eq.1)   write(*,*)  ntyp, nat, ibrav, celldm
   if(nrank.eq.1)   write(*,*)  ntyp_,nat_,ibrav_,celldm_
     call warning('read mat3R','something wrong',1)
  endif
  do ii=1,6
     if(celldm_(ii).ne.celldm(ii)) then 
      if(nrank.eq.1)   write(*,*) celldm(ii), celldm_(ii)
        call warning('read mat3R','wrong celldm',ii)
     end if
  end do
  if (ibrav_.eq.0) then
     read(iun,*) ((at_(ic,jc),ic=1,3),jc=1,3)
     do jc = 1, 3
     do ic = 1, 3
        if ( abs(at_(ic,jc)-at(ic,jc)).gt.eps8 ) then
         if(nrank.eq.1)   write(*,*) ic, jc
         if(nrank.eq.1)   write(*,*) at
          if(nrank.eq.1)  write(*,*) at_
           call warning('read mat3R','wrong at',1)
        end if
     enddo
     enddo
  endif
  do it = 1, ntyp_
     read(iun,*) it_, atm_(it), amass_(it)
     if (it.ne.it_) call errore('read mat3R','wrong type read',it)
  end do
  do ia = 1, nat_
     read(iun,*) ia_, ityp_(ia),(tau_(ic,ia),ic=1,3)
     if (ia.ne.ia_)             call errore('read mat3R','wrong ia read',ia)
     if (ANY(ABS(tau_(:,ia)-tau(:,ia))>eps8)) &
                                call errore('read mat3R','wrong tau read',ia)
     if (ityp_(ia).ne.ityp(ia)) call errore('read mat3R','wrong ityp read',ia)
     do ic = 1, 3
        if ( abs(tau_(ic,ia)-tau(ic,ia)).gt.eps8 ) then
          if(nrank.eq.1)   write(*,*) ic, ia
          if(nrank.eq.1)  write(*,*) tau
          if(nrank.eq.1)  write(*,*) tau_
          call warning('read mat3R','wrong tau',1)
        endif
     enddo
  enddo

  read (iun,*) nr1_, nr2_, nr3_ 

  first = .true.
  do na1 = 1, nat
  do na2 = 1, nat 
  do na3 = 1, nat
     do j1 = 1, 3
     do j2 = 1, 3 
     do j3 = 1, 3
        read(iun,*) j1_,j2_,j3_,na1_,na2_,na3_
        if ( j1.ne.j1_   .or. j2.ne.j2_   .or. j3_.ne.j3   .or. &
             na1.ne.na1_ .or. na2.ne.na2_ .or. na3_.ne.na3      ) then
         if(nrank.eq.1)   write(*,*) j1, j2, j3, na1, na2, na3
         if(nrank.eq.1)   write(*,*) j1_, j2_,j3_, na1_, na2_, na3_
           call warning('reading matR3','something wrong',1)
        endif
        jn1 = j1_ + 3*(na1_-1)
        jn2 = j2_ + 3*(na2_-1)
        jn3 = j3_ + 3*(na3_-1)
        read(iun,*) nRbig3_in
        if (first) then
          allocate( iRbig3_in ( 3, 2, nRbig3_in )             )
          allocate(  mat3_in  ( nat3, nat3, nat3, nRbig3_in ) )
          first = .false.
        endif
        do iR = 1, nRbig3_in
           read(iun,*) (iRbig3_in(ic,1,iR),ic=1,3), &
                       (iRbig3_in(ic,2,iR),ic=1,3), raux
           mat3_in(jn1,jn2,jn3,iR) = raux
        end do
     end do
     end do
     end do
  end do
  end do
  end do
  close(iun)

  ! I redefine the meaning of the iR vactors in order to simplifiy
  ! the Fourier transform in setupmat3; if you comment the following lines
  ! you should change setupmat3 accordingly
  do iR=1,nRbig3_in
    do ic = 1, 3
       iRbig3_in(ic,1,iR) = iRbig3_in(ic,2,iR) - iRbig3_in(ic,1,iR)
    enddo
  enddo

return
end subroutine read_mat3R
!-----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
subroutine asr2_s(nRbig2t,nat,mat2,iR2,asr2)
  use constants, only : DP
  implicit none 
  logical :: asr2
  integer, intent(in):: nat,nRbig2t,iR2(3,nRbig2t)
  integer :: j1,j2,na1,na2,iR, iR0, im1, im2, im2_
  complex(DP) :: sum1
  complex(DP) :: mat2(3*nat,3*nat,nRbig2t)

  if (.not.asr2) return
  iR0 = -1
  do ir = 1, nRbig2t
     if ( iR2(1,ir)==0 .and. iR2(2,ir)==0 .and. iR2(3,ir)==0 ) iR0 = ir
  enddo
  if (iR0.eq.-1) call errore('asr2',' something wrong 123',1)

  do j1=1,3
  do j2=1,3
  do na1=1,nat
     sum1 = 0.0d0
     im1 = j1 + (na1-1)*3
     do na2=1,nat
     im2 = j2 + (na2-1)*3
     do iR=1,nRbig2t
        sum1 = sum1 + mat2(im1,im2,iR)
     enddo
     enddo
     im2_ = j2 + (na1-1)*3
     mat2(im1,im2_,iR0) = mat2(im1,im2_,iR0) - sum1
  enddo
  enddo
  enddo
  return
end subroutine asr2_s

!-----------------------------------------------------------------------
subroutine change_m2r_index
  !-----------------------------------------------------------------------
  use common_variables
  implicit none
  logical, allocatable :: ldone(:)
  integer :: ir, jr
  integer imin(3), imax(3), ii, ic, i1, i2, i3

  if (ispeed2.eq.1 .or. ispeed2.eq.2 ) then
     nRbig2t   = nRbig2_in
     allocate ( iRbig2 (3, nRbig2t)        )
     allocate (  Rbig2 (3, nRbig2t)        )
     allocate (   mat2 (nat3,nat3,nRbig2t) )
     iRbig2     = iRbig2_in
     Rbig2(:,:) = dfloat(iRbig2(:,:))
     mat2 = mat2_in
     deallocate ( iRbig2_in )
     deallocate (   mat2_in )
     return
  endif
     
     do ic = 1, 3
        imin(ic) = iRbig2_in(ic,1)
        imax(ic) = iRbig2_in(ic,1)
        do ii = 1, nRbig2_in
           imin(ic) = min ( iRbig2_in(ic,ii),imin(ic) )
           imax(ic) = max ( iRbig2_in(ic,ii),imax(ic) )
        enddo
        nRbig2(ic) = imax(ic) - imin(ic) + 1
     enddo
     nRbig2t = nRbig2(1)*nRbig2(2)*nRbig2(3)
     allocate ( iRbig2 (3, nRbig2t)        )
     allocate (  Rbig2 (3, nRbig2t)        )
     allocate (   mat2 (nat3,nat3,nRbig2t) )
     ii = 0
     do i1 = imin(ic_out), imax(ic_out)
        do i2 = imin(ic_med), imax(ic_med)
           do i3 = imin(ic_inn), imax(ic_inn)
              ii = ii + 1
              iRbig2(ic_out,ii) = i1
              iRbig2(ic_med,ii) = i2
              iRbig2(ic_inn,ii) = i3
           enddo
        enddo
     enddo
     
     allocate (ldone(nRbig2_in))
     mat2 = 0.d0
     ldone = .false.
     do ir = 1, nRbig2t
        do jr = 1, nRbig2_in
           if ( iRbig2(1,ir)==iRbig2_in(1,jr) .and. &
                iRbig2(2,ir)==iRbig2_in(2,jr) .and. &
                iRbig2(3,ir)==iRbig2_in(3,jr) ) then
              mat2(:,:,ir) =  mat2_in(:,:,jr)
              if (ldone(jr)) call errore('change_mr2_index','something wrong 1',1)
              ldone(jr) = .true.
              goto 120
           endif
        enddo
        mat2(:,:,ir) = 0.d0
120     continue
     enddo
     do jr = 1, nRbig2_in
        if (.not.ldone(jr)) call errore('change_nr2_index','something wrong 2',1)
     enddo
     
     Rbig2(:,:) = dfloat(iRbig2(:,:))
     
     deallocate ( iRbig2_in )
     deallocate (   mat2_in )
     deallocate ( ldone     )
     
     return
   end subroutine change_m2r_index
   !-----------------------------------------------------------------------
   
!
!-----------------------------------------------------------------------
subroutine change_m3r_index_asr3 (mat3_input)
  !-----------------------------------------------------------------------
  use constants, only : DP
  use common_variables
  IMPLICIT NONE
  complex(DP) :: mat3_input(nat33,nRbig3_in)
  complex(DP), ALLOCATABLE :: mat3_rr_in(:,:,:)
  integer, parameter :: nrmax=10000
  integer :: iRou_(3,nrmax), ir, jr, il, iir
  LOGICAL, ALLOCATABLE :: lmatin(:,:)
  logical :: already_there, ldone(nRbig3_in)

  !write(*,*) ' lab 001'

  nRbig32 = 1
  iRou_(:,nRbig32) = iRbig3_in(:,1,1)
  IF ( iRbig3_in(1,2,1).NE.iRou_(1,1) .OR.   &
       iRbig3_in(2,2,1).NE.iRou_(2,1) .OR.   &
       iRbig3_in(3,2,1).NE.iRou_(3,1) ) THEN
     nRbig32 = 2
     iRou_(:,nRbig32) = iRbig3_in(:,2,1)
  END IF

  DO il = 1, 2
  DO ir = 2, nRbig3_in
     already_there = .FALSE.
     DO jr = 1, nRbig32
        already_there = already_there .or.     &
        ( iRbig3_in(1,il,ir).eq.iRou_(1,jr) .and.   &
          iRbig3_in(2,il,ir).eq.iRou_(2,jr) .and.   &
          iRbig3_in(3,il,ir).eq.iRou_(3,jr) )
     END DO
     if (.not.already_there) then
        nRbig32 = nRbig32 + 1
        if (nRbig32.gt.nrmax) call errore('change_m3r_index','nrmax exceeded',1)
        iRou_(:,nRbig32) = iRbig3_in(:,il,ir)
     endif
  END DO
  END DO

  allocate (iRbig32(3,nRbig32))
  allocate ( Rbig32(3,nRbig32))
  ALLOCATE ( lmatin(nRbig32,nRbig32) )
  do ir = 1, nRbig32
     iRbig32(:,ir) = iRou_(:,ir)
     Rbig32(:,ir) = dfloat(iRbig32(:,ir))
  enddo

  allocate (mat3_rr(nat33,nRbig32,nRbig32))
  mat3_rr = 0.d0
  ldone(:) = .false.
  lmatin(:,:) = .FALSE.
  do ir = 1, nRbig32
    do jr = 1, nRbig32
       do iir = 1, nRbig3_in
          if ( iRbig32(1,ir)==iRbig3_in(1,1,iir) .and. &
               iRbig32(2,ir)==iRbig3_in(2,1,iir) .and. &
               iRbig32(3,ir)==iRbig3_in(3,1,iir) .and. &
               iRbig32(1,jr)==iRbig3_in(1,2,iir) .and. &
               iRbig32(2,jr)==iRbig3_in(2,2,iir) .and. &
               iRbig32(3,jr)==iRbig3_in(3,2,iir) ) then
             mat3_rr(:,ir,jr) = mat3_input(:,iir)
             lmatin(ir,jr) = .TRUE.
             if (ldone(iir)) call errore('change_m3r_index','something wrong 1',1)
             ldone(iir) = .true.
             goto 130
          endif
       enddo
130 continue
     enddo
  enddo

  do ir = 1, nRbig3_in
     if (.not.ldone(ir)) call errore('change_m3r_index','something wrong 2',1)
  enddo

!  IF (asr3) CALL set_asr3 ( nat3, nRbig32, iRbig32, mat3_rr  )
  IF (asr3) CALL set_asr3_ml ( nat3, nRbig32, iRbig32, mat3_rr, lmatin  )

! NOTA: DOVREBBERO FUNGERE ENTRAMBI
!  IF ( ispeed3 .LE. 2 ) THEN
  IF ( ispeed3 .EQ. 1 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     nRbig31t = nRbig32
     allocate (iRbig31(3,nRbig31t))
     allocate ( Rbig31(3,nRbig31t))
     iRbig31(:,:) = iRbig32(:,:)
      Rbig31(:,:) =  Rbig32(:,:)
  ELSE
     do ir = 1, nRbig32
        iRou_(:,ir) = iRbig32(:,ir)
     enddo
     nRbig31t = nRbig32
     call reorder_r_index ( nrmax, nRbig31t, nRbig31, iRou_)
     allocate (iRbig31(3,nRbig31t))
     allocate ( Rbig31(3,nRbig31t))
     do ir = 1, nRbig31t
        iRbig31(:,ir) = iRou_(:,ir)
        Rbig31(:,ir) = dfloat(iRbig31(:,ir))
     enddo

     ALLOCATE ( mat3_rr_in(nat33,nRbig32,nRbig32) )
     mat3_rr_in(:,:,:) = mat3_rr(:,:,:)
     DEALLOCATE (mat3_rr)
     allocate (mat3_rr(nat33,nRbig31t,nRbig32))
     mat3_rr = 0.d0
     do jr = 1, nRbig32
        ldone(:) = .false.
        do ir = 1, nRbig31t
           do iir = 1, nRbig32
              if ( iRbig31(1,ir)==iRbig32(1,iir) .and. &
                   iRbig31(2,ir)==iRbig32(2,iir) .and. &
                   iRbig31(3,ir)==iRbig32(3,iir) ) then
                 mat3_rr(:,ir,jr) = mat3_rr_in(:,iir,jr)
                 if (ldone(iir)) call errore('change_m3r_index','something wrong 3',1)
                 ldone(iir) = .true.
                 goto 110
              endif
           enddo
110 continue
        enddo
        do ir = 1, nRbig32
           if (.not.ldone(ir)) call errore('change_m3r_index','something wrong 4',1)
        enddo
     enddo
  END IF

  deallocate ( mat3_in   )
  deallocate ( iRbig3_in )
!  write(*,*) ' lab 002'

  return
end subroutine change_m3r_index_asr3
!-----------------------------------------------------------------------
!

!--------------------------------------------------------------------------------
subroutine set_asr3_ml ( nat3, nRbig, iRbig, mat, lmatin  )
  !------------------------------------------------------------------------------
  use constants, only : DP
  use common_variables, only : asr3iters
  use mpi_base
  implicit none

  INTEGER :: nat3, nRbig, iRbig(3,nRbig)
  COMPLEX(DP) :: mat(nat3,nat3,nat3,nRbig,nRbig)
  LOGICAL :: lmatin(nRbig,nRbig)

  INTEGER :: ir0, ii
  INTEGER, ALLOCATABLE :: irm(:), irmm(:,:)
  REAL(DP) :: chi

  ALLOCATE ( irm  (nRbig)       )
  ALLOCATE ( irmm (nRbig,nRbig) )

  write(*,'(5x,a,i5)') "ASR3: ml style symmetrized with asr3iters = ", asr3iters

  CALL calc_iR_tables ( nRbig, iRbig, lmatin, ir0, irm, irmm )

  CALL calc_chi_ortog ( nat3, nRbig, mat, chi )
   if(nrank.eq.1)  write(*,'(7x,a,e12.3)') 'chi0:', chi

  DO ii = 1, asr3iters
     CALL ortog_1st_index ( nat3, nRbig, mat )
     CALL symm_m3_indexes ( nat3, nRbig, mat, lmatin, irm, irmm )
     CALL calc_chi_ortog ( nat3, nRbig, mat, chi )
    if(nrank.eq.1)   write(*,'(7x,a,i4,a,e12.3)') 'chi',ii,':', chi
  END DO

!  WRITE(*,*) 'STOPPING....'
!  STOP

  DEALLOCATE ( irm  )
  DEALLOCATE ( irmm )

  RETURN
END SUBROUTINE set_asr3_ml
!--------------------------------------------------------------------------------!--------------------------------------------------------------------------------
SUBROUTINE ortog_1st_index ( nat3, nRbig, mat )
  !------------------------------------------------------------------------------
  USE constants, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nat3, nRbig
  COMPLEX(DP) :: mat(nat3,nat3,nat3,nRbig,nRbig)

  COMPLEX(DP), ALLOCATABLE :: mat_(:,:,:,:,:)
  COMPLEX(DP) :: sum1, sum2
  INTEGER :: nat, ir1, ir2, im1, im2, im3, ic, ia
  ALLOCATE ( mat_(nat3,nat3,nat3,nRbig,nRbig) )

  nat = nat3/3
  DO ir2 = 1, nRbig
  DO im3 = 1, nat3
  DO im2 = 1, nat3
  DO  ic = 1, 3
     sum1 = 0.d0
     sum2 = 0.d0
     DO ir1 = 1, nRbig
     DO ia = 1, nat
        im1 = ic + (ia-1)*3
        sum1 = sum1 + mat ( im1, im2, im3, ir2 ,ir1 )
        sum2 = sum2 + mat ( im1, im2, im3, ir2 ,ir1 ) **2
     END DO
     END DO
     IF ( sqrt(abs(real(sum2))) .GT. 1.d-10 ) THEN
        DO ir1 = 1, nRbig
        DO ia = 1, nat
           im1 = ic + (ia-1)*3
           mat_(im1,im2,im3,ir2,ir1) = mat(im1,im2,im3,ir2,ir1) &
                                     - mat(im1,im2,im3,ir2,ir1)**2 * sum1 / sum2
        END DO
        END DO
     END IF
  END DO
  END DO
  END DO
  END DO

  mat(:,:,:,:,:) = mat_(:,:,:,:,:)

  DEALLOCATE (mat_)
  RETURN
END SUBROUTINE ortog_1st_index
!--------------------------------------------------------------------------------!--------------------------------------------------------------------------------
SUBROUTINE calc_chi_ortog ( nat3, nRbig, mat, chi )
  !------------------------------------------------------------------------------
  USE constants, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nat3, nRbig
  REAL(DP) :: chi
  COMPLEX(DP) :: mat(nat3,nat3,nat3,nRbig,nRbig)

  COMPLEX(DP) :: sum1
  INTEGER :: nat, ir1, ir2, im1, im2, im3, ic, ia

  nat = nat3/3
  chi = 0.d0
  DO ir2 = 1, nRbig
  DO im3 = 1, nat3
  DO im2 = 1, nat3
  DO  ic = 1, 3
     sum1 = 0.d0
     DO ir1 = 1, nRbig
     DO ia = 1, nat
        im1 = ic + (ia-1)*3
        sum1 = sum1 + mat ( im1, im2, im3, ir2 ,ir1 )
     END DO
     END DO
     chi = chi + sqrt(sum1 * conjg(sum1))
  END DO
  END DO
  END DO
  END DO

  RETURN
END SUBROUTINE calc_chi_ortog
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
SUBROUTINE symm_m3_indexes ( nat3, nRbig, mat, lmatin, irm, irmm )
  !------------------------------------------------------------------------------
  USE constants, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nat3, nRbig, irm(nRbig), irmm(nRbig,nRbig)
  LOGICAL :: lmatin(nRbig,nRbig)
  COMPLEX(DP) :: mat(nat3,nat3,nat3,nRbig,nRbig)

  COMPLEX(DP), ALLOCATABLE :: mat_(:,:,:,:,:)
  INTEGER :: ir1, ir2, irm1, irm2, ir1m2, ir2m1, im1, im2, im3
  ALLOCATE ( mat_(nat3,nat3,nat3,nRbig,nRbig) )

  DO ir2 = 1, nRbig
  DO ir1 = 1, nRbig
     IF ( lmatin(ir1,ir2) ) THEN
        irm1  = irm(ir1)
        irm2  = irm(ir2)
        ir1m2 = irmm(ir1,ir2)
        ir2m1 = irmm(ir2,ir1)
        DO im1 = 1, nat3
        DO im2 = 1, nat3
        DO im3 = 1, nat3
! Devi usare questa se la ridefinizione non viene fatta
!        mat_ (im1,im2,im3,ir1,ir2) = 1.d0/6.d0 * (     &
!                mat ( im1, im2, im3, ir1   , ir2   ) + &
!                mat ( im2, im3, im1, ir2m1 , irm1  ) + &
!                mat ( im3, im1, im2, irm2  , ir1m2 ) + &
!                mat ( im1, im3, im2, ir2   , ir1   ) + &
!                mat ( im2, im1, im3, irm1  , ir2m1 ) + &
!                mat ( im3, im2, im1, ir1m2 , irm2  ) )
! Devi usare questa a causa della ridefinizione di R1 ed R2 fatta in read mat
           mat_ (im1,im2,im3,ir1,ir2) = 1.d0/6.d0 * (     &
                   mat ( im1, im2, im3, ir1   , ir2   ) + &
                   mat ( im2, im3, im1, irm2  , ir1m2 ) + &
                   mat ( im3, im1, im2, ir2m1 , irm1  ) + &
                   mat ( im1, im3, im2, irm1  , ir2m1 ) + &
                   mat ( im2, im1, im3, ir2   , ir1   ) + &
                   mat ( im3, im2, im1, ir1m2 , irm2  ) )
!           WRITE(*,'(5i6)') ir2, ir1, im1, im2, im3
!           WRITE(*,'(6f12.6)') &
!                   real( mat ( im1, im2, im3, ir1   , ir2   )) ,&
!                   real( mat ( im2, im3, im1, irm2  , ir1m2 )) ,&
!                   real( mat ( im3, im1, im2, ir2m1 , irm1  )) ,&
!                   real( mat ( im1, im3, im2, irm1  , ir2m1 )) ,&
!                   real( mat ( im2, im1, im3, ir2   , ir1   )) ,&
!                   real( mat ( im3, im2, im1, ir1m2 , irm2  ))
        END DO
        END DO
        END DO
     END IF
  END DO
  END DO
  mat(:,:,:,:,:) = mat_(:,:,:,:,:)

  DEALLOCATE (mat_)

  RETURN
END SUBROUTINE symm_m3_indexes
!--------------------------------------------------------------------------------!--------------------------------------------------------------------------------
SUBROUTINE calc_iR_tables ( nRbig, iRbig, lmatin, ir0, irm, irmm )
  !------------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nRbig, iRbig(3,nRbig), ir0, irm(nRbig), irmm(nRbig,nRbig)
  LOGICAL :: lmatin(nRbig,nRbig)

  INTEGER :: ir, jr, ir1, ir2, irmm_


  ir0 = -1
  DO ir = 1, nRbig
     IF ( iRbig(1,ir)==0 .AND. iRbig(2,ir)==0 .AND. iRbig(3,ir)==0 ) ir0 = ir
  END DO
  if (ir0.eq.-1) CALL errore('set asr3','unexpected 0022',ir)

  DO ir = 1, nRbig
     DO jr = 1, nRbig
        IF ( iRbig(1,ir) .eq. -iRbig(1,jr) .AND. &
             iRbig(2,ir) .eq. -iRbig(2,jr) .AND. &
             iRbig(3,ir) .eq. -iRbig(3,jr) ) THEN
           irm(ir) = jr
           GOTO 140
        END IF
     END DO
     CALL errore('set asr3','unexpected 0023',ir)
140 CONTINUE
  END DO

  DO ir1 = 1, nRbig
  DO ir2 = 1, nRbig
     IF ( lmatin(ir1,ir2) ) THEN
        DO irmm_ = 1, nRbig
           IF ( iRbig(1,irmm_) .eq. iRbig(1,ir1) - iRbig(1,ir2) .AND. &
                iRbig(2,irmm_) .eq. iRbig(2,ir1) - iRbig(2,ir2) .AND. &
                iRbig(3,irmm_) .eq. iRbig(3,ir1) - iRbig(3,ir2) ) THEN
              irmm(ir1,ir2) = irmm_
              GOTO 150
           END IF
        END DO
        write(*,*) ir1, iRbig(:,ir1)
        write(*,*) ir2, iRbig(:,ir2)
        CALL errore('set asr3','unexpected 0024',ir)
     END IF
150 CONTINUE
  END DO
  END DO

  RETURN
END SUBROUTINE calc_iR_tables
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine reorder_r_index ( nrx, nr, nrc, ir)
  !-----------------------------------------------------------------------
  ! Rorder the R vector indeces to mqke them more suitqble for the Fourier
  ! interpolation
  !
  use common_variables, only : ic_inn, ic_med, ic_out
  implicit none
  integer :: nrx, nr, nrc(3), ir(3,nrx)
  integer :: imin(3), imax(3), ii, ic, i1, i2, i3

  do ic = 1, 3
     imin(ic) = ir(ic,1)
     imax(ic) = ir(ic,1)
     do ii = 1, nr
        imin(ic) = min ( ir(ic,ii),imin(ic) )
        imax(ic) = max ( ir(ic,ii),imax(ic) )
     enddo
  enddo

  ii = 0
  do i1 = imin(ic_out), imax(ic_out)
  do i2 = imin(ic_med), imax(ic_med)
  do i3 = imin(ic_inn), imax(ic_inn)
     ii = ii + 1
     if (ii.gt.nrx) call errore('reorder_r_index','nrmax exceeded',1)
     ir(ic_out,ii) = i1
     ir(ic_med,ii) = i2
     ir(ic_inn,ii) = i3
  enddo
  enddo
  enddo

  nr = ii
  do ic = 1, 3
     nrc(ic) = imax(ic) - imin(ic) + 1
  enddo
  if (nr.ne.nrc(1)*nrc(2)*nrc(3)) call errore('reorder_r_index','something wrong',1)

  return
end subroutine reorder_r_index
!-----------------------------------------------------------------------


subroutine div_mass_mat2 (nat3,nRbig2t,mat,sqrtm1)
  use constants, only : DP
  implicit none
  integer :: nat3, nRbig2t
  complex(DP) :: mat(nat3,nat3,nRbig2t)
  real(DP) :: sqrtm1(nat3)
  integer :: im, jm, ir

  do ir = 1, nRbig2t
  do jm = 1, nat3
  do im = 1, nat3
     mat(im, jm, ir) = mat(im, jm, ir) * sqrtm1(im) * sqrtm1(jm)
  enddo
  enddo
  enddo

end subroutine div_mass_mat2
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine div_mass_mat3 (nat3,nRbig31t,nRbig32,mat,sqrtm1)
  use constants, only : DP
  implicit none
  integer     :: nat3, nRbig31t, nRbig32
  complex(DP) :: mat(nat3,nat3,nat3,nRbig31t,nRbig32)
  real(DP)    :: sqrtm1(nat3)
  integer     :: im, jm, km, ir, jr

  do jr = 1, nRbig32
  do ir = 1, nRbig31t
  do km = 1, nat3
  do jm = 1, nat3
  do im = 1, nat3
     mat(im, jm, km, ir, jr) = mat(im, jm, km, ir, jr) * sqrtm1(im) * sqrtm1(jm) * sqrtm1(km)
  enddo
  enddo
  enddo
  enddo
  enddo

end subroutine div_mass_mat3
!-----------------------------------------------------------------------------------------

END MODULE mat_input
