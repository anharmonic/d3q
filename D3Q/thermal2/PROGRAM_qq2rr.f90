!
! Written by Lorenzo Paulatto, 2015.
!
MODULE program_qq2rr
  USE kinds, ONLY : DP
  !
  TYPE d3_list
    COMPLEX(DP),ALLOCATABLE :: d(:,:,:, :,:,:)
    REAL(DP) :: xq1(3), xq2(3), xq3(3)
    !REAL(DP) :: xq1_true(3), xq2_true(3), xq3_true(3)
  END TYPE
  
  CONTAINS
  ! \/o\________\\\________________\\/\_________________________/^>
  FUNCTION check_in_grid(nq, xq, at, skip) RESULT(iq)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3)
    LOGICAL,INTENT(out) :: skip
    REAL(DP),INTENT(in) :: xq(3), at(3,3)
    REAL(DP) :: aq(3)
    INTEGER  :: iq(3), i
    aq = xq
    CALL cryst_to_cart (1,aq,at,-1)
    aq = aq*nq
    iq = NINT(aq)
    IF(ANY(ABS(iq-aq)>1.d-4)) THEN
      !WRITE(*,'(a,3f10.4,3x,3f10.4,3i3)') "skippa", xq, aq, iq
      skip = .true.
      iq = -1
      RETURN
    ELSE
      !iq = iq+1
      skip = .false.
      DO i =1,3
        iq(i) = MODULO(iq(i), nq(i))+1
      ENDDO
    ENDIF
    
  END FUNCTION
  !
  ! \/o\________\\\________________\\/\_________________________/^>
  SUBROUTINE read_d3_matrices(nq, nq_trip, S, d3grid)
    USE kinds,       ONLY : DP
    USE d3matrix_io, ONLY : read_d3dyn_xml
    USE d3_shuffle,  ONLY : nperms, d3perms_order, d3_shuffle_global, d3_shuffle_equiv
    USE d3_basis,    ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx
    USE parameters,  ONLY : ntypx
    USE input_fc,    ONLY : ph_system_info, same_system, aux_system
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3), nq_trip
    TYPE(ph_system_info),INTENT(out) :: S
    TYPE(d3_list),INTENT(inout) ::  d3grid(nq_trip)

    REAL(DP) :: xq(3,3)
    COMPLEX(DP),ALLOCATABLE   :: d3(:,:,:, :,:,:), p3(:,:,:)
    INTEGER :: i,j,k, a,b,c, ios, ntyp, iperm, iq_trip, nxr_list, nxr_list_old, iq_aux
    LOGICAL :: first, skip
    INTEGER,ALLOCATABLE :: found(:,:,:,:,:,:)
    INTEGER :: countq(nq_trip)
    CHARACTER(len=512) :: filename
    INTEGER :: iqa(3), iqb(3), iqc(3)
    REAL(DP) :: xq_shift(3,3)
    TYPE(ph_system_info) :: Sx
    !
    ALLOCATE(found(nq(1),nq(2),nq(3),nq(1),nq(2),nq(3)))
    found = -1
    countq = 0
    first = .true.
    iq_trip = 0
    !
    DO 
      READ(*, '(a512)',iostat=ios) filename
      IF(ios/=0) EXIT
!       WRITE(*,*) "Reading '", TRIM(filename),"'..."

      IF(first)THEN
        first=.false.
        CALL read_d3dyn_xml(filename, xq(:,1), xq(:,2), xq(:,3), d3=d3, nat=S%nat, atm=S%atm, ntyp=S%ntyp, &
                            ityp=S%ityp, ibrav=S%ibrav, celldm=S%celldm, at=S%at, amass=S%amass,&
                            tau=S%tau, seek=.false.)
        CALL aux_system(S)
        CALL latgen( S%ibrav, S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega )
        S%at=S%at/S%celldm(1)
        CALL recips(S%at(:,1), S%at(:,2), S%at(:,3), S%bg(:,1), S%bg(:,2), S%bg(:,3))
        ALLOCATE(p3(S%nat3, S%nat3, S%nat3))
        !
          !WRITE(*,*), 100, TRIM(filename)
      ELSE
        CALL read_d3dyn_xml(filename, xq(:,1), xq(:,2), xq(:,3), d3=d3, nat=Sx%nat, atm=Sx%atm, ntyp=Sx%ntyp, &
                            ityp=Sx%ityp, ibrav=Sx%ibrav, celldm=Sx%celldm, at=Sx%at, amass=Sx%amass,&
                            tau=Sx%tau, seek=.false.)
        CALL latgen( Sx%ibrav, Sx%celldm, Sx%at(:,1), Sx%at(:,2), Sx%at(:,3), Sx%omega )
        Sx%at=Sx%at/Sx%celldm(1)
        CALL recips(Sx%at(:,1), Sx%at(:,2), Sx%at(:,3), Sx%bg(:,1), Sx%bg(:,2), Sx%bg(:,3))
        IF(.not.same_system(S,Sx)) CALL errore("qq2rr", "not same system "//TRIM(filename), 1)
          !WRITE(*,*), 200, TRIM(filename)
      ENDIF
      PERM_LOOP : &
      DO iperm = 1,nperms
        a = d3perms_order(1,iperm)
        b = d3perms_order(2,iperm)
        c = d3perms_order(3,iperm)
          !WRITE(*,*), 300+iperm, TRIM(filename)
        !
        ! Check if the points are in the grid:
        iqb = check_in_grid(nq, xq(:,b), S%at, skip)
          !WRITE(*,*), 310+iperm, TRIM(filename)
        IF(skip) CYCLE PERM_LOOP
        iqc = check_in_grid(nq, xq(:,c), S%at, skip)
          !WRITE(*,*), 320+iperm, TRIM(filename)
        IF(skip) CYCLE PERM_LOOP
        !
        ! Add this point to the grid, if an equivalent was not already added
        IF(found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3)) < 0) THEN
          iq_trip = iq_trip + 1
          WRITE(*,'(i6,3x,a,3(3f10.4,3x),2(3i3,3x),3i2,3x,a,i6  )') iq_trip, "xq ", xq(:,a), xq(:,b), xq(:,c), iqb, iqc, &
            a,b,c, TRIM(filename),iperm
          found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3)) = iq_trip
          countq(iq_trip) = 1
          !
          ! Repack the matrix, shuffle it and repack it again
          CALL d3_6idx_2_3idx(S%nat, d3, p3)
          !CALL d3_shuffle_simple( a,b,c, .false., p3)
          CALL d3_shuffle_global( 1,2,3, a,b,c, .false., p3)
          !CALL d3_shuffle_global( 1,2,3, a,b,c, d3perms_conjg(iperm), p3)
          !CALL d3_shuffle_equiv( 1,2,3, a,b,c, .not.d3perms_conjg(iperm), p3)
          CALL d3_3idx_2_6idx(S%nat, p3, d3)
          !
          ALLOCATE(d3grid(iq_trip)%d(3,3,3, S%nat,S%nat,S%nat))
          d3grid(iq_trip)%d   = d3
          !
          !  Recenter xq2 qnd xq3 in the first unit cell (it makes no difference really)
!           xq_shift(:,2) = (iqb-1)/DBLE(nq)
!           xq_shift(:,3) = (iqc-1)/DBLE(nq)
!           iqa = check_in_grid(nq, xq(:,a), S%at, skip)
!           xq_shift(:,1) = (iqa-1)/DBLE(nq)
!           
!           !xq_shift(:,1) = -(xq_shift(:,2)+xq_shift(:,3))
!           CALL cryst_to_cart(3, xq_shift, S%bg, +1)
!           d3grid(iq_trip)%xq1(:) = xq_shift(:,1)
!           d3grid(iq_trip)%xq2(:) = xq_shift(:,2)
!           d3grid(iq_trip)%xq3(:) = xq_shift(:,3)
          !
          d3grid(iq_trip)%xq1(:) = xq(:,a)
          d3grid(iq_trip)%xq2(:) = xq(:,b)
          d3grid(iq_trip)%xq3(:) = xq(:,c)
        ELSE
!           iq_aux = found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3))
!           countq(iq_aux) = countq(iq_aux) +1
!           CALL d3_6idx_2_3idx(S%nat, d3, p3)
!           CALL d3_shuffle_global( 1,2,3, a,b,c, .false., p3)
!           CALL d3_3idx_2_6idx(S%nat, p3, d3)
!           IF(ANY( ABS(d3grid(iq_aux)%d - d3)>1.d-6) )THEN
!             WRITE(*,'(i6,3x,a,3(3f10.4,3x),2(3i3,3x),3i2,3x,a)') iq_aux, "xqR", xq(:,a), xq(:,b), xq(:,c), iqb, iqc, &
!             a,b,c, TRIM(filename)
!             print*, iq_aux, iq_trip
!             print*, d3grid(iq_aux)%xq1(:),d3grid(iq_aux)%xq2(:) ,d3grid(iq_aux)%xq3(:) 
!             print*, xq
! !             WRITE(*,'(3(2f12.6,2x))') d3grid(iq_aux)%d - d3
!             !stop 666
!          ENDIF
          !d3grid(iq_aux)%d   = d3 + d3grid(iq_aux)%d
          
        ENDIF
      ENDDO &
      PERM_LOOP
      
    ENDDO
    
!     DO iq_aux = 1,nq_trip
!       print*, iq_aux, countq(iq_aux)
!       d3grid(iq_aux)%d = d3grid(iq_aux)%d/DBLE(countq(iq_aux))
!     ENDDO
    
    IF(first) CALL errore("read_d3_matrices","I found nothing to read",1)
    DEALLOCATE(p3)
    
    IF(ANY(found<0)) THEN
      WRITE(*,*) "Expecting:", nq_trip, "triplets, found:", iq_trip
      WRITE(*,*) "List of missing ones follow:"
      DO c = 1, nq(3)
      DO b = 1, nq(2)
      DO a = 1, nq(1)
        DO k = 1, nq(3)
        DO j = 1, nq(2)
        DO i = 1, nq(1)
          IF(found(i,j,k,a,b,c)<0) print*, i,j,k,a,b,c
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL errore('read_d3_matrices', "missing triplets!!", 1)
    ENDIF
    
  END SUBROUTINE
  !
  ! \/o\________\\\________________\\/\_________________________/^>
  SUBROUTINE bwfft_d3_interp(nq, nq_trip, nat, tau, at, bg, d3grid, fc3)
    USE constants,        ONLY : tpi
    USE fc3_interpolate,  ONLY : grid
    USE d3_basis,         ONLY : d3_6idx_2_3idx
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3), nq_trip, nat
    REAL(DP),INTENT(in) :: tau(3,nat), at(3,3), bg(3,3)
    TYPE(d3_list),INTENT(in) ::  d3grid(nq_trip)
    TYPE(grid),INTENT(inout) :: fc3
    !
    INTEGER :: i,j,k, iq, ifar, jfar, nxr, ixr, jxr, iat, jat, kat, irx
    REAL(DP),ALLOCATABLE :: xr(:,:), sc_xr(:,:)
    REAL(DP) :: d1(3), d2(3), d3(3), &
                xtau1(3), xtau2(3), xtau3(3)
    ! recentering of the Fcs 
    INTEGER,PARAMETER :: far = 2, nfar = (2*far+1)**3, nperix=8**2
    ! probably nperix = 64 is sufficient (i.e. both points in a WS cell corner, cubic)
    REAL(DP) :: peri_min, perix, perinorm
    REAL(DP),ALLOCATABLE :: farx_list(:,:,:), rx_list(:,:,:)
    INTEGER :: nxr_list, nxr_list_old, rx_idx(nperix), iperi, nperi
    REAL(DP),PARAMETER :: eps_peri = 1.d-5, eps_imag = 1.d-5
    ! fft aux variable
    REAL(DP) :: arg, norm
    COMPLEX(DP) :: phase, fc(3,3,3)
    COMPLEX(DP),ALLOCATABLE :: mat(:,:,:, :,:,:, :), matx(:,:,:, :,:,:, :), fcx(:,:,:)
    ! test variables
    REAL(DP) :: sum_imag, max_imag, avg_imag, max_imag_frac
    COMPLEX(DP) :: max_imag_mag
    INTEGER  :: count_imag
    !
    ! Generate a super-grid composed of nq(1) x nq(2) x nq(3) supercells of the unit cell
    ALLOCATE(sc_xr(3,nfar))
    ifar = 0
    DO k = -far,far
    DO j = -far,far
    DO i = -far,far
      ifar = ifar+1
      sc_xr(:,ifar) = at(:,1)*nq(1)*i + at(:,2)*nq(2)*j + at(:,3)*nq(3)*k
    ENDDO
    ENDDO
    ENDDO
    IF(ifar/=nfar) CALL errore("bwfft_d3_interp", "something wrong nfar", 1)
    !
    ! Compose the base of cells inside the nq(1) x nq(2) x nq(3) supercell
    nxr = nq(1)*nq(2)*nq(3)
    ALLOCATE(xr(3,nxr))
    ixr = 0 
    DO k = 0, nq(3)-1
    DO j = 0, nq(2)-1
    DO i = 0, nq(1)-1
      ixr = ixr+1
      xr(:,ixr) = at(:,1)*i + at(:,2)*j + at(:,3)*k
    ENDDO
    ENDDO
    ENDDO
    IF(ixr/=nxr) CALL errore("bwfft_d3_interp", "something wrong nxr", 1)
    ! xr + sc_xr is a super-cell equivalent copy of xr, for any sc_xr
    !
    ! For every xr2 and xr3 picked from xr we form the triangle 0-xr2-xr3 and 
    ! for any of the super-cell equivalent copies we take the one with the 
    ! shortest perimeter. 
    ALLOCATE(farx_list(3,2,nperix)) 
    nxr_list = 0
    norm = 1._dp/DBLE(nq_trip)
    !
    DO kat = 1,nat
    DO jat = 1,nat
    DO iat = 1,nat
      !
      xtau1 = tau(:,iat)
      xtau2 = tau(:,jat)
      xtau3 = tau(:,kat)
      !
      DO jxr = 1, nxr
      DO ixr = 1, nxr
        !
        fc = 0._dp
        DO iq = 1, nq_trip
          arg = tpi * ( SUM( xr(:,ixr)*d3grid(iq)%xq2) + SUM(xr(:,jxr)*d3grid(iq)%xq3) )
          phase = CMPLX( Cos(arg), Sin(arg), kind=DP)
          fc = fc + phase*norm* d3grid(iq)%d(:,:,:,iat,jat,kat)
          !write(*,'(3(2f10.4,3x))') d3grid(iq)%xq2,d3grid(iq)%xq3
        ENDDO
        !if(iat/=jat) then
          write(999,*)
          write(999,'(99i6)') iat,jat,kat,ixr,jxr
          write(999,'(3(2f10.4,3x))') fc
        !  stop 10
        !endif
        !stop 10
        !
        ! Look among the super-cell equivalent tripets of R points for the one(s)
        ! with the shortest perimeter
        nperi = 0
        peri_min = 0._dp ! it is set on first loop
        DO jfar = 1,nfar
        DO ifar = 1,nfar
          !
          d1 =                           xtau1 - (xtau2+xr(:,ixr)+sc_xr(:,ifar))
          d2 = (xtau2+xr(:,ixr)+sc_xr(:,ifar)) - (xtau3+xr(:,jxr)+sc_xr(:,jfar))
          d3 = (xtau3+xr(:,jxr)+sc_xr(:,jfar)) -  xtau1 
          !perix = SUM(ABS(d1)) + SUM(ABS(d2)) + SUM(ABS(d3))
          perix = DSQRT(SUM(d1**2)) + DSQRT(SUM(d2**2)) + DSQRT(SUM(d3**2))
          !perix = DSQRT(NORM2(d1) + NORM2(d2) + NORM2(d3))
          !
          IF (perix < peri_min-eps_peri .or. nperi==0 ) THEN
            nperi = 1
            farx_list = 0._dp
            farx_list(:,1,nperi) = xr(:,ixr)+sc_xr(:,ifar)
            farx_list(:,2,nperi) = xr(:,jxr)+sc_xr(:,jfar)
            peri_min = perix
          ELSE IF ( ABS(perix-peri_min) <= eps_peri ) THEN
            nperi = nperi + 1
            IF(nperi > nperix) CALL errore("bwfft_d3_interp", "nperix is too small, please report to developers", 1)
            farx_list(:,1,nperi) = xr(:,ixr)+sc_xr(:,ifar)
            farx_list(:,2,nperi) = xr(:,jxr)+sc_xr(:,jfar)
            !peri_min = perix
          ENDIF
          !
          !
        ENDDO
        ENDDO
        !
        ! Add the 2*nperi vectors from farx list to rx_list, if they are not already in the list
        ! in any case, return the indeces of the vectors in the list
        nxr_list_old = nxr_list
        CALL update_rlist(nperi, farx_list, nxr_list, rx_list, rx_idx)
        !IF(nperi>1) WRITE(*,'(3i3,3x,2i3,3x,i6,999i8)') iat,jat,kat, ixr,jxr, nperi, rx_idx(1:nperi)
        !
        ! Make room for new force constants, if necessary
        IF(nxr_list>nxr_list_old)THEN
          !PRINT*, nxr_list, nxr_list_old
          IF(allocated(mat)) THEN
            ALLOCATE(matx(3,3,3, nat,nat,nat, nxr_list_old))
            matx = mat
            DEALLOCATE(mat)
            ALLOCATE(mat(3,3,3, nat,nat,nat, nxr_list))
            mat(:,:,:, :,:,:, nxr_list_old+1:nxr_list) = 0._dp
            mat(:,:,:, :,:,:, 1:nxr_list_old) = matx(:,:,:, :,:,:, 1:nxr_list_old)
            DEALLOCATE(matx)
          ELSE
            ALLOCATE(mat(3,3,3, nat,nat,nat, nxr_list))
            mat = 0._dp
          ENDIF
        ENDIF
        !
        IF(nperi>0)THEN
          perinorm = 1._dp/REAL(nperi,kind=DP)
          !print*, perinorm * fc(:,:,:)
          DO iperi = 1, nperi
            !
            mat(:,:,:, iat,jat,kat, rx_idx(iperi)) = perinorm * fc(:,:,:)
            !PRINT*, SUM(ABS(mat(:,:,:, iat,jat,kat, rx_idx(iperi)))), iat,jat,kat, rx_idx(iperi)
          ENDDO
        ELSE
          CALL errore("bwfft_d3_interp", "found no perimeters", 1)
        ENDIF
        !
      ENDDO
      ENDDO
      !
    ENDDO
    ENDDO
    ENDDO
    !
    !WRITE(998,*) mat
    !
    sum_imag = 0._dp
    max_imag = 0._dp
    max_imag_frac = 0._dp
    max_imag_mag = 0._dp
    count_imag = 0
    !
    fc3%n_R = nxr_list
    ALLOCATE(fcx(3*nat,3*nat,3*nat))
    ALLOCATE(fc3%fc(3*nat,3*nat,3*nat, nxr_list))
    ALLOCATE(fc3%ifc(3*nat,3*nat,3*nat, nxr_list))
    ALLOCATE(fc3%xR2(3,nxr_list), fc3%xR3(3,nxr_list))
    DO irx = 1, nxr_list
      !
      ! Transfer the matrices to force constant type
      CALL d3_6idx_2_3idx(nat, mat(:,:,:, :,:,:, irx), fcx)
      fc3%fc(:,:,:,irx)  = DBLE(fcx)
      fc3%ifc(:,:,:,irx) = DIMAG(fcx)
      fc3%xR2(:,irx) = rx_list(:,1,irx)
      fc3%xR3(:,irx) = rx_list(:,2,irx)
      !
      !
      ! Run some tests of imaginary parts:
      sum_imag = sum_imag + SUM(ABS(DIMAG(fcx)))
      max_imag = MAX(max_imag, MAXVAL(ABS(DIMAG(fcx))))
      DO k = 1, 3*nat
      DO j = 1, 3*nat
      DO i = 1, 3*nat
        IF(ABS(REAL(fcx(i,j,k)))>eps_imag)THEN
          count_imag = count_imag + 1
          IF(ABS(DIMAG(fcx(i,j,k))/REAL(fcx(i,j,k))) > max_imag_frac) THEN
            max_imag_frac = ABS(DIMAG(fcx(i,j,k))/REAL(fcx(i,j,k)))
            max_imag_mag  = fcx(i,j,k)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      !
    ENDDO
    DEALLOCATE(fcx,mat)
    avg_imag = sum_imag / count_imag
    !
    WRITE(*,'(4x,a,es10.3)') "Sum of imaginary parts:    ", sum_imag
    WRITE(*,'(4x,a,es10.3)') "Average imaginary part:    ", avg_imag
    WRITE(*,'(4x,a,es10.3,a,2es10.3)') "Largest imaginary fraction:", max_imag_frac, " for ", max_imag_mag
    WRITE(*,'(4x,a,es10.3)') "Largest imaginary part:    ", max_imag
    !
    WHERE(ABS(fc3%fc) < 1.d-30)  fc3%fc = 0._dp
    WHERE(ABS(fc3%ifc) < 1.d-30) fc3%ifc = 0._dp
    !
    !
    ALLOCATE(fc3%yR2(3,nxr_list), fc3%yR3(3,nxr_list))
    CALL cryst_to_cart(2*nxr_list, rx_list, bg, -1)
    DO irx = 1, nxr_list
      fc3%yR2(:,irx) = NINT(rx_list(:,1,irx))
      fc3%yR3(:,irx) = NINT(rx_list(:,2,irx))
    ENDDO
    DEALLOCATE(rx_list, farx_list)
    !
  END SUBROUTINE bwfft_d3_interp
  !
  ! \/o\________\\\________________\\/\_________________________/^>
  SUBROUTINE update_rlist(nperi, farx_list, nxr_list, rx_list, rx_idx)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nperi
    REAL(DP),INTENT(in) :: farx_list(3,2,nperi)
    INTEGER, INTENT(inout) :: nxr_list
    REAL(DP),INTENT(inout),ALLOCATABLE :: rx_list(:,:,:) !(3,2,nxr_list)
    INTEGER,INTENT(inout) :: rx_idx(nperi)         !(nxr_list)
    !
    INTEGER :: iperi, irx, nxr_list_new, nxr_list_out
    REAL(DP),ALLOCATABLE :: rx_aux(:,:,:) !(3,2,nxr_list)
    REAL(DP),PARAMETER :: eps_rx = 1.d-5
    !
    rx_idx = -1
    !
    IF(nxr_list>0)THEN
      !
      DO iperi = 1, nperi
        DO irx = 1, nxr_list
          !
          IF( ALL(ABS(farx_list(:,:,iperi)-rx_list(:,:,irx))<eps_rx) ) THEN
            rx_idx(iperi) = irx 
            EXIT
          ENDIF
          !
        ENDDO
      ENDDO
      !
      ! If I found all the couple of Rs, return here
      IF (ALL(rx_idx>0)) RETURN
      !
    ENDIF
    !
    ! Add enough space for the new vectors
    nxr_list_new = COUNT(rx_idx < 0)
    !print*, "adding ", nxr_list_new
    IF(nxr_list>0)THEN
      ALLOCATE(rx_aux(3,2,nxr_list))
      rx_aux(:,:,1:nxr_list) = rx_list(:,:,1:nxr_list)
      DEALLOCATE(rx_list)
      !
      nxr_list_out = nxr_list+nxr_list_new
      ALLOCATE(rx_list(3,2,nxr_list_out))
      rx_list(:,:,1:nxr_list) = rx_aux(:,:,1:nxr_list)
      DEALLOCATE(rx_aux)
    ELSE
      nxr_list_out = nxr_list_new
      ALLOCATE(rx_list(3,2,nxr_list_new))
    ENDIF
    !
    ! Set the new indeces
    DO iperi = 1, nperi
      IF(rx_idx(iperi)<0)THEN
        nxr_list = nxr_list+1
        rx_list(:,:,nxr_list) = farx_list(:,:,iperi)
        rx_idx(iperi) = nxr_list
      ENDIF
    ENDDO
    IF(nxr_list /= nxr_list_out) CALL errore("update_rlist", "something wrong",1)
    IF (ANY(rx_idx<0)) CALL errore('update_rlist', 'something wrong', 2)
    
  END SUBROUTINE update_rlist
  ! \/o\________\\\________________\\/\_________________________/^>
  SUBROUTINE test_fwfft_d3(nq_trip, S, d3grid, fc3)
    USE input_fc,         ONLY : ph_system_info
    USE d3_basis,         ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx
    USE fc3_interpolate,  ONLY : grid
    IMPLICIT NONE
    INTEGER :: nq_trip
    TYPE(d3_list),INTENT(in) :: d3grid(nq_trip)
    TYPE(grid),INTENT(in)    :: fc3
    TYPE(ph_system_info),INTENT(in) :: S
    
    COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
    REAL(DP) :: rmaxi, imaxi
    INTEGER :: iq
    LOGICAL :: found
    
    found = .false.
    ALLOCATE(D3(S%nat3,S%nat3,S%nat3))
    ALLOCATE(P3(S%nat3,S%nat3,S%nat3))
    DO iq = 1, nq_trip
      !CALL fc3%interpolate(d3grid(iq)%xq2_true, d3grid(iq)%xq3_true, S%nat3, D3)
      CALL fc3%interpolate(d3grid(iq)%xq2, d3grid(iq)%xq3, S%nat3, D3)
      CALL d3_6idx_2_3idx(S%nat, d3grid(iq)%d, P3)
      rmaxi = MAXVAL(ABS(DBLE(D3)-DBLE(P3)))
      imaxi = MAXVAL(ABS(DIMAG(D3)-DIMAG(P3)))
      IF(rmaxi > 1.d-6) THEN
        found = .true.
        WRITE(*,'("    Real:",i8,3(3f12.4,3x),es20.8)') iq, d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3, rmaxi
      ENDIF
      IF(imaxi > 1.d-6) THEN
        found = .true.
        WRITE(*,'("    Imag:",i8,3(3f12.4,3x),es20.8)') iq, d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3, imaxi
      ENDIF
    ENDDO
    !
    IF(found) THEN
      WRITE(*,*) "  WARNING! BW/FW FFT discrepancy was found! "
    ELSE
      WRITE(*,*) "  BW/FW FFT fine."
    ENDIF

  END SUBROUTINE test_fwfft_d3

END MODULE program_qq2rr


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM qq2rr
  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : ph_system_info
  USE iso_c_binding,   ONLY : c_int
  USE fc3_interpolate, ONLY : grid
  USE program_qq2rr
  
  IMPLICIT NONE
  TYPE(d3_list),ALLOCATABLE ::  d3grid(:)
  TYPE(grid)                :: fc3
  INTEGER :: nq(3)
  INTEGER :: nq_trip, nq_grid
  INTEGER(kind=c_int)    :: kb
  !
  TYPE(ph_system_info) :: S
  COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
  !
  INTEGER :: count, i, j
  CHARACTER(len=512) :: argv

  count = command_argument_count()
  IF(count>=3)THEN
    DO i = 1,3
      CALL get_command_argument(i,argv)
      READ(argv,*) nq(i)
    ENDDO
  ELSE
    nq = (/ 0,0,0 /)
  ENDIF
  WRITE(*,*) "Reading grid", nq
  
  nq_grid = nq(1)*nq(2)*nq(3)
  nq_trip = nq_grid**2
  ALLOCATE(d3grid(nq_trip))
  !
  WRITE(*,*) "Reading D3 matrices..."
  CALL read_d3_matrices(nq, nq_trip, S, d3grid)
  WRITE(*,*) "Reading D3 matrices done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  WRITE(*,*) "Doing Backward FFT..."
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fc3)
  fc3%nq = nq
  WRITE(*,*) "Backward FFT done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  WRITE(*,*) "Writing FCs to file..."
  CALL fc3%write("mat3R.xxx", S)
  !
  WRITE(*,*) "Testing Forward FFT, with imaginary part..."
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing Forward FFT, without imaginary part..."
  DEALLOCATE(fc3%ifc)
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing forward FFT done"
  !
END PROGRAM qq2rr


