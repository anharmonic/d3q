!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Recentering of the two bodies force constants based on a previous
! implementation by G. Fugallo
!
MODULE f3_bwfft
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
    USE d3matrix_io2,ONLY : read_d3dyn_xml2
    USE d3_shuffle !,  ONLY : nperms, d3perms_order2, d3_shuffle_global, d3_shuffle_equiv
    USE d3_basis,    ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx
    USE parameters,  ONLY : ntypx
    USE input_fc,    ONLY : ph_system_info, same_system, aux_system
    !USE write_d3dyn_ascii, ONLY : write_d3dyn_XXX, zero_d3dyn_XXX
    USE functions,   ONLY : refold_bz
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3), nq_trip
    TYPE(ph_system_info),INTENT(out) :: S
    TYPE(d3_list),INTENT(inout) ::  d3grid(nq_trip)

    REAL(DP) :: xq(3,3)
    COMPLEX(DP),ALLOCATABLE   :: d3(:,:,:, :,:,:), &
                                 d3_shuffled(:,:,:,:,:,:), &
                                 p3(:,:,:)
    INTEGER :: i,j,k, a,b,c, ios, ntyp, iperm, iq_trip, &
               nxr_list, nxr_list_old, iq_aux
    LOGICAL :: first, skip
    INTEGER,ALLOCATABLE :: found(:,:,:,:,:,:)
    INTEGER :: countq(nq_trip)
    CHARACTER(len=512) :: filename
    CHARACTER(len=5) :: format_version
    INTEGER :: iqa(3), iqb(3), iqc(3)
    REAL(DP) :: xq_shift(3,3), thresh
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
      WRITE(*,*) "Reading '", TRIM(filename),"'..."

      IF(first)THEN
        first=.false.
        CALL read_d3dyn_xml2(filename, xq(:,1), xq(:,2), xq(:,3), d3=d3, &
                            nat=S%nat, atm=S%atm, ntyp=S%ntyp, &
                            ityp=S%ityp, ibrav=S%ibrav, celldm=S%celldm, at=S%at,&
                            amass=S%amass, tau=S%tau, seek=.false.,&
                            file_format_version=format_version)
        CALL aux_system(S)
        CALL latgen( S%ibrav, S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega )
        S%at=S%at/S%celldm(1)
        CALL recips(S%at(:,1), S%at(:,2), S%at(:,3), S%bg(:,1), S%bg(:,2), S%bg(:,3))
        ALLOCATE(p3(S%nat3, S%nat3, S%nat3))
        ALLOCATE(d3_shuffled(3,3,3, S%nat, S%nat, S%nat))
        !
      ELSE
        CALL read_d3dyn_xml2(filename, xq(:,1), xq(:,2), xq(:,3), d3=d3, &
                            nat=Sx%nat, atm=Sx%atm, ntyp=Sx%ntyp, &
                            ityp=Sx%ityp, ibrav=Sx%ibrav, celldm=Sx%celldm, at=Sx%at,&
                            amass=Sx%amass,tau=Sx%tau, seek=.false.,&
                            file_format_version=format_version)
        CALL latgen( Sx%ibrav, Sx%celldm, Sx%at(:,1), Sx%at(:,2), Sx%at(:,3), Sx%omega )
        Sx%at=Sx%at/Sx%celldm(1)
        CALL recips(Sx%at(:,1), Sx%at(:,2), Sx%at(:,3), Sx%bg(:,1), Sx%bg(:,2), Sx%bg(:,3))
        IF(.not.same_system(S,Sx)) CALL errore("qq2rr", "not same system "//TRIM(filename), 1)
      ENDIF
      !
      IF( format_version == "1.0.0") THEN
        WRITE(*,*) "Workaround for d3q =< 1.7b: taking complex conjugate of D3 matrix"
        d3 = DCONJG(d3)
      ELSEIF ( format_version /= "1.1.0")THEN
        CALL errore("read_d3_matrices","unknow file format version: "//TRIM(format_version),1)
      ENDIF
      !
      xq(:,1) = refold_bz(xq(:,1), S%bg)
      xq(:,2) = refold_bz(xq(:,2), S%bg)
      xq(:,3) = refold_bz(xq(:,3), S%bg)
      
      PERM_LOOP : &
      DO iperm = 1,nperms
        a = d3perms_order2(1,iperm)
        b = d3perms_order2(2,iperm)
        c = d3perms_order2(3,iperm)
        !
        ! Check if the points are in the grid (this allows one to select a sub-grid):
        iqb = check_in_grid(nq, xq(:,b), S%at, skip)
        IF(skip) CYCLE PERM_LOOP
        !
        iqc = check_in_grid(nq, xq(:,c), S%at, skip)
        IF(skip) CYCLE PERM_LOOP
        !
        ! Add this point to the grid, if an equivalent was not already added
        IF(found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3)) < 0) THEN
          iq_trip = iq_trip + 1
          found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3)) = iq_trip
          countq(iq_trip) = 1
          !
          WRITE(*,'(i6,3x,a,3(3f10.4,3x),2(3i3,3x),3i2,3x,a,i6  )') &
            iq_trip, "xq ", xq(:,a), xq(:,b), xq(:,c), iqb, iqc, &
            a,b,c, TRIM(filename),iperm
          !
          ! Repack the matrix, shuffle it and repack it again
          CALL d3_6idx_2_3idx(S%nat, d3, p3)
          CALL d3_shuffle_global(S%nat, 1,2,3, a,b,c, .false., p3)
          !
!           write(4999,*) iperm
!           write(4999,'(i6,3(3f8.4,3x))') iq_trip, xq(:,a), xq(:,b), xq(:,c)
!           write(4999,'(3(2f10.4,3x))') 1.d+4*p3
          !
          CALL d3_3idx_2_6idx(S%nat, p3, d3_shuffled)
          !
          ALLOCATE(d3grid(iq_trip)%d(3,3,3, S%nat,S%nat,S%nat))
          d3grid(iq_trip)%d   = d3_shuffled
          !
          d3grid(iq_trip)%xq1(:) = xq(:,a)
          d3grid(iq_trip)%xq2(:) = xq(:,b)
          d3grid(iq_trip)%xq3(:) = xq(:,c)
        ELSE !IF(.false.) THEN
          ! Average out different dyn mat files that can be transformed to the
          ! same by permutation, this has the effect to apply those symmetries
          ! which have the effect of changing the order of the q points in the 
          ! triplet (these symmetries are not used in the d3 code)
          iq_aux = found(iqb(1),iqb(2),iqb(3), iqc(1),iqc(2),iqc(3))
          countq(iq_aux) = countq(iq_aux) +1
          CALL d3_6idx_2_3idx(S%nat, d3, p3)
          CALL d3_shuffle_global(S%nat, 1,2,3, a,b,c, .false., p3)
          CALL d3_3idx_2_6idx(S%nat, p3, d3_shuffled)
          d3grid(iq_aux)%d   = d3grid(iq_aux)%d + d3_shuffled
          WRITE(*,'(i6,3x,a,3(3f10.4,3x),2(3i3,3x),3i2,3x,a,2i6  )') &
            iq_trip, "xq'", xq(:,a), xq(:,b), xq(:,c), iqb, iqc, &
            a,b,c, TRIM(filename),iperm, countq(iq_aux)
        ENDIF
      ENDDO &
      PERM_LOOP
      
    ENDDO
    
    DO iq_aux = 1,nq_trip
      IF(countq(iq_aux)>1)THEN
        d3grid(iq_aux)%d = d3grid(iq_aux)%d/countq(iq_aux)
      ENDIF
    ENDDO
    
    IF(first) CALL errore("read_d3_matrices","I found nothing to read",1)
    DEALLOCATE(p3)
    
    IF(ANY(found<0)) THEN
      WRITE(*,*) "Expecting:", nq_trip, "triplets, found:", iq_trip
      WRITE(*,*) "List of missing ones follows:"
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
  SUBROUTINE bwfft_d3_interp(nq, nq_trip, nat, tau, at, bg, d3grid, fc3, far)
    USE constants,        ONLY : tpi
    USE fc3_interpolate,  ONLY : grid
    USE d3_basis,         ONLY : d3_6idx_2_3idx
    USE functions,        ONLY : norm, cross
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3), nq_trip, nat
    REAL(DP),INTENT(in) :: tau(3,nat), at(3,3), bg(3,3)
    TYPE(d3_list),INTENT(in) ::  d3grid(nq_trip)
    TYPE(grid),INTENT(inout) :: fc3
    INTEGER,INTENT(in) :: far
    !
    INTEGER :: i,j,k, iq, ifar, jfar, nxr, ixr, jxr, iat, jat, kat, irx
    REAL(DP),ALLOCATABLE :: xr(:,:), sc_xr(:,:)
    REAL(DP) :: d1(3), d2(3), d3(3), p1(3), p2(3), p3(3), p0(3), &
                xtau1(3), xtau2(3), xtau3(3)
    ! recentering of the Fcs 
    INTEGER :: nfar
    INTEGER, PARAMETER :: nperix=512 !8**2
    ! probably nperix = 64 is sufficient (i.e. both points in a WS cell corner, cubic)
    REAL(DP) :: peri_min, perix, perinorm
    REAL(DP),ALLOCATABLE :: farx_list(:,:,:), rx_list(:,:,:)
    INTEGER :: nxr_list, nxr_list_old, rx_idx(nperix), iperi, nperi, PASS
    REAL(DP),PARAMETER :: eps_peri = 1.d-5, eps_imag = 1.d-5
    ! fft aux variable
    REAL(DP) :: arg, pref
    COMPLEX(DP) :: phase, fc(3,3,3)
    COMPLEX(DP),ALLOCATABLE :: mat(:,:,:, :,:,:, :), matx(:,:,:, :,:,:, :), fcx(:,:,:)
    ! test variables
    REAL(DP) :: sum_imag, max_imag, avg_imag, max_imag_frac
    COMPLEX(DP) :: max_imag_mag
    INTEGER  :: count_imag
    !
    nfar = (2*far+1)**3
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
    pref = 1._dp/DBLE(nq_trip)
    !
    TWO_PASS : &
    DO PASS = 1,2
      IF(PASS==2)THEN
        WRITE(*,*) "Found", nxr_list, "possible couples of R2,R3"
        nxr_list_old = nxr_list ! keep track of this to be sure it does not change anymore
        ALLOCATE(mat(3,3,3, nat,nat,nat, nxr_list))
      ENDIF
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
          ! On first pass, only count the perimeters
          IF (PASS==2) THEN
            fc = (0._dp, 0._dp)
            DO iq = 1, nq_trip
              arg = tpi * ( SUM( xr(:,ixr)*d3grid(iq)%xq2) + SUM(xr(:,jxr)*d3grid(iq)%xq3) )
              phase = CMPLX( Cos(arg), Sin(arg), kind=DP)
              fc = fc + phase*pref* d3grid(iq)%d(:,:,:,iat,jat,kat)
            ENDDO
          ENDIF
          !
          ! Look among the super-cell equivalent triplets of R points for the one(s)
          ! with the shortest perimeter
          nperi = 0
          peri_min = 0._dp ! it is set on first loop
          DO jfar = 1,nfar
          DO ifar = 1,nfar
            !
            p1 = xtau1
            p2 = xtau2+xr(:,ixr)+sc_xr(:,ifar)
            p3 = xtau3+xr(:,jxr)+sc_xr(:,jfar)
            d1 = p1 - p2
            d2 = p2 - p3
            d3 = p3 - p1
            ! Perimeter of the triangle
            perix = DSQRT(SUM(d1**2)) + DSQRT(SUM(d2**2)) + DSQRT(SUM(d3**2))
            !
!            ! Radius of the smallest sphere containing the three atoms
!            ! (or radius of the circle circumscribing the triangle)
!            perix = 1.d+6
!            ! If the triangle is acute, then the smallest sphere has the circumcircle at its equator,
!            ! this is the circumcircle radius according to Wikipedia:
!            IF(norm(cross(d1,d2)) > 0._dp) perix = MIN(perix,0.5_dp * norm(d1)*norm(d2)*norm(d3) / norm(cross(d1,d2)))
!            ! In obtuse triangles, the circumcenter falls outside the triangle, I can find a smaller
!            ! sphere which has its center in the mid point of the longest side, two points will be on the surface
!            ! of the sphere, one inside.
!            ! I check this by checking if the third point is inside the sphere with the center in the mid point of
!            ! the other two and passing through them.
!            ! These three condition also work in degenerate case, i.e. three points aligned, or some points coincide
!            IF(norm((p1+p2)/2 - p3) <= norm(d1)/2) perix=MIN(perix,norm(d1)/2)
!            IF(norm((p1+p3)/2 - p2) <= norm(d3)/2) perix=MIN(perix,norm(d3)/2)
!            IF(norm((p2+p3)/2 - p1) <= norm(d2)/2) perix=MIN(perix,norm(d2)/2)
!            ! Note that in the case of a right triangle, the circumcenter is the middle point of the
!            ! hypotenuse, the two methods give the same value!
            !
!             p0 = (p1+p2+p3)/3._dp
!             d1 = p0 - p1
!             d2 = p0 - p2
!             d3 = p0 - p3
!             perix = SUM((/ SUM(d1**2), SUM(d2**2), SUM(d3**2) /))
!             perix = DSQRT(SUM(d1**2)) + DSQRT(SUM(d2**2)) + DSQRT(SUM(d3**2))
            !
            IF (perix < peri_min-eps_peri .or. nperi==0 ) THEN
              nperi = 1
              farx_list = 0._dp
              farx_list(:,1,nperi) = xr(:,ixr)+sc_xr(:,ifar)
              farx_list(:,2,nperi) = xr(:,jxr)+sc_xr(:,jfar)
              peri_min = perix
            ELSE IF ( ABS(perix-peri_min) <= eps_peri ) THEN
              nperi = nperi + 1
              IF(nperi > nperix) CALL errore("bwfft_d3_interp", "nperix is too small", 1)
              farx_list(:,1,nperi) = xr(:,ixr)+sc_xr(:,ifar)
              farx_list(:,2,nperi) = xr(:,jxr)+sc_xr(:,jfar)
              peri_min = (peri_min*(nperi-1)+perix)/DBLE(nperi)
            ENDIF
            !
          ENDDO
          ENDDO
          !
          IF (nperi==0) CALL errore("bwfft_d3_interp", "found no perimeters", 1)
!           IF (nperi>1 .and. PASS==1)THEN
!             WRITE(*,*)
!             d1 = xr(:,ixr)
!             CALL cryst_to_cart(1, d1, bg, -1)
!             d2 = xr(:,jxr)
!             CALL cryst_to_cart(1, d2, bg, -1)
!             WRITE(*,'(2(3i3,3x))') NINT(d1), NINT(d2)
!             DO iperi = 1, nperi
!               d1 = farx_list(:,1,iperi)
!               CALL cryst_to_cart(1, d1, bg, -1)
!               d2 = farx_list(:,2,iperi)
!               CALL cryst_to_cart(1, d2, bg, -1)
!               WRITE(*,'(2(3i3,3x))') NINT(d1), NINT(d2)
!             ENDDO
!           ENDIF
          ! Add the 2*nperi vectors from farx list to rx_list, if they are 
          ! not already in the list in any case, return the indexes 
          ! of the vectors in the list
          CALL update_rlist(nperi, farx_list, nxr_list, rx_list, rx_idx)
          !
          IF(PASS==2)THEN
            IF(nxr_list > nxr_list_old) CALL errore("bwfft_d3_interp", &
                                              "unexpected new triangle", 1)
            perinorm = 1._dp/REAL(nperi,kind=DP)
            !print*, perinorm * fc(:,:,:)
            DO iperi = 1, nperi
              !
              mat(:,:,:, iat,jat,kat, rx_idx(iperi)) = perinorm * fc(:,:,:)
            ENDDO
            !
          ENDIF
          !
        ENDDO
        ENDDO
        !
      ENDDO
      ENDDO
      ENDDO
    ENDDO TWO_PASS
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
        IF(ABS(DBLE(fcx(i,j,k)))>eps_imag)THEN
          count_imag = count_imag + 1
          IF(ABS(DIMAG(fcx(i,j,k))/DBLE(fcx(i,j,k))) > max_imag_frac) THEN
            max_imag_frac = ABS(DIMAG(fcx(i,j,k))/DBLE(fcx(i,j,k)))
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
    WRITE(*,'(4x,a,es10.3,a,2(es10.3,2x))') "Largest imaginary fraction:", max_imag_frac, " for ", max_imag_mag
    WRITE(*,'(4x,a,es10.3)') "Largest imaginary part:    ", max_imag
    !
!     WHERE(ABS(fc3%fc) < 1.d-30)  fc3%fc = 0._dp
!     WHERE(ABS(fc3%ifc) < 1.d-30) fc3%ifc = 0._dp
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
    INTEGER,INTENT(out) :: rx_idx(nperi)
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
  SUBROUTINE test_fwfft_d3(nq_trip, S, d3grid, fc3, &
                diff_stop, write_diff, write_diff_prefix)
    USE input_fc,         ONLY : ph_system_info
    USE d3_basis,         ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx
    USE fc3_interpolate,  ONLY : grid
    USE d3matrix_io2,     ONLY : write_d3dyn_xml2
    IMPLICIT NONE
    INTEGER :: nq_trip
    TYPE(d3_list),INTENT(in) :: d3grid(nq_trip)
    TYPE(grid),INTENT(in)    :: fc3
    TYPE(ph_system_info),INTENT(in) :: S
    LOGICAL,INTENT(in) :: diff_stop, write_diff
    CHARACTER(len=*),INTENT(in) :: write_diff_prefix
    
    COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:), D3_6idx(:,:,:,:,:,:)
    REAL(DP) :: rmaxi, imaxi
    INTEGER :: iq
    LOGICAL :: found
    
    found = .false.
    ALLOCATE(D3(S%nat3,S%nat3,S%nat3))
    ALLOCATE(D3_6idx(3,3,3, S%nat,S%nat,S%nat))
    ALLOCATE(P3(S%nat3,S%nat3,S%nat3))
    DO iq = 1, nq_trip
      !CALL fc3%interpolate(d3grid(iq)%xq2_true, d3grid(iq)%xq3_true, S%nat3, D3)
      CALL fc3%interpolate(d3grid(iq)%xq2, d3grid(iq)%xq3, S%nat3, D3)
      CALL d3_6idx_2_3idx(S%nat, d3grid(iq)%d, P3)
      rmaxi = MAXVAL(ABS(DBLE(D3)-DBLE(P3)))
      imaxi = MAXVAL(ABS(DIMAG(D3)-DIMAG(P3)))
      IF(rmaxi > 1.d-5) THEN
        found = .true.
        WRITE(*,'("    Real:",i8,3(3f12.4,3x),es20.8)') iq, d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3, rmaxi
      ENDIF
      IF(imaxi > 1.d-5) THEN
        found = .true.
        WRITE(*,'("    Imag:",i8,3(3f12.4,3x),es20.8)') iq, d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3, imaxi
      ENDIF
      IF(write_diff .and. (rmaxi>1.d-5 .or. imaxi>1.d-5) )THEN
        CALL d3_3idx_2_6idx(S%nat, D3, D3_6idx)
        CALL write_d3dyn_xml2(trim(write_diff_prefix), &
                            d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3,&
                            D3_6idx, S%ntyp, S%nat, S%ibrav, S%celldm, S%at, S%ityp, &
                            S%tau, S%atm, S%amass)
      ENDIF
    ENDDO
    !
    IF(found) THEN
      IF(diff_stop)THEN
        WRITE(*,*) "  ERROR!! BW/FW FFT discrepancy was found! "
        WRITE(*,*) "  This indicates a serious problem with D3 matrices "
        STOP 1
      ELSE
        WRITE(*,*) "  BEWARE! Possibly harmless BW/FW FFT discrepancy was found! "
        WRITE(*,*) "  This can be caused by low cutoff energy and/or translational symmetry"
        WRITE(*,*) "  violation by GGA/PBE functionals. "
      ENDIF
    ELSE
      WRITE(*,*) "  BW/FW FFT fine."
    ENDIF
    !
    DEALLOCATE(D3,P3,D3_6idx)

  END SUBROUTINE test_fwfft_d3
  ! \/o\________\\\________________\\/\_________________________/^>
  ! This subroutine does a forward FFt from Force constants to a specified 
  ! grid of \vec(nq)^2 triplets
  !
  SUBROUTINE regen_fwfft_d3(nq, nq_trip, S, d3grid, fc3, writed3)
    USE input_fc,         ONLY : ph_system_info
    USE d3_basis,         ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx
    USE fc3_interpolate,  ONLY : grid
    USE d3matrix_io2,     ONLY : write_d3dyn_xml2
    USE functions,        ONLY : refold_bz
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nq(3)
    INTEGER,INTENT(in) :: nq_trip
    TYPE(d3_list),INTENT(out) :: d3grid(nq_trip)
    TYPE(grid),INTENT(in)    :: fc3
    TYPE(ph_system_info),INTENT(in) :: S
    LOGICAL,INTENT(in) :: writed3
    
    COMPLEX(DP),ALLOCATABLE :: D3(:,:,:)
    INTEGER :: iq, i1,j1,k1, i2,j2,k2
    REAL(DP) :: bgi(3,3)
    
    bgi(:,1) = S%bg(:,1)/DBLE(nq(1))
    bgi(:,2) = S%bg(:,2)/DBLE(nq(2))
    bgi(:,3) = S%bg(:,3)/DBLE(nq(3))
    
    iq = 0
    DO i1 = 0,nq(1)-1
    DO j1 = 0,nq(2)-1
    DO k1 = 0,nq(3)-1
      DO i2 = 0,nq(1)-1
      DO j2 = 0,nq(2)-1
      DO k2 = 0,nq(3)-1
        iq = iq+1
        IF(iq>nq_trip) CALL errore("regen_fwfft_d3","wrong nq_trip",1)
        ALLOCATE(d3grid(iq)%d(3,3,3, S%nat,S%nat,S%nat))
        d3grid(iq)%xq1 = i1*bgi(:,1) + j1*bgi(:,2) + k1*bgi(:,3)
        d3grid(iq)%xq2 = i2*bgi(:,1) + j2*bgi(:,2) + k2*bgi(:,3)
        d3grid(iq)%xq3 = -d3grid(iq)%xq1-d3grid(iq)%xq2
        ! refold in the Brillouin zone
        d3grid(iq)%xq1 = refold_bz(d3grid(iq)%xq1, S%bg)
        d3grid(iq)%xq2 = refold_bz(d3grid(iq)%xq2, S%bg)
        d3grid(iq)%xq3 = refold_bz(d3grid(iq)%xq3, S%bg)
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO

    ALLOCATE(D3(S%nat3,S%nat3,S%nat3))
    DO iq = 1, nq_trip
      !CALL fc3%interpolate(d3grid(iq)%xq2_true, d3grid(iq)%xq3_true, S%nat3, D3)
      CALL fc3%interpolate(d3grid(iq)%xq2, d3grid(iq)%xq3, S%nat3, D3)
      CALL d3_3idx_2_6idx(S%nat, D3, d3grid(iq)%d)
      IF(writed3) THEN
        CALL write_d3dyn_xml2("atmp", d3grid(iq)%xq1, d3grid(iq)%xq2, d3grid(iq)%xq3, &
                              d3grid(iq)%d, S%ntyp, S%nat, S%ibrav, S%celldm, S%at,  &
                              S%ityp, S%tau, S%atm, S%amass)
      ENDIF
    ENDDO

  END SUBROUTINE regen_fwfft_d3
END MODULE f3_bwfft

