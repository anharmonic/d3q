!
!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Uses PW/wsweights.f90 from Quantum-ESPRESSO
MODULE quter_module

  USE kinds, ONLY : DP

  PUBLIC  :: quter
  PRIVATE :: R_list_idx, expand_matR

  CONTAINS

  SUBROUTINE quter(nq1, nq2, nq3, nat,tau,at,bg, matq, gridq, fc, far)
    USE input_fc,  ONLY : forceconst2_grid, allocate_fc2_grid, write_fc2
    USE constants, ONLY : tpi
    USE functions, ONLY : default_if_not_present
    IMPLICIT NONE
    ! Dummy arguments
    INTEGER,INTENT(in)     :: nq1, nq2, nq3 ! dimensions of the q-points grid
    INTEGER,INTENT(in)     :: nat ! number of atoms
    REAL(DP),INTENT(in)    :: tau(3,nat) ! atom positions (alat units)
    REAL(DP),INTENT(in)    :: at(3,3), bg(3,3) ! real, reciprocal lattice
    COMPLEX(DP),INTENT(in) :: matq(3,3,nat,nat,nq1*nq2*nq3)
    REAL(DP),INTENT(in)    :: gridq(3,nq1*nq2*nq3)
    TYPE(forceconst2_grid),INTENT(out) :: fc
    INTEGER,OPTIONAL,INTENT(in) :: far
    !
    REAL(DP),PARAMETER :: eps = 1.d-8
    INTEGER :: i, na1, na2, j1,j2, jn1,jn2,  iiq, iR, l1,l2,l3, far_
    INTEGER :: nRbig, nqt
    REAL(DP),ALLOCATABLE :: Rbig(:,:)
    REAL(DP) :: aus, arg, dist(3), totalweight, Rx(3)
    REAL(DP) :: imag_tot, imag_max
    !
    ! For the list of used R vectors:
    REAL(DP),POINTER :: Rout(:,:) => null()
    INTEGER :: iout, nRout
    !
    ! Stuff used to compute Wigner-Seitz weights:
    INTEGER, PARAMETER:: nrwsx=5000
    INTEGER :: nrws
    REAL(DP) :: atws(3,3) ! supercell of size nq1 x nq2 x nq3
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    !
    COMPLEX(DP),POINTER :: matR(:,:,:,:,:) => null()
    !
    far_ = default_if_not_present(2, far)

    nqt = nq1*nq2*nq3
    atws(:,1) = nq1*at(:,1)
    atws(:,2) = nq2*at(:,2)
    atws(:,3) = nq3*at(:,3)
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,atws)
    !
    ! Construct a big enough lattice of Supercell vectors
    IF(far_>0)THEN
        nRbig = (2*far_*nq1+1)*(2*far_*nq2+1)*(2*far_*nq3+1)
        ALLOCATE(Rbig(3,nRbig))
        nRbig=0
        DO l1=-far_*nq1,far_*nq1
        DO l2=-far_*nq2,far_*nq2
        DO l3=-far_*nq3,far_*nq3
        nRbig=nRbig+1
        Rbig(:, nRbig) = at(:,1)*l1 +at(:,2)*l2 +at(:,3)*l3
        END DO
        END DO
        END DO
        IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)  
 !       WRITE(*,*) "seeking over ", nRbig," vectors"
    ELSE
        nRbig = nq1*nq2*nq3
        ALLOCATE(Rbig(3,nRbig))
        nRbig=0
        DO l1=0,nq1-1
        DO l2=0,nq2-1
        DO l3=0,nq3-1
        nRbig=nRbig+1
        Rbig(:, nRbig) = at(:,1)*l1 +at(:,2)*l2 +at(:,3)*l3
        END DO
        END DO
        END DO
        IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)  
!        WRITE(*,*) "nfar==0 => standard FT over ", nRbig," vectors"
    ENDIF
    !
    ! dyn.mat. FFT
    !
    nRout=0
    DO na1=1,nat
    DO na2=1,nat
      DO j1=1,3
      DO j2=1,3
        totalweight = 0._dp
        DO iR= 1,nRbig
          !
          aus= 0._dp
          IF(far_>0)THEN
            dist(:) = Rbig(:,iR) + tau(:,na1) - tau(:,na2)
            wg = wsweight(dist,rws,nrws)
            wg = wg/DFLOAT(nqt)
          ELSE
            wg = 1._dp/DFLOAT(nqt)
          ENDIF
          IF(wg /= 0) THEN 
            !
            DO iiq=1,nqt
              !   
              arg=tpi*SUM(gridq(:,iiq)*Rbig(:,iR)) 
              aus = aus + CMPLX(cos(arg),sin(arg),kind=DP)*matq(j1,j2,na1,na2,iiq)
              !
            END DO
            !
            Rx = Rbig(:,iR)
            CALL cryst_to_cart(1,Rx,bg,-1)
            ! find if we have already used this R point ...
            iout = R_list_idx(nRout,Rout,Rx)
            ! ... if not, increase storage
            CALL expand_matR(iout,nRout,nat,matR)
            matR(j1,j2,na1,na2,iout) = aus * wg  
            !
            totalweight=totalweight+wg
            !
          END IF
        END DO
        IF(ABS(totalweight-1._dp)>eps) CALL errore('main','wrong totalweight',1)  
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    !
    imag_tot = 0._dp
    imag_max = 0._dp
    !
    CALL allocate_fc2_grid(nRout, nat, fc)
    fc%nq = (/ nq1, nq2, nq3 /)
    DO na1=1,nat
    DO na2=1,nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3            
          !
          DO i = 1, fc%n_R
            fc%FC(jn1,jn2,i) = DBLE(matR(j1,j2,na1,na2,i))
            fc%yR(:,i) = NINT(Rout(:,i))
            IF( ALL(fc%yR(:,i)==0) ) fc%i_0 = i
            !
            imag_tot = imag_tot + ABS(IMAG(matR(j1,j2,na1,na2,i)))
            imag_max = MAX( imag_max, ABS(IMAG(matR(j1,j2,na1,na2,i))) )
          ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    fc%xR = DBLE(fc%yR)
    CALL cryst_to_cart(fc%n_R,fc%xR,at,+1)
    !
    IF(imag_tot>eps) WRITE(*,*) "CHECK! Sum of Imaginary parts:",imag_tot
    IF(imag_max>eps) WRITE(*,*) "CHECK! Maximum imaginary part:",imag_max
    !
  END SUBROUTINE quter
  !
  SUBROUTINE expand_matR(iR,nR,nat,matR)
    IMPLICIT NONE
    COMPLEX(DP),POINTER,INTENT(inout) :: matR(:,:,:,:,:)
    INTEGER,INTENT(in)    :: iR,nat
    INTEGER,INTENT(inout) :: nR
    COMPLEX(DP),POINTER   :: auxR(:,:,:,:,:)
    !
    IF(.not.associated(matR))THEN
!       print*, "new matR space created"
      ALLOCATE(matR(3,3,nat,nat,1))
      nR = 1
      matR = 0._dp
      RETURN
    ENDIF
    !
    IF(iR<=nR) RETURN
    !
!     print*, "expanded matR to hold", nR, "R vectors"
    auxR => matR
    ALLOCATE(matR(3,3,nat,nat,iR))
    matR(1:3,1:3,1:nat,1:nat,1:nR) = auxR(1:3,1:3,1:nat,1:nat,1:nR) 
    matR(1:3,1:3,1:nat,1:nat,nR+1:iR) = 0._dp
    DEALLOCATE(auxR)
    nR = iR
  END SUBROUTINE
  !
  ! given a list of R vectors, find the index of R in the list, 
  ! if R is not found it is appended to the list. NR is changed accordingly
  INTEGER FUNCTION R_list_idx(NR, listR, R) RESULT(idx)
    IMPLICIT NONE
    INTEGER,INTENT(inout)  :: NR
    REAL(DP),INTENT(inout),POINTER :: listR(:,:)
    REAL(DP),INTENT(in)    :: R(3)
    !
    INTEGER :: i
    REAL(DP),PARAMETER :: eps = 1.d-12
    REAL(DP),POINTER :: tmpR(:,:)
    !
    IF(.not.associated(listR))THEN
!       print*, "new list of R initialized"
      ALLOCATE(listR(3,1))
      idx=1
      listR(:,1) = R
    ENDIF
    !
    DO i = 1,NR
      IF( SUM((R-listR(1:3,i))**2)<eps )THEN
        idx = i
        RETURN
      ENDIF
    ENDDO
    !
    tmpR => listR
    ALLOCATE(listR(3,NR+1))
    listr(1:3,1:NR) = tmpR(1:3,1:NR)
    DEALLOCATE(tmpR)
    idx = NR+1
    listR(1:3,idx)   = R
    ! NR=NR+1 ==> this is done later in expand_matR
  END FUNCTION R_list_idx
  !
END MODULE quter_module
