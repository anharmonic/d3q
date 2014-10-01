!
! Copyright (C) Lorenzo Paulatto, 2014
! Algorithm from Giorgia Fugallo and Michele Lazzeri
! Uses wsweight from PW/src
! 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quter_standalone

  USE kinds, ONLY : DP
  !
  PUBLIC  :: quter           ! computer FCs from dynamical matrices
  PUBLIC  :: fftinterp_mat2  ! computes one dynamical matrix from FCs interpolation
  PUBLIC  :: impose_asr2     ! applies "simple" acoustic sum rules
  !
  ! MANUAL:
  !
  ! In order to compute the force constants this code requires that the 
  ! dynamical matrices of the regular grid are stored in 5 variables:
  !  nq1, nq2, nq3 (aka nr1, nr2, nr3) : the number of points in the grid
  !  gridq(3,nq1*nq2*nq3) : the list of q points, any order
  !  matq(3,3,nat,nat,nq1*nq2*nq3) : the dynamical matrices for each point
  !
  ! The force constants are computed and stored in 3 variables:
  !   nR : the number of real space lattice vectors required to store the FCs
  !   gridR(3,nR) : the list of aforementioned vectors
  !   matR(3*nat,3*nat,nR) : the force constants
  !
  ! Optionally you can use impose_asr2 to enforce simple sum rules. 
  !   IMPORTANT: FCs must NOT be divided by the sqrt of mass when enforcing sum rules!
  !
  ! Interfaces:
  !  CALL quter(nq1, nq2, nq3, nat,tau,at,bg, matq, gridq, nR, gridR, matR)
  !  nat,tau,at,bg -> as everywhere in quantum espresso 
  !                  (cartesian coordinates, units of alat or 2pi/alat)
  !
  !  CALL fftinterp_mat2(xq, nat, nR, gridR, matR, D)
  !  D(3*nat,3*nat) -> the complex dynamical matrix at q (no division by mass is performed)
  !
  !  CALL impose_asr2(nat,nR,gridR,matR)
  !
  !  CALL mat2_diag(nat, D, w2)
  !  D -> on input the dyn.mat., on output the eigenvectors
  !  w2 -> phonon frequencies squared, \omega^2 
  !
  PRIVATE :: R_list_idx, expand_cmatR

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE impose_asr2(nat,nR,gridR,matR)
    USE input_fc,       ONLY : forceconst2_grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)     :: nat, nR
    REAL(DP),INTENT(in)    :: gridR(3,nR)
    REAL(DP),INTENT(inout) :: matR(3*nat,3*nat,nR)
    !
    INTEGER :: iR, a,b, i,j, nu,mu, i0
    REAL(DP):: delta
    !
    i0 = -1
    DO iR = 1,nR
      IF(ALL( gridR(:,iR)==0._dp ))THEN
        i0=iR
        EXIT
      ENDIF
    ENDDO
    IF(i0<0) CALL errore("impose asr2","cannot find R==0",1)
    !
    DO a = 1,3
    DO b = 1,3
      DO i = 1,nat
      nu = 3*(i-1)+a
        !
        delta = 0._dp
        !
        DO j = 1,nat
        mu = 3*(j-1)+b
        DO iR = 1,nR
          delta = delta+matR(mu,nu,iR)
        ENDDO
        ENDDO

        mu = 3*(i-1)+b
        matR(mu,nu,i0) = matR(mu,nu,i0) - delta
        !
      ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE impose_asr2
  ! \/o\________\\\_________________________________________/^>
  ! IN PLACE diagonalization of D
  SUBROUTINE mat2_diag(nat, D, w2)
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    COMPLEX(DP),INTENT(inout) :: D(3*nat, 3*nat)
    REAL(DP),INTENT(out)      :: w2(3*nat)
    !
    INTEGER  :: nb    ! block size
    ! BEWARE the next two variables are saved:
    INTEGER,save :: nat3=-1, lwork=-1 ! aux. var

    INTEGER :: info      ! flag saying if the exec. of libr. routines was ok

    INTEGER,EXTERNAL ::ILAENV ! function which gives block size
    !
    REAL(DP), ALLOCATABLE    :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    ! Compute optimal storage space, only do it once unless the number of 
    ! atoms has misteriously changed (if it changes often, make lwork a parameter)
    IF ( nat3 /= 3*nat .or. lwork < 0 ) THEN
      !
      IF(nat3 > 0 ) PRINT*, "WARNING! Redoing ILAENV"
      nat3 = 3*nat
      !     check for the block size
      nb = ILAENV( 1, 'ZHETRD', 'U', nat3, -1, -1, -1 )
      IF (nb<1) nb = MAX(1,nat3)
      IF (nb==1 .or. nb>=nat3) then
        lwork=2*nat3-1
      ELSE
        lwork = (nb+1)*nat3
      ENDIF
      !
    ENDIF
    !
    ALLOCATE(work (lwork))
    ALLOCATE(rwork (3*nat3-2))
    !
    CALL ZHEEV('V','U',nat3,D,nat3,w2,work,lwork,rwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(rwork)
    DEALLOCATE(work)
    !
  END SUBROUTINE mat2_diag  
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE fftinterp_mat2(xq, nat, nR, gridR, matR, D)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq(3)
    !
    INTEGER,INTENT(in)  :: nat, nR
    REAL(DP),INTENT(in) :: matR(3*nat,3*nat,nR)
    REAL(DP),INTENT(in) :: gridR(3,nR)
    !
    COMPLEX(DP),INTENT(out) :: D(3*nat, 3*nat)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
    DO i = 1, nR
      arg = tpi * SUM(xq(:)*gridR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:, :) = D(:, :) + phase * matR(:, :, i)
    END DO
    !
  END SUBROUTINE fftinterp_mat2
  ! \/o\________\\\_________________________________________/^>
  !
  SUBROUTINE quter(nq1, nq2, nq3, nat,tau,at,bg, matq, gridq, nR, gridR, matR)
    USE constants, ONLY : tpi
    IMPLICIT NONE
    ! Dummy arguments
    INTEGER,INTENT(in)     :: nq1, nq2, nq3 ! dimensions of the q-points grid
    INTEGER,INTENT(in)     :: nat ! number of atoms
    REAL(DP),INTENT(in)    :: tau(3,nat) ! atom positions (alat units)
    REAL(DP),INTENT(in)    :: at(3,3), bg(3,3) ! real, reciprocal lattice
    COMPLEX(DP),INTENT(in) :: matq(3,3,nat,nat,nq1*nq2*nq3)
    REAL(DP),INTENT(in)    :: gridq(3,nq1*nq2*nq3)
    INTEGER,INTENT(out)              :: nR
    REAL(DP),ALLOCATABLE,INTENT(out) :: matR(:,:,:)
    REAL(DP),ALLOCATABLE,INTENT(out) :: gridR(:,:)
    !
    INTEGER,PARAMETER :: far = 2
    REAL(DP),PARAMETER :: eps = 1.d-8
    INTEGER :: i, na1, na2, j1,j2, jn1,jn2,  iiq, iR, l1,l2,l3
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
    INTEGER, PARAMETER:: nrwsx=200
    INTEGER :: nrws
    REAL(DP) :: atws(3,3) ! supercell of size nq1 x nq2 x nq3
    REAL(DP) :: wg, rws(0:3,nrwsx)
    REAL(DP),EXTERNAL :: wsweight
    !
    COMPLEX(DP),POINTER :: cmatR(:,:,:,:,:) => null()
    !
    nqt = nq1*nq2*nq3
    atws(:,1) = nq1*at(:,1)
    atws(:,2) = nq2*at(:,2)
    atws(:,3) = nq3*at(:,3)
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,atws)
    !
    ! Construct a big enough lattice of Supercell vectors
    nRbig = (2*far*nq1+1)*(2*far*nq2+1)*(2*far*nq3+1)
    ALLOCATE(Rbig(3,nRbig))
    nRbig=0
    DO l1=-far*nq1,far*nq1
    DO l2=-far*nq2,far*nq2
    DO l3=-far*nq3,far*nq3
      nRbig=nRbig+1
      Rbig(:, nRbig) = at(:,1)*l1 +at(:,2)*l2 +at(:,3)*l3
    END DO
    END DO
    END DO
    IF(nRbig/=size(Rbig)/3) call errore('main','wrong nRbig',1)  
    print*, "seeking over ", nRbig," vectors"
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
          dist(:) = Rbig(:,iR) + tau(:,na1) - tau(:,na2)
          wg = wsweight(dist,rws,nrws)
          wg = wg/DFLOAT(nqt)
          IF(wg /= 0) THEN 
            !
            DO iiq=1,nqt
              !   
              arg=tpi*SUM(gridq(:,iiq)*Rbig(:,iR)) 
              aus = aus + DCMPLX(cos(arg),sin(arg))*matq(j1,j2,na1,na2,iiq)
              !
            END DO
            !
            Rx = Rbig(:,iR)
            CALL cryst_to_cart(1,Rx,bg,-1)
            ! find if we have already used this R point ...
            iout = R_list_idx(nRout,Rout,Rx)
            ! ... if not, increase storage
            CALL expand_cmatR(iout,nRout,nat,cmatR)
            cmatR(j1,j2,na1,na2,iout) = aus * wg  
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
    nR = nRout
    IF(allocated(matR)) DEALLOCATE(matR)
    ALLOCATE(matR(3*nat,3*nat,nR))
    IF(allocated(gridR)) DEALLOCATE(gridR)
    ALLOCATE(gridR(3,nR))
    !
    imag_tot = 0._dp
    imag_max = 0._dp
    !
    DO na1=1,nat
    DO na2=1,nat 
      DO j1=1,3
      jn1 = j1 + (na1-1)*3
      DO j2=1,3     
      jn2 = j2 + (na2-1)*3            
          !
          DO i = 1, nR
            matR(jn1,jn2,i) = DBLE(cmatR(j1,j2,na1,na2,i))
            gridR(:,i) = NINT(Rout(:,i))
            !
            imag_tot = imag_tot + ABS(IMAG(cmatR(j1,j2,na1,na2,i)))
            imag_max = MAX(imag_max, IMAG(cmatR(j1,j2,na1,na2,i)))
          ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    CALL cryst_to_cart(nR,gridR,at,+1)
    !
    IF(imag_tot>eps) WRITE(*,*) "CHECK! Sum of Imaginary parts:",imag_tot
    IF(imag_max>eps) WRITE(*,*) "CHECK! Maximum imaginary part:",imag_max
    !
  END SUBROUTINE quter
  !
  ! -------------------- Auxiliary subroutines --------------------
  !
  SUBROUTINE expand_cmatR(iR,nR,nat,cmatR)
    IMPLICIT NONE
    COMPLEX(DP),POINTER,INTENT(inout) :: cmatR(:,:,:,:,:)
    INTEGER,INTENT(in)    :: iR,nat
    INTEGER,INTENT(inout) :: nR
    COMPLEX(DP),POINTER   :: auxR(:,:,:,:,:)
    !
    IF(.not.associated(cmatR))THEN
!       print*, "new cmatR space created"
      ALLOCATE(cmatR(3,3,nat,nat,1))
      nR = 1
      RETURN
    ENDIF
    !
    IF(iR<=nR) RETURN
    !
!     print*, "expanded cmatR to hold", nR, "R vectors"
    auxR => cmatR
    ALLOCATE(cmatR(3,3,nat,nat,nR+1))
    cmatR(1:3,1:3,1:nat,1:nat,1:nR) = auxR(1:3,1:3,1:nat,1:nat,1:nR) 
    cmatR(1:3,1:3,1:nat,1:nat,nR+1) = 0._dp
    DEALLOCATE(auxR)
    nR = nR+1
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
      IF( NORM2(R-listR(1:3,i))<eps )THEN
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
    ! NR=NR+1 ==> this is done later in expand_cmatR
  END FUNCTION R_list_idx
  !
END MODULE quter_standalone