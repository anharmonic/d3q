!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! This program applies acousting sum rules to 3rd order dynamical matrices. 
! It is currently quite sketchy, but it works. It reads a force constant file in grid format
! (no sparse) and applies ASR iteratively.
! This file also contains tentative subroutines that applies asr in one go, but 
! they are not working properly.
! This code only reads a file called "mat3R" and writes to a file names "mat3R_asr".
!
!
MODULE asr3_module
  USE kinds,    ONLY : DP
  TYPE forceconst3_ofRR
    REAL(DP),ALLOCATABLE :: F(:,:,:, :,:,:) !dimension: (3,3,3, nat,nat,nat)
  ENDTYPE forceconst3_ofRR
  TYPE todo3_ofRR
    LOGICAL,ALLOCATABLE :: F(:,:,:, :,:,:) !dimension: (3,3,3, nat,nat,nat)
  ENDTYPE todo3_ofRR
  TYPE index_r_type
    ! list of R vectors, used only for acoustic sum rule
    INTEGER :: nR                    ! number of distinc R
    INTEGER,ALLOCATABLE :: nRi(:)    ! how many times each R appears
    INTEGER :: nRx                   ! maxval(nRi); nRi <= nRx
    INTEGER,ALLOCATABLE :: yR(:,:)   ! the list of distinct R
    INTEGER,ALLOCATABLE :: idR(:,:)  ! for each of the nR distinct R, the index
                                     ! of its nRi-th occurence in the global list
                                     ! yR(:, idR(1:nRi(1:nR2), 1:nR2) (padded with zeroes)
    INTEGER :: iRe0                  ! index of R = (0,0,0)                                    
    INTEGER,ALLOCATABLE :: idmR(:)   ! for each R, the index of -R
    INTEGER,ALLOCATABLE :: idRmR(:,:)! for each couple i,j=1,..nR
                                     ! index of R_i-R_j
  END TYPE
  
  INTEGER,SAVE :: iter = 1
  REAL(DP),PARAMETER :: eps  = 1.e-12_dp, eps2 = 1.e-24_dp
  REAL(DP),PARAMETER :: eps0 = 1.e-20_dp
  
  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! In the original form (as written by Q2R3) you have a matrix of 3rd order
  ! force constants for each couple of (R2, R3). The long list of (R2, R3) couples
  ! repeats many times R2 and R3 in no specific order.
  ! This subroutine takes the long list of R2 *or* R3 and make a list of the used ones.
  ! The ordering is stable as long as the limits are the same.
  ! Furthermore, it finds the posion of -R for each R 
  ! and of R-R' for each R and R' (if it is in the list)
  SUBROUTINE stable_index_R(nR0, R0, nR,nRx, nRi, R, idR, iRe0, idmR, idRmR, nq, use_modulo)
    USE input_fc, ONLY : ph_system_info
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nR0, R0(3,nR0)
    !
    INTEGER,INTENT(out) :: nR                    ! number of distinc R
    INTEGER,ALLOCATABLE,INTENT(out) :: nRi(:)    ! how many times each R appears
    INTEGER,INTENT(out) :: nRx                   ! maxval(nRi); nRi <= nRx
    INTEGER,ALLOCATABLE,INTENT(out) :: R(:,:)    ! the list of distinct R
    INTEGER,ALLOCATABLE,INTENT(out) :: idR(:,:)  ! for each of the nR distinct R, the index
                                ! of its nRi-th occurence in the global list
                                ! yR(:, idR(1:nRi(1:nR2), 1:nR2) (padded with zeroes)
    INTEGER,ALLOCATABLE,INTENT(out) :: idmR(:)   ! for each of the nR, index of -R
    INTEGER,ALLOCATABLE,INTENT(out) :: idRmR(:,:)! for each couple i,j=1,..nR
                                                 ! index of R_i-R_j
    INTEGER,INTENT(out) :: iRe0                  ! index of R = (0,0,0)      
    INTEGER,INTENT(in)  :: nq(3)      ! size of the initial grid, only used for use_modulo
    LOGICAL,INTENT(in)  :: use_modulo ! use periodic boundary conditions to find R -> -R correspondences
    !
    INTEGER :: i, j, k, kk, n
    INTEGER,ALLOCATABLE :: nRi_(:) ! how many times each R2,R3 appears
    INTEGER,ALLOCATABLE :: R_(:,:)
    !
    INTEGER :: minR(3), maxR(3), Rgr(3)
    !
    ALLOCATE(R_(3,nR0),nRi_(nR0))
    !
    ! Find lowest and highest possible value of x, y z among all R's
    DO i = 1,3
      minR(i) = MINVAL(R0(i,:))
      maxR(i) = MAXVAL(R0(i,:))
    ENDDO
    !
    nR = 0
    R_(:,:) = 0
    nRi_ = 0
    !
    DO i = minR(1), maxR(1)
    DO j = minR(2), maxR(2)
    DO k = minR(3), maxR(3)
      Rgr = (/ i,j,k /)
      ! Find if this R is used
      IF(find_index(nR0, R0, Rgr)>0)THEN
        nR = nR+1
        R_(:,nR) = Rgr
        ! Count them
        DO n = 1,nR0
          IF(ALL(R0(:,n)==Rgr)) nRi_(nR) = nRi_(nR)+1
        ENDDO
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    !
    ALLOCATE(R(3,nR),nRi(nR))
    R(:,1:nR) = R_(:,1:nR)
    nRi(1:nR) = nRi_(1:nR)
    DEALLOCATE(R_, nRi_)
    
    nRx = MAXVAL(nRi)
    ALLOCATE(idR(nRx,nR))
    idR = 0
    !
    kk = 0
    DO j = 1,nR
      !
      k = 0
      DO i = 1, nR0
        IF ( ALL(R(:,j) == R0(:,i)) ) THEN
          k = k+1
          IF(k>nRi(j)) CALL errore("index_R", "too many r",1)
          idR(k,j) = i
        ENDIF
        IF( ALL(R(:,j) == 0)) iRe0 = j
      ENDDO
      kk = kk +k
      !
    ENDDO
    !
    IF(kk/=nR0) CALL errore('index_R', "some R not found",2)
    !
    ! Find -R for each R
    ALLOCATE(idmR(nR))
    idmR = -1
    DO j = 1,nR
    DO i = 1,nR
      IF(use_modulo)THEN
        IF((MODULO(-R(1,i),nq(1)) == R(1,j)) .and. &
           (MODULO(-R(2,i),nq(2)) == R(2,j)) .and. &
           (MODULO(-R(3,i),nq(3)) == R(3,j)) ) idmR(i) = j
      ELSE
        IF ( ALL(R(:,i) == -R(:,j)) ) idmR(i) = j
      ENDIF
    ENDDO
    ENDDO
    !
    ! Find R_k = R_i-R_j for each i and j
    ALLOCATE(idRmR(nR,nR))
    idRmR = -1
    DO k = 1,nR
    DO j = 1,nR
    DO i = 1,nR
      IF(use_modulo)THEN
        IF((MODULO(R(1,i)-R(1,j),nq(1)) == R(1,k)) .and. &
           (MODULO(R(2,i)-R(2,j),nq(2)) == R(2,k)) .and. &
           (MODULO(R(3,i)-R(3,j),nq(3)) == R(3,k)) ) idRmR(i,j) = k
      ELSE
        IF ( ALL(R(:,k) == R(:,i)-R(:,j)) ) idRmR(i,j) = k
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    !
    !
    !print*, "R ===================", nR, iRe0
    !DO i = 1, nR
    !  WRITE(*,'(2i4,2x,3i3,2x,i3,2x,1000i5)') i,idmR(i), R(:,i), nRi(i)!, idR(1:nRi(i),i)
    !ENDDO
    
  END SUBROUTINE stable_index_R
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE index_2R(idx2,idx3,fc,idR23)
  ! Find the correspondence between two lists of R vectors (list2 and list3) 
  ! a single list of couples (R2,R3) 
  ! i.e. for each R_2 in list2 and R3 in list3 it find the index of (R2,R3) in the big list
    USE fc3_interpolate,       ONLY : grid
    IMPLICIT NONE
    !
    TYPE(index_r_type),INTENT(in) :: idx2, idx3
    TYPE(grid),INTENT(inout) :: fc
    INTEGER,INTENT(out),ALLOCATABLE :: idR23(:,:)
    INTEGER :: i2,i3, j
    
    ALLOCATE(idR23(idx2%nR,idx3%nR))
    idR23 = -1
    
    DO i2 = 1,idx2%nR
    DO i3 = 1,idx3%nR
      DO j = 1,fc%n_R
        IF(     ALL(idx2%yR(:,i2) == fc%yR2(:,j)) &
          .and. ALL(idx3%yR(:,i3) == fc%yR3(:,j)) ) THEN
          idR23(i2,i3) = j
          !WRITE(*,'(2i4,2x,2i4,2x,i4)') i2,i3, idR23(i2,i3)
        ENDIF
      ENDDO
      !WRITE(333,'(2i4,2x,2i4,2x,i4)') i2,i3, idR23(i2,i3)
    ENDDO
    ENDDO
    
  END SUBROUTINE index_2R
  ! \/o\________\\\_________________________________________/^>
  ! Find the position of an R in a list of Rs
  INTEGER FUNCTION find_index(n, list, R) RESULT(idR)
    IMPLICIT NONE
    INTEGER,INTENT(in) ::  n, list(3,n), R(3)
    DO idR = 1,n
      IF( ALL(list(:,idR) == R(:)) ) RETURN
    ENDDO
    idR = -1
  END FUNCTION
  ! \/o\________\\\_________________________________________/^>
  ! Reshape 3rd order force constant matrix from 3 to 6 indexes (dir>0), or back (dir<0)
  SUBROUTINE fc_3idx_2_6idx(nat, fc3, fc6, dir)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP),INTENT(inout) :: fc3(3*nat, 3*nat, 3*nat)   ! d3 in cartesian basis
    REAL(DP),INTENT(inout) :: fc6(3,3,3, nat,nat,nat)    ! d3 in pattern basis
    INTEGER,INTENT(in) :: nat
    INTEGER,INTENT(in) :: dir
    !
    INTEGER :: i, j, k
    INTEGER :: at_i, at_j, at_k
    INTEGER :: pol_i, pol_j, pol_k
    !
    DO k = 1, 3*nat
      at_k  = 1 + (k-1)/3
      pol_k = k - 3*(at_k-1)
      !
      DO j = 1, 3*nat
        at_j  = 1 + (j-1)/3
        pol_j = j - 3*(at_j-1)
        !
        DO i = 1, 3*nat
          at_i  = 1 + (i-1)/3
          pol_i = i - 3*(at_i-1)
          !
          IF(dir > 0) THEN
            fc6(pol_i, pol_j, pol_k, at_i, at_j, at_k) &
              = fc3(i,j,k)
          ELSE IF( dir<0) THEN
            fc3(i,j,k) = &
                fc6(pol_i, pol_j, pol_k, at_i, at_j, at_k)
          ELSE
            CALL errore("fc_3idx_2_6idx","wrong direction",1)
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE fc_3idx_2_6idx
  ! \/o\________\\\_________________________________________/^>
  ! Initially the force constants are indexed on coupe of (R2,R3)
  ! this subroutine reshuffle the elements and assign them to two different
  ! indexes: one for R2 and one for R3. 
  SUBROUTINE reindex_fc3(nat,fc,idR23,idx2,idx3,fx,dir)
    USE fc3_interpolate,       ONLY : grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat, dir
    TYPE(grid),INTENT(inout) :: fc
    TYPE(index_r_type) :: idx2, idx3
    INTEGER,INTENT(in) :: idR23(idx2%nR,idx3%nR)
    TYPE(forceconst3_ofRR),INTENT(inout) :: fx(idx2%nR,idx3%nR)
    !
    INTEGER :: iR2,iR3

    IF(dir>0) THEN
        DO iR2 = 1,idx2%nR
        DO iR3 = 1,idx3%nR
          IF(.not.ALLOCATED(fx(iR2,iR3)%F)) &
            ALLOCATE(fx(iR2,iR3)%F(3,3,3,nat,nat,nat))
          fx(iR2,iR3)%F = 0._dp
        ENDDO
        ENDDO
    ELSE IF(dir<0) THEN
      IF(.not.ALLOCATED(fc%fc)) &
        ALLOCATE(fc%fc(3*nat,3*nat,3*nat, fc%n_R))
      fc%fc = 0._dp
    ELSE
          CALL errore("reindex_fc3","bad direction (+1/-1)",1)
    ENDIF
    !
    DO iR2 = 1,idx2%nR
    DO iR3 = 1,idx3%nR
      IF(idR23(iR2,iR3)>0) THEN
        CALL fc_3idx_2_6idx(nat, &
                            fc%fc(:,:,:,  idR23(iR2,iR3)), &
                            fx(iR2,iR3)%F, dir )
      
      ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE reindex_fc3
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE upindex_fcx(nat,idx,fx,up)
    USE fc3_interpolate,       ONLY : grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(index_r_type) :: idx
    TYPE(forceconst3_ofRR),INTENT(in) :: fx(idx%nR,idx%nR)
    REAL(DP),ALLOCATABLE,INTENT(out) :: up(:,:,:,:,:,:)
    !
    INTEGER :: iR2,iR3, a,b,c,i,j,k, ii,jj,kk, iR,iRp,iRpp
    INTEGER :: natR
    
    natR = nat*idx%nR
    ALLOCATE(up(3,3,3, natR,natR,natR))
    !
    DO iRpp = 1,idx%nR
    DO iRp = 1,idx%nR
    DO iR = 1,idx%nR
      iR2 = idx%idRmR(iRp,iR)   ! R2 = R'-R
      iR3 = idx%idRmR(iRpp,iR)  ! R3 = R''-R
      IF(iR2>0 .and. iR3>0)THEN
        DO k = 1,nat
        kk = k + (iRpp-1)*nat
          DO j = 1,nat
          jj = j + (iRp-1)*nat
            DO i = 1,nat
            ii = i + (iR-1)*nat
                WRITE(9999,'(100i3)') ii,jj,kk,iR,iRp,iRpp,i,j,k, idx%nR, natR
                up(:,:,:,ii,jj,kk)=fx(iR2,iR3)%F(:,:,:, i,j,k)
                WHERE(ABS(up(:,:,:,ii,jj,kk))<1.d-15) up(:,:,:,ii,jj,kk) = 0._dp
                WRITE(9999,'(3e25.15)') up(:,:,:,ii,jj,kk)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        print*, "missing one"
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    
    CALL upsum_fcx(natR, up)
    !
  END SUBROUTINE upindex_fcx
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE upsum_fcx(natR,up)
    USE fc3_interpolate,       ONLY : grid
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: natR
    REAL(DP),INTENT(inout) :: up(3,3,3,natR,natR,natR)
    REAL(DP) :: dw(3,3,3,natR,natR,natR)
    !
    INTEGER :: i,j,k
    INTEGER :: i_,j_,k_
    REAL(DP) :: q1,q2,q3, fact
    !
    dw = 0._dp
    fact = 1./16.
    DO k = 1, natR
    DO k_ = 1, natR
      q3 = d(k,k_) - fact
      DO j = 1, natR
      DO j_ = 1, natR
        q2 = d(j,j_) - fact
        DO i = 1, natR
        DO i_ = 1, natR
          q1 = d(i,i_) - fact
            dw(:,:,:,i,j,k) = dw(:,:,:,i,j,k) &
                            + q1*q2*q3*up(:,:,:, i_,j_,k_)
        ENDDO
        ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    
    DO k = 1, natR
      DO j = 1, natR
        DO i = 1, natR
          WRITE(1001,'(100i3)') i,j,k
          WRITE(1001,'(3e25.15)') up(:,:,:,i,j,k)
          WRITE(1002,'(100i3)') i,j,k
          WRITE(1002,'(3e25.15)') dw(:,:,:,i,j,k)
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE upsum_fcx
  !
  ! \/o\________\\\_________________________________________/^>
  ! symmetrize 3-body force constants wrt permutation of the indexes
  ! Note that the first R of the FC is always zero, hence we have to
  ! use translational symmetry, i.e.:
  ! F^3(R2,0,R3| b,a,c | j,i,k )  --> F^3(0,-R2,R3-R2| b,a,c | j,i,k )
  FUNCTION perm_symmetrize_fc3(nat,idx,fx) RESULT(delta)
    USE timers, ONLY : t_asr3s
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(index_r_type) :: idx
    TYPE(forceconst3_ofRR),INTENT(inout) :: fx(idx%nR,idx%nR)
    !TYPE(forceconst3_ofRR),ALLOCATABLE,SAVE :: fsym(:,:)
    TYPE(todo3_ofRR),ALLOCATABLE,SAVE       :: todo(:,:)
    !
    INTEGER :: iR2,iR3, a,b,c, i,j,k
    INTEGER :: mR2, mR3, iR2mR3, iR3mR2
    REAL(DP) :: avg, delta
    REAL(DP),PARAMETER :: onesixth = 1._dp/6._dp !0.1666666666666667_dp

    CALL t_asr3s%start()
    delta = 0._dp

    IF(.not.allocated(todo)) THEN
      !ALLOCATE(fsym(idx%nR,idx%nR))
      ALLOCATE(todo(idx%nR,idx%nR))
      DO iR2 = 1,idx%nR
      DO iR3 = 1,idx%nR
        ALLOCATE(todo(iR2,iR3)%F(3,3,3, nat,nat,nat))
        todo(iR2,iR3)%F = .true.
      ENDDO
      ENDDO
    ELSE
      DO iR2 = 1,idx%nR
      DO iR3 = 1,idx%nR
        todo(iR2,iR3)%F = .true.
      ENDDO
      ENDDO
    ENDIF
    !
    DO iR2 = 1,idx%nR
    R3_LOOP : &
    DO iR3 = 1,idx%nR
      mR2 = idx%idmR(iR2)
      mR3 = idx%idmR(iR3)
      !
      iR2mR3 = idx%idRmR(iR2,iR3)
      iR3mR2 = idx%idRmR(iR3,iR2)
      !
      IF(iR2mR3<0 .or. iR3mR2<0.or.mR3<0 .or. mR2<0) THEN
        IF( ANY(ABS(fx(iR2,iR3)%F) > eps0) ) &
           CALL errore('impose_asr3', "A matrix element that should be zero isn't. If using non-center mat3R use option '-m' ", 1)
        fx(iR2,iR3)%F = 0._dp
        CYCLE R3_LOOP
      ENDIF
      !
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,c,i,j,k,avg)  &
!$OMP REDUCTION(+:delta) COLLAPSE(6) 
      DO a = 1,3
      DO b = 1,3
      DO c = 1,3
        DO i = 1,nat
        DO j = 1,nat
        DO k = 1,nat
        
          IF( todo(iR2,iR3)%F(a,b,c, i,j,k) ) THEN
            avg = onesixth * ( &
                      fx(iR2,   iR3   )%F(a,b,c, i,j,k) &
                    + fx(iR3,   iR2   )%F(a,c,b, i,k,j) &
                    + fx(mR2,   iR3mR2)%F(b,a,c, j,i,k) &
                    + fx(iR3mR2,mR2   )%F(b,c,a, j,k,i) &
                    + fx(mR3,   iR2mR3)%F(c,a,b, k,i,j) &
                    + fx(iR2mR3,mR3   )%F(c,b,a, k,j,i) &
                    )

            delta = delta + (fx(iR2,iR3)%F(a,b,c, i,j,k) - avg)**2

            fx(iR2,iR3)%F(a,b,c, i,j,k) = avg
            fx(iR3,   iR2   )%F(a,c,b, i,k,j) = avg
            fx(mR2,   iR3mR2)%F(b,a,c, j,i,k) = avg
            fx(iR3mR2,mR2   )%F(b,c,a, j,k,i) = avg
            fx(mR3,   iR2mR3)%F(c,a,b, k,i,j) = avg
            fx(iR2mR3,mR3   )%F(c,b,a, k,j,i) = avg

            todo(iR2,   iR3   )%F(a,b,c, i,j,k) = .false.
            todo(iR3,   iR2   )%F(a,c,b, i,k,j) = .false.
            todo(mR2,   iR3mR2)%F(b,a,c, j,i,k) = .false.
            todo(iR3mR2,mR2   )%F(b,c,a, j,k,i) = .false.
            todo(mR3,   iR2mR3)%F(c,a,b, k,i,j) = .false.
            todo(iR2mR3,mR3   )%F(c,b,a, k,j,i) = .false.
          ENDIF
        
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO R3_LOOP
    ENDDO
    !
    delta = DSQRT(delta)
    !
    CALL t_asr3s%stop()
    !
  END FUNCTION perm_symmetrize_fc3
  !
  ! \/o\________\\\_________________________________________/^>
  ! Imposes sum rule on the first index of the FCs
  REAL(DP) FUNCTION impose_asr3_1idx(nat,idx,fx,pow) RESULT(delta)
    USE timers, ONLY : t_asr3a
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(index_r_type) :: idx
    TYPE(forceconst3_ofRR),INTENT(inout) :: fx(idx%nR,idx%nR)
    REAL(DP),INTENT(in) :: pow
    !
    INTEGER :: iR2,iR3, a,b,c, i,j,k
    REAL(DP):: deltot, delsig, delperm
    REAL(DP):: d1, q1, r1
    REAL(DP) :: invpow
    !
    CALL t_asr3a%start()
    !
    invpow = 1._dp/pow

    deltot = 0._dp
    delsig = 0._dp
    DO iR2 = 1,idx%nR
      !
      DO j = 1,nat
      DO i = 1,nat
        DO c = 1,3
        DO b = 1,3
        DO a = 1,3
          !
          ! The sum is on one R (it does not matter which)
          ! and the first atom index
          d1 = 0._dp
          q1 = 0._dp
!$OMP PARALLEL DO PRIVATE(k,iR3) REDUCTION(+:d1,q1) COLLAPSE(2)
          R3_LOOP : &
          DO iR3 = 1,idx%nR
          DO k = 1,nat
              d1 = d1+fx(iR2, iR3 )%F(a,b,c, i,j,k)
              q1 = q1+ABS(fx(iR2, iR3 )%F(a,b,c, i,j,k))**pow
          ENDDO
          ENDDO R3_LOOP
!$OMP END PARALLEL DO
          !
          deltot = deltot + ABS(d1)**pow
          delsig = delsig + d1
          !
          IF(q1>0._dp) THEN
            r1 = d1/q1
            !
!$OMP PARALLEL DO PRIVATE(k,iR3) COLLAPSE(2)
            R3_LOOPb : &
            DO iR3 = 1,idx%nR
            DO k = 1,nat
                fx(iR2,iR3)%F(a,b,c, i,j,k) = fx(iR2,iR3)%F(a,b,c, i,j,k) &
                                          - r1 * ABS(fx(iR2,iR3)%F(a,b,c, i,j,k))**pow
            ENDDO
            ENDDO R3_LOOPb
!$OMP END PARALLEL DO
            !
          ENDIF
          !
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
    ENDDO
    !
    CALL t_asr3a%stop()
    !
    ! Re-symmetrize the matrix
    delperm = perm_symmetrize_fc3(nat,idx,fx)
    !
    IF(iter==1) THEN
    WRITE(*,'(2x,a)') "Minimization started: create a file named 'STOP' to stop."
    WRITE(*,'(2x,a,a10,3a20)')   "asr3", "iter", "SQRT(deltot)", "SQRT(ABS(deltot-delsig**2))", "delperm"
    ENDIF
    WRITE(*,'(2x,a,i10,3e20.6)') "asr3", iter, &
      (deltot)**invpow, (ABS(deltot-delsig**2))**invpow, delperm
    iter = iter+1
    delta = ABS(deltot)**invpow
    !
    !
  END FUNCTION impose_asr3_1idx
  ! \/o\________\\\_________________________________________/^>

  ! \/o\________\\\_________________________________________/^>
  ! Imposes sum rule on the first index of the FCs
  REAL(DP) FUNCTION impose_asr3_mauri(nat,idx,fx) RESULT(delta)
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: nat
    TYPE(index_r_type) :: idx
    TYPE(forceconst3_ofRR),INTENT(inout) :: fx(idx%nR,idx%nR)
    !
    INTEGER :: iR2,iR3, a,b,c, i,j,k, iR, iRp, iRpp
    INTEGER :: iR2_,iR3_,i_,j_,k_, iR_, iRp_, iRpp_
    INTEGER :: dropped, kept
    REAL(DP) :: q1,q2,q3, delperm
    LOGICAL,ALLOCATABLE :: todo(:,:), todo_(:,:)
    TYPE(forceconst3_ofRR),ALLOCATABLE :: faux(:,:)
    REAL(DP) :: fact
    !
    ALLOCATE(todo(idx%nR,idx%nR), todo_(idx%nR,idx%nR))
    todo  = .true.
    todo_ = .true.

    ALLOCATE(faux(idx%nR,idx%nR))
    DO iR3 = 1,idx%nR
    DO iR2 = 1,idx%nR
      ALLOCATE(faux(iR2,iR3)%F(3,3,3, nat,nat,nat))
      faux(iR2,iR3)%F = fx(iR2,iR3)%F
      fx(iR2,iR3)%F = 0._dp
    ENDDO
    ENDDO
    !
    fact = 1._dp/DBLE(2*2*2 *2) !DBLE(idx%nR*nat)
    !
    DO iR2 = 1,idx%nR
    print*, 100*iR2/DBLE(idx%nR), "%"
    DO iR3 = 1,idx%nR
    iR   = idx%iRe0
    iRp  = iR2
    iRpp = iR3
!     DO iRpp = 1,idx%nR
!     DO iRp = 1,idx%nR
!     DO iR = 1,idx%nR
! 
!     
!     iR2 = idx%idRmR(iRp,iR)   ! R2 = R'-R
!     iR3 = idx%idRmR(iRpp,iR)  ! R3 = R''-R
!     ROUTER_IF : IF(iR2>0 .and. iR3>0 .and. todo(iR2,iR3))THEN
!     todo(iR2,iR3) = .false.
    
      todo_ = .true.
      DO iRpp_ = 1,idx%nR
      DO iRp_  = 1,idx%nR
      DO iR_   = 1,idx%nR
      
      iR2_ = idx%idRmR(iRp_,iR_)   ! R2 = R'-R
      iR3_ = idx%idRmR(iRpp_,iR_)  ! R3 = R''-R
      
      RINNER_IF : IF(iR2_>0 .and. iR3_>0 .and. todo_(iR2_,iR3_))THEN
!       todo_(iR2_,iR3_) = .false.
!       DO iR2_ = 1,idx%nR
!       DO iR3_ = 1,idx%nR
        DO c = 1,3
        DO b = 1,3
        DO a = 1,3
          DO k_ = 1,nat
          DO k = 1,nat
!             q3 = d(iR3,iR3_)*d(k,k_)-fact
            q3 = d(iRpp,iRpp_)*d(k,k_)-fact

            DO j_ = 1,nat
            DO j = 1,nat
!               q2 = d(iR2,iR2_)*d(j,j_)-fact
              q2 = d(iRp,iRp_)*d(j,j_)-fact

              DO i_ = 1,nat
              DO i = 1,nat
!                 q1 = d(i,i_)-fact
                q1 = d(iR,iR_)*d(i,i_)-fact
                !
                fx(iR2, iR3)%F(a,b,c, i,j,k) = &
                fx(iR2, iR3)%F(a,b,c, i,j,k) + &
                    q1*q2*q3*faux(iR2_, iR3_)%F(a,b,c, i_,j_,k_)
              ENDDO
              ENDDO
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        ENDDO

      ELSE
!         print*, iR2_, iR3_, todo(iR2_,iR3_)
!         STOP 888
      ENDIF RINNER_IF
      ENDDO
      ENDDO
      ENDDO
!     ENDIF ROUTER_IF
!     ENDDO
    ENDDO
    ENDDO
    !
!     delperm = perm_symmetrize_fc3(nat,idx,fx)
    !print*, "-->", delperm
    delta = 0._dp
    !
  END FUNCTION impose_asr3_mauri  
  ! \/o\________\\\_________________________________________/^>
  !
  INTEGER PURE FUNCTION d(i,j)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: i,j
    IF(i==j) THEN
      d=1
    ELSE
      d=0
    ENDIF
  END FUNCTION
  
  
END MODULE asr3_module


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM asr3

  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : aux_system, ph_system_info
  USE fc3_interpolate, ONLY : grid
!  USE io_global,       ONLY : *
  USE asr3_module
  USE more_constants,  ONLY : print_citations_linewidth
  USE clib_wrappers,        ONLY : memstat
  USE cmdline_param_module
  USE timers
  IMPLICIT NONE
  !
  TYPE(grid)             :: fc
  TYPE(ph_system_info)   :: s
  TYPE(index_r_type)     :: idx2,idx3,idx_wrap
  INTEGER                :: kb
  INTEGER,ALLOCATABLE :: idR23(:,:)
  TYPE(forceconst3_ofRR),ALLOCATABLE :: fx(:,:)
  !
  INTEGER :: j, ios
  REAL(DP) :: delta, threshold
  ! 
  CHARACTER(len=256) :: self, filein, fileout, aux
  INTEGER :: niter_max
  REAL(DP) :: pow
  !
  TYPE(grid) :: fcb
  INTEGER :: nq(3), nq_trip
  CHARACTER(len=:),ALLOCATABLE :: cmdline
  LOGICAL :: use_modulo

  filein   = cmdline_param_char("i", "mat3R")
  fileout  = cmdline_param_char("o", TRIM(filein)//".asr")
  threshold    = cmdline_param_dble("t", 1.d-12)
  niter_max   = cmdline_param_int("n", 1000)
  pow         = cmdline_param_dble("p", 2._dp)
  use_modulo  = cmdline_param_logical("m")
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: d3_asr3.x [-i FILEIN] [-o FILEOUT] [-t THR] [-n NITER] [-p POWER] [-m]"
      WRITE(*,*) ""
      WRITE(*,'(a)') " FILEIN  : input 3rd order force constants in grid form (default: mat3R)"
      WRITE(*,'(a)') " FILEOUT : output file with with sum rule applied (default: 'FILEIN.asr3')"
      WRITE(*,'(a)') "           "
      WRITE(*,'(a)') " THR     : stop when residual violation of ASR is smaller than THR (default: 1.d-12)"
      WRITE(*,'(a)') " NITER   : maximum number of iterations (default: 1000)"
      WRITE(*,'(a)') " POWER   : apply the correction to each matrix element proportionally to itself "
      WRITE(*,'(a)') "           to power POWER (default: 2), the value 1 can be more effective some times,"
      WRITE(*,'(a)') "           fractional numbers (e.g. 1.5) are also accepted."
      WRITE(*,'(a)') " -m : use periodic boundary conditions when comparing supercell vectors,"
      WRITE(*,'(a)') "      only use for non-centered mat3R files (wrong results in any other case)"
      STOP 1
  ENDIF
  
  CALL cmdline_check_exausted()
  !
  WRITE(*,*) "Note: create a file called 'STOP' to stop the code at next iteration and write out the result."
  WRITE(*,*) "Use option '-h' for more help."
  WRITE(*,*)
  !
  WRITE(*,*) " --- PARAMETERS : ---"
  WRITE(*,*) " power to distribute ASR correction term: ", pow
  WRITE(*,*) " threshold of ASR violation: ", threshold
  WRITE(*,*) " number of iterations max: ", niter_max
  WRITE(*,*) " use periodic bounday conditions: ", use_modulo
  WRITE(*,*) " --- ------------ ---"
  ! ----------------------------------------------------------------------------
    CALL t_asr3io%start()
  CALL fc%read(filein, S)
  CALL aux_system(s)
    CALL t_asr3io%stop()
  !
  CALL memstat(kb)
  WRITE(*,*) " grid size", fc%nq 
  WRITE(*,*) "Reading : done. //  Mem used:", DBLE(kb)/1000._dp, "Mb"
  !
  ! ----------------------------------------------------------------------------
    CALL t_asr3idx%start()
  CALL stable_index_R(fc%n_R, fc%yR2, idx2%nR, idx2%nRx, idx2%nRi, &
                idx2%yR, idx2%idR, idx2%iRe0, idx2%idmR, idx2%idRmR, fc%nq, use_modulo)
  !
  CALL stable_index_R(fc%n_R, fc%yR3, idx3%nR, idx3%nRx, idx3%nRi, &
                idx3%yR, idx3%idR, idx3%iRe0, idx3%idmR, idx3%idRmR, fc%nq, use_modulo)

!   CALL stable_index_R(fc%n_R, fc%yR3, idx_wrap%nR, idx_wrap%nRx, idx_wrap%nRi, &
!           idx_wrap%yR, idx_wrap%idR, idx_wrap%iRe0, idx_wrap%idmR, idx_wrap%idRmR, .true.)
  !
  ! Check that indexes are identical, this is not strictly necessary, 
  ! but it makes the rest easier
  IF(idx2%nR/=idx3%nR .or. ANY(idx2%yR/=idx3%yR) )&
    CALL errore("asr3", "problem with R",1)
  IF(ANY(idx2%idmR/=idx3%idmR))&
    CALL errore("asr3", "problem with -R",2)
  IF(ANY(idx2%idRmR/=idx3%idRmR))&
    CALL errore("asr3", "problem with R-R",3)
  !
  ! Map couple of in
  CALL index_2R(idx2,idx3, fc, idR23)
  
  CALL memstat(kb)
  WRITE(*,*) "R indexing : done. //  Mem used:", kb/1000, "Mb"
  ! ----------------------------------------------------------------------------
  !
  ALLOCATE(fx(idx2%nR, idx3%nR))

  CALL reindex_fc3(S%nat,fc,idR23,idx2,idx3,fx,+1)
    CALL t_asr3idx%stop()
  CALL memstat(kb)
  WRITE(*,*) "FC3 reindex : done. //  Mem used:", kb/1000, "Mb"

!   CALL upindex_fcx(S%nat, idx_wrap, fx, up)
!   WRITE(*,*) "Upscale reindex : done. //  Mem used:", kb/1000, "Mb"
! 
!!  threshold = impose_asr3_mauri(S%nat,idx2,fx)


     WRITE(*,*) "Pre-symmetrization:", perm_symmetrize_fc3(S%nat,idx2,fx)
     !
     APPLY_ASR : &
     DO j = 1,niter_max
       IF( impose_asr3_1idx(S%nat,idx2,fx,pow) < threshold) EXIT
       OPEN(unit=100, file="STOP", status='OLD', iostat=ios)
       IF(ios==0) THEN
         CLOSE(100,status="DELETE")
         EXIT APPLY_ASR 
       ENDIF
     ENDDO APPLY_ASR 
     !
   CALL memstat(kb)
   WRITE(*,*) "Impose asr3 : done. //  Mem used:", kb/1000, "Mb"
  ! ----------------------------------------------------------------------------
  !
    CALL t_asr3idx%start()
  CALL reindex_fc3(S%nat,fc,idR23,idx2,idx3,fx,-1)
    CALL t_asr3idx%stop()

    CALL t_asr3io%start()
  CALL fc%write(fileout, S)
    CALL t_asr3io%stop()
  !
  CALL print_citations_linewidth()
  !
  WRITE(*,'("   * WALL : ",f12.4," s")') get_wall()
  CALL print_timers_header()
  CALL t_asr3a%print()
  CALL t_asr3s%print()
  CALL t_asr3io%print()
  CALL t_asr3idx%print()
  !
END PROGRAM asr3
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
