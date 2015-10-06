!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE fc2_interpolate
  !
  USE kinds,    ONLY : DP
  USE nanoclock
  !
  ! I take forceconst2_grid globally in this module so I can USE it from here 
  ! elsewhere, which makes more sense
  USE input_fc, ONLY : forceconst2_grid
  
  ! \/o\________\\\_________________________________________/^>
  ! Moved to input_fc to avoid circular dependencies
!   TYPE forceconst2_grid
!     ! q points
!     INTEGER :: n_R = 0, i_0 = -1
!     INTEGER,ALLOCATABLE  :: yR(:,:) ! crystalline coords  3*n_R
!     REAL(DP),ALLOCATABLE :: xR(:,:) ! cartesian coords    3*n_R
!     REAL(DP),ALLOCATABLE :: FC(:,:,:) ! 3*nat,3*nat, n_R
!     INTEGER :: nq(3) ! initial grid size, only kept for reconstruction purposes
!   END TYPE forceconst2_grid
  !
  ! We can use two versions of the diagonalisation: one that saves the temporary work array dimension,
  ! the other that recomputes it each time. Saving it can be faster, but is not thread save
  ! i.e. it does not work with OMP unless we pass the workspace upstream, which is annoying. Here we use the
  ! pure version, which redoes everything locally and is safe.
  INTERFACE mat2_diag
!     MODULE PROCEDURE mat2_diag_save     ! temporary space is reused (save)
    MODULE PROCEDURE mat2_diag_pure       ! temporary space is reallocated every time
!     MODULE PROCEDURE mat2_diag_pure_dac ! same, using Divide & Conquer algorithm
  END INTERFACE

  ! I have three versions of the interpolation routine, they only differ by the 
  ! way OMP is used to parallelise the operations. In principle they are all three valid,
  ! but only the _flat one seems to work. Probably some bug with reduce operations on
  ! arrays, or something I do not understand about OMP.
  INTERFACE fftinterp_mat2
#define __PRECOMPUTE_PHASES
#ifdef __PRECOMPUTE_PHASES
!dir$ message "----------------------------------------------------------------------------------------------" 
!dir$ message "Using MKL vectorized Sin and Cos implementation, this can use more memory but should be faster" 
!dir$ message "----------------------------------------------------------------------------------------------" 
    MODULE PROCEDURE fftinterp_mat2_flat_mkl   ! use standard OMP but paralelize on the innermost cycle
#else
    MODULE PROCEDURE fftinterp_mat2_flat   ! use standard OMP but paralelize on the innermost cycle
#endif
!    MODULE PROCEDURE fftinterp_mat2_reduce  ! use standard OMP reduce on dummy variables
!    MODULE PROCEDURE fftinterp_mat2_safe    ! do not use dummy variable in OMP clauses
  END INTERFACE

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_reduce(xq, nat3, fc, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)   :: nat3
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,arg,phase) &
!$OMP REDUCTION(+: D)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D(:, :) = D(:, :) + phase * fc%fc(:, :, i)
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE fftinterp_mat2_reduce
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_flat(xq, nat3, fc, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)   :: nat3
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    INTEGER :: i, mu,nu
    !
    D = (0._dp, 0._dp)
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mu,nu)
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(2)
      DO mu= 1, nat3
      DO nu= 1, nat3
        D(nu,mu) = D(nu,mu) + phase * fc%fc(nu,mu, i)
      ENDDO
      ENDDO
!$OMP END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE fftinterp_mat2_flat
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_flat_mkl(xq, nat3, fc, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq(3)
    INTEGER,INTENT(in)   :: nat3
    TYPE(forceconst2_grid),INTENT(in) :: fc
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3)
    !
    INTEGER :: i, mu,nu
    REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
    COMPLEX(DP) :: vphase(fc%n_R)
    !
    D = (0._dp, 0._dp)
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
    FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq(:)*fc%xR(:,i))
#if defined(__INTEL)
    CALL vdCos(fc%n_R, varg, vcos)
    CALL vdSin(fc%n_R, varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mu,nu)
    DO i = 1, fc%n_R
!       arg = tpi * SUM(xq(:)*fc%xR(:,i))
!       phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(2)
      DO mu= 1, nat3
      DO nu= 1, nat3
        D(nu,mu) = D(nu,mu) + vphase(i) * fc%fc(nu,mu, i)
      ENDDO
      ENDDO
!$OMP END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE fftinterp_mat2_flat_mkl
  !
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc (safe version, for OpenMP)
  SUBROUTINE fftinterp_mat2_safe(xq, nat3, fc, D)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    INTEGER,INTENT(in)   :: nat3
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(nat3, nat3)
    !
    REAL(DP) :: arg
    COMPLEX(DP) :: phase
    COMPLEX(DP),ALLOCATABLE :: D_aux(:,:)
    INTEGER :: i
    !
    ! I'm using an aux variable because using a dummy variable in an OMP reduction
    ! clause is often miscompiled by ifort 14 (Feb 2015)
    !
    D = (0._dp, 0._dp)
!$OMP PARALLEL DEFAULT(none) SHARED(D,nat3,fc,xq) PRIVATE(i,arg,phase,D_aux)
    ALLOCATE(D_aux(nat3,nat3))
    D_aux = (0._dp, 0._dp)
!$OMP DO
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D_aux(:, :) = D_aux(:, :) + phase * fc%fc(:, :, i)
    END DO
!$OMP END DO
!$OMP CRITICAL(fc_interpolate_fftinterp_mat2_10)
    D = D + D_aux
!$OMP FLUSH(D)
!$OMP END CRITICAL(fc_interpolate_fftinterp_mat2_10)
    DEALLOCATE(D_aux)
!$OMP END PARALLEL
    !
  END SUBROUTINE fftinterp_mat2_safe
  
  ! \/o\________\\\_________________________________________/^>
  ! IN PLACE diagonalization of D
  SUBROUTINE mat2_diag_save(n, D, w2)
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: n
    COMPLEX(DP),INTENT(inout) :: D(n, n)
    REAL(DP),INTENT(out)      :: w2(n)
    !
    INTEGER  :: nb    ! block size
    INTEGER,save :: nat3=-1, lwork=-1 ! aux. var

    INTEGER :: info      ! flag saying if the exec. of libr. routines was ok

    INTEGER,EXTERNAL ::ILAENV ! function which gives block size
    !
    REAL(DP), ALLOCATABLE    :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    IF ( nat3 /= n .or. lwork < 0 ) THEN
      !     check for the block size
      nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
      IF (nb<1) nb = MAX(1,n)
      IF (nb==1 .or. nb>=n) then
        lwork=2*n-1
      ELSE
        lwork = (nb+1)*n
      ENDIF
      !
      IF(nat3 > 0 ) PRINT*, "WARNING! Redoing ILAENV"
      nat3 = n
      !
    ENDIF
    !
    ALLOCATE(work (lwork))
    ALLOCATE(rwork (3*n-2))
    !
    CALL ZHEEV('V','U',n,D,n,w2,work,lwork,rwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(rwork)
    DEALLOCATE(work)
    !
  END SUBROUTINE mat2_diag_save
  ! \/o\________\\\_________________________________________/^>
  ! IN PLACE diagonalization of D
  ! Cannot declare explicitly pure, because ZHEEV isn't pure
  SUBROUTINE mat2_diag_pure(n, D, w2)
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: n
    COMPLEX(DP),INTENT(inout) :: D(n, n)
    REAL(DP),INTENT(out)      :: w2(n)
    !
    INTEGER :: lwork, info
    REAL(DP),ALLOCATABLE    :: rwork(:)
    COMPLEX(DP),ALLOCATABLE :: work(:)
    !
    lwork = 2*n-1
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*n-2))
    !
    CALL ZHEEV('V','U',n,D,n,w2,work,lwork,rwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(rwork)
    DEALLOCATE(work)
    !
  END SUBROUTINE mat2_diag_pure
  !
  ! IN PLACE diagonalization of D using divide and conquer, 
  ! should be faster, but it'sactually slower for small number of atoms
  ! Cannot declare explicitly pure, because ILAENZ and ZHEEV aren't pure
  SUBROUTINE mat2_diag_pure_dac(n, D, w2)
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: n
    COMPLEX(DP),INTENT(inout) :: D(n, n)
    REAL(DP),INTENT(out)      :: w2(n)
    !
    INTEGER  :: lwork, lrwork, liwork ! aux
    INTEGER :: info
    COMPLEX(DP),ALLOCATABLE :: work(:)
    REAL(DP),ALLOCATABLE    :: rwork(:)
    INTEGER,ALLOCATABLE     :: iwork(:)
    !
    ! n is always greater than one:
    lwork  = 2*n+n**2
    lrwork = 1 + 5*n + 2*n**2
    liwork = 3+5*n
    
    
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(lrwork))
    ALLOCATE(iwork(liwork))
    !
    CALL ZHEEVD('V','U',n,D,n,w2,work,lwork,rwork,lrwork,iwork,liwork,info)
    CALL errore ('mat2_diag','ZHEEV info =/= 0',ABS(info))
    !
    DEALLOCATE(work,rwork,iwork)
    !
  END SUBROUTINE mat2_diag_pure_dac  !
  !
  ! Auxiliary subroutines follow:
  !
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq_safe(xq, S, fc2, freq, U)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: epsq = 1.e-6_dp
    REAL(DP) :: cq(3)
    !
    CALL fftinterp_mat2(xq, S%nat3, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    U = CONJG(U)
    
    ! Set patterns and frequency to exactly zero for Gamma (and Gamma+G)
    cq = xq
    CALL cryst_to_cart(1,cq,S%at,-1)
    IF(ALL( ABS(cq-INT(cq))<epsq ) )THEN
      freq(1:3) = 0._dp
!         U(:,1:3) = (0._dp, 0._dp)
    ENDIF
    
    WHERE    (freq >=  0._dp)
      freq = SQRT(freq)
    ELSEWHERE(freq < 0._dp)
      freq = -SQRT(-freq)
    ENDWHERE
    !
  END SUBROUTINE freq_phq_safe
  !
  ! Return 1 if q-point is not Gamma, 4 if it is (3 points in input)
  ! avoid divergence because of 1/\omega at Gamam for acoustic bands
  FUNCTION set_nu0(xq, at) RESULT(nu0)
    IMPLICIT NONE
    REAL(DP), INTENT(in) :: xq(3), at(3,3)
    INTEGER :: nu0
    !
    REAL(DP),PARAMETER :: epsq = 1.d-6
    REAL(DP) :: cq(3)
    cq = xq
    CALL cryst_to_cart(1,cq,at,-1)
    cq = cq-NINT(cq)
    IF(ALL( ABS(cq)<epsq ) )THEN
      nu0 = 4
    ELSE
      nu0 = 1
    ENDIF
  END FUNCTION
  !
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq(xq, S, fc2, freq, U)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: eps = 1.e-16_dp
    !
    CALL fftinterp_mat2(xq, S%nat3, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    U = CONJG(U)
    WHERE    (freq >  eps)
      freq = SQRT(freq)
    ELSEWHERE(freq < -eps)
      freq = -SQRT(-freq)
    ELSEWHERE ! i.e. freq=0
      freq = 0._dp
    ENDWHERE
      !
  END SUBROUTINE freq_phq
  !
  ! Compute Bose-Einstein distribution of freq
  PURE SUBROUTINE bose_phq(T, nat3, freq, bose)
    USE functions, ONLY : f_bose
    IMPLICIT NONE
    REAL(DP),INTENT(in)  :: T
    INTEGER,INTENT(in)   :: nat3
    REAL(DP),INTENT(in)  :: freq(nat3)
    REAL(DP),INTENT(out) :: bose(nat3)
    !
    ! Is the following mess really necessary? (3 days later: it is)
    !WHERE    (freq*T >  eps)
    !WHERE    (freq /=  0._dp)
    WHERE    (freq*T > 0._dp)
      bose = f_bose(freq, T)
    ELSEWHERE
      bose = 0._dp
    ENDWHERE
      !
  END SUBROUTINE bose_phq
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE dyn_cart2pat(d2in, nat3, u, dir, d2out)
    !   Rotates third derivative of the dynamical basis from cartesian axis
    !   to the basis of the modes. Rotation is not really in place
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! d2 matrix, input: in cartesian basis, output: on the patterns basis
    COMPLEX(DP),INTENT(inout) :: d2in(nat3, nat3)
    COMPLEX(DP),OPTIONAL,INTENT(out) :: d2out(nat3, nat3)
    INTEGER,INTENT(in)        :: nat3
    ! patterns (transposed, with respect to what we use in the d3q.x code)
    INTEGER,INTENT(IN)        :: dir  ! +1 -> cart2pat, -1 -> pat2car
    COMPLEX(DP),INTENT(in)    :: u(nat3, nat3)
    !
    COMPLEX(DP),ALLOCATABLE  :: d2tmp(:,:)
    COMPLEX(DP),ALLOCATABLE  :: u_fw(:,:), u_bw(:,:)
    !
    INTEGER :: a, b, i, j
    COMPLEX(DP) :: AUX
    REAL(DP),PARAMETER :: EPS = 1.e-8_dp
    !
    ALLOCATE(d2tmp(nat3, nat3))
    ALLOCATE(u_fw(nat3, nat3))
    ALLOCATE(u_bw(nat3, nat3))
    d2tmp = (0._dp, 0._dp)
    !
    !STOP "untested"
    !
    IF(dir==+1)THEN
      u_bw = u
      u_fw = CONJG(u)
    ELSEIF(dir==-1)THEN
      u_bw = TRANSPOSE(CONJG(u))
      u_fw = TRANSPOSE(u)
!       u_fw = CONJG(u)
!       u_bw = u
    ELSE
      CALL errore("dyn_cart2pat", "wrong action",1)
    ENDIF

    DO j = 1,nat3
    DO b = 1,nat3
      !
      DO i = 1,nat3
      DO a = 1,nat3
          d2tmp(i, j) = d2tmp(i, j) &
                          + u_fw(a,i) * d2in(a, b) * u_bw(b,j)
      ENDDO
      ENDDO
      !
    ENDDO
    ENDDO    
    !
    IF(present(d2out)) THEN
      d2out = d2tmp
    ELSE
      d2in  = d2tmp
    ENDIF
    DEALLOCATE(d2tmp, u_fw, u_bw)
    RETURN
  END SUBROUTINE dyn_cart2pat
  
  !
  ! \/o\________\\\_________________________________________/^>
  SUBROUTINE ip_cart2pat(d3in, nat3, u1, u2, u3, d3out)
    !   Rotates third derivative of the dynamical basis from cartesian axis
    !   to the basis of the modes. Rotation is not really in place
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! d3 matrix, input: in cartesian basis, output: on the patterns basis
    COMPLEX(DP),INTENT(inout) :: d3in(nat3, nat3, nat3)
    COMPLEX(DP),OPTIONAL,INTENT(out) :: d3out(nat3, nat3, nat3)
    INTEGER,INTENT(in)        :: nat3
    ! patterns (transposed, with respect to what we use in the d3q.x code)
    COMPLEX(DP),INTENT(in)    :: u1(nat3, nat3), u2(nat3, nat3), u3(nat3, nat3) 
    !
    COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
    !
    INTEGER :: a, b, c, i, j, k
    COMPLEX(DP) :: AUX
    REAL(DP),PARAMETER :: EPS = 1.e-8_dp
    !
    ALLOCATE(d3tmp(nat3, nat3, nat3))
    d3tmp = (0._dp, 0._dp)
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,a,b,c,AUX) REDUCTION(+: d3tmp) COLLAPSE(2)
    DO k = 1,nat3
    DO c = 1,nat3
    IF(ABS(u3(c,k))>EPS)THEN
      !
      DO j = 1,nat3
      DO b = 1,nat3
      IF(ABS(u2(b,j))>EPS)THEN
        !
        AUX = u2(b,j) * u3(c,k)
        DO i = 1,nat3
        DO a = 1,nat3
            d3tmp(i, j, k) = d3tmp(i, j, k) &
                            + u1(a,i) * AUX * d3in(a, b, c) 
        ENDDO
        ENDDO
        !
      ENDIF
      ENDDO
      ENDDO
      !
    ENDIF
    ENDDO
    ENDDO
!$OMP END PARALLEL DO
! !!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,a,b,c) REDUCTION(+: d3tmp) COLLAPSE(2)
!     DO k = 1,nat3
!     DO j = 1,nat3
!     DO i = 1,nat3
!       !
!       DO c = 1,nat3
!         IF(ABS(u3(c,k))>EPS)THEN
!         DO b = 1,nat3
!           IF(ABS(u2(b,j))>EPS)THEN
!           AUX = u2(b,j) * u3(c,k)
!           DO a = 1,nat3
!             !
!             d3tmp(i, j, k) = d3tmp(i, j, k) &
!                             + u1(a,i) * AUX * d3in(a, b, c) 
!           ENDDO
!           ENDIF
!         ENDDO
!         ENDIF
!       ENDDO
!         !
!     ENDDO
!     ENDDO
!     ENDDO
! !!!$OMP END PARALLEL DO

    !
    IF(present(d3out)) THEN
      d3out = d3tmp
    ELSE
      d3in  = d3tmp
    ENDIF
    DEALLOCATE(d3tmp)
    !
    RETURN
  END SUBROUTINE ip_cart2pat
  !
END MODULE fc2_interpolate


