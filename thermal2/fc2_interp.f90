
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE fc2_interpolate
  !
  USE kinds,    ONLY : DP
  USE nanoclock
  USE rigid, ONLY : rgd_blk
#include "mpi_thermal.h"  
  !
  ! I take forceconst2_grid globally in this module so I can USE it from here 
  ! elsewhere, which makes more sense
  USE input_fc, ONLY : forceconst2_grid, ph_system_info
  
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
    MODULE PROCEDURE fftinterp_mat2_flat_mkl   ! use standard OMP but paralelize on the innermost cycle./fc2_interp.f90:329:
#else
    MODULE PROCEDURE fftinterp_mat2_flat   ! use standard OMP but paralelize on the innermost cycle
#endif
!    MODULE PROCEDURE fftinterp_mat2_reduce  ! use standard OMP reduce on dummy variables
!    MODULE PROCEDURE fftinterp_mat2_safe    ! do not use dummy variable in OMP clauses
  END INTERFACE

  ! Only one version of the routine to interpolate the derivative of D
  ! is available, if you need the other one, write it..
  ! (i.e. the code will not work if __PRECOMPUTE_PHASES is undefined)
  INTERFACE fftinterp_dmat2
    MODULE PROCEDURE fftinterp_dmat2_flat_mkl   
  END INTERFACE

  CONTAINS
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_reduce(xq, S, fc, D, xq_hat)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in) :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(in),OPTIONAL :: xq_hat(3)
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
    !
    IF(S%lrigid) CALL add_rgd_blk(xq, S, fc, D, xq_hat)

  END SUBROUTINE fftinterp_mat2_reduce
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_flat(xq, S, fc, D, xq_hat)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in) :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(in),OPTIONAL :: xq_hat(3)
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
      DO mu= 1, S%nat3
      DO nu= 1, S%nat3
        D(nu,mu) = D(nu,mu) + phase * fc%fc(nu,mu, i)
      ENDDO
      ENDDO
!$OMP END DO
    END DO
!$OMP END PARALLEL
    !
    IF(S%lrigid) CALL add_rgd_blk(xq, S, fc, D, xq_hat)
    !
  END SUBROUTINE fftinterp_mat2_flat
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc
  SUBROUTINE fftinterp_mat2_flat_mkl(xq, S, fc, D, xq_hat)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq(3)
    TYPE(ph_system_info),INTENT(in) :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(in),OPTIONAL :: xq_hat(3)
    !
    INTEGER :: i, mu,nu
    REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
    COMPLEX(DP) :: vphase(fc%n_R)
    !
    D = (0._dp, 0._dp)
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
    FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq(:)*fc%xR(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R, varg, vcos)
    CALL vdSin(fc%n_R, varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    ! NOTE: do not use these OMP directive, as it is MUCH more effective to
    ! parallelize outside this subroutine!
!/!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,mu,nu) REDUCTION(+:D)
!/!$OMP DO
    DO i = 1, fc%n_R
      DO mu= 1, S%nat3
      DO nu= 1, S%nat3
        D(nu,mu) = D(nu,mu) + vphase(i) * fc%fc(nu,mu, i)
      ENDDO
      ENDDO
    END DO
!/!$OMP END DO
!/!$OMP END PARALLEL
    !
    IF(S%lrigid) CALL add_rgd_blk(xq, S, fc, D, xq_hat)
    !
  END SUBROUTINE fftinterp_mat2_flat_mkl
  !
  ! Compute the Dynamical matrix D by Fourier-interpolation of the
  ! force constants fc (safe version, for OpenMP)
  SUBROUTINE fftinterp_mat2_safe(xq, S, fc, D, xq_hat)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(out) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(in),OPTIONAL :: xq_hat(3)
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
!/!$OMP PARALLEL DEFAULT(none) SHARED(D,S,fc,xq) PRIVATE(i,arg,phase,D_aux)
    ALLOCATE(D_aux(S%nat3,S%nat3))
    D_aux = (0._dp, 0._dp)
!/!$OMP DO
    DO i = 1, fc%n_R
      arg = tpi * SUM(xq(:)*fc%xR(:,i))
      phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
      D_aux(:, :) = D_aux(:, :) + phase * fc%fc(:, :, i)
    END DO
!/!$OMP END DO
!/!$OMP CRITICAL(fc_interpolate_fftinterp_mat2_10)
    D = D + D_aux
!/!$OMP FLUSH(D)
!/!$OMP END CRITICAL(fc_interpolate_fftinterp_mat2_10)
    DEALLOCATE(D_aux)
!/!$OMP END PARALLEL
    !
    IF(S%lrigid) CALL add_rgd_blk(xq, S, fc, D, xq_hat)
    !
  END SUBROUTINE fftinterp_mat2_safe
  !
  SUBROUTINE add_rgd_blk(xq, S, fc, D, xq_hat)
    USE random_numbers,   ONLY : randy
    USE rigid,            ONLY : nonanal
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    REAL(DP),INTENT(in) :: xq(3)
    COMPLEX(DP),INTENT(inout) :: D(S%nat3, S%nat3)
    REAL(DP),INTENT(in),OPTIONAL :: xq_hat(3)
    !
    COMPLEX(DP) :: phi(3,3,S%nat,S%nat)
    INTEGER :: nu,mu,i,j,ia,ja
    !
    INTEGER :: itau(S%nat)
    REAL(DP) :: qhat(3), qnorm
    !
    !
    phi = 0._dp
    ! Non-analitical is only for q=0 and depends on the direction.
    ! If it is not present in input, it will be set to a random direction
    ! this can break symmetry in crystals where there is a non-analytical
    ! LT-TO splitting at Gamma
    IF(ALL(ABS(xq)<1.d-12))THEN
      IF(present(xq_hat)) THEN
        qhat = xq_hat
      ELSE
        qhat = 0._dp
        ! choose the direction x randomly, but we use always x
        ! qhat matters only in gamma.
        ! the group velocity is ill-defined in Gamma (since there is the nonanalytic term)
        ! the nonanalytic term is real, direction-dependent and is added only in Gamma.
        ! When we compute the group velocity, the nonanalytic term is not added because 
        ! with the finite differences we copute D(h)-D(-h) [with h=1E-7]
        ! thus the nonanalytic term plays a role only in the choice of the 
        ! eigenvectors at Gamma and is added when all the components of q are smaller than 1E-12 
        qhat(1) = 1.0 !randy()
        qhat(2) = 0.0 !randy()
        qhat(3) = 0.0 !randy()
      ENDIF
      qnorm = DSQRT(SUM(qhat**2))
      IF(qnorm>0._dp) qhat = qhat/qnorm
      FORALL(i=1:S%nat) itau(i) = i
      !
      CALL nonanal (S%nat, S%nat, itau, S%epsil, qhat, S%zeu, S%omega, phi)
    ENDIF
    !
    ! Add the long-range term rom effective charges
    CALL rgd_blk(fc%nq(1),fc%nq(2),fc%nq(3),S%nat,phi,xq, &
                  S%tau,S%epsil,S%zeu,S%bg,S%omega,S%alat,.false.,+1.d0)
    !
    DO ja = 1,S%nat
    DO ia = 1,S%nat
      DO j = 1,3
      mu = j + 3*(ja-1)
      DO i = 1,3
        nu = i + 3*(ia-1)
        !print*, mu,nu,i,j,ia,ja
        D(nu,mu) = D(nu,mu) + phi(i,j,ia,ja)&
                              *S%sqrtmm1(nu)*S%sqrtmm1(mu)
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    
    RETURN
  END SUBROUTINE
  !
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
      nb = ILAENV( 1, 'ZHEEV', 'U', n, -1, -1, -1 )
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
  ! Interpolate dynamical matrice at q and diagonalize it, 
  ! put acoustic gamma frequencies at exactly zero, 
  ! STOP if we get any negative frequencies
  SUBROUTINE freq_phq_safe(xq, S, fc2, freq, U)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants,          ONLY : RY_TO_CMM1
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: epsq = 1.e-8_dp
    REAL(DP) :: cq(3)
    INTEGER :: i
    LOGICAL :: gamma
    !
    ! RAF
    !U = CONJG(U)
    CALL fftinterp_mat2(xq, S, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    cq = xq
    CALL cryst_to_cart(1,cq,S%at,-1)
    gamma = ALL( ABS(cq-NINT(cq))<epsq)
    IF( gamma )THEN
      freq(1:3) = 0._dp
      U(:,1:3) = (0._dp, 0._dp)
    ENDIF
    
    IF(ANY(freq<0._dp)) THEN
      WRITE(*,*) gamma
      WRITE(*,"('cq = ',3f12.6)") cq
      WRITE(*,"('xq = ',3f12.6)") xq
      WHERE    (freq >  0.0)
        freq = DSQRT(freq)
      ELSEWHERE(freq < 0.0)
        freq = -DSQRT(-freq)
      ENDWHERE 
      WRITE(*,"('negative freq = ',12e12.4)")freq*RY_TO_CMM1
      CALL errore("freq_phq_safe", "cannot continue with negative frequencies",1)
    ELSE
    !
       freq = DSQRT(freq)
    END IF
    !
  END SUBROUTINE freq_phq_safe
  !
  ! Return 1 if q-point is not Gamma, 4 if it is (3 points in input)
  ! used to avoid divergence because of 1/freq at Gamma for acoustic bands
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
    REAL(DP),PARAMETER :: eps = 0._dp 
    !
    CALL fftinterp_mat2(xq, S, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    !U = CONJG(U)
    WHERE    (freq >  eps)
      freq = DSQRT(freq)
    ELSEWHERE(freq < -eps)
      freq = -DSQRT(-freq)
    ELSEWHERE ! i.e. freq=0
      freq = 0._dp
    ENDWHERE
      !
  END SUBROUTINE freq_phq
  !
  ! Interpolate dynamical matrice at q and diagonalize it
  SUBROUTINE freq_phq_positive(xq, S, fc2, freq, U)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    REAL(DP),INTENT(in)               :: xq(3)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: eps = 0._dp 
    !
    CALL fftinterp_mat2(xq, S, fc2, U)
    CALL mat2_diag(S%nat3, U, freq)
    !U = CONJG(U)
    WHERE    (freq >  eps)
      freq = DSQRT(freq)
    ELSEWHERE(freq < -eps)
      freq = DSQRT(-freq)
    ELSEWHERE ! i.e. freq=0
      freq = 0._dp
    ENDWHERE
      !
  END SUBROUTINE freq_phq_positive
  !
  !
  !
  ! Interpolate a point in a path, keep in account LO-TO at Gamma
  SUBROUTINE freq_phq_path(nq, iq, xq, S, fc2, freq, U)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    IMPLICIT NONE
    INTEGER,INTENT(in)                :: nq, iq
    REAL(DP),INTENT(in)               :: xq(3,nq)
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    REAL(DP),INTENT(out)              :: freq(S%nat3)
    COMPLEX(DP),INTENT(out)           :: U(S%nat3,S%nat3)
    REAL(DP),PARAMETER :: eps = 0._dp 
    REAL(DP) :: xq_hat(3), yq(3)
    !
    yq = xq(:,iq)
    CALL cryst_to_cart(1, yq, S%at, -1)
    yq = yq-DBLE(NINT(yq))
    CALL cryst_to_cart(1, yq, S%bg, +1)
    ! Check for exactly Gamma, even a tiny displacement works already
    IF(ALL(ABS(yq)< 1.d-8))THEN
      ! We set xq to EXACTLY zero in this case, to avoid double lo-to counting
      yq=0._dp
      IF(nq == 0)THEN
        ! this should never happen: zero length path?
        CALL errore("freq_phq_path","no points in this path",1)
      ELSE IF(nq == 1)THEN
        xq_hat = xq(:,1)
      ELSE IF(iq==1) THEN
        xq_hat = xq(:,2)-xq(:,1)
      ELSE
        xq_hat = xq(:,iq) - xq(:,iq-1) ! <- the most common case
      ENDIF
    ENDIF
    !
    ! WARNING: dirty hack!
    ! Instead of computing the non-analitical contribution (which seemed to be missing some
    ! piece that I did not manage to track down) I'm computing the frequencies at 
    ! a tiny displacement out of Gamma
    !CALL fftinterp_mat2(xq(:,iq)+xq_hat*1.d-6, S, fc2, U, xq_hat)
    ! attention! Here this routine is always called with the xq_hat argument!
    CALL fftinterp_mat2(yq, S, fc2, U, xq_hat)
    !
    CALL mat2_diag(S%nat3, U, freq)
    !U = CONJG(U)
    WHERE    (freq >  eps)
      freq = DSQRT(freq)
    ELSEWHERE(freq < -eps)
      freq = -DSQRT(-freq)
    ELSEWHERE ! i.e. freq=0
      freq = 0._dp
    ENDWHERE
      !
  END SUBROUTINE freq_phq_path
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
    IF(T==0._dp)THEN
     bose = 0._dp
     RETURN
    ENDIF

    WHERE    (freq > 0._dp)
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
  SUBROUTINE ip_cart2pat(d3in, nat3, u1, u2, u3)
    ! Rotates D3 matrix from cartesian axis to the basis
    ! of the modes. Rotation is not really in place
    USE kinds, ONLY : DP
    IMPLICIT NONE
    ! d3 matrix, input: in cartesian basis, output: on the patterns basis
    COMPLEX(DP),INTENT(inout) :: d3in(nat3, nat3, nat3)
    INTEGER,INTENT(in)        :: nat3
    ! patterns (transposed, with respect to what we use in the d3q.x code)
    COMPLEX(DP),INTENT(in)    :: u1(nat3, nat3), u2(nat3, nat3), u3(nat3, nat3) 
    !
    INTEGER :: a, b, c, i, j, k
    COMPLEX(DP),ALLOCATABLE  :: d3tmp(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: AUX(:,:)
    COMPLEX(DP),PARAMETER :: Z1 = (1._dp, 0._dp)
    COMPLEX(DP) :: u1t(nat3,nat3)
    !
    ALLOCATE(AUX(nat3,nat3))
    ALLOCATE(d3tmp(nat3, nat3, nat3))
    d3tmp = 0._dp
    !
    u1t = TRANSPOSE(CONJG(u1))
    !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,a,b,c,AUX) REDUCTION(+: d3tmp) COLLAPSE(2)
    DO c = 1,nat3
    DO b = 1,nat3
      ! Precompute u2*u3 to save some FLOPS without 
      ! compromising memory access order
      DO k = 1,nat3
      DO j = 1,nat3
        AUX(j,k) = CONJG(u2(b,j) * u3(c,k))
      ENDDO
      ENDDO
      !
      DO a = 1,nat3
        DO k = 1,nat3
        DO j = 1,nat3
        DO i = 1,nat3
              d3tmp(i, j, k) = d3tmp(i, j, k) &
                              + u1t(i,a) * AUX(j,k) * d3in(a, b, c) 
        ENDDO
        ENDDO
        ENDDO
      ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !
    d3in = d3tmp
    DEALLOCATE(d3tmp)
    !
    RETURN
    !
    RETURN
  END SUBROUTINE ip_cart2pat
  !
  ! \/o\________\\\_________________________________________/^>
  ! Compute the derivative of the Dynamical matrix D by 
  ! FT of force constants fc * iR
  SUBROUTINE fftinterp_dmat2_flat_mkl(xq, S, fc, dD)
    USE input_fc, ONLY : ph_system_info, forceconst2_grid
    USE constants, ONLY : tpi
    IMPLICIT NONE
    !
    REAL(DP),INTENT(in) :: xq(3)
    TYPE(ph_system_info),INTENT(in) :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    COMPLEX(DP),INTENT(out) :: dD(S%nat3, S%nat3,3)
    !
    INTEGER :: i, mu,nu, alpha
    REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
    COMPLEX(DP) :: vphase(fc%n_R), fac
    !
    dD = (0._dp, 0._dp)
    !
    ! Pre-compute phase to use the vectorized MKL subroutines
    FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq(:)*fc%xR(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
    CALL vdCos(fc%n_R, varg, vcos)
    CALL vdSin(fc%n_R, varg, vsin)
#else
    vcos = DCOS(varg)
    vsin = DSIN(varg)
#endif
    vphase =  CMPLX( vcos, -vsin, kind=DP  )
    !
    ! NOTE: do not use these OMP directive, as it is MUCH more effective to
    ! parallelize outside this subroutine!
    DO i = 1, fc%n_R
      DO alpha = 1,3
        fac = vphase(i)*CMPLX(0._dp, -fc%xR(alpha,i)*S%alat,kind=DP)
        DO mu= 1, S%nat3
        DO nu= 1, S%nat3
          dD(nu,mu,alpha) = dD(nu,mu,alpha)+ fac * fc%fc(nu,mu, i)
        ENDDO
        ENDDO
      ENDDO
    END DO
    !
    !IF(S%lrigid) CALL add_rgd_blk(xq, S, fc, D)
    !
  END SUBROUTINE fftinterp_dmat2_flat_mkl
  !
END MODULE fc2_interpolate
