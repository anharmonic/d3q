!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Module that uses object oriented features of Fortran 2003 to deal with both
! regular-grid and sparse representations of Force constants in a transparent way.
! May not compile with some less maintained compilers, should work at least with
! gfortran, xlf, ifort and pgi compiler. Does not compile with g95.
!
! It is not written defensively, e.g. writing and interpolating FCs do not check
! if the data has been read. May be improved in the future.
!
#define __PRECOMPUTE_PHASES
#ifdef __PRECOMPUTE_PHASES
!dir$ message "----------------------------------------------------------------------------------------------"
!dir$ message "Using MKL vectorized Sin and Cos implementation, this can use more memory but should be faster"
!dir$ message "----------------------------------------------------------------------------------------------"
#endif
!
MODULE fc3_interpolate
   USE kinds,            ONLY : DP
   USE parameters,       ONLY : ntypx
   USE mpi_thermal,      ONLY : ionode
#include "mpi_thermal.h"
   !USE input_fc,              ONLY : ph_system_info
   ! \/o\________\\\______________________//\/___________________/~^>>
  INTERFACE ip_cart2pat
!!#define __IP_WO_ZGEMM
#ifdef __IP_WO_ZGEMM
!dir$ message "----------------------------------------------------------------------------------------------"
!dir$ message "D3 matrix rotation: using fortran loop (slow but safe)"
!dir$ message "----------------------------------------------------------------------------------------------"
    MODULE PROCEDURE ip_cart2pat_byhand
#else
!dir$ message "----------------------------------------------------------------------------------------------"
!dir$ message "D3 matrix rotation: using zgemm (fast but can cause problem with buggy lapack)"
!dir$ message "----------------------------------------------------------------------------------------------"
    MODULE PROCEDURE ip_cart2pat_gemm
#endif
!    MODULE PROCEDURE fftinterp_mat2_reduce  ! use standard OMP reduce on dummy variables
!    MODULE PROCEDURE fftinterp_mat2_safe    ! do not use dummy variable in OMP clauses
  END INTERFACE ip_cart2pat

#if defined(__INTEL)
#define ZGEMM zgemm3m
#else
!dir$ message "----------------------------------------------------------------------------------------------"
!dir$ message "Not using ZGEMM3M: if you have this libray you can enable it in d3q/thermal2/fc3_interp.f90
!dir$ message "----------------------------------------------------------------------------------------------"
#endif



   ! Abstract implementation: all methods are abstract
   TYPE,ABSTRACT :: forceconst3
      ! q points
      INTEGER :: n_R = 0
      INTEGER :: i_00 = -1 ! set to the index where both R2 and R3 are zero
      INTEGER :: nq(3) ! initial q-point grid size (unused)
      INTEGER,ALLOCATABLE  :: yR2(:,:), yR3(:,:) ! crystalline coords  3*n_R
      REAL(DP),ALLOCATABLE :: xR2(:,:), xR3(:,:) ! cartesian coords    3*n_R
   CONTAINS
      procedure(fft_interp_sum_R2), DEFERRED :: sum_R2
      procedure(fftinterp_mat3_error),deferred :: interpolate
      !> gradient obtained in reciprocal space, only for sparse, NOT TESTED
      procedure(fftinterp_mat3_grad), deferred :: interpolate_grad
      procedure(fft_doubleinterp_mat3_error),deferred :: double_interpolate
      procedure :: destroy      => destroy_fc3_
      procedure(read_fc3_error),deferred :: read
      procedure(write_fc3_error),deferred :: write
      procedure(div_mass_fc3_error),deferred :: div_mass
   END TYPE forceconst3
   !
   TYPE :: d3_mixed
      LOGICAL, ALLOCATABLE :: R3(:)
      COMPLEX(DP), ALLOCATABLE :: DR3(:,:,:,:)
      INTEGER :: minr(3), maxr(3)
   END TYPE
   ! Interfaces to the deferred subroutines:
   ABSTRACT INTERFACE
      SUBROUTINE fft_interp_sum_R2(fc, xq, nat3, Dqr)
         USE kinds,     ONLY : DP
         IMPORT forceconst3
         IMPORT d3_mixed
         CLASS(forceconst3), INTENT(IN) :: fc
         REAL(DP), INTENT(IN) :: xq(3)
         INTEGER, INTENT(IN) :: nat3
         TYPE(d3_mixed), INTENT(OUT)  :: Dqr
      END SUBROUTINE
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE fftinterp_mat3_error(fc, xq2,xq3, nat3, D)
         USE kinds,     ONLY : DP
         IMPORT forceconst3
         CLASS(forceconst3),INTENT(in) :: fc
         INTEGER,INTENT(in)   :: nat3
         REAL(DP),INTENT(in) :: xq2(3), xq3(3)
         COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      END SUBROUTINE fftinterp_mat3_error
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE fftinterp_mat3_grad(fc, xq2,xq3, nat3, D, Dgrad)
         USE kinds,     ONLY : DP
         IMPORT forceconst3
         CLASS(forceconst3),INTENT(in) :: fc
         INTEGER,INTENT(in)   :: nat3
         REAL(DP),INTENT(in) :: xq2(3), xq3(3)
         COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
         COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)
      END SUBROUTINE fftinterp_mat3_grad
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE fft_doubleinterp_mat3_error(fc, xq2,xq3,xq3b, nat3, D, Db)
         USE kinds,     ONLY : DP
         IMPORT forceconst3
         CLASS(forceconst3),INTENT(in) :: fc
         INTEGER,INTENT(in)   :: nat3
         REAL(DP),INTENT(in) :: xq2(3), xq3(3), xq3b(3)
         COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
         COMPLEX(DP),INTENT(out) :: Db(nat3, nat3, nat3)
      END SUBROUTINE fft_doubleinterp_mat3_error
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE read_fc3_error(fc, filename, S)
         USE input_fc,              ONLY : ph_system_info
         IMPORT forceconst3
         CHARACTER(*),INTENT(in)        :: filename
         TYPE(ph_system_info),INTENT(inout) :: S ! = System
         CLASS(forceconst3),INTENT(inout)   :: fc
      END SUBROUTINE read_fc3_error
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE write_fc3_error(fc, filename, S)
         USE input_fc,              ONLY : ph_system_info
         IMPORT forceconst3
         CHARACTER(*),INTENT(in)     :: filename
         TYPE(ph_system_info),INTENT(in) :: S ! = System
         CLASS(forceconst3),INTENT(in)   :: fc
      END SUBROUTINE write_fc3_error
   END INTERFACE
   !
   ABSTRACT INTERFACE
      SUBROUTINE div_mass_fc3_error(fc, S)
         USE input_fc,              ONLY : ph_system_info
         IMPORT forceconst3
         CLASS(forceconst3),INTENT(inout) :: fc
         TYPE(ph_system_info),INTENT(in)  :: S ! = System
      END SUBROUTINE div_mass_fc3_error
   END INTERFACE
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   !
   ! First implementation: not regular double grid
   TYPE,EXTENDS(forceconst3) :: grid
      REAL(DP),ALLOCATABLE :: FC(:,:,:,:) ! 3*nat,3*nat,3*nat, n_R
      ! Optional, imaginary part of the force constants, for testing only (should be zero)
      REAL(DP),ALLOCATABLE :: iFC(:,:,:,:) ! 3*nat,3*nat,3*nat, n_R
      !
   CONTAINS
      ! There are 2 version of the interpolation routine, that differ on the
      ! parallelism model via OMP, in principle they are both correct but only
      ! _flat works, because of some unclear problem with the reduction of arrays
      !procedure :: interpolate  => fftinterp_mat3_grid_reduce
#ifdef __PRECOMPUTE_PHASES
      procedure :: interpolate  => fftinterp_mat3_grid_flat_prec
      procedure :: interpolate_grad  => fftinterp_mat3_grid_flat_prec_grad
#else
!DEC$ SUPPRESS ALL
      procedure :: interpolate  => fftinterp_mat3_grid_flat
      procedure :: interpolate_grad  => fftinterp_mat3_grid_flat_grad
!DEC$ SUPPRESS NONE
#endif
      procedure :: double_interpolate  => fft_doubleinterp_mat3_grid_flat
      !
      procedure :: sum_R2       => todo_sum_R2
      procedure :: destroy      => destroy_fc3_grid
      procedure :: read         => read_fc3_grid
      procedure :: write        => write_fc3_grid
      procedure :: div_mass     => div_mass_fc3_grid
   END TYPE grid
   INTERFACE grid
      MODULE PROCEDURE create_fc3_grid
   END INTERFACE
   !
   ! \/o\________\\\_________________________________________/^>
   !
   ! Second implementation: sparse grid
   TYPE forceconst3_sparse_helper
      REAL(DP),ALLOCATABLE :: fc(:)    ! the matrix elements
      INTEGER,ALLOCATABLE  :: idx(:,:) ! the index 3*nat,3*nat,3*nat
   END TYPE forceconst3_sparse_helper
   !
   TYPE,EXTENDS(forceconst3) :: sparse
      INTEGER,ALLOCATABLE  :: n_terms(:) ! number of non-zero elements for each R2,R3
      TYPE(forceconst3_sparse_helper),ALLOCATABLE :: dat(:)
      !
   CONTAINS
#ifdef __PRECOMPUTE_PHASES
      procedure :: interpolate  => fftinterp_mat3_sparse_prec
      procedure :: interpolate_grad  => fftinterp_mat3_sparse_prec_grad

#else
      procedure :: interpolate  => fftinterp_mat3_sparse
      procedure :: interpolate_grad  => fftinterp_mat3_sparse_grad

#endif
      procedure :: sum_R2 => sum_R2_sparse
      procedure :: double_interpolate  => fft_doubleinterp_mat3_sparse_prec
      procedure :: destroy      => destroy_fc3_sparse
      procedure :: read         => read_fc3_sparse
      procedure :: write        => write_fc3_sparse
      procedure :: div_mass     => div_mass_fc3_sparse
   END TYPE sparse
   INTERFACE sparse
      MODULE PROCEDURE create_fc3_sparse
   END INTERFACE
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   !
   ! Third implementation: always return a constant (divided by the masses?)
   TYPE,EXTENDS(forceconst3) :: constant
      REAL(DP) :: constant = 0._dp
      COMPLEX(DP),ALLOCATABLE :: D(:,:,:)
      !
   CONTAINS
      !
      procedure :: sum_R2       => todo_sum_R2_const
      procedure :: interpolate  => fftinterp_mat3_constant
      procedure :: interpolate_grad => fftinterp_mat3_constant_grad
      procedure :: double_interpolate  => fft_doubleinterp_mat3_constant
      procedure :: destroy      => destroy_fc3_constant
      procedure :: read         => read_fc3_constant
      procedure :: write        => write_fc3_constant
      procedure :: div_mass     => div_mass_fc3_constant
   END TYPE constant
   INTERFACE constant
      MODULE PROCEDURE create_fc3_constant
   END INTERFACE
   ! \/o\________\\\_________________________________________/^>
CONTAINS
   ! \/o\________\\\_________________________________________/^>
   !
   ! Returns a pointer to FC matrices, call as:
   !    TYPE(ph_system_info),INTENT(in)   :: S
   !    CLASS(forceconst3),POINTER :: fc
   !     ...
   !    fc => read_fc3("file_mat3", S)
   ! This will read both sparse and grid files, there is room for improvement but not much.
   !
   SUBROUTINE todo_sum_R2(fc, xq, nat3, Dqr)
      ! USE kinds,     ONLY : DP
      ! IMPORT forceconst3
      CLASS(grid), INTENT(IN)               :: fc
      REAL(DP), INTENT(IN)                  :: xq(3)
      INTEGER, INTENT(IN)                   :: nat3
      TYPE(d3_mixed), INTENT(OUT)           :: Dqr
      ! TODO
      CALL errore('fc3','not implemented case',1)
   END SUBROUTINE

   SUBROUTINE todo_sum_R2_const(fc, xq, nat3, Dqr)
      ! USE kinds,     ONLY : DP
      ! IMPORT forceconst3
      CLASS(constant), INTENT(IN)               :: fc
      REAL(DP), INTENT(IN)                  :: xq(3)
      INTEGER, INTENT(IN)                   :: nat3
      TYPE(d3_mixed), INTENT(OUT)           :: Dqr
      CALL errore('fc3','not implemented case',2)
      ! TODO
   END SUBROUTINE

   FUNCTION read_fc3(filename, S) RESULT(fc)
      USE input_fc, ONLY : read_system, ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(out)   :: S ! = System
      CLASS(forceconst3),POINTER :: fc
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: unit, ios
      CHARACTER(32) :: buf
      !
      ! Scan the type of file:
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
      IF(ios/=0) CALL errore("read_fc3","opening '"//TRIM(filename)//"'", 1)
      READ(unit, *) buf
      CLOSE(unit)
      !
      IF(buf=="sparse") THEN
         ALLOCATE(sparse:: fc)
      ELSE IF(buf=="constant") THEN
         ALLOCATE(constant::   fc)
      ELSE
         ALLOCATE(grid::   fc)
      ENDIF
      !
      ! Use class specific read method:
      CALL fc%read(filename, S)
      !
      RETURN
      !
   END FUNCTION read_fc3
   ! \/o\________\\\______________________//\/___________________/~^>>
   ! Placeholder creators:
   TYPE(grid) FUNCTION create_fc3_grid()
      RETURN
   END FUNCTION
   TYPE(sparse) FUNCTION create_fc3_sparse()
      RETURN
   END FUNCTION
   TYPE(constant) FUNCTION create_fc3_constant()
      RETURN
   END FUNCTION
   !
   ! Convert regular-grid FCs to sparse
   ! \/o\________\\\_________________________________________/^>
   SUBROUTINE fc3_grid_to_sparse(nat,fc, sfc, thr)
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(in) :: nat
      TYPE(grid),INTENT(in)    :: fc
      TYPE(sparse),INTENT(out) :: sfc
      REAL(DP),OPTIONAL,INTENT(in) :: thr
      !
      REAL(DP),PARAMETER :: eps12 = 1.d-12
      INTEGER :: iR, n_found
      INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
      REAL(DP) :: eps
      !
      eps = eps12
      IF(present(thr)) eps = thr
      !
      sfc%n_R  = fc%n_R
      sfc%i_00 = fc%i_00
      sfc%nq   = fc%nq
      !
      ALLOCATE(sfc%yR2(3,sfc%n_R), sfc%yR3(3,sfc%n_R))
      ALLOCATE(sfc%xR2(3,sfc%n_R), sfc%xR3(3,sfc%n_R))
      sfc%yR2 = fc%yR2;    sfc%yR3 = fc%yR3
      sfc%xR2 = fc%xR2;    sfc%xR3 = fc%xR3
      ALLOCATE(sfc%dat(sfc%n_R))
      ALLOCATE(sfc%n_terms(sfc%n_R))

      DO iR = 1, sfc%n_R
         sfc%n_terms(iR) = COUNT(  ABS(fc%FC(:,:,:,iR))>eps )
         ALLOCATE(sfc%dat(iR)%idx(3,sfc%n_terms(iR)))
         ALLOCATE(sfc%dat(iR)%fc(sfc%n_terms(iR)))
      ENDDO

      DO iR = 1, sfc%n_R
         !
         n_found = 0
         !
         DO na3=1,nat
            DO na2=1,nat
               DO na1=1,nat
                  DO j3=1,3
                     jn3 = j3 + (na3-1)*3
                     DO j2=1,3
                        jn2 = j2 + (na2-1)*3
                        DO j1=1,3
                           jn1 = j1 + (na1-1)*3
                           !
                           IF(ABS(fc%FC(jn1,jn2,jn3,iR))>eps)THEN
                              n_found = n_found + 1
                              IF(n_found > sfc%n_terms(iR)) CALL errore("fc3_grid_to_sparse", "too many terms", iR)
                              sfc%dat(iR)%idx(1,n_found) = jn1
                              sfc%dat(iR)%idx(2,n_found) = jn2
                              sfc%dat(iR)%idx(3,n_found) = jn3
                              sfc%dat(iR)%FC(n_found) = fc%FC(jn1,jn2,jn3,iR)
                           ENDIF
                           !
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !
         IF(n_found /= sfc%n_terms(iR)) CALL errore("fc3_grid_to_sparse", "wrong number of terms", iR)
         !
      ENDDO
      !
   END SUBROUTINE fc3_grid_to_sparse
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_grid_reduce(fc, xq2,xq3, nat3, D)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      REAL(DP) :: arg
      COMPLEX(DP) :: phase
      INTEGER :: i, a,b,c
      !
      D = (0._dp, 0._dp)
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,arg,phase) REDUCTION(+: D)
      DO i = 1, fc%n_R
         arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
         D(:,:,:) = D(:,:,:) + phase * fc%fc(:,:,:, i)
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE fftinterp_mat3_grid_reduce
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_grid_flat(fc, xq2,xq3, nat3, D)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      REAL(DP) :: arg
      COMPLEX(DP) :: phase
      INTEGER :: i, a,b,c
      !
      D = (0._dp, 0._dp)
      !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
      DO i = 1, fc%n_R
         arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(3)
         DO c = 1,nat3
            DO b = 1,nat3
               DO a = 1,nat3
                  D(a,b,c) = D(a,b,c) + phase * fc%fc(a,b,c, i)
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO
      END DO
!$OMP END PARALLEL
   END SUBROUTINE fftinterp_mat3_grid_flat
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_grid_flat_prec(fc, xq2,xq3, nat3, D)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      INTEGER :: i, a,b,c
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R)
      !
      D = (0._dp, 0._dp)
      !
      ! Pre-compute phase, use vectorized MKL subroutines if available
      FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
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
! ==================== with imaginary part, only for testing during qq2rr
      IMAGINARY: &
         IF(allocated(fc%ifc))THEN
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
         DO i = 1, fc%n_R
            !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(3)
            DO c = 1,nat3
               DO b = 1,nat3
                  DO a = 1,nat3
                     D(a,b,c) = D(a,b,c) + vphase(i) * CMPLX(fc%fc(a,b,c, i), fc%ifc(a,b,c, i), kind=DP)
                  ENDDO
               ENDDO
            ENDDO
!$OMP END DO
         END DO
!$OMP END PARALLEL
! ==================== without imaginary part
      ELSE IMAGINARY
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
         DO i = 1, fc%n_R
            !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(3)
            DO c = 1,nat3
               DO b = 1,nat3
                  DO a = 1,nat3
                     D(a,b,c) = D(a,b,c) + vphase(i) * fc%fc(a,b,c, i)
                  ENDDO
               ENDDO
            ENDDO
!$OMP END DO
         END DO
!$OMP END PARALLEL
      ENDIF &
         IMAGINARY
! ==================== without imaginary part
   END SUBROUTINE fftinterp_mat3_grid_flat_prec

   SUBROUTINE fftinterp_mat3_grid_flat_grad(fc, xq2,xq3, nat3, D, Dgrad)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)

      !
      REAL(DP) :: arg, cosine, sine
      INTEGER :: i, a,b,c
      !
      D = (0._dp, 0._dp)
      !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
      DO i = 1, fc%n_R
         arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         cosine = Cos(arg)
         sine = Sin(arg)
!$OMP DO COLLAPSE(3)
         DO c = 1,nat3
            DO b = 1,nat3
               DO a = 1,nat3
                  D(a,b,c) = D(a,b,c) + CMPLX(cosine, -sine, kind=DP) * fc%fc(a,b,c, i)
                  Dgrad(a,b,c,:) = Dgrad(a,b,c,:) + (fc%xR2(:,i) - fc%xR3(:,i)) * &
                     CMPLX(sine, cosine, kind=DP) * fc%fc(a,b,c, i)
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO
      END DO
!$OMP END PARALLEL
   END SUBROUTINE fftinterp_mat3_grid_flat_grad
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_grid_flat_prec_grad(fc, xq2,xq3, nat3, D, Dgrad)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)

      !
      INTEGER :: i, a,b,c
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R), vphaseR(fc%n_R)
      !
      D = (0._dp, 0._dp)
      !
      ! Pre-compute phase, use vectorized MKL subroutines if available
      FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R, varg, vcos)
      CALL vdSin(fc%n_R, varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      vphaseR = CMPLX( vsin, vcos, kind=DP )
      !
! ==================== with imaginary part, only for testing during qq2rr
      IMAGINARY: &
         IF(allocated(fc%ifc))THEN
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
         DO i = 1, fc%n_R
            !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(3)
            DO c = 1,nat3
               DO b = 1,nat3
                  DO a = 1,nat3
                     D(a,b,c) = D(a,b,c) + vphase(i) * CMPLX(fc%fc(a,b,c, i), fc%ifc(a,b,c, i), kind=DP)
                  ENDDO
               ENDDO
            ENDDO
!$OMP END DO
         END DO
!$OMP END PARALLEL
! ==================== without imaginary part
      ELSE IMAGINARY
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
         DO i = 1, fc%n_R
            !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
!$OMP DO COLLAPSE(3)
            DO c = 1,nat3
               DO b = 1,nat3
                  DO a = 1,nat3
                     D(a,b,c) = D(a,b,c) + vphase(i) * fc%fc(a,b,c, i)
                     Dgrad(a,b,c,:) = Dgrad(a,b,c,:) + (fc%xR2(:,i) - fc%xR3(:,i)) * vphaseR(i) * fc%fc(a,b,c, i)
                  ENDDO
               ENDDO
            ENDDO
!$OMP END DO
         END DO
!$OMP END PARALLEL
      ENDIF &
         IMAGINARY
! ==================== without imaginary part
   END SUBROUTINE fftinterp_mat3_grid_flat_prec_grad
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fft_doubleinterp_mat3_grid_flat(fc, xq2,xq3,xq3b, nat3, D, Db)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(grid),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3), xq3b(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Db(nat3, nat3, nat3)
      !
      REAL(DP) :: arg1, arg2
      COMPLEX(DP) :: phase1, phase2
      INTEGER :: i, a,b,c
      !
      D  = (0._dp, 0._dp)
      Db = (0._dp, 0._dp)
      !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a,b,c)
      DO i = 1, fc%n_R
         arg1 = tpi *  SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         phase1 = CMPLX(Cos(arg1),-Sin(arg1), kind=DP)
         arg2 = tpi * SUM(-xq2(:)*fc%xR2(:,i) + xq3b(:)*fc%xR3(:,i))
         phase2 = CMPLX(Cos(arg2),-Sin(arg2), kind=DP)
!$OMP DO COLLAPSE(3)
         DO c = 1,nat3
            DO b = 1,nat3
               DO a = 1,nat3
                  D(a,b,c)  = D(a,b,c)  + phase1 * fc%fc(a,b,c, i)
                  Db(a,b,c) = Db(a,b,c) + phase2 * fc%fc(a,b,c, i)
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO
      END DO
!$OMP END PARALLEL
   END SUBROUTINE fft_doubleinterp_mat3_grid_flat
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_sparse(fc, xq2,xq3, nat3, D)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(sparse),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      REAL(DP) :: arg
      COMPLEX(DP) :: phase
      INTEGER :: i,j
      !
      D = (0._dp, 0._dp)
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,arg,phase) REDUCTION(+: D)
      DO i = 1, fc%n_R
         arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
         DO j = 1, fc%n_terms(i)
            D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + phase * fc%dat(i)%fc(j)
         ENDDO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE fftinterp_mat3_sparse
   !

   SUBROUTINE sum_R2_sparse(fc, xq, nat3, Dqr)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      CLASS(sparse), INTENT(IN)             :: fc
      REAL(DP), INTENT(IN)                  :: xq(3)
      INTEGER, INTENT(IN)                   :: nat3
      TYPE(d3_mixed), INTENT(OUT)           :: Dqr
      !
      INTEGER :: i,j, index3(3), iR3, nr, limits(3)
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R)
      COMPLEX(DP), ALLOCATABLE :: TMP(:,:,:,:)

      !
      DO i = 1, 3
         Dqr%minr(i) = MINVAL(fc%yR2(i,:))
         Dqr%maxr(i) = MAXVAL(fc%yR2(i,:))
      ENDDO
      limits = Dqr%maxr - Dqr%minr + 1
      nr = PRODUCT(limits)
      ALLOCATE(Dqr%R3(nr))
      ALLOCATE(TMP(nat3, nat3, nat3, nr))
      TMP = (0._dp, 0._dp)
      Dqr%R3 = .false.
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
      FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq(:)*fc%xR2(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
      !dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R, varg, vcos)
      CALL vdSin(fc%n_R, varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) REDUCTION(+: TMP)
      !/!$ACC DATA COPYIN(fc%dat, fc%idx, vphase) COPY(TMP)
      DO i = 1, fc%n_R
         !arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
         DO j = 1, fc%n_terms(i)
            index3 = fc%yR3(:,i) - Dqr%minr
            iR3 = index3(1)*limits(2)*limits(3) + index3(2)*limits(3) + index3(3) + 1
            Dqr%R3(iR3) = .true.
            TMP(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),iR3) &
               = TMP(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),iR3) &
               + vphase(i) * fc%dat(i)%fc(j)
         ENDDO
      END DO
      !/!$ACC END DATA
      !$OMP END PARALLEL DO

      ALLOCATE(Dqr%DR3(nat3, nat3, nat3, nr))
      Dqr%DR3 = TMP

      DEALLOCATE(TMP)
   END SUBROUTINE

   SUBROUTINE sum_R3(S, xq, Dqr, D)
      use constants, only : tpi
      use ph_system, only : ph_system_info
      TYPE(ph_system_info) :: S
      REAL(DP), INTENT(IN) :: xq(3)
      TYPE(d3_mixed), INTENT(IN)           :: Dqr
      COMPLEX(DP),INTENT(OUT) :: D(S%nat3, S%nat3, S%nat3)
      !
      INTEGER :: irx, iry, irz, ir
      REAL(DP) :: xR3(3), arg
      EXTERNAL :: cryst_to_cart
      !
      D = (0._dp, 0._dp)

      ir = 0
      DO irx = Dqr%minr(1), Dqr%maxr(1)
         DO iry = Dqr%minr(2), Dqr%maxr(2)
            DO irz = Dqr%minr(3), Dqr%maxr(3)
               ir = ir + 1
               IF( .not. Dqr%R3(ir)) CYCLE
               xR3 = REAL((/irx, iry, irz/), kind=DP)
               CALL cryst_to_cart(1, xR3, S%at, 1)
               arg = tpi * DOT_PRODUCT(xq, xR3)
               D = D + Dqr%DR3(:,:,:,ir) * CMPLX( COS(arg), -SIN(arg), kind=DP  )
            END DO
         END DO
      END DO
   END SUBROUTINE

   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_sparse_prec(fc, xq2,xq3, nat3, D)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(sparse),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      INTEGER :: i,j
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R)
      !
      D = (0._dp, 0._dp)
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
      FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R, varg, vcos)
      CALL vdSin(fc%n_R, varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      !    !
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) REDUCTION(+: D)
      !/!$ACC DATA COPYIN(fc%dat, fc%idx, vphase) COPY(D)
      DO i = 1, fc%n_R
         !arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
         DO j = 1, fc%n_terms(i)
            D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + vphase(i) * fc%dat(i)%fc(j)
         ENDDO
      END DO
      !/!$ACC END DATA
      !$OMP END PARALLEL DO
   END SUBROUTINE fftinterp_mat3_sparse_prec

   SUBROUTINE fftinterp_mat3_sparse_grad(fc, xq2,xq3, nat3, D, Dgrad)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(sparse),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)

      !
      REAL(DP) :: arg, cosine, sine
      INTEGER :: i,j
      !
      D = (0._dp, 0._dp)
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,arg) REDUCTION(+: D)
      DO i = 1, fc%n_R
         arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         cosine = Cos(arg)
         sine = Sin(arg)
         DO j = 1, fc%n_terms(i)
            D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + CMPLX(cosine, -sine, kind=DP) * fc%dat(i)%fc(j)
            Dgrad(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),:) &
               = Dgrad(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),:) &
               + (fc%xR2(:,i) - fc%xR3(:,i)) * CMPLX(sine,  cosine, kind=DP) * fc%dat(i)%fc(j)
         ENDDO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE fftinterp_mat3_sparse_grad
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_sparse_prec_grad(fc, xq2,xq3, nat3, D, Dgrad)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(sparse),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)


      !
      INTEGER :: i,j
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R), vphaseR(fc%n_R)
      !
      D = (0._dp, 0._dp)
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
      FORALL(i=1:fc%n_R) varg(i) =  tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R, varg, vcos)
      CALL vdSin(fc%n_R, varg, vsin)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      vphaseR = CMPLX( vsin, vcos, kind=DP  )
      !    !
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) REDUCTION(+: D)
      !/!$ACC DATA COPYIN(fc%dat, fc%idx, vphase) COPY(D)
      DO i = 1, fc%n_R
         !arg = tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
         !phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
         DO j = 1, fc%n_terms(i)
            D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + vphase(i) * fc%dat(i)%fc(j)
            Dgrad(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),:) &
               = Dgrad(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j),:) &
               + (fc%xR2(:,i) - fc%xR3(:,i)) * vphaseR(i) * fc%dat(i)%fc(j)
         ENDDO
      END DO
      !/!$ACC END DATA
      !$OMP END PARALLEL DO
   END SUBROUTINE fftinterp_mat3_sparse_prec_grad
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fft_doubleinterp_mat3_sparse_prec(fc, xq2,xq3,xq3b, nat3, D, Db)
      USE constants, ONLY : tpi
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(sparse),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3), xq3b(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Db(nat3, nat3, nat3)
      !
      INTEGER :: i,j
      REAL(DP) :: varg(fc%n_R), vcos(fc%n_R), vsin(fc%n_R)
      REAL(DP) :: vargb(fc%n_R), vcosb(fc%n_R), vsinb(fc%n_R)
      COMPLEX(DP) :: vphase(fc%n_R)
      COMPLEX(DP) :: vphaseb(fc%n_R)
      !
      D  = (0._dp, 0._dp)
      Db = (0._dp, 0._dp)
      !
      ! Pre-compute phase to use the vectorized MKL subroutines
      FORALL(i=1:fc%n_R) varg(i)  =  tpi * SUM(xq2(:)*fc%xR2(:,i) + xq3(:)*fc%xR3(:,i))
      FORALL(i=1:fc%n_R) vargb(i) =  tpi * SUM(-xq2(:)*fc%xR2(:,i) + xq3b(:)*fc%xR3(:,i))
#if defined(__INTEL) && defined(__HASVTRIG)
!dir$ message "Using MKL vectorized Sin and Cos implementation, if this does not compile, remove -D__HASVTRIG from Makefile"
      CALL vdCos(fc%n_R, varg, vcos)
      CALL vdSin(fc%n_R, varg, vsin)
      CALL vdCos(fc%n_R, vargb, vcosb)
      CALL vdSin(fc%n_R, vargb, vsinb)
#else
      vcos = DCOS(varg)
      vsin = DSIN(varg)
      vcosb = DCOS(vargb)
      vsinb = DSIN(vargb)
#endif
      vphase =  CMPLX( vcos, -vsin, kind=DP  )
      vphaseb =  CMPLX( vcos, -vsin, kind=DP  )
      !    !
!
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) REDUCTION(+: D)
      DO i = 1, fc%n_R
         DO j = 1, fc%n_terms(i)
            D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = D(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + vphase(i) * fc%dat(i)%fc(j)

            Db(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               = Db(fc%dat(i)%idx(1,j),fc%dat(i)%idx(2,j),fc%dat(i)%idx(3,j)) &
               + vphaseb(i) * fc%dat(i)%fc(j)
         ENDDO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE fft_doubleinterp_mat3_sparse_prec
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE read_fc3_sparse(fc, filename, S)
      USE input_fc, ONLY : read_system, ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(inout)   :: S ! = System
      CLASS(sparse),INTENT(inout) :: fc
      !
      CHARACTER(13),PARAMETER :: sub = "read_fc3_sparse"
      !
      INTEGER :: unit, ios
      INTEGER, EXTERNAL :: find_free_unit
      CHARACTER(32) :: buf
      CHARACTER(6) :: dummy
      REAL(DP) :: factor
      !
      INTEGER :: i, j
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      READ(unit, '(a32)') buf
      IF(buf(1:6)/="sparse") CALL errore(sub, "cannot read this format", 1)
      READ(buf, *, iostat=ios) dummy, factor
      IF(ios/=0) factor = 1._dp
      IF(factor/=1._dp .and. ionode) WRITE(stdout, *) "WARNING! rescaling D3", factor
      !
      ioWRITE(stdout,*) "** Reading sparse FC3 file ", TRIM(filename)
      CALL read_system(unit, S)
      !
      READ(unit, *) fc%nq
      READ(unit, *) fc%n_R
      ioWRITE(stdout,*) "   Original FC3 grid:", fc%nq
      ioWRITE(stdout,*) "   Number of R:      ", fc%n_R
      !
      ALLOCATE(fc%n_terms(fc%n_R))
      ALLOCATE(fc%yR2(3,fc%n_R), fc%yR3(3,fc%n_R))
      DO i = 1,fc%n_R
         READ(unit, *) fc%n_terms(i),  fc%yR2(:,i), fc%yR3(:,i)
      ENDDO
      !
      ALLOCATE(fc%dat(fc%n_R))
      !
      DO i = 1, fc%n_R
         !
         ALLOCATE(fc%dat(i)%idx(3,fc%n_terms(i)))
         ALLOCATE(fc%dat(i)%fc(fc%n_terms(i)))
         !
         DO j = 1, fc%n_terms(i)
            READ(unit,*) fc%dat(i)%idx(:,j), fc%dat(i)%fc(j)
         ENDDO
         fc%dat(i)%fc = fc%dat(i)%fc*factor
      ENDDO
      !
      CLOSE(unit)
      ! Prepare the lists of R in cartesian coords
      ALLOCATE(fc%xR2(3,fc%n_R), fc%xR3(3,fc%n_R))
      fc%xR2 = REAL(fc%yR2, kind=DP)
      CALL cryst_to_cart(fc%n_R,fc%xR2,S%at,1)
      fc%xR3 = REAL(fc%yR3, kind=DP)
      CALL cryst_to_cart(fc%n_R,fc%xR3,S%at,1)
      !
   END SUBROUTINE read_fc3_sparse
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE write_fc3_sparse(fc, filename, S)
      USE input_fc, ONLY : write_system, ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(in)   :: S ! = System
      CLASS(sparse),INTENT(in) :: fc
      !
      CHARACTER(14),PARAMETER :: sub = "write_fc3_sparse"
      !
      INTEGER :: unit, ios
      INTEGER, EXTERNAL :: find_free_unit
      CHARACTER (6), EXTERNAL :: int_to_char
      !
      INTEGER :: i, j, n_digits(3)
      CHARACTER(64) :: cformat
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      WRITE(unit, '(a)') "sparse representation"
      !
      CALL write_system(unit, S)
      WRITE(unit, '(3i9)') fc%nq
      WRITE(unit, '(3i9)') fc%n_R
      !
      ! Writing with the minimum number of white spaces saves loads of disk space (unless you zip)
      ! still keeping everything aligned for readability
      IF(MAXVAL(fc%n_terms)>0)THEN
         n_digits(1) = CEILING(LOG10(DBLE(MAXVAL(ABS(fc%n_terms))))) ! no need to leave space in front
      ELSE
         n_digits(1) = 1
      ENDIF
      n_digits(2) = CEILING(LOG10(DBLE(MAXVAL(ABS(fc%yR2)))+0.1_dp))+2 ! +2 to account minus signs
      n_digits(3) = CEILING(LOG10(DBLE(MAXVAL(ABS(fc%yR3)))+0.1_dp))+2
      cformat = '(i'//int_to_char(n_digits(1))// &
         ',3i'//int_to_char(n_digits(2))// &
         ',3i'//int_to_char(n_digits(3))//')'
      !cformat = '(i9,3i6,3i6)' ! this also works
      DO i = 1,fc%n_R
         WRITE(unit, cformat) fc%n_terms(i),  fc%yR2(:,i), fc%yR3(:,i)
      ENDDO
      !
      n_digits(1) = CEILING(LOG10(DBLE(3*S%nat)+0.1_dp))
      n_digits(2) = n_digits(1)+1
      !
      cformat = '(i'//int_to_char(n_digits(1))// &
         ',2i'//int_to_char(n_digits(2))//',1pe25.15)'
      !cformat = '(i9,2i6,1pe25.14)' ! this also works
      !
      DO i = 1, fc%n_R
         DO j = 1, fc%n_terms(i)
            !
            WRITE(unit,cformat) fc%dat(i)%idx(1:3,j), fc%dat(i)%fc(j)
            !
         ENDDO
         WRITE(unit,*)
      ENDDO
      !
      CLOSE(unit)
      !
   END SUBROUTINE write_fc3_sparse
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE destroy_fc3_sparse(fc)
      USE input_fc, ONLY : write_system, ph_system_info
      IMPLICIT NONE
      CLASS(sparse),INTENT(inout) :: fc
      INTEGER :: i
      !
      DEALLOCATE(fc%yR2)
      DEALLOCATE(fc%yR3)
      DEALLOCATE(fc%xR2)
      DEALLOCATE(fc%xR3)
      DO i = 1, fc%n_R
         !
         DEALLOCATE(fc%dat(i)%idx)
         DEALLOCATE(fc%dat(i)%fc)
         !
      ENDDO
      DEALLOCATE(fc%dat)
      !
   END SUBROUTINE destroy_fc3_sparse
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE div_mass_fc3_sparse (fc,S)
      USE kinds, only : DP
      USE input_fc, ONLY : ph_system_info
      IMPLICIT NONE
      CLASS(sparse),INTENT(inout) :: fc
      TYPE(ph_system_info),INTENT(in) :: S
      !
      INTEGER :: i, j
      !
      IF(.not.ALLOCATED(S%sqrtmm1)) &
         call errore('div_mass_fc3_sparse', 'missing sqrtmm1, call aux_system first', 1)

      DO i = 1, fc%n_R
         DO j = 1, fc%n_terms(i)
            fc%dat(i)%fc(j) = fc%dat(i)%fc(j) &
               *S%sqrtmm1( fc%dat(i)%idx(1,j) ) &
               *S%sqrtmm1( fc%dat(i)%idx(2,j) ) &
               *S%sqrtmm1( fc%dat(i)%idx(3,j) )
         ENDDO
      ENDDO
      !
   END SUBROUTINE div_mass_fc3_sparse
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE read_fc3_grid(fc, filename, S)
      USE input_fc, ONLY : read_system
      USE input_fc, ONLY : ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(inout)   :: S ! = System
      CLASS(grid),INTENT(inout) :: fc
      !
      CHARACTER(13),PARAMETER :: sub = "read_fc3_grid"
      !
      INTEGER :: unit, ios
      INTEGER, EXTERNAL :: find_free_unit
      !
      INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
      INTEGER :: na1_, na2_, na3_, j1_, j2_, j3_
      INTEGER :: n_R, i
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      ioWRITE(stdout,*) "** Reading full grid FC3 file ", TRIM(filename)
      !
      CALL read_system(unit, S)
      !
      READ(unit, *) fc%nq
      ioWRITE(stdout,*) "   Original FC3 grid:", fc%nq
      !
      DO na1=1,S%nat
         DO na2=1,S%nat
            DO na3=1,S%nat
               DO j1=1,3
                  jn1 = j1 + (na1-1)*3
                  DO j2=1,3
                     jn2 = j2 + (na2-1)*3
                     DO j3=1,3
                        jn3 = j3 + (na3-1)*3
                        !
                        READ(unit,*) j1_, j2_, j3_, na1_, na2_, na3_
                        IF ( ANY((/na1,na2,na3,j1,j2,j3/) /= (/na1_,na2_,na3_,j1_,j2_,j3_/)) ) THEN
                           print*, (/na1,na2,na3,j1,j2,j3/)
                           print*, (/na1_,na2_,na3_,j1_,j2_,j3_/)
                           CALL errore(sub,'not matching na1,na2,na3,j1,j2,j3 in file "'//TRIM(filename)//"'",1)
                        ENDIF
                        !
                        READ(unit,*) n_R
                        IF ( fc%n_R == 0) THEN
                           IF(allocated(fc%yR2) .or. allocated(fc%xR2) .or. &
                              allocated(fc%yR3) .or. allocated(fc%xR3) ) &
                              CALL errore(sub, 'some element are already allocated', 1)
                           !
                           ioWRITE(stdout,*) "   Number of R:      ", n_R
                           ALLOCATE(fc%yR2(3,n_R), fc%yR3(3,n_R))
                           ALLOCATE(fc%xR2(3,n_R), fc%xR3(3,n_R))
                           ALLOCATE(fc%FC(3*S%nat,3*S%nat,3*S%nat,n_R))
                           fc%n_R = n_R
                        ELSE
                           IF(n_R/=fc%n_R) CALL errore(sub, "cannot read this fc format",1)
                        ENDIF
                        !
                        DO i = 1, fc%n_R
                           READ(unit,*) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i)
                           ! also, find the index of R=0
                           IF( ALL(fc%yR2(:,i)==0) .and. ALL(fc%yR3(:,i)==0) ) fc%i_00 = i
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      CLOSE(unit)
      ! Prepare the lists of R in cartesian coords
      fc%xR2 = REAL(fc%yR2, kind=DP)
      CALL cryst_to_cart(n_R,fc%xR2,S%at,1)
      fc%xR3 = REAL(fc%yR3, kind=DP)
      CALL cryst_to_cart(n_R,fc%xR3,S%at,1)
      !
   END SUBROUTINE read_fc3_grid
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE write_fc3_grid(fc, filename, S)
      USE input_fc, ONLY : write_system, ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(in)   :: S ! = System
      CLASS(grid),INTENT(in) :: fc
      !
      CHARACTER(14),PARAMETER :: sub = "write_fc3_grid"
      !
      INTEGER :: unit, ios
      INTEGER, EXTERNAL :: find_free_unit
      ! format:
      CHARACTER (6), EXTERNAL :: int_to_char
      INTEGER :: n_digits(3)
      CHARACTER(64) :: cformat, cformat2
      !
      INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
      INTEGER :: i
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      CALL write_system(unit, S)
      !
      WRITE(unit, '(3i9)') fc%nq
      !
      n_digits(1) = CEILING(LOG10(DBLE(S%nat)+0.1_dp))+1
      cformat='(i1,2i2,x,3i'//int_to_char(n_digits(1))//')'
      !
      n_digits(2) = CEILING(LOG10(DBLE(MAXVAL(ABS(fc%yR2)))+0.1_dp))+2 ! +2 to account minus signs
      n_digits(1) = n_digits(2)-1 ! no space at lines beginning
      n_digits(3) = CEILING(LOG10(DBLE(MAXVAL(ABS(fc%yR3)))+0.1_dp))+2
      !
      cformat2 = '( i'//int_to_char(n_digits(1))// &
         ',2i'//int_to_char(n_digits(2))// &
         ',3i'//int_to_char(n_digits(3))//',1pe25.17,1pe13.3)'
      !cformat2 = '(i6,2i6,3i6,1pe25.17,1pe13.3)'
      !
      DO na1=1,S%nat
         DO na2=1,S%nat
            DO na3=1,S%nat
               DO j1=1,3
                  jn1 = j1 + (na1-1)*3
                  DO j2=1,3
                     jn2 = j2 + (na2-1)*3
                     DO j3=1,3
                        jn3 = j3 + (na3-1)*3
                        !
                        WRITE(unit,cformat) j1, j2, j3, na1, na2, na3
                        !
                        WRITE(unit,'(i9)') fc%n_R
                        !
                        DO i = 1, fc%n_R
                           IF(allocated(fc%iFc))THEN
                              WRITE(unit,cformat2) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i), fc%iFC(jn1,jn2,jn3,i)
                           ELSE
                              WRITE(unit,cformat2) fc%yR2(:,i), fc%yR3(:,i), fc%FC(jn1,jn2,jn3,i)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      CLOSE(unit)
      !
   END SUBROUTINE write_fc3_grid
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE destroy_fc3_grid(fc)
      USE input_fc, ONLY : write_system
      IMPLICIT NONE
      CLASS(grid),INTENT(inout) :: fc
      !
      DEALLOCATE(fc%yR2)
      DEALLOCATE(fc%yR3)
      DEALLOCATE(fc%xR2)
      DEALLOCATE(fc%xR3)
      DEALLOCATE(fc%FC)
      !
   END SUBROUTINE destroy_fc3_grid
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE div_mass_fc3_grid(fc,S)
      USE kinds,    ONLY : DP
      USE input_fc, ONLY : ph_system_info
      IMPLICIT NONE
      CLASS(grid),INTENT(inout) :: fc
      TYPE(ph_system_info),INTENT(in) :: S
      !
      INTEGER :: i, j, k, i_R
      !
      IF(.not.ALLOCATED(S%sqrtmm1)) &
         call errore('div_mass_fc3_grid', 'missing sqrtmm1, call aux_system first', 1)

      DO i_R = 1, fc%n_R
         DO k = 1, S%nat3
            DO j = 1, S%nat3
               DO i = 1, S%nat3
                  fc%FC(i, j, k, i_R) = fc%FC(i, j, k, i_R) &
                     * S%sqrtmm1(i)*S%sqrtmm1(j)*S%sqrtmm1(k)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
   END SUBROUTINE div_mass_fc3_grid
   !
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE destroy_fc3_(fc)
      IMPLICIT NONE
      CLASS(forceconst3),INTENT(inout) :: fc
      !
      DEALLOCATE(fc%yR2)
      DEALLOCATE(fc%yR3)
      DEALLOCATE(fc%xR2)
      DEALLOCATE(fc%xR3)
      !
   END SUBROUTINE destroy_fc3_
   !
!       procedure :: interpolate  => fftinterp_mat3_constant
!       procedure :: destroy      => destroy_fc3_constant
!       procedure :: read         => read_fc3_constant
!       procedure :: write        => write_fc3_constant
!       procedure :: div_mass     => div_mass_fc3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_constant(fc, xq2,xq3, nat3, D)
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(constant),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      !
      D = fc%D
      !
   END SUBROUTINE fftinterp_mat3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fftinterp_mat3_constant_grad(fc, xq2,xq3, nat3, D, Dgrad)
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(constant),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Dgrad(nat3, nat3, nat3,3)

      !
      D = fc%D
      Dgrad = 0._dp
      !
   END SUBROUTINE fftinterp_mat3_constant_grad
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE fft_doubleinterp_mat3_constant(fc, xq2,xq3,xq3b, nat3, D, Db)
      IMPLICIT NONE
      !
      INTEGER,INTENT(in)   :: nat3
      CLASS(constant),INTENT(in) :: fc
      REAL(DP),INTENT(in) :: xq2(3), xq3(3), xq3b(3)
      COMPLEX(DP),INTENT(out) :: D(nat3, nat3, nat3)
      COMPLEX(DP),INTENT(out) :: Db(nat3, nat3, nat3)
      D  = fc%D
      Db = fc%D
   END SUBROUTINE fft_doubleinterp_mat3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE div_mass_fc3_constant(fc,S)
      USE kinds,    ONLY : DP
      USE input_fc, ONLY : ph_system_info
      IMPLICIT NONE
      CLASS(constant),INTENT(inout) :: fc
      TYPE(ph_system_info),INTENT(in) :: S
      !
      INTEGER :: i, j, k, i_R
      !
      IF(.not.ALLOCATED(S%sqrtmm1)) &
         call errore('div_mass_fc3_grid', 'missing sqrtmm1, call aux_system first', 1)

      ALLOCATE(fc%D(S%nat3, S%nat3, S%nat3))
      fc%D = fc%constant
      !
      DO k = 1, S%nat3
         DO j = 1, S%nat3
            DO i = 1, S%nat3
               fc%D(i, j, k) = fc%D(i, j, k) &
                  * S%sqrtmm1(i)*S%sqrtmm1(j)*S%sqrtmm1(k)
            ENDDO
         ENDDO
      ENDDO
      !
   END SUBROUTINE div_mass_fc3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE destroy_fc3_constant(fc)
      IMPLICIT NONE
      CLASS(constant),INTENT(inout) :: fc
      !
      DEALLOCATE(fc%D)
      fc%constant = 0._dp
      !
   END SUBROUTINE destroy_fc3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE read_fc3_constant(fc, filename, S)
      USE input_fc, ONLY : read_system
      USE input_fc, ONLY : ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(inout)   :: S ! = System
      CLASS(constant),INTENT(inout) :: fc
      INTEGER :: unit, ios
      CHARACTER(17),PARAMETER :: sub = "read_fc3_constant"
      CHARACTER(32) :: buf
      INTEGER,EXTERNAL :: find_free_unit
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='read',status='old',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      READ(unit, '(a32)') buf
      IF(buf/="constant") CALL errore(sub, "cannot read this format", 1)
      !
      ioWRITE(stdout,*) "** Reading constant FC3 file ", TRIM(filename)
      CALL read_system(unit, S)
      READ(unit, *) fc%constant
      WRITE(stdout,*) "Constant FC", fc%constant
      RETURN
      !
   END SUBROUTINE read_fc3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>
   SUBROUTINE write_fc3_constant(fc, filename, S)
      USE input_fc, ONLY : write_system, ph_system_info
      IMPLICIT NONE
      CHARACTER(*),INTENT(in)          :: filename
      TYPE(ph_system_info),INTENT(in)   :: S ! = System
      CLASS(constant),INTENT(in) :: fc
      CHARACTER(14),PARAMETER :: sub = "write_fc3_constant"
      INTEGER :: unit, ios
      INTEGER, EXTERNAL :: find_free_unit
      !
      unit = find_free_unit()
      OPEN(unit=unit,file=filename,action='write',status='unknown',iostat=ios)
      IF(ios/=0) CALL errore(sub,"opening '"//TRIM(filename)//"'", 1)
      !
      WRITE(unit, '(a)') "constant"
      CALL write_system(unit, S)
      WRITE(unit,*) FC%constant
   END SUBROUTINE write_fc3_constant
   ! \/o\________\\\______________________//\/___________________/~^>>

   SUBROUTINE ip_cart2pat_byhand(d3in, nat3, u1, u2, u3)
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
      ! u1t = TRANSPOSE(CONJG(u1))
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,a,b,c,AUX) REDUCTION(+: d3tmp) COLLAPSE(2)
      DO c = 1,nat3
         DO b = 1,nat3
            ! Precompute u2*u3 to save some FLOPS without
            ! compromising memory access order
            DO k = 1,nat3
               DO j = 1,nat3
                  AUX(j,k) = u2(b,j) * u3(c,k)
               ENDDO
            ENDDO
            !
            DO a = 1,nat3
               DO k = 1,nat3
                  DO j = 1,nat3
                     DO i = 1,nat3
                        d3tmp(i, j, k) = d3tmp(i, j, k) &
                           + CONJG(u1(a,i) * AUX(j,k)) * d3in(a, b, c)
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
   END SUBROUTINE ip_cart2pat_byhand
   !
   SUBROUTINE ip_cart2pat_gemm(d3in, nat3, u1, u2, u3)
      USE kinds, ONLY : DP
      IMPLICIT NONE
      ! Inputs/outputs
      COMPLEX(DP), INTENT(INOUT) :: d3in(nat3, nat3, nat3)
      INTEGER, INTENT(IN)        :: nat3
      COMPLEX(DP), INTENT(IN)    :: u1(nat3, nat3), u2(nat3, nat3), u3(nat3, nat3)
   
      ! Temporary arrays
      COMPLEX(DP), ALLOCATABLE :: tmp(:,:,:)
      INTEGER :: i, k, nat32
      COMPLEX(DP) :: u3c(nat3, nat3)
   
      ! Allocate temporary arrays
      ALLOCATE(tmp(nat3, nat3, nat3))
      ! ALLOCATE(tmp2(nat3, nat3, nat3))
   
      ! Initialize temporary arrays
      tmp = 0._dp
      ! tmp2 = 0._dp
      nat32 = nat3**2
   
      ! d3(i,j,k) = u(a,i) * u(b,j) * u(c,k) * d3in(a,b,c)
   
      u3c = CONJG(u3)
      ! Step 1: tmp(a,b,k) = d3in(a,b,c) * CONJG(u3(c,k))
      CALL ZGEMM('N', 'N', nat32, nat3, nat3, 1.0_DP, d3in, nat32, u3c, nat3, 0.0_DP, tmp, nat32)
   
      ! ! Step 2: tmp2(i,b,k) = TRANSCONJ(u(a,i)) * tmp(a,b,k)
      CALL ZGEMM('C', 'N', nat3, nat32, nat3, 1.0_DP, u1, nat3, tmp, nat3, 0.0_DP, d3in, nat3)
   
      ! Step 3: tmp(b,i,k) = RESHAPE(tmp2(i,b,k))
      do k = 1, nat3
         tmp(:, :, k) = TRANSPOSE(d3in(:, :, k))
      enddo
   
      ! Step 3: d3in(i,j,k) = TRANSCONJ(u2(b,j)) * tmp(b,i,k)
      CALL ZGEMM('C', 'N', nat3, nat32, nat3, 1.0_DP, u2, nat3, tmp, nat3, 0.0_DP, d3in, nat3)
   
      do k = 1, nat3
         d3in(:, :, k) = TRANSPOSE(d3in(:, :, k))
      enddo
   
      ! Deallocate temporary arrays
      DEALLOCATE(tmp)
      ! DEALLOCATE(tmp2)
   
      RETURN
   END SUBROUTINE

END MODULE
