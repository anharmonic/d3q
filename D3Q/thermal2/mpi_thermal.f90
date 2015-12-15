
MODULE mpi_thermal
  USE kinds,  ONLY : DP
  !USE timers, ONLY : t_mpicom
  include "mpif.h"
#define __MPI_THERMAL
#include "mpi_thermal.h"
  
  INTEGER :: my_id=0, num_procs=1, ierr
  LOGICAL :: ionode = .TRUE. ! everyone is ionode before I start MPI
  LOGICAL :: mpi_started = .FALSE.
  INTEGER :: omp_tot_thr=1

  INTERFACE mpi_bsum
     MODULE PROCEDURE mpi_bsum_int

     MODULE PROCEDURE mpi_bsum_scl
     MODULE PROCEDURE mpi_bsum_vec
     MODULE PROCEDURE mpi_bsum_mat
     MODULE PROCEDURE mpi_bsum_tns

     MODULE PROCEDURE mpi_bsum_zscl
     MODULE PROCEDURE mpi_bsum_zvec
     MODULE PROCEDURE mpi_bsum_zmat
     MODULE PROCEDURE mpi_bsum_ztns
     MODULE PROCEDURE mpi_bsum_ztns4
  END INTERFACE


  CONTAINS

  SUBROUTINE start_mpi()
    IMPLICIT NONE
#ifdef __OPENMP
    INTEGER,EXTERNAL ::  omp_get_max_threads
#endif
    call MPI_INIT ( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
    ionode = (my_id == 0)
    IF(ionode .and. num_procs>1) WRITE(*,"(2x,a,i6,a)") "Using ", num_procs, " MPI processes"
#ifdef __OPENMP
    omp_tot_thr =  omp_get_max_threads()
#else
   omp_tot_thr = 1
#endif
    CALL mpi_bsum(omp_tot_thr)
    IF(ionode .and. omp_tot_thr>num_procs) WRITE(*,"(2x,a,i6,a)") "Using",  omp_tot_thr, " total MPI+OpenMP threads"
    mpi_started = .true.
  END SUBROUTINE

  SUBROUTINE stop_mpi()
     call MPI_FINALIZE ( ierr )
  END SUBROUTINE

  SUBROUTINE abort_mpi(errorcode)
        INTEGER :: ierr
        INTEGER, INTENT(IN):: errorcode
        CALL mpi_abort(mpi_comm_world, errorcode, ierr)
  END SUBROUTINE

  SUBROUTINE mpi_any(lgc)
    IMPLICIT NONE
    LOGICAL,INTENT(inout) :: lgc

      !timer_CALL t_mpicom%start()
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, lgc, 1, MPI_LOGICAL, MPI_LOR,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_all(lgc)
    IMPLICIT NONE
    LOGICAL,INTENT(inout) :: lgc

      !timer_CALL t_mpicom%start()
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, lgc, 1, MPI_LOGICAL, MPI_LAND,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
 
  ! In-place MPI sum of integer, scalar, vector and matrix
  SUBROUTINE mpi_bsum_int(scl)
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: scl

      !timer_CALL t_mpicom%start()
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_INTEGER, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_scl(scl)
    IMPLICIT NONE
    REAL(DP),INTENT(inout) :: scl

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_vec(nn, vec)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nn
    REAL(DP),INTENT(inout) :: vec(nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vec, nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_mat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    REAL(DP),INTENT(inout) :: mat(mm,nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mat, mm*nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_tns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    REAL(DP),INTENT(inout) :: tns(ll, mm,nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns, ll*mm*nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
!!  ! --------- ------------- --- -- -- -- - - - complex numbers follow
  SUBROUTINE mpi_bsum_zscl(scl)
    IMPLICIT NONE
    COMPLEX(DP),INTENT(inout) :: scl

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_zvec(nn, vec)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nn
    COMPLEX(DP),INTENT(inout) :: vec(nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vec, nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_zmat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    COMPLEX(DP),INTENT(inout) :: mat(mm,nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mat, mm*nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_ztns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    COMPLEX(DP),INTENT(inout) :: tns(ll, mm,nn)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns, ll*mm*nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE
  SUBROUTINE mpi_bsum_ztns4(ll, mm, nn, oo, tns4)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn, oo
    COMPLEX(DP),INTENT(inout) :: tns4(ll, mm,nn, oo)

      !timer_CALL t_mpicom%start() 
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns4, ll*mm*nn*oo, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
      !timer_CALL t_mpicom%stop()
  END SUBROUTINE


  ! Scatter in-place a vector
  SUBROUTINE scatteri_vec(nn, vec, ii)
    IMPLICIT NONE
    INTEGER,INTENT(inout)  :: nn
    REAL(DP),INTENT(inout),ALLOCATABLE :: vec(:)
    INTEGER,OPTIONAL,INTENT(out)  :: ii
    !
    INTEGER  :: nn_send, nn_recv
    REAL(DP),ALLOCATABLE :: vec_send(:), vec_recv(:)
    !
    IF(.not.allocated(vec)) CALL errore('scatteri_vec', 'input vector must be allocated', 1)
    IF(size(vec)/=nn)      CALL errore('scatteri_vec', 'input vector must be of size nn', 2)
    !
    nn_send = nn
    ALLOCATE(vec_send(nn_send))
    vec_send(1:nn_send) = vec(1:nn_Send)
    CALL scatter_vec(nn_send, vec_send, nn_recv, vec_recv, ii)
    DEALLOCATE(vec)
    ALLOCATE(vec(nn_recv))
    vec(1:nn_recv) = vec_recv(1:nn_recv)
    nn = nn_recv
  END SUBROUTINE

   ! Scatter in-place a matrix, along the second dimension
  SUBROUTINE scatteri_mat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: mm
    INTEGER,INTENT(inout)  :: nn
    REAL(DP),INTENT(inout),ALLOCATABLE :: mat(:,:)
    !
    INTEGER  :: nn_send, nn_recv
    REAL(DP),ALLOCATABLE :: mat_send(:,:), mat_recv(:,:)
    !
    IF(.not.allocated(mat)) CALL errore('scatteri_mat', 'input matrix must be allocated', 1)
    IF(size(mat,1)/=mm)      CALL errore('scatteri_mat', 'input matrix must be of size mm*nn', 2)
    IF(size(mat,2)/=nn)      CALL errore('scatteri_mat', 'input matrix must be of size mm*nn', 3)
    !
    nn_send = nn
    ALLOCATE(mat_send(mm,nn_send))
    mat_send(:,1:nn_send) = mat(:,1:nn_Send)
    CALL scatter_mat(mm,nn_send, mat_send, nn_recv, mat_recv)
    DEALLOCATE(mat)
    ALLOCATE(mat(mm,nn_recv))
    mat(:,1:nn_recv) = mat_recv(:,1:nn_recv)
    nn = nn_recv
  END SUBROUTINE
 
  ! Divide a vector among all the CPUs
  SUBROUTINE scatter_vec(nn_send, vec_send, nn_recv, vec_recv, ii_recv)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nn_send
    REAL(DP),INTENT(in) :: vec_send(nn_send)
    INTEGER,INTENT(out)  :: nn_recv
    REAL(DP),ALLOCATABLE,INTENT(out) :: vec_recv(:) 
    INTEGER,OPTIONAL,INTENT(out)  :: ii_recv
    !
    INTEGER :: nn_residual, i
    INTEGER,ALLOCATABLE  :: nn_scatt(:), ii_scatt(:)
    CHARACTER(len=11),PARAMETER :: sub='scatter_vec'
    !
    IF(.not.mpi_started) CALL errore(sub, 'MPI not started', 1)
    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

    ALLOCATE(nn_scatt(num_procs))
    ALLOCATE(ii_scatt(num_procs))

    nn_recv = INT(nn_send/num_procs)
    nn_residual = nn_send - nn_recv*num_procs
    nn_scatt=nn_recv
    IF(nn_residual>0) nn_scatt(1:nn_residual) = nn_scatt(1:nn_residual) +1 
    nn_recv = nn_scatt(my_id+1)
    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    ALLOCATE(vec_recv(nn_recv))

    ii_scatt = 0
    DO i = 2,num_procs
     ii_scatt(i) = ii_scatt(i-1)+nn_scatt(i-1)
    ENDDO
    IF(present(ii_recv)) ii_recv = ii_scatt(my_id+1)

    CALL MPI_scatterv(vec_send, nn_scatt, ii_scatt, MPI_DOUBLE_PRECISION, &
                      vec_recv, nn_recv, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr )
    
  END SUBROUTINE

  ! Divide a matrix, along the last dimension, among all the CPUs
  SUBROUTINE scatter_mat(mm, nn_send, mat_send, nn_recv, mat_recv)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: mm, nn_send
    REAL(DP),INTENT(in) :: mat_send(mm,nn_send)
    INTEGER,INTENT(out)  :: nn_recv
    REAL(DP),ALLOCATABLE,INTENT(out) :: mat_recv(:,:) 
    !
    INTEGER :: nn_residual, i
    INTEGER,ALLOCATABLE  :: nn_scatt(:), ii_scatt(:)
    CHARACTER(len=11),PARAMETER :: sub='scatter_mat'
    !
    IF(.not.mpi_started) CALL errore(sub, 'MPI not started', 1)
    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

    ALLOCATE(nn_scatt(num_procs))
    ALLOCATE(ii_scatt(num_procs))

    nn_recv = INT(nn_send/num_procs)
    nn_residual = nn_send - nn_recv*num_procs
    nn_scatt=nn_recv
    IF(nn_residual>0) nn_scatt(1:nn_residual) = nn_scatt(1:nn_residual) +1 
    nn_recv = nn_scatt(my_id+1)
    IF(allocated(mat_recv)) DEALLOCATE(mat_recv)
    ALLOCATE(mat_recv(mm,nn_recv))
    !
    ! account for the first dimension:
    nn_scatt=nn_scatt*mm
    !
    ii_scatt = 0
    DO i = 2,num_procs
     ii_scatt(i) = ii_scatt(i-1)+nn_scatt(i-1)
    ENDDO
    !
    CALL MPI_scatterv(mat_send, nn_scatt, ii_scatt, MPI_DOUBLE_PRECISION, &
                      mat_recv, nn_scatt(my_id+1), MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr )
    
  END SUBROUTINE



END MODULE mpi_thermal



