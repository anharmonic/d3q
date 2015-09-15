
MODULE mpi_thermal
  USE kinds, ONLY : DP
  include "mpif.h"

  INTEGER :: my_id, num_procs, ierr
  LOGICAL :: ionode

  CONTAINS

  SUBROUTINE start_mpi()
    IMPLICIT NONE
    call MPI_INIT ( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
    ionode = (my_id == 0)
    IF(ionode) WRITE(*,*) "Running on ", num_procs, "cpus"
  END SUBROUTINE

  ! Divide a vector among all the CPUs
  SUBROUTINE scatter_v(nn_send, v_send, nn_recv, v_recv)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nn_send
    REAL(DP),INTENT(in) :: v_send(nn_send)
    !
    INTEGER,INTENT(out)  :: nn_recv
    REAL(DP),ALLOCATABLE,INTENT(out) :: v_recv(:)

    !
    INTEGER :: nn_residual, i
    INTEGER,ALLOCATABLE  :: nn_scatt(:), ii_scatt(:)

    nn_recv = INT(nn_send/num_procs)
    nn_residual = nn_send - nn_recv*num_procs
    ALLOCATE(nn_scatt(num_procs))
    ALLOCATE(ii_scatt(num_procs))
    IF(nn_residual>0) nn_scatt(1:nn_residual) = nn_scatt(1:nn_residual) +1 
    IF(SUM(nn_scatt)/=nn_send) STOP 100
    ii_scatt = 0
    DO i = 2,num_procs
     ii_scatt(i) = ii_scatt(i-1)+nn_scatt(i-1)
    ENDDO
    nn_recv = nn_scatt(my_id+1)
    ALLOCATE(v_recv(nn_recv))

    CALL MPI_scatterv(v_send, nn_scatt, ii_scatt, MPI_DOUBLE_PRECISION, &
                      v_recv, nn_recv, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr )
    

  END SUBROUTINE


END MODULE mpi_thermal



