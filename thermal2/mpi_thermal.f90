!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! hack:
#if defined(__XLF)
SUBROUTINE FLUSH(un)
 INTEGER :: un
 RETURN
END SUBROUTINE
#endif

MODULE mpi_thermal
  USE kinds,  ONLY : DP
  !USE timers, ONLY : t_mpicom
#ifdef __MPI
  include "mpif.h"
#endif
#define __MPI_THERMAL
#include "mpi_thermal.h"

!dir$ message "----------------------------------------------------------------------------------------------" 
#ifdef __MPI
!dir$ message "Compiling _with_ MPI support" 
#else
!dir$ message "Compiling _without_ MPI support (mpi subroutines will do nothing)" 
#endif
!
#ifdef _OPENMP
!dir$ message "Compiling _with_ OpenMP directives" 
#else
!dir$ message "Compiling _without_ OpenMP directives" 
#endif
!dir$ message "----------------------------------------------------------------------------------------------" 

  INTEGER :: my_id=0, num_procs=1, ierr
  LOGICAL :: ionode = .TRUE. ! everyone is ionode before I start MPI
  LOGICAL :: mpi_started = .FALSE.
  INTEGER :: omp_tot_thr=1

  INTERFACE mpi_broadcast
     MODULE PROCEDURE mpi_bcast_logical
     !
     MODULE PROCEDURE mpi_bcast_scl
     MODULE PROCEDURE mpi_bcast_vec
     MODULE PROCEDURE mpi_bcast_mat
     MODULE PROCEDURE mpi_bcast_tns
     MODULE PROCEDURE mpi_bcast_tns4
     !
     MODULE PROCEDURE mpi_bcast_zscl
     MODULE PROCEDURE mpi_bcast_zvec
     MODULE PROCEDURE mpi_bcast_zmat
     MODULE PROCEDURE mpi_bcast_ztns
     MODULE PROCEDURE mpi_bcast_ztns4
     !
     MODULE PROCEDURE mpi_bcast_integer
     MODULE PROCEDURE mpi_bcast_integer_vec
     !
     MODULE PROCEDURE mpi_bcast_character
  END INTERFACE
  !
  INTERFACE mpi_bsum
     MODULE PROCEDURE mpi_bsum_int
     MODULE PROCEDURE mpi_bsum_ivec

     MODULE PROCEDURE mpi_bsum_scl
     MODULE PROCEDURE mpi_bsum_vec
     MODULE PROCEDURE mpi_bsum_mat
     MODULE PROCEDURE mpi_bsum_tns
     MODULE PROCEDURE mpi_bsum_tns4

     MODULE PROCEDURE mpi_bsum_zscl
     MODULE PROCEDURE mpi_bsum_zvec
     MODULE PROCEDURE mpi_bsum_zmat
     MODULE PROCEDURE mpi_bsum_ztns
     MODULE PROCEDURE mpi_bsum_ztns4
  END INTERFACE


  CONTAINS

  SUBROUTINE start_mpi()
    IMPLICIT NONE
#ifdef _OPENMP
    INTEGER,EXTERNAL ::  omp_get_max_threads
#endif
#ifdef __MPI
    call MPI_INIT ( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
    ionode = (my_id == 0)
    IF(ionode) WRITE(*,"(2x,a,i6,a)") "Using ", num_procs, " MPI processes"
#else
    my_id = 0
    num_procs = 1
    ionode = .true.
    IF(ionode) WRITE(*,"(2x,a,i6,a)") "Running without MPI support"
#endif
!
#ifdef _OPENMP
    omp_tot_thr =  omp_get_max_threads()
#else
   omp_tot_thr = 1
#endif
    CALL mpi_bsum(omp_tot_thr)
    IF(ionode .and. omp_tot_thr>num_procs) WRITE(*,"(2x,a,i6,a)") &
      "Using",  omp_tot_thr, " total MPI+OpenMP threads"
    mpi_started = .true.
  END SUBROUTINE

  SUBROUTINE stop_mpi()
#ifdef __MPI
     call MPI_FINALIZE ( ierr )
#endif
  END SUBROUTINE

  SUBROUTINE abort_mpi(errorcode)
        IMPLICIT NONE
        INTEGER :: ierr
        INTEGER, INTENT(IN):: errorcode
#ifdef __MPI
        CALL mpi_abort(mpi_comm_world, errorcode, ierr)
#endif
  END SUBROUTINE

  SUBROUTINE mpi_wbarrier()
        IMPLICIT NONE
        INTEGER :: ierr
#ifdef __MPI
        CALL mpi_barrier(mpi_comm_world, ierr)
#endif
  END SUBROUTINE
  
  SUBROUTINE mpi_any(lgc)
    IMPLICIT NONE
    LOGICAL,INTENT(inout) :: lgc

#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, lgc, 1, MPI_LOGICAL, MPI_LOR,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  SUBROUTINE mpi_all(lgc)
    IMPLICIT NONE
    LOGICAL,INTENT(inout) :: lgc

#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, lgc, 1, MPI_LOGICAL, MPI_LAND,&
                       MPI_COMM_WORLD, ierr)
#endif
   END SUBROUTINE
 
  ! In-place MPI sum of integer, scalar, vector and matrix
  SUBROUTINE mpi_bsum_int(scl)
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: scl
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_INTEGER, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_ivec(nn, vec)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nn
    INTEGER,INTENT(inout) :: vec(nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vec, nn, MPI_INTEGER, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE

  !
  SUBROUTINE mpi_bsum_scl(scl)
    IMPLICIT NONE
    REAL(DP),INTENT(inout) :: scl
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_vec(nn, vec)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nn
    REAL(DP),INTENT(inout) :: vec(nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vec, nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_mat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    REAL(DP),INTENT(inout) :: mat(mm,nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mat, mm*nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_tns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    REAL(DP),INTENT(inout) :: tns(ll, mm,nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns, ll*mm*nn, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_tns4(ll, mm, nn, oo, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn, oo
    REAL(DP),INTENT(inout) :: tns(ll, mm,nn, oo)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns, ll*mm*nn*oo, MPI_DOUBLE_PRECISION, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  
  !!  ! --------- ------------- --- -- -- -- - - - complex numbers follow
  SUBROUTINE mpi_bsum_zscl(scl)
    IMPLICIT NONE
    COMPLEX(DP),INTENT(inout) :: scl

#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, scl, 1, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_zvec(nn, vec)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: nn
    COMPLEX(DP),INTENT(inout) :: vec(nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vec, nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_zmat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    COMPLEX(DP),INTENT(inout) :: mat(mm,nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mat, mm*nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_ztns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    COMPLEX(DP),INTENT(inout) :: tns(ll, mm,nn)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns, ll*mm*nn, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bsum_ztns4(ll, mm, nn, oo, tns4)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn, oo
    COMPLEX(DP),INTENT(inout) :: tns4(ll, mm,nn, oo)
#ifdef __MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tns4, ll*mm*nn*oo, MPI_DOUBLE_COMPLEX, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  ! --------- ------------- --- -- -- -- - - - complex numbers follow
  !
  SUBROUTINE mpi_bcast_logical(logi)
    IMPLICIT NONE
    LOGICAL,INTENT(inout) :: logi
#ifdef __MPI
    CALL MPI_BCAST(logi, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_integer(inte)
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: inte
#ifdef __MPI
    CALL MPI_BCAST(inte, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_integer_vec(nn,inte)
    IMPLICIT NONE
    INTEGER,INTENT(in)    :: nn
    INTEGER,INTENT(inout) :: inte(nn)
#ifdef __MPI
    CALL MPI_BCAST(inte, nn, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_character(schar)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(inout) :: schar
#ifdef __MPI
    CALL MPI_BCAST(schar, LEN(schar), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_scl(scl)
    IMPLICIT NONE
    REAL(DP),INTENT(inout) :: scl
#ifdef __MPI
    CALL MPI_BCAST(scl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_vec(nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     ::  nn
    REAL(DP),INTENT(inout) :: mat(nn)
#ifdef __MPI
    CALL MPI_BCAST(mat, nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_mat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    REAL(DP),INTENT(inout) :: mat(mm,nn)
#ifdef __MPI
    CALL MPI_BCAST(mat, mm*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_tns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    REAL(DP),INTENT(inout) :: tns(ll, mm,nn)
#ifdef __MPI
    CALL MPI_BCAST(tns, ll*mm*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_tns4(ll, mm, nn, oo, tns4)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn, oo
    REAL(DP),INTENT(inout) :: tns4(ll, mm,nn, oo)
#ifdef __MPI
    CALL MPI_BCAST(tns4, ll*mm*nn*oo, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  ! Now broadcast for complex numbers
  !
  SUBROUTINE mpi_bcast_zscl(scl)
    IMPLICIT NONE
    COMPLEX(DP),INTENT(inout) :: scl
#ifdef __MPI
    CALL MPI_BCAST(scl, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_zvec(nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     ::  nn
    COMPLEX(DP),INTENT(inout) :: mat(nn)
#ifdef __MPI
    CALL MPI_BCAST(mat, nn, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_zmat(mm, nn, mat)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: mm, nn
    COMPLEX(DP),INTENT(inout) :: mat(mm,nn)
#ifdef __MPI
    CALL MPI_BCAST(mat, mm*nn, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_ztns(ll, mm, nn, tns)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn
    COMPLEX(DP),INTENT(inout) :: tns(ll, mm,nn)
#ifdef __MPI
    CALL MPI_BCAST(tns, ll*mm*nn, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  SUBROUTINE mpi_bcast_ztns4(ll, mm, nn, oo, tns4)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: ll, mm, nn, oo
    COMPLEX(DP),INTENT(inout) :: tns4(ll, mm,nn, oo)
#ifdef __MPI
    CALL MPI_BCAST(tns4, ll*mm*nn*oo, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE
  !
  !!
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
#ifdef __MPI
    nn_send = nn
    ALLOCATE(vec_send(nn_send))
    vec_send(1:nn_send) = vec(1:nn_send)
    CALL scatter_vec(nn_send, vec_send, nn_recv, vec_recv, ii)
    DEALLOCATE(vec)
    ALLOCATE(vec(nn_recv))
    vec(1:nn_recv) = vec_recv(1:nn_recv)
    nn = nn_recv
    DEALLOCATE(vec_send, vec_recv)
#else
    ! do nothing
#endif
  END SUBROUTINE

!  ! Scatter in-place a vector
!  SUBROUTINE gatheri_vec(nn, vec, ii)
!    IMPLICIT NONE
!    INTEGER,INTENT(inout)  :: nn
!    REAL(DP),INTENT(inout),ALLOCATABLE :: vec(:)
!    INTEGER,OPTIONAL,INTENT(out)  :: ii
!    !
!    INTEGER  :: nn_send, nn_recv
!    REAL(DP),ALLOCATABLE :: vec_send(:), vec_recv(:)
!    !
!    IF(.not.allocated(vec)) CALL errore('scatteri_vec', 'input vector must be allocated', 1)
!    IF(size(vec)/=nn)      CALL errore('scatteri_vec', 'input vector must be of size nn', 2)
!    !
!#ifdef __MPI
!    nn_send = nn
!    ALLOCATE(vec_send(nn_send))
!    vec_send(1:nn_send) = vec(1:nn_send)
!    CALL gather_vec(nn_send, vec_send, nn_recv, vec_recv, ii)
!    DEALLOCATE(vec)
!    ALLOCATE(vec(nn_recv))
!    vec(1:nn_recv) = vec_recv(1:nn_recv)
!    nn = nn_recv
!#else
!    ! do nothing
!#endif
!  END SUBROUTINE
  
  
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
#ifdef __MPI
    nn_send = nn
    ALLOCATE(mat_send(mm,nn_send))
    mat_send(:,1:nn_send) = mat(:,1:nn_Send)
    CALL scatter_mat(mm,nn_send, mat_send, nn_recv, mat_recv)
    DEALLOCATE(mat)
    ALLOCATE(mat(mm,nn_recv))
    mat(:,1:nn_recv) = mat_recv(:,1:nn_recv)
    nn = nn_recv
    DEALLOCATE(mat_send, mat_recv)
#else
    ! do nothing
#endif
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
!    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

#ifdef __MPI
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
    DEALLOCATE(nn_scatt, ii_scatt)
#else
    nn_recv = nn_send
    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    ALLOCATE(vec_recv(nn_recv))
    vec_recv = vec_send
    IF(present(ii_recv)) ii_recv = 0
#endif
  END SUBROUTINE

  ! Divide a vector among all the CPUs
  SUBROUTINE gather_vec(nn_send, vec_send, vec_recv) !, ii_recv)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nn_send
    REAL(DP),INTENT(in) :: vec_send(nn_send)
    REAL(DP),ALLOCATABLE,INTENT(out) :: vec_recv(:) 
    !INTEGER,OPTIONAL,INTENT(out)  :: ii_recv
    !
    !INTEGER :: nn_summed, i
    INTEGER :: i, nn_recv
    INTEGER,ALLOCATABLE  :: nn_scatt(:), ii_scatt(:)
    CHARACTER(len=11),PARAMETER :: sub='scatter_vec'
    !
    IF(.not.mpi_started) CALL errore(sub, 'MPI not started', 1)
!    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

#ifdef __MPI
    ALLOCATE(nn_scatt(num_procs))
    ALLOCATE(ii_scatt(num_procs))

    nn_scatt = 0
    nn_scatt(my_id+1) = nn_send
    CALL mpi_bsum(num_procs, nn_scatt)
    ii_scatt = 0
    ii_scatt(my_id+1) = SUM(nn_scatt(1:my_id))
    CALL mpi_bsum(num_procs, ii_scatt)
    nn_recv = nn_send
    CALL mpi_bsum(nn_recv)

    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    IF(my_id==0) THEN
      ALLOCATE(vec_recv(nn_recv))
    ELSE
      ALLOCATE(vec_recv(0))
    ENDIF
    !WRITE(*,'(99i6)'), my_id, nn_scatt, ii_scatt
    !WRITE(*,'(99i6)'), my_id, nn_send, size(vec_send)
    

    CALL MPI_gatherv(vec_send, nn_send, MPI_DOUBLE_PRECISION, &
                      vec_recv, nn_scatt, ii_scatt, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr )
#else
    nn_recv = nn_send
    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    ALLOCATE(vec_recv(nn_recv))
    vec_recv = vec_send
#endif
  END SUBROUTINE
  !
  ! Divide a vector among all the CPUs
  SUBROUTINE gather_mat(mm, nn_send, vec_send, vec_recv) !, ii_recv)
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: mm, nn_send
    REAL(DP),INTENT(in) :: vec_send(nn_send)
    REAL(DP),ALLOCATABLE,INTENT(out) :: vec_recv(:,:) 
    !INTEGER,OPTIONAL,INTENT(out)  :: ii_recv
    !
    !INTEGER :: nn_summed, i
    INTEGER :: i, nn_recv, mmnn_send
    INTEGER,ALLOCATABLE  :: nn_scatt(:), ii_scatt(:)
    CHARACTER(len=11),PARAMETER :: sub='scatter_vec'
    !
    IF(.not.mpi_started) CALL errore(sub, 'MPI not started', 1)
!    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

#ifdef __MPI
    ALLOCATE(nn_scatt(num_procs))
    ALLOCATE(ii_scatt(num_procs))

    nn_scatt = 0
    nn_scatt(my_id+1) = nn_send
    CALL mpi_bsum(num_procs, nn_scatt)
    ii_scatt = 0
    ii_scatt(my_id+1) = SUM(nn_scatt(1:my_id))
    CALL mpi_bsum(num_procs, ii_scatt)
    nn_recv = nn_send
    CALL mpi_bsum(nn_recv)
    

    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    IF(my_id==0) THEN
      ALLOCATE(vec_recv(mm,nn_recv))
    ELSE
      ALLOCATE(vec_recv(0,0))
    ENDIF
    !WRITE(*,'(99i6)'), my_id, nn_scatt, ii_scatt
    !WRITE(*,'(99i6)'), my_id, nn_send, size(vec_send)
    
    nn_scatt = nn_scatt*mm
    ii_scatt = ii_scatt*mm
    mmnn_send = nn_send*mm

    CALL MPI_gatherv(vec_send, mmnn_send, MPI_DOUBLE_PRECISION, &
                     vec_recv, nn_scatt, ii_scatt, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr )
#else
    nn_recv = nn_send
    IF(allocated(vec_recv)) DEALLOCATE(vec_recv)
    ALLOCATE(vec_recv(mm,nn_recv))
    vec_recv = vec_send
#endif
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
!    IF(num_procs>nn_send) CALL errore(sub, 'num_procs > nn_send, this can work but makes no sense', 1)

#ifdef __MPI
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
#else
    nn_recv = nn_send
    IF(allocated(mat_recv)) DEALLOCATE(mat_recv)
    ALLOCATE(mat_recv(mm,nn_recv))
    mat_recv = mat_send
#endif
    
  END SUBROUTINE



END MODULE mpi_thermal



