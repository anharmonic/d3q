
#ifndef __MPI_THERMAL
USE mpi_thermal, ONLY : ionode
#endif

#define ioWRITE IF(ionode) WRITE
#define ioREAD  IF(ionode) READ
#define ioFLUSH IF(ionode) FLUSH
#define timer_CALL CALL
#define stdin  5
#define stdout 6
#define stderr 0

!IMPLICIT NONE
!#define i2c(a) = TRIM(int_to_char(a))
!CHARACTER (LEN=6), EXTERNAL :: int_to_char

