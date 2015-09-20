
#ifndef __MPI_THERMAL
USE mpi_thermal, ONLY : ionode
#endif

#define ioWRITE IF(ionode) WRITE
#define timer_CALL CALL
#define stdin  5
#define stdout 6
#define stderr 0

