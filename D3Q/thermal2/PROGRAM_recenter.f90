!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE recenter_module
  USE kinds, ONLY : DP
#include "mpi_thermal.h"
  IMPLICIT NONE
  !
END MODULE


PROGRAM recenter
  USE fc3_interpolate, ONLY : grid
  USE input_fc,        ONLY : read_system, ph_system_info
  USE ph_system,       ONLY : aux_system
  USE recenter_module
  USE f3_bwfft
  IMPLICIT NONE
  INTEGER :: narg, i, far
  CHARACTER(len=512) :: argv
  CHARACTER(len=256) :: filein, fileout
  TYPE(grid) :: fc, fcb
  TYPE(ph_system_info) :: S
  !
  INTEGER :: nq(3), nq_trip
  TYPE(d3_list),ALLOCATABLE :: d3grid(:)
  !
  !fileout  = "mat3R.centered"

  !
  narg = command_argument_count()
  IF(narg>=3)THEN
    DO i = 1,3
      CALL get_command_argument(i,argv)
      READ(argv,*) nq(i)
    ENDDO
    WRITE(*,*) "Output supercell size:", nq
  ELSE
      WRITE(*,*) "Syntax: d3_qq2rr.x NQX NQY NQZ [NFAR] [mat3R.INPUT] [mat3R.OUTPUT]"
      WRITE(*,*) ""
      WRITE(*,*) "Reads force constants from mat3R.INPUT, interpolate them on a grid"
      WRITE(*,*) "of NQX x NQY x NQZ points, recenter them on a Wigner-Seitz cell "
      WRITE(*,*) "constructed up to NFAR unit cells and save the result in mat3R.OUTPUT"
      WRITE(*,*) ""
      STOP 1
  ENDIF
  !
  IF(narg>=4)THEN
    CALL get_command_argument(4,argv)
    READ(argv,*) far
  ELSE
    far = 2
  ENDIF
  !
  IF(narg>=5)THEN
    CALL get_command_argument(5,filein)
  ELSE
    filein = 'mat3R'
  ENDIF
  !
  IF(narg>=6)THEN
    CALL get_command_argument(6,fileout)
  ELSE
    fileout = TRIM(filein)//"_recentered"
  ENDIF
  !
  IF(TRIM(fileout)==TRIM(filein)) &
    CALL errore("recenter","filein and fileout are the same, I refuse to do that",1)
  !
  WRITE(*,*) "Number of neighbours to check for BZ", far
  WRITE(*,*) "Input file", TRIM(filein)
  WRITE(*,*) "Output file", TRIM(fileout)
  
  CALL fc%read(filein, S)
  CALL aux_system(s)
  !
  nq_trip = (nq(1)*nq(2)*nq(3))**2
  ALLOCATE(d3grid(nq_trip))
  CALL regen_fwfft_d3(nq, nq_trip, S, d3grid, fc)
  
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fcb, far)
  CALL fcb%write(fileout, S)
  !
  !
END PROGRAM
!
