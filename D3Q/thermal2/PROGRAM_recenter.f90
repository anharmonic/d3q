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
  USE cmdline_param_module
  IMPLICIT NONE
  INTEGER :: far, ios
  CHARACTER(len=256) :: filein, fileout, dummy
  TYPE(grid) :: fc, fcb
  TYPE(ph_system_info) :: S
  !
  INTEGER :: nq(3), nq_trip
  TYPE(d3_list),ALLOCATABLE :: d3grid(:)
  !
  CHARACTER(len=:),ALLOCATABLE :: cmdline
  LOGICAL :: writed3

  filein  = cmdline_param_char("i", "mat3R")
  fileout = cmdline_param_char("o", TRIM(filein)//".recentered")
  far      = cmdline_param_int("f", 2)
  writed3  = cmdline_param_logical("w")
  !
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: d3_recenter.x NQX NQY NQZ [-i FILEIN] [-o FILEOUT] [-f NFAR] [-w]"
      WRITE(*,*) ""
      WRITE(*,*) "Reads force constants from FILEIN (default: mat3R), interpolate them on a grid"
      WRITE(*,*) "of NQX x NQY x NQZ points, recenter them on a Wigner-Seitz cell constructed up"
      WRITE(*,*) "to NFAR unit cells and save the result in FILEOUT (default: FILEIN_recenter)"
      WRITE(*,*) "If '-w' is specified, write the intermediate D3 matrices to files called "
      WRITE(*,*) "   atmp_Q1x_Q2x_Q3x (default: don't write, lot of output!)"
     
      STOP 1
  ENDIF
  !
  cmdline = cmdline_residual()
  READ(cmdline, *, iostat=ios) nq, dummy
  IF(ios==0 .and. len(dummy)>0) CALL errore("qq2rr", "too many argument use command '-h' for help",1)
  READ(cmdline, *, iostat=ios) nq
  IF(ios/=0) CALL errore("recenter", "missing argument use command '-h' for help",1)
  !
  IF(TRIM(fileout)==TRIM(filein)) &
    CALL errore("recenter","filein and fileout are the same, I refuse to do that",1)
  !
  WRITE(*,*) "Number of neighbours to check for BZ", far
  WRITE(*,*) "Input file ", TRIM(filein)
  WRITE(*,*) "Output file ", TRIM(fileout)
  
  CALL fc%read(filein, S)
  CALL aux_system(s)
  !
  nq_trip = (nq(1)*nq(2)*nq(3))**2
  ALLOCATE(d3grid(nq_trip))
  CALL regen_fwfft_d3(nq, nq_trip, S, d3grid, fc, writed3)
  
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fcb, far)
  fcb%nq = nq
  CALL fcb%write(fileout, S)
  !
  !
END PROGRAM
!
