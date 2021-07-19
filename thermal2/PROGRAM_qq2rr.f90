!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! Code contributions from Giorgia Fugallo and Michele Lazzeri
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
PROGRAM qq2rr
  USE kinds,           ONLY : DP
  USE input_fc,        ONLY : ph_system_info
  USE fc3_interpolate, ONLY : grid
  USE f3_bwfft,        ONLY : d3_list, read_d3_matrices, bwfft_d3_interp, test_fwfft_d3
  USE clib_wrappers,        ONLY : memstat
  USE cmdline_param_module
  
  IMPLICIT NONE
  TYPE(d3_list),ALLOCATABLE ::  d3grid(:)
  TYPE(grid)                :: fc3
  INTEGER :: nq(3)
  INTEGER :: nq_trip, nq_grid
  INTEGER :: kb
  !
  TYPE(ph_system_info) :: S
  COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
  !
  INTEGER :: far, ios
  CHARACTER(len=512) :: argv, filename, dummy
  CHARACTER(len=:),ALLOCATABLE :: cmdline
  LOGICAL :: write_diff, skip_test

  filename = cmdline_param_char("o", "mat3R")
  far      = cmdline_param_int("f", 2)
  write_diff = cmdline_param_logical("w")
  skip_test = cmdline_param_logical("s")
  !
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: ls anh*| d3_qq2rr.x NQX NQY NQZ [-o FILEOUT] [-f NFAR] [-w] [-s]"
      WRITE(*,*) ""
      WRITE(*,*) "Selects a grid of (NQX x NQY x NQZ) points from the anh* files"
      WRITE(*,*) "Apply the inverse Fourier transform, and saves it to FILEOUT (default: mat3R)."
      WRITE(*,*) "Check for shortes perimeter up to NFAR unit cells away (default: 2)."
      WRITE(*,*)
      WRITE(*,*) "-w : when performing the FFT test, if the re-computed D3 matrix differs"
      WRITE(*,*) "     significantly from the initial one it will be printed to a file."
      WRITE(*,*) "     The file name will start with prefix 'anh_cmplx' if the test was"
      WRITE(*,*) "     performed with complex force constant and 'anh_real' if the test"
      WRITE(*,*) "     was performed with just the real part of the FCs (they should be real)"
      WRITE(*,*) ""
      WRITE(*,*) "-s : skip the test"
      WRITE(*,*) ""
      STOP 1
  ENDIF
  
  cmdline = cmdline_residual()
  READ(cmdline, *, iostat=ios) nq, dummy
  IF(ios==0 .and. len(dummy)>0) CALL errore("qq2rr", "too many argument use command '-h' for help",1)
  READ(cmdline, *, iostat=ios) nq
  IF(ios/=0) CALL errore("qq2rr", "missing argument use command '-h' for help",1)
  !
  WRITE(*,*) "Reading grid", nq
  WRITE(*,*) "Number of neighbours to check for BZ", far
  
  nq_grid = nq(1)*nq(2)*nq(3)
  nq_trip = nq_grid**2
  ALLOCATE(d3grid(nq_trip))
  !
  WRITE(*,*) "Reading D3 matrices..."
  CALL read_d3_matrices(nq, nq_trip, S, d3grid)
  WRITE(*,*) "Reading D3 matrices done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  WRITE(*,*) "Doing Backward FFT..."
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fc3, far)
  fc3%nq = nq
  WRITE(*,*) "Backward FFT done"
  CALL memstat(kb)
  WRITE(*,*) "Total memory used : ", kb/1000, "Mb"
  !
  IF(filename /="none")THEN
    WRITE(*,*) "Writing FCs to file..."
    CALL fc3%write(filename, S)
  ENDIF
  !
  IF(.not.skip_test)THEN
    WRITE(*,*) "Testing Forward FFT, with imaginary part..."
    WRITE(*,*) "(you can stop the code with CTRL-C to avoid running tests)  "
    CALL test_fwfft_d3(nq_trip, S, d3grid, fc3, .true., write_diff, "anh_cmplx")
    WRITE(*,*) "Testing Forward FFT, without imaginary part..."
    DEALLOCATE(fc3%ifc)
    CALL test_fwfft_d3(nq_trip, S, d3grid, fc3, .false., write_diff, "anh_real")
    WRITE(*,*) "Testing forward FFT done"
  ENDIF
  !
END PROGRAM qq2rr


