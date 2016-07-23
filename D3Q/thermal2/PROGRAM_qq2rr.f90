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
  USE iso_c_binding,   ONLY : c_int
  USE fc3_interpolate, ONLY : grid
  USE f3_bwfft,        ONLY : d3_list, read_d3_matrices, bwfft_d3_interp, test_fwfft_d3
  
  IMPLICIT NONE
  TYPE(d3_list),ALLOCATABLE ::  d3grid(:)
  TYPE(grid)                :: fc3
  INTEGER :: nq(3)
  INTEGER :: nq_trip, nq_grid
  INTEGER(kind=c_int)    :: kb
  !
  TYPE(ph_system_info) :: S
  COMPLEX(DP),ALLOCATABLE :: D3(:,:,:), P3(:,:,:)
  !
  INTEGER :: narg, i, j
  CHARACTER(len=512) :: argv, filename

  narg = command_argument_count()
  IF(narg>=3)THEN
    DO i = 1,3
      CALL get_command_argument(i,argv)
      READ(argv,*) nq(i)
    ENDDO
    IF(narg>=4)THEN
      CALL get_command_argument(i,filename)
    ELSE
      filename="mat3R"
    ENDIF
  ELSE
      WRITE(*,*) "Syntax: ls anh*| d3_qq2rr.x NQX NQY NQZ [FILENAME]"
      WRITE(*,*) ""
      WRITE(*,*) "Selects a grid of (NQX x NQY x NQZ) points from the anh* files"
      WRITE(*,*) "apply the inverse Fourier transform, and saves it to FILENAME."
      STOP 1
  ENDIF
  WRITE(*,*) "Reading grid", nq
  
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
!   DO i = 1,nq_trip
!     WRITE(2999,'(2(3f8.4,3x))') d3grid(i)%xq2, d3grid(i)%xq3
!     !WRITE(2999,*)
!   ENDDO
  !
  WRITE(*,*) "Doing Backward FFT..."
  CALL bwfft_d3_interp(nq, nq_trip, S%nat, S%tau, S%at, S%bg, d3grid, fc3)
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
  WRITE(*,*) "Testing Forward FFT, with imaginary part..."
  WRITE(*,*) "(you can stop the code with CTRL-C to avoid running tests)  "
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing Forward FFT, without imaginary part..."
  DEALLOCATE(fc3%ifc)
  CALL test_fwfft_d3(nq_trip, S, d3grid, fc3)
  WRITE(*,*) "Testing forward FFT done"
  !
END PROGRAM qq2rr


