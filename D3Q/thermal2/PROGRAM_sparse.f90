!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
MODULE gen_sparse_program
  USE kinds, ONLY : DP
  
  CONTAINS
  !
  SUBROUTINE test_classes(fc)
    USE fc3_interpolate
    IMPLICIT NONE
    !
    CLASS(forceconst3),INTENT(in) :: fc
    !
    SELECT TYPE (fc)
    CLASS IS  ( forceconst3 )
      print*, "base"
    TYPE IS ( grid )
      print*,"grid"
    TYPE IS ( sparse )
      print*,"sparse"
!     DEFAULT TYPE
!       print*,4
    END SELECT
    !
  END SUBROUTINE
END MODULE gen_sparse_program

PROGRAM gen_sparse

    USE kinds,          ONLY : DP
    USE iso_c_binding,  ONLY : c_int, c_long
    USE fc3_interpolate,      ONLY : grid, sparse, forceconst3, fc3_grid_to_sparse, read_fc3
    USE input_fc,       ONLY : aux_system,  ph_system_info
    USE io_global,      ONLY : stdout
    USE nanoclock,      ONLY : nanotimer, print_timers_header
    USE gen_sparse_program
    USE random_numbers, ONLY : randy
    USE wrappers,       ONLY : f_get_vm_size
    USE mpi_thermal,    ONLY : start_mpi, stop_mpi
    IMPLICIT NONE
    !
    TYPE(grid)   :: fc
    TYPE(sparse) :: sfc
    TYPE(ph_system_info)   :: S
    INTEGER(kind=c_int)    :: kb
    COMPLEX(DP),ALLOCATABLE :: D1(:,:,:), D2(:,:,:)
    !
    REAL(DP) :: xq1(3), xq2(3),x 
    
    CLASS(forceconst3),POINTER :: afc
    
    INTEGER :: i
    !
    TYPE(nanotimer) :: t_fc  = nanotimer("Full matrix form")
    TYPE(nanotimer) :: t_sfc  = nanotimer("Sparse matrix form")
    !
    INTEGER,INTRINSIC :: iargc
    CHARACTER(len=256) :: argx, filein, fileout
    INTEGER :: nargs, ntest
    REAL(DP) :: thr, delta, deltasum, deltamax
    !
    nargs = iargc()
    !
    CALL getarg(0,argx)
    IF(nargs>0) THEN
      CALL getarg(1, filein)
    ELSE 
      WRITE(*,'(a)') "Syntax: "//TRIM(argx)//" FILE_IN [FILE_OUT] [threshold] [n_test]"
      WRITE(*,'(a)') "        FILE_IN  : input 3rd order force constants in grid form"
      WRITE(*,'(a)') "        FILE_OUT : output force constants in sparse form; use 'none' to avoid writing FILE_OUT"
      WRITE(*,'(a)') "                   (default: 'FILE_IN_sparse')"
      WRITE(*,'(a)') "        threshold: force constants smaller than this times the largest FC will be ignored; "
      WRITE(*,'(a)') "                   can make the result faster, at the expense of some precision (default: 0.d0)"
      WRITE(*,'(a)') "        n_test   : put an integer number to test speed of grid vs. sparse over"
      WRITE(*,'(a)') "                   interpolation of n_test random q-points (default: no test)"
      STOP
      !CALL errore("asr3", "missing arguments", 1)
    ENDIF
    IF(nargs>1)THEN
      CALL getarg(2, fileout)
    ELSE
      fileout = TRIM(filein)//"_sparse"
    ENDIF
    
!     afc => read_fc3("mat3R_sparse", S)

    CALL fc%read(filein, S)
    CALL aux_system(S)
    CALL memstat(kb)
    WRITE(stdout,*) "FC Memory used : ", kb/1000, "Mb"

    IF(nargs>2)THEN
      CALL getarg(3, argx)
      READ(argx, *,iostat=i) thr
      IF(i/=0) CALL errore("gen_sparse","bad threshold, run without arguments for help",1)
    ELSE
      thr = 0._dp
    ENDIF
    !thr = thr*MAXVAL(ABS(fc%FC))
    WRITE(stdout,*) "Cutting off FCs smaller than : ", thr, "Ry/bohr^3"

    CALL fc3_grid_to_sparse(S%nat, fc, sfc, thr)
    WRITE(stdout,*) "FC+Sparse Memory used : ", kb/1000, "Mb"
    IF(fileout/="none") CALL sfc%write(fileout, S)

    IF(nargs>3)THEN
      !CALL start_mpi()

      CALL getarg(4, argx)
      READ(argx, *,iostat=i) ntest
      IF(i/=0) CALL errore("gen_sparse","bad test number, run without arguments for help",1)
      
      WRITE(*,*) "Running", ntest, "test configurations"
      
      ALLOCATE(D1(S%nat3, S%nat3, S%nat3))
      ALLOCATE(D2(S%nat3, S%nat3, S%nat3))

      deltasum=0._dp
      deltamax=0._dp

      x = randy(ntest)
      DO i = 1,ntest
        xq1 = (/ randy(), randy(), randy() /)
        xq2 = (/ randy(), randy(), randy() /)
        
        CALL t_fc%start()
          CALL fc%interpolate(xq1,xq2,S%nat3,D1)
        CALL t_fc%stop()

        CALL t_sfc%start()
          CALL sfc%interpolate(xq1,xq2,S%nat3,D2)
        CALL t_sfc%stop()

        delta=MAXVAL(ABS(D1-D2))
        deltasum=deltasum+delta
        deltamax=MAX(delta,deltamax)

      ENDDO

      DEALLOCATE(D1,D2)

      WRITE(*,*) "Max delta: ", deltamax, deltamax/DSQRT(DBLE(ntest))
      WRITE(*,*) "Avg delta: ", deltasum/DBLE(ntest)
      WRITE(*,*) "Speedup:   ", t_fc%tot/t_sfc%tot

      CALL t_fc%print()
      CALL t_sfc%print()

      !CALL stop_mpi()
    ENDIF

    CALL fc%destroy()
    WRITE(stdout,*) "Sparse Memory used : ", kb/1000, "Mb"
    
    CALL sfc%destroy()
    
END PROGRAM gen_sparse
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
