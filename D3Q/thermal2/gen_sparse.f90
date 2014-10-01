!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
MODULE gen_sparse_program
  USE kinds, ONLY : DP
  
  CONTAINS
  !
  SUBROUTINE test_classes(fc)
    USE sparse_fc
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
    USE sparse_fc,      ONLY : grid, sparse, forceconst3, fc3_grid_to_sparse, read_fc3
    USE input_fc,       ONLY : aux_system,  ph_system_info
    USE io_global,      ONLY : stdout
    USE nanoclock
    USE gen_sparse_program
    USE random_numbers, ONLY : randy
    USE wrappers,       ONLY : f_get_vm_size
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
    TYPE(nanotimer) :: t_fc  = nanotimer("t_fc")
    TYPE(nanotimer) :: t_sfc  = nanotimer("t_sfc")
    !
    INTEGER,INTRINSIC :: iargc
    CHARACTER(len=256) :: argx, filein, fileout
    INTEGER :: nargs, ntest
    !
    nargs = iargc()
    !
    CALL getarg(0,argx)
    IF(nargs>0) THEN
      CALL getarg(1, filein)
    ELSE 
      WRITE(*,*) "Syntax: "//TRIM(argx)//" FILE_IN [FILE_OUT] [n_test]"
      WRITE(*,*) "        FILE_IN  : input 3rd order force constants in grid form"
      WRITE(*,*) "        FILE_OUT : output force constants in sparse form; use 'none' to avoid writing FILE_OUT"
      WRITE(*,*) "                   (default: 'FILE_IN_sparse')"
      WRITE(*,*) "        n_test   : put an integer number to test speed of grid vs. sparse over"
      WRITE(*,*) "                   interpolation of n_test random q-points (default: no test)"
      
      CALL errore("asr3", "missing arguments")
    ENDIF
    IF(nargs>1)THEN
      CALL getarg(2, fileout)
    ELSE
      fileout = TRIM(filein)//"_sparse"
    ENDIF
    
!     afc => read_fc3("mat3R_sparse", S)

    CALL fc%read(filein, S)
    CALL memstat(kb)
    WRITE(stdout,*) "FC Memory used : ", kb/1000, "Mb"

    CALL fc3_grid_to_sparse(S%nat, fc, sfc, 1.d-8)
    WRITE(stdout,*) "FC+Sparse Memory used : ", kb/1000, "Mb"
    IF(fileout/="none") CALL sfc%write(fileout, S)

    IF(nargs>2)THEN
      CALL getarg(3, argx)
      READ(argx, *,iostat=i) ntest
      IF(i/=0) CALL errore("gen_sparse","bad test number, run without arguments for help",1)
      
      WRITE(*,*) "Running", ntest, "test configurations"
      
      ALLOCATE(D1(S%nat3, S%nat3, S%nat3))
      ALLOCATE(D2(S%nat3, S%nat3, S%nat3))

      x = randy(ntest)
      CALL start_nanoclock(t_fc)
      DO i = 1,ntest
        xq1 = (/ randy(), randy(), randy() /)
        xq2 = (/ randy(), randy(), randy() /)
        CALL fc%interpolate(xq1,xq2,S%nat3,D1)
      ENDDO
      CALL stop_nanoclock(t_fc)
      CALL print_nanoclock(t_fc)

      x = randy(ntest)
      CALL start_nanoclock(t_sfc)
      DO i = 1,ntest
        xq1 = (/ randy(), randy(), randy() /)
        xq2 = (/ randy(), randy(), randy() /)
        CALL sfc%interpolate(xq1,xq2,S%nat3,D2)
      ENDDO
      CALL stop_nanoclock(t_sfc)
      CALL print_nanoclock(t_sfc)

      DEALLOCATE(D1,D2)

    ENDIF

    CALL fc%destroy()
    WRITE(stdout,*) "Sparse Memory used : ", kb/1000, "Mb"
    
    CALL sfc%destroy()
    
END PROGRAM gen_sparse
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!