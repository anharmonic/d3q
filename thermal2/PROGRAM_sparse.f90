!
! Written by Lorenzo Paulatto (2013-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! \/o\________\\\________________\\/\_________________________/^>
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
! \/o\________\\\________________\\/\_________________________/^>
PROGRAM gen_sparse

    USE kinds,          ONLY : DP
    USE fc3_interpolate,      ONLY : grid, sparse, forceconst3, fc3_grid_to_sparse, read_fc3
    USE input_fc,       ONLY : aux_system,  ph_system_info
    USE io_global,      ONLY : stdout
    USE nanoclock,      ONLY : nanotimer
    USE gen_sparse_program
    USE random_numbers, ONLY : randy
    USE mpi_thermal,    ONLY : start_mpi, stop_mpi, ionode
    USE more_constants,  ONLY : print_citations_linewidth
    USE clib_wrappers,        ONLY : memstat
    USE cmdline_param_module
    IMPLICIT NONE
    !
    TYPE(grid)   :: fc
    TYPE(sparse) :: sfc
    TYPE(ph_system_info)   :: S
    INTEGER      :: kb
    COMPLEX(DP),ALLOCATABLE :: D1(:,:,:), D2(:,:,:)
    !
    REAL(DP) :: xq1(3), xq2(3),x 
    
    CLASS(forceconst3),POINTER :: afc
    
    INTEGER :: i
    !
    TYPE(nanotimer) :: t_fc  = nanotimer("Full matrix form")
    TYPE(nanotimer) :: t_sfc  = nanotimer("Sparse matrix form")
    !
    CHARACTER(len=256) :: filein, fileout
    INTEGER ::  ntest
    REAL(DP) :: thr, delta, deltasum, deltamax
    !
    !
    CHARACTER(len=:),ALLOCATABLE :: cmdline

    filein   = cmdline_param_char("i", "mat3R")
    fileout  = cmdline_param_char("o", TRIM(filein)//".sparse")
    thr      = cmdline_param_dble("t", 0.d0)
    ntest    = cmdline_param_int("n", -1)
    
    IF (cmdline_param_logical('h')) THEN
        WRITE(*,*) "Syntax: d3_sparse.x [-i FILEIN] [-o FILEOUT] [-t THRESH] [-n NTESTS]"
        WRITE(*,*) ""
        WRITE(*,'(a)') "        FILEIN  : input 3rd order force constants in grid form (default: mat3R)"
        WRITE(*,'(a)') "        FILEOUT : output force constants in sparse form; use 'none' to avoid writing FILE_OUT"
        WRITE(*,'(a)') "                  (default: 'FILEIN.sparse')"
        WRITE(*,'(a)') "        THRESH  : force constants smaller than this times the largest FC will be ignored; "
        WRITE(*,'(a)') "                  can make the result faster, at the expense of some precision (default: 0.d0)"
        WRITE(*,'(a)') "        NTESTS  : put an integer number to test speed of grid vs. sparse over"
        WRITE(*,'(a)') "                  interpolation of n_test random q-points (default: no test)"
        STOP 1
    ENDIF
    CALL cmdline_check_exausted()

    IF(TRIM(fileout)==TRIM(filein)) &
      CALL errore("gen_sparse","filein and fileout are the same, I refuse to do that",1)

    CALL fc%read(filein, S)
    CALL aux_system(S)
    CALL memstat(kb)
    WRITE(stdout,*) "FC Memory used : ", kb/1000, "Mb"

    WRITE(stdout,*) "Cutting off FCs smaller than : ", thr, "Ry/bohr^3"

    CALL fc3_grid_to_sparse(S%nat, fc, sfc, thr)
    WRITE(stdout,*) "FC+Sparse Memory used : ", kb/1000, "Mb"
    IF(fileout/="none") CALL sfc%write(fileout, S)

    IF(ntest>0)THEN
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

    ENDIF

    CALL fc%destroy()
    WRITE(stdout,*) "Sparse Memory used : ", kb/1000, "Mb"
    
    CALL sfc%destroy()
    
    CALL print_citations_linewidth()
    
END PROGRAM gen_sparse
! \/o\________\\\________________\\/\_________________________/^>
