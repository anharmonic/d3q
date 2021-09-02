!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE interpolate2_module
  USE kinds, ONLY : DP
#include "mpi_thermal.h"
  IMPLICIT NONE

  CONTAINS
  PURE FUNCTION rotate_d2(nat3, D, U)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    COMPLEX(DP),INTENT(in) :: D(nat3,nat3), U(nat3,nat3)
    COMPLEX(DP) :: rotate_d2(nat3,nat3)
    rotate_d2 = matmul(transpose(conjg(U)), matmul(D,U))
  END FUNCTION
  PURE FUNCTION backrotate_d2(nat3, D, U)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nat3
    COMPLEX(DP),INTENT(in) :: D(nat3,nat3), U(nat3,nat3)
    COMPLEX(DP) :: backrotate_d2(nat3,nat3)
    backrotate_d2 = matmul(U, matmul(D,transpose(conjg(U))))
  END FUNCTION


  !
END MODULE


PROGRAM interpolate2
  USE constants, ONLY : RY_TO_CMM1
  USE input_fc,        ONLY : read_system, ph_system_info, read_fc2, write_dyn
  USE ph_system,       ONLY : aux_system
  USE interpolate2_module
  USE cmdline_param_module
  USE quter_module,       ONLY : quter
  USE input_fc,    ONLY : ph_system_info, forceconst2_grid, write_fc2, multiply_mass_dyn, div_mass_fc2
  USE fc2_interpolate,     ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE asr2_module,        ONLY : impose_asr2
  USE rigid, ONLY : rgd_blk

  IMPLICIT NONE
  INTEGER :: far, ios
  CHARACTER(len=256) :: filein, fileout
  TYPE(forceconst2_grid) :: fcin, fcout
  TYPE(ph_system_info) :: S
  INTEGER :: i,j,k, nqi, nqj, nqk, nq, na,nb,a,b, nua,nub, numax
  REAL(DP),ALLOCATABLE :: gridq(:,:), w2ref(:), w2mld(:), w2in(:), w2out(:)
  COMPLEX(DP),ALLOCATABLE :: matq(:,:,:,:,:), &
                             Din(:,:), Dref(:,:), Dmld(:,:), Dout(:,:), Dout2(:,:), &
                             U(:,:), Uref(:,:), Umld(:,:), Umld0(:,:), WWout(:,:)
  REAL(DP) :: xq(3), xqh(3), olap, maxolap
  LOGICAL :: lrigid_save
  !

  filein  = cmdline_param_char("i", "mat2R")
  fileout = cmdline_param_char("o", TRIM(filein)//".out")
  !
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: d3_interpolate2.x [-i FILEIN] [-o FILEOUT] [-r FILE.reference] [-m FILE.mould]"
      WRITE(*,*) ""
      WRITE(*,*) "Reads:"
      WRITE(*,*) " * initial force constants FCi from FILEIN (default: mat2R)"
      WRITE(*,*) " * reference FCr from FILE.reference (default mat2R.reference)"
      WRITE(*,*) " *  mould FCm from FILE.mould (default: mat2R.mould)"
      WRITE(*,*) "Interpolates the difference of FCi-FCr on the grid of FCm and writes the interpolated"
      WRITE(*,*) "(FCi-FCr)+FCm to the output file FILEOUT"
     
      STOP 1
  ENDIF
  !
!  cmdline = cmdline_residual()
  CALL cmdline_check_exausted()
  !
  IF(TRIM(fileout)==TRIM(filein)) &
    CALL errore("interpolate2","filein and fileout are the same, I refuse to do that",1)
  !
  WRITE(*,*) "Input file", TRIM(filein)
  WRITE(*,*) "Output file", TRIM(fileout)
 
  CALL read_fc2(filein,  S,  fcin)


  CALL aux_system(S)


  CALL impose_asr2("simple", S%nat, fcin)


  CALL div_mass_fc2(S,fcin)

  !
  nqi = fcin%nq(1)
  nqj = fcin%nq(2)
  nqk = fcin%nq(3)
  ALLOCATE(matq(3,3,S%nat,S%nat,nqi*nqj*nqk))
  ALLOCATE(gridq(3,nqi*nqj*nqk))
  ALLOCATE(Din(S%nat3,S%nat3))
  ALLOCATE(Dout(S%nat3,S%nat3))
  ALLOCATE(U(S%nat3,S%nat3))

  nq = 0
  DO i = 0,nqi-1
  DO j = 0,nqj-1
  DO k = 0,nqk-1
    nq = nq+1
    xq = s%bg(:,1)*i/DBLE(nqi) + s%bg(:,2)*j/DBLE(nqj) + s%bg(:,3)*k/DBLE(nqk)
    gridq(:,nq) = xq
    
    IF(nq==1) THEN
       xq(1) = 1.d-8
    ENDIF

    CALL fftinterp_mat2(xq, S, fcin, Din)
    Din = multiply_mass_dyn(S, Din)


    DO nb = 1,S%nat
    DO na = 1,S%nat
    DO b = 1,3
     nub = (nb-1)*3+b
     DO a = 1,3
      nua = (na-1)*3+a
      matq(a,b,na,nb, nq) = Din(nua,nub)
     ENDDO
    ENDDO
    ENDDO
    ENDDO
    !CALL rgd_blk (nqi,nqj,nqk,S%nat,matq(:,:,:,:,nq),gridq(:,nq), &
    !              S%tau,Smld%epsil,Smld%zeu,S%bg,S%omega,S%celldm(1), .false.,-1._dp)
  ENDDO
  ENDDO
  ENDDO

  S%lrigid    = .false.
  S%zeu = 0._dp
  S%epsil = 0._dp

  CALL quter(nqi, nqj, nqk, S%nat,S%tau,S%at,S%bg, matq, gridq, fcout, far=2)
  CALL write_fc2(fileout, S, fcout)
  !
  CONTAINS

END PROGRAM
!
