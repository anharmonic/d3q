!
! Written by Lorenzo Paulatto (2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE interpolate2_module
  USE kinds, ONLY : DP
  USE functions,        ONLY : rotate_d2, backrotate_d2
#include "mpi_thermal.h"
  IMPLICIT NONE

  CONTAINS

  !
END MODULE


PROGRAM interpolate2
  !USE constants,       ONLY : RY_TO_CMM1
  !USE input_fc,        ONLY : read_system, ph_system_info, read_fc2, write_dyn, 
  USE ph_system,       ONLY : aux_system
  USE interpolate2_module
  USE cmdline_param_module
  USE quter_module,       ONLY : quter
  USE input_fc,           ONLY : read_fc2, ph_system_info, forceconst2_grid, write_fc2, allocate_fc2_grid
  !div_mass_fc2, multiply_mass_dyn
  USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat, fc2_recenter
  !USE asr2_module,        ONLY : impose_asr2
  !USE rigid_d3,           ONLY : rgd_blk_d3

  IMPLICIT NONE
  INTEGER :: far, ios
  REAL(DP),PARAMETER :: eps_r = 1.d-4 ! tolerance on atomic positions, in bohr 
  CHARACTER(len=256) :: filein_uc, filein_sc, fileout
  TYPE(forceconst2_grid) :: fc_uc, fc_sc, fc_out, fc_out_c
  TYPE(ph_system_info) :: S_uc, S_sc
  LOGICAL :: lrigid_save
  INTEGER :: nfar, method, sc_size(3), sc_factor, nat_uc
  INTEGER :: i, j, k, ir_uc, jr_sc, i_copy, j_copy, i_uc, j_uc, i_sc, j_sc, found
  INTEGER :: a,b,nu,mu,nu2,mu2
  INTEGER,ALLOCATABLE :: map_sc(:,:)
  REAL(DP) :: at_b(3,3), sc_b(3,3), bg_b_uc(3,3), bg_b_sc(3,3), R(3), dist(3), dist_sc(3), delta(3), weight, pref
  REAL(DP),ALLOCATABLE :: tau_b(:,:), Rtau_b(:,:)
  !

  filein_uc  = cmdline_param_char("u", "mat2R_uc")
  filein_sc  = cmdline_param_char("p", "mat2R_sc")
  fileout    = cmdline_param_char("o", "mat2R_out")
  weight     = cmdline_param_dble("w", -1._dp)
  nfar    = cmdline_param_int ("f", 2)
  method  = cmdline_param_int ("m", 3)
  !
  IF (cmdline_param_logical('h')) THEN
      WRITE(*,*) "Syntax: d3_interpolate2.x [-u FILEIN_unitcell] [-p FILEIN_perturbed_supercell] [-o FILEOUT]"
      WRITE(*,*) "                          [-f NFAR] [-m METHOD]"
      WRITE(*,*) ""
      STOP 1
  ENDIF
  !
!  cmdline = cmdline_residual()
  CALL cmdline_check_exausted()
  !
  IF(TRIM(fileout)==TRIM(filein_sc) .or. TRIM(fileout)==TRIM(filein_uc)) &
    CALL errore("interpolate2","filein and fileout are the same, I refuse to do that",1)
  !
  WRITE(*,*) "Input file", TRIM(filein_uc)
  WRITE(*,*) "Perturbed supercell file", TRIM(filein_sc)
  WRITE(*,*) "Output file", TRIM(fileout)
 
  CALL read_fc2(filein_uc,  S_uc,  fc_uc)
  CALL aux_system(S_uc)
  !CALL impose_asr2("simple", S_uc%nat, fc_uc)
  !CALL div_mass_fc2(S,fc_uc)
  !
  CALL read_fc2(filein_sc,  S_sc,  fc_sc)
  CALL aux_system(S_sc)
  !CALL impose_asr2("simple", S_uc%nat, fc_sc)
  !CALL div_mass_fc2(S,fc_sc)
  !
  CALL allocate_fc2_grid(fc_uc%n_R, S_uc%nat, fc_out)
  !
  at_b = S_uc%at * S_uc%alat
  CALL recips(at_b(:,1), at_b(:,2), at_b(:,3), bg_b_uc(:,1), bg_b_uc(:,2), bg_b_uc(:,3))
  sc_b = S_sc%at * S_sc%alat
  CALL recips(sc_b(:,1), sc_b(:,2), sc_b(:,3), bg_b_sc(:,1), bg_b_sc(:,2), bg_b_sc(:,3))
  !
  FORALL(i=1:3) sc_size(i) = NINT(NORM2(sc_b(:,i))/NORM2(at_b(:,i)))
  sc_factor = sc_size(1)*sc_size(2)*sc_size(3)
  print*, "Supercell factor", sc_factor, sc_size
  !
  IF(fc_sc%n_R>1) CALL errore('impurity','supercell must be at Gamma and FC periodic',1)

  ! lrigid_save = S_uc%lrigid
  ! S_uc%lrigid    = .false.
  ! S_sc%lrigid    = .false.
  !
  ! check that supercell is consistent
  nat_uc =  S_sc%nat/sc_factor
  IF(nat_uc /= S_uc%nat) CALL errore('impurity', 'supercell is inconsistent', 1)
  IF(nat_uc*sc_factor /= S_sc%nat) CALL errore('impurity', 'supercell is inconsistent', 2)

  ! unit cell atoms positions in bohr
  ALLOCATE(tau_b(3,S_uc%nat))
  tau_b = S_uc%tau * S_uc%alat
  ! supercell atoms positions in bohr
  ALLOCATE(Rtau_b(3,S_sc%nat))
  Rtau_b = S_sc%tau * S_sc%alat
  !
  ! build a map that for each atom in the UC (SC?) associates it's sc_factor copies in the SC
  ALLOCATE(map_sc(sc_factor,S_uc%nat))
  !
  ! for any atom of the SC find which other atoms of the SC are connected by a R vector of the UC lattice
  ! do not care about the atomic type (we are looking for an impurity)
  DO i_uc = 1, S_uc%nat
    i_copy = 0
    DO j_sc = 1, S_sc%nat
        delta = Rtau_b(:,i_uc)-Rtau_b(:,j_sc)
        CALL cryst_to_cart(1, delta, bg_b_uc, -1)
        !WRITE(*,'(2i4,(3f7.2,4x),3i3)') i_uc, j_sc, delta, NINT(delta)
        !
        IF( SUM(ABS(delta-NINT(delta))) < eps_r) THEN
          i_copy = i_copy+1
          IF(i_copy>sc_factor) CALL errore('impurity', 'found too many copies',1)
          map_sc(i_copy,i_uc) = j_sc
        ENDIF
    ENDDO
    IF(i_copy<sc_factor) CALL errore('impurity', 'did not find a copy',1)
  ENDDO
  !
  fc_out%FC = 0._dp
  !
  ! For every couple tau, tau+R in the unit-cell find its equivalent tau,tau in the supercell. 
  ! using periodic boudary condition can be required
  DO i_uc = 1, S_uc%nat
  DO j_uc = 1, S_uc%nat
  DO ir_uc = 1,fc_uc%n_R
    ! compute their distance vector
    dist = tau_b(:,j_uc)-fc_uc%xR(:,ir_uc)*S_uc%alat-tau_b(:,i_uc) !check
    !
    ! now go get the sc_factor atoms in the SC equivalent to the 1st one
    DO i_copy = 1, sc_factor
      found = 0
      i_sc = map_sc(i_copy,i_uc)
      ! and for each one, look the the ONLY other atom that respects the distance i-(j+R) of the UC+grid
      DO j_copy = 1, sc_factor
        j_sc = map_sc(j_copy,j_uc)
        dist_sc = Rtau_b(:,j_sc) - Rtau_b(:,i_sc)
        !WRITE(*,'(3i4,a,2i4, 2(3f7.2,4x))') i_uc, j_uc, ir_uc, "->", i_sc,j_sc, dist, dist_sc
        delta = dist-dist_sc
        CALL cryst_to_cart(1, delta, bg_b_sc, -1)
        !WRITE(*,'(3i4,a,2i4, 3(3f7.2,4x),3i3)') i_uc, j_uc, ir_uc, "->", i_sc,j_sc, dist, dist_sc, delta, NINT(delta)
        IF( SUM(ABS(delta-NINT(delta))) < eps_r) THEN
          !WRITE(*,'(3i4,a,3i4)') i_uc, j_uc, ir_uc, "->", i_sc,j_sc
          found = found+1
          DO a = 1,3
          DO b = 1,3
            nu = a + 3*(i_uc-1)
            mu = b + 3*(j_uc-1)
            nu2 = a + 3*(i_sc-1)
            mu2 = b + 3*(j_sc-1)
            fc_out%FC(nu,mu,ir_uc) = fc_out%FC(nu,mu,ir_uc) &
                                    +fc_sc%FC(nu2,mu2,1)
          ENDDO
          ENDDO
          ! EXIT ! when confident that the algorithm is ok, one could exit here and remove the check below
        ENDIF
      ENDDO
      IF(found>1) CALL errore('impurity', 'found too many', found)
      IF(found<1) CALL errore('impurity', 'not found enough', 1)
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  !
  print*, "Supercell force contants refolded to unit cell"
  fc_out%FC = fc_out%FC/sc_factor
  IF(weight/=0._dp)THEN
     fc_out%FC = fc_out%FC + weight*fc_uc%FC
     print*, "UC force constants added back with weight", weight
  ENDIF
  fc_out%xR = fc_uc%xR
  fc_out%yR = fc_uc%yR
  fc_out%nq = fc_uc%nq
  fc_out%periodic = fc_uc%periodic
  fc_out%centered = fc_uc%centered
  CALL write_fc2(fileout, S_uc, fc_out)
  
  CALL fc2_recenter(S_uc, fc_out, fc_out_c, 2)
  CALL write_fc2(TRIM(fileout)//'_c', S_uc, fc_out_c)
  !
  CONTAINS

END PROGRAM
!
