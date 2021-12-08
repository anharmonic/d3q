 ! This module contains subroutines and functions for compute_force.90 code
! 10th Feb., 2021
!
!-------------------------------------------------------------------------
MODULE harmonic_module
#include "mpi_thermal.h"
  !
  INTEGER, ALLOCATABLE   :: idx_R_map(:,:) 
  !
 CONTAINS
  !
  !
  !
 SUBROUTINE harmonic_force_md(n_steps,nat_tot,S,fc2,u_disp,h_force,h_energy)
  !----------------------------------------------------------------------
  ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
  ! at finite temperature md.out. 
  !
  USE Kinds, ONLY     : DP
  USE input_fc, ONLY  : read_fc2, forceconst2_grid, ph_system_info, & 
                          read_system, aux_system ! use only: .. to know precisely what comes from where, a good practice
  !USE harmonic_subroutines
  !
  IMPLICIT NONE
  !
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc2
  INTEGER,INTENT(in)    :: nat_tot, n_steps
  REAL(DP),INTENT(in) :: u_disp(3,nat_tot, n_steps)

  INTEGER :: i_step
  REAL(DP),ALLOCATABLE,INTENT(out)   :: h_force(:,:,:), h_energy(:)
  IF(nat_tot /= S%nat*fc2%n_R) CALL errore("harmonic_force", "nat_tot doers not match supercell",1)

  IF(.NOT.ALLOCATED(h_force)) ALLOCATE(h_force(3,nat_tot,n_steps))
  IF(.NOT.ALLOCATED(h_energy)) ALLOCATE(h_energy(n_steps))
  !
  ! compute force
  ! OPEN(331, file="F_harm.dat", status="unknown")
  h_force(:,:,:) = 0._dp
  DO i_step = 1, n_steps
    CALL harmonic_force(S,fc2,nat_tot,u_disp(:,:,i_step),h_force(:,:,i_step),h_energy(i_step))
  ENDDO

 ! 
 END SUBROUTINE 
 !
 SUBROUTINE harmonic_force(S,fc2,nat_tot,u_disp,h_force,h_energy)
  !----------------------------------------------------------------------
  ! Compute the harmonic force from DFPT FCs (mat2R) and molecular dynamics displacement
  ! at finite temperature md.out. 
  !
  USE Kinds, ONLY     : DP
  USE input_fc, ONLY  : read_fc2, forceconst2_grid, ph_system_info, & 
                          read_system, aux_system ! use only: .. to know precisely what comes from where, a good practice
  !USE harmonic_subroutines
  !
  IMPLICIT NONE
  !
  TYPE(ph_system_info)   :: S
  TYPE(forceconst2_grid) :: fc2
  INTEGER,INTENT(in)     :: nat_tot
  REAL(DP),INTENT(in)    :: u_disp(3,nat_tot)
  REAL(DP),INTENT(out)   :: h_force(3,nat_tot), h_energy

  INTEGER          :: i, j, k, jj, kk, nu, mu, beta, nat_sc, jat, kat, cR(3), alpha
  !
  ! Only allocate and compute the index map on first call
  IF(.not.ALLOCATED(idx_R_map)) THEN
    ALLOCATE(idx_R_map(fc2%n_R,fc2%n_R))
    ! Build map (i,j)->i such that 
    !    (R_j-R_k) = (R_i + RR),
    ! for R_i, R_j and R_k unit-cell vectors in the super-cell and RR a super-lattice vector
    idx_R_map = -1
    !fc2%i_0 
    DO j = 1, fc2%n_R
    DO k = 1, fc2%n_R
      cR = fc2%yR(:,j) - fc2%yR(:,k) ! R~ = R-R'
      cR(1) = MODULO( cR(1), fc2%nq(1))
      cR(2) = MODULO( cR(2), fc2%nq(2))
      cR(3) = MODULO( cR(3), fc2%nq(3))   ! R^ is R~ but taken inside the grid
      DO i = 1, fc2%n_R
        IF(ALL(cR==fc2%yR(:,i))) THEN
          ! is the index that we are looking for
          IF(idx_R_map(j,k)>0) CALL errore('map','found twice',1)
          idx_R_map(j,k) = i
          !EXIT ! <- stop the loop and exit
        ENDIF
      ENDDO
      IF(idx_R_map(j,k)==-1) CALL errore("harm_force", "could not find some R,R'", 1)
    ENDDO
    ENDDO
  ENDIF

  ! compute force
  ! OPEN(331, file="F_harm.dat", status="unknown")
  h_force(:,:) = 0._dp
  h_energy     = 0._dp
  jj=0 
  DO j = 1, fc2%n_R  ! loop over all atoms in the super cell
  DO jat = 1, S%nat
    jj=jj+1
    DO alpha = 1,3      
      nu = (jat-1)*3 + alpha
      ! ((j,jat)->jj, alpha) -> nu
      kk=0
      DO k = 1, fc2%n_R   !loop over index of vector of the cell in the supercell
      DO kat = 1, S%nat ! ... over all atoms in the unit cell
      kk=kk+1
      DO beta = 1,3 ! ... over cartesian axes x,y,z for each atom
        ! ((k,kat)->kk, beta) -> mu
        mu = (kat-1)*3 + beta
        i = idx_R_map(j,k)
        h_force(alpha,jj) = h_force(alpha,jj) &
                  - fc2%FC(nu,mu,i) * u_disp(beta,kk)
        h_energy = h_energy &
                  + 0.5_dp*fc2%FC(nu,mu,i) * u_disp(alpha,kk)*u_disp(beta,jj)

      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDDO
  ENDDO
 ! 
 END SUBROUTINE 


END MODULE
