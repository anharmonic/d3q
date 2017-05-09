!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Common for d3toten
!
MODULE d3_parameteres
  USE kinds, ONLY : DP
  REAL(DP) :: eps_delta = 1.e-5_dp
END MODULE d3_parameteres
! 
! Terms that are composed as ion-wise sum, radial FFT with form factor
!
MODULE d3_onecenter
  USE kinds, only: DP
  !
  TYPE d3_local_potential
    REAL(DP),POINTER    :: loc(:,:) => null()
  END TYPE d3_local_potential
  !TYPE(d3_local_potential),SAVE :: d3v(-3:3)
  TYPE(d3_local_potential),ALLOCATABLE :: d3v(:)
  !
  TYPE d3_nlcc_variation
    COMPLEX(DP),POINTER    :: drc(:,:) => null()
  END TYPE d3_nlcc_variation
  !TYPE(d3_nlcc_variation),SAVE :: d3c(-3:3)
  TYPE(d3_nlcc_variation),ALLOCATABLE :: d3c(:)

END MODULE d3_onecenter
!
!   the units of the files and the record lengths
!
MODULE d3_kgrid
  USE kinds, ONLY : DP
  INTEGER :: d3_nk1, d3_nk2, d3_nk3
  INTEGER :: d3_k1,  d3_k2,  d3_k3
  REAL(DP) :: d3_degauss
END MODULE d3_kgrid
!
!    third order dynamical matrix
!
MODULE d3_control
  USE kinds, only: DP
  !
  CHARACTER(len=8),PARAMETER :: code = 'D3_toten'
  CHARACTER(len=256) :: d3dir = '' ! scratch directory for d3 calculation (heavy I/O)
  !
  REAL(DP) :: ethr_ph ! eigenvalues convergence threshold
  INTEGER  :: istop
  LOGICAL  :: restart, safe_io
  LOGICAL  :: print_star, print_perm, print_trev
  character(len=64) :: fild3dyn
  INTEGER :: max_seconds = -1
END MODULE d3_control
!
! In the parallel version of the program some loop on perturbations
! may be split betweem pools. npert_i and npert_f are the initial
! and final value for a counter on the modes to be split among pools
!
MODULE npert_mod
  INTEGER :: &
       npert_i,    &! starting value for the mode counter
       npert_f      ! final value for the mode counter
END MODULE npert_mod
!
! Variables used for computing and writing only selected modes at q=0
! --the first index of the dthird matrix--
!
MODULE d3com
  use d3_parameteres
  use d3_onecenter
  use d3_control
  use npert_mod
  character(len=6),parameter :: version = "2.0.0b"
  character(len=64) :: d3_mode
END MODULE d3com
