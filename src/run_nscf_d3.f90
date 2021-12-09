!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE nscf_d3
!----------------------------------------------------------------------------

  PUBLIC :: run_nscf_d3

  PRIVATE

 CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE reload_nscf_d3()
  !-----------------------------------------------------------------------
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvecs,            ONLY : doublegrid
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : v, vltot, vrs, kedtau
  !
  WRITE(stdout, '(/,5x,a)') "Loading non-scf wavefunctions from interrupted D3 calculation."
  CALL setup_nscf_d3()
  CALL clean_pw( .true. )
  CALL close_files( .true. )
  CALL read_file()
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE run_nscf_d3(do_band)
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the pwscf program called from the
  ! ... phonon code. 
  !
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : at
  USE control_flags,   ONLY : conv_ions, pw_restart=> restart
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir
  USE io_global,       ONLY : stdout, ionode
  USE input_parameters,ONLY : pseudo_dir
  USE save_ph,         ONLY : tmp_dir_save
  USE d3_restart,      ONLY : done_nscf, d3_check_restart
  USE d3matrix_io2,    ONLY : d3matrix_filename2
  USE kplus3q,         ONLY : kplusq
  USE d3_control,      ONLY : d3dir
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE fft_base,        ONLY : dffts, dfftp
  USE fft_types,       ONLY : fft_type_allocate
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft

  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  CHARACTER(len=512) :: fileout
  !
  INTEGER :: stdout_tmp = -1
  !
  ! ---------> restart begin
  ! If we are restarting and wfcs had been computed, just reload them
  IF(done_nscf)THEN
    CALL reload_nscf_d3()
    RETURN
  ENDIF
  ! <--------- restart end
  !
  ! DO NSCF CALCULATION OF WAVEFUNCTIONS AT K, K+Q1, K+Q2, K+Q3, K-Q1, K-Q2, K-Q3
  !
  WRITE(stdout, '(/,5x,a)') "Starting non-scf calculation of ground state wavefunctions."
  !
  CALL start_clock( 'nscf' )
  !
  CALL clean_pw( .false. )
  !
  CALL close_files( .true. )
  !
  ! ... Setting the values for the nscf run
  !
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  !pw_restart        = .false.
  pseudo_dir= TRIM(tmp_dir_save) // TRIM(prefix) // '.save'
  !CALL restart_from_file()
  conv_ions=.true.
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  CALL setup_nscf_d3()
  !
  IF(ionode)THEN
    fileout = TRIM(d3dir)//"/"//&
      TRIM(d3matrix_filename2(kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq,&
                             at, 'nscf'))//'.out'
    WRITE(stdout, '(/,5x,a)') "--> output from 'electrons' written to "//TRIM(fileout)
    FLUSH( stdout )
    stdout_tmp = stdout
    stdout = 900
    OPEN(unit= stdout, file=TRIM(fileout), access='sequential') !, POSITION='append')
  ENDIF
  !
  ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  CALL init_run()
  CALL non_scf()
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  !
  IF(ionode)THEN
    CLOSE(stdout)
    stdout=stdout_tmp
  ENDIF
  !
  !twfcollect=.false.
  CALL punch( 'all' )
  !
  CALL close_files( .true. )
  !
  ! Store restart information
  done_nscf = .true.
  CALL d3_check_restart('write')
  !
  CALL stop_clock( 'nscf' )
  !
  RETURN
  !----------------------------------------------------------------------------
END SUBROUTINE run_nscf_d3
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE setup_nscf_d3()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes variables for the non-scf calculations at k 
  ! ... and k+q required by the linear response calculation at finite q.
  ! ... In particular: finds the symmetry group of the crystal that leaves
  ! ... the phonon q-vector (xq) or the single atomic displacement (modenum)
  ! ... unchanged; determines the k- and k+q points in the irreducible BZ
  ! ... Needed on input (read from data file):
  ! ... "nsym" crystal symmetries s, ftau, t_rev, "nrot" lattice symetries "s"
  ! ... "nkstot" k-points in the irreducible BZ wrt lattice symmetry
  ! ... Produced on output:
  ! ... symmetries ordered with the "nsymq" phonon symmetries first
  ! ... "nkstot" k- and k+q-points in the IBZ calculated for the phonon sym.)
  ! ... Misc. data needed for running thte non-scf calculation
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg
  USE ions_base,          ONLY : nat, ityp
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symm_base,          ONLY : s, t_rev, nrot, time_reversal, &
                                 copy_sym, s_axis_to_cart
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, max_cg_iter, &
                                 noinv !, use_para_diag
  USE mp_pools,           ONLY : kunit, inter_pool_comm
  !USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin,domag
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3, reset_grid
  USE klist,              ONLY : degauss
  USE lr_symm_base,       ONLY : nsymq, invsymq, minus_q !, gi, gimq, irgq, irotmq, minus_q
  USE kplus3q,            ONLY : kplus3q_grids, kplusq
  USE mp,                 ONLY : mp_sum
  USE upf_ions,           ONLY : n_atom_wfc
  USE d3_kgrid,           ONLY : d3_nk1, d3_nk2, d3_nk3, d3_k1, d3_k2, d3_k3, d3_degauss
  USE check_stop,         ONLY : check_stop_init
  USE control_lr,      ONLY : ethr_nscf
  !
  IMPLICIT NONE
  !
  LOGICAL  ::  magnetic_sym, skip_equivalence, newgrid
  INTEGER :: iq=0
!  LOGICAL, EXTERNAL  :: check_para_diag
  !
  CALL check_stop_init()
  !
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr - should be good for all cases
  !
  ethr= 1.0D-12 / nelec
  ethr_nscf = 1.0D-12 / nelec
  !ethr = ethr_ph
  !
  ! ... variables for iterative diagonalization (Davidson is assumed)
  !
  isolve = 0
  david = 4
  max_cg_iter=100
  nbndx = david*nbnd
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  
  !lkpoint_dir=.false. ! to avoid very large number of directories
  IF(d3_degauss>0._dp) THEN
    degauss = d3_degauss
    WRITE(stdout, '(7x,a,f10.4)') "* new degauss set from d3 input:", degauss
  ENDIF
  !
!  use_para_diag = check_para_diag( nbnd )
  !
  ! ... Symmetry and k-point section
  !
  ! ... time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  IF(time_reversal) WRITE(stdout, "(7x,a)") "* has time reversal symmetry."
  !
  ! check if inversion (I) is a symmetry. If so, there should be nsymq/2
  ! symmetries without inversion, followed by nsymq/2 with inversion
  ! Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  !
  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
  !
  !  Since the order of the s matrices is changed we need to recalculate:
  !
!  call s_axis_to_cart () 
!   IF (okpaw) CALL d_matrix(d1,d2,d3)
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  !
  newgrid = reset_grid (d3_nk1, d3_nk2, d3_nk3, d3_k1, d3_k2, d3_k3)
  IF( nks_start > 0 .and. .not. newgrid ) then
     !
     !  In this case I keep the same points of the Charge density
     !  calculations
     !
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nkstot)
     wk(1:nkstot)   = wk_start(1:nkstot)
  ELSE
     !
     ! In this case I generate a new set of k-points
     skip_equivalence = .false.
     CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                        bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
  ENDIF
  !
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  ! FIXME: minus_q should work now.. but it does not
  CALL irreducible_BZ(nrot, s, nsymq, minus_q, magnetic_sym, at, bg, npk, &
                      nkstot, xk, wk, t_rev)
  !
  ! ... add k+q to the list of k
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = MAXVAL((/ SQRT(SUM(kplusq(1)%xq**2)), &
                    SQRT(SUM(kplusq(2)%xq**2)), &
                    SQRT(SUM(kplusq(3)%xq**2)) /) )
  !
!#ifdef __MPI
  !
  ! Divide original symmetry-opened kpoints among pools
  kunit = 1
  !CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  WRITE(stdout,'(5x,a)') "Total number of kpoints:"
  WRITE(stdout,'(7x,a,i6)') "--> after reducing symmetry:", nkstot
  !
  ! Expand kpoints to k,k+q1,k+q2,k+q3,k-q1,k-q2,k-q3 and write a igk file
  ! for each k+qX sub-grid
  CALL kplus3q_grids( xk, wk, nks, npk, nat, kunit )
  !
  nkstot=nks
  CALL mp_sum(nkstot, inter_pool_comm )
  !
  WRITE(stdout,'(7x,a,i6," (",i1,")")') "--> after including k+/-q_i (kunit):", nkstot, kunit
  !
!#else
!  !
!  nks = nkstot
!  !
!#endif
!   IF ( lsda ) THEN
!      !
!      ! ... LSDA case: two different spin polarizations,
!      ! ...            each with its own kpoints
!      !
!      if (nspin /= 2) call errore ('setup','nspin should be 2; check iosys',1)
!      !
!      CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
!      !
!   ELSE IF ( noncolin ) THEN
!      !
!      ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
!      !
!      if (nspin /= 4) call errore ('setup','nspin should be 4; check iosys',1)
!      current_spin = 1
!      !
!   ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk = wk * degspin
     DO iq = -3,3
       !print*, allocated(kplusq(iq)%wk)
       kplusq(iq)%wk = kplusq(iq)%wk *degspin
     ENDDO
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'd3_setup', 'nspin should be 1; check iosys', 1 )
     !
!   END IF
  !
  !
  IF ( nkstot > npk ) CALL errore( 'setup', 'too many k points', nkstot )
  !
  !
  RETURN
  !
  !----------------------------------------------------------------------------
END SUBROUTINE setup_nscf_d3
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
END MODULE nscf_d3
!----------------------------------------------------------------------------
