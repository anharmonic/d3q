!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE linter_d3q

  USE kinds, ONLY : DP

  ! This variables need to stay here because they are used by ch_psi_all2,
  ! which is called by cgsolve_all. ch_psi_all2 is passed to cgsolve_all
  ! as an interface argument and is invoked as h_psi inside there.
  ! Actually it would be sufficient to keep igk_prj in the module, igk_wfc
  ! is only kept for concistency.
  ! TODO: use the same trick for vkb which is currently picked by
  ! ch_psi_all2 from the uspp module.
  !
  USE d3_h_psi, ONLY :  igk_wfc, igk_prj, vkb_wfc, vkb_prj, psi_prj, psi_wfc, g2kin_prj
  !
! #define PRECONDITIONING_FACTOR 3._dp 1.5_dp
#define PRECONDITIONING_FACTOR 1.35_dp
 !original values, e.g. in phonon code: 1.35_dp

CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE generate_dwfc2()
  !-----------------------------------------------------------------------
  !
  ! Driver to calculate and writes  | P_c d/du(q_x) psi(k+q_y) >
  !
  USE io_global,   ONLY : stdout !, ionode
  USE pwcom,       ONLY : degauss
  USE d3_basis,    ONLY : patq
  USE d3_symmetry, ONLY : symq
  USE d3_open,     ONLY : reopn_d3
  USE d3_control,  ONLY : safe_io
  USE kplus3q,     ONLY : kplusq, q_names, q_names2, q_sum_rule
  USE d3_iofiles,  ONLY : iu_dwfc, iu_psi_dH_psi
  USE d3_efermi_shift,ONLY : write_efsh
  USE mp,          ONLY : mp_barrier
  USE mp_world,    ONLY : world_comm
  !
  USE d3_restart, ONLY : done_dwfc, done_pdvp, done_lmetq0, d3_check_restart
  !
  IMPLICIT NONE
  !
  INTEGER irr, irr1, imode0
  ! switch
  ! the number of irreducible representation
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  INTEGER :: iq_wfc_x, iq_prt_x, iq_prj_x, iq_prj, iq_prt, iq_wfc, idwfc
  INTEGER :: unit_psidqvpsi, unit_dpsi
  LOGICAL :: lmetq0, ldwfc
  INTEGER,PARAMETER :: max_dwfc_todo = 12
  INTEGER,PARAMETER :: max_dH_todo = 6
  INTEGER,PARAMETER :: dwfc_todo(2,max_dwfc_todo+max_dH_todo) = &
        RESHAPE( (/ &
          1,0,  1,-1, -1,0, -1,1, &  !
          2,0,  2,-2, -2,0, -2,2, &  ! <-- do dpsi and <psi|dV|psi>
          3,0,  3,-3, -3,0, -3,3, &  !
          1,2, 1,3, 2,1, 2,3, 3,1, 3,2 & ! <-- only do <psi|dV|psi>
         /), (/2,max_dwfc_todo+max_dH_todo/) )
  ! Note: the <psi|dV|psi> term for iq_prt < 0 is never actually used, we keep
  ! it as it is almost free (may be fixed in the future)
  !
  CALL start_clock('generate_dwfc')
  !
  DO idwfc = 1, max_dwfc_todo+max_dH_todo
    iq_prt = dwfc_todo(1,idwfc)
    iq_prt_x = kplusq(iq_prt)%copy_of
    !
    iq_wfc = dwfc_todo(2,idwfc)
    iq_wfc_x = kplusq(iq_wfc)%copy_of
    !
    iq_prj = q_sum_rule(iq_wfc,iq_prt)
    iq_prj_x = kplusq(iq_prj)%copy_of
    !
    ldwfc = (idwfc <= max_dwfc_todo)
    ! Special rules to compute only 3 dwfc (instead of 5) when one of the q-points is Gamma
    IF(kplusq(1)%lgamma .and. (iq_prt == 3 .or. iq_prt_x==3)) ldwfc = .false.
    IF(kplusq(2)%lgamma .and. (iq_prt == 3 .or. iq_prt_x==3)) ldwfc = .false.
    IF(kplusq(3)%lgamma .and. (iq_prt == 2 .or. iq_prt_x==2)) ldwfc = .false.
    !
    IF( .not.  done_pdvp(iq_wfc_x, iq_prt_x) ) THEN !.or. &
!         .not. (done_dwfc(iq_wfc_x, iq_prt_x).and.ldwfc) ) THEN
      !
      lmetq0 = kplusq(iq_wfc)%lgamma .and. kplusq(iq_prt)%lgamma &
                .and. (degauss /= 0._dp) .and. (.not. done_lmetq0)
      !
      IF (ldwfc) THEN
        WRITE(stdout, '(/,5x,"Computing P_c^",a," |d^",a," psi_k",a,">")') &
           TRIM(q_names(iq_prj)), TRIM(q_names(iq_prt)), TRIM(q_names2(iq_wfc))
       ELSE
        WRITE(stdout, '(/,5x,"Computing <psi_k",a," |d^",a,"V| psi_k",a,">")') &
           TRIM(q_names2(iq_prj)), TRIM(q_names(iq_prt)), TRIM(q_names2(iq_wfc))
       ENDIF
      !
      unit_dpsi      = iu_dwfc(iq_wfc, iq_prt)
      unit_psidqvpsi = iu_psi_dH_psi(iq_wfc, iq_prt)
      !
      DO irr = 1, symq(ABS(iq_prt))%nirr
        imode0 = 0
        DO irr1 = 1, irr - 1
            imode0 = imode0 +  symq(ABS(iq_prt))%npert(irr1)
        ENDDO
        !
        CALL solve_linter_d3q (irr, imode0, symq(ABS(iq_prt))%npert(irr), &
                                iq_wfc, iq_prj, iq_prt, patq(iq_prt)%u, kplusq(iq_prt)%xq, &
                                unit_psidqvpsi, unit_dpsi, lmetq0, ldwfc)
        !
        CALL mp_barrier(world_comm)
        !
      ENDDO
      !
      ! close and reopen dpsi and psidHpsi files, this is necessary for restart
      IF(safe_io)THEN
        CALL reopn_d3(unit_dpsi)
        CALL reopn_d3(unit_psidqvpsi)
      ENDIF
      !
      IF ( lmetq0 ) THEN 
        WRITE(stdout,'(7x,"Gamma-Gamma perturbation and for metal: E_Fermi shift has been computed.")')
        CALL write_efsh()
       !
      ENDIF
      !
      ! Restart info:
      done_pdvp(iq_wfc, iq_prt)     = .true.
      done_pdvp(iq_wfc_x, iq_prt_x) = .true.
      done_dwfc(iq_wfc, iq_prt)     = ldwfc
      done_dwfc(iq_wfc_x, iq_prt_x) = ldwfc
      IF(lmetq0) done_lmetq0=.true.
      !
    ELSE
      ! Restart info:
      done_pdvp(iq_wfc, iq_prt) = .true.
      done_dwfc(iq_wfc, iq_prt) = ldwfc
      !
      IF (ldwfc) THEN
        WRITE(stdout, '(/,5x,"Skipping P_c^",a," |d^",a," psi_k",a,">")') &
              TRIM(q_names(iq_prj)), TRIM(q_names(iq_prt)), TRIM(q_names2(iq_wfc))
        WRITE(stdout, '(7x,"--> already done as P_c^",a," |d^",a," psi_k",a,">")') &
              TRIM(q_names(iq_prj_x)), TRIM(q_names(iq_prt_x)), TRIM(q_names2(iq_wfc_x))
      ELSE
        WRITE(stdout, '(/,5x,"Skipping <psi_k",a," |d^",a,"V| psi_k",a,">")') &
           TRIM(q_names2(iq_prj)), TRIM(q_names(iq_prt)), TRIM(q_names2(iq_wfc))
        WRITE(stdout, '(7x,"--> already done as <psi_k",a," |d^",a,"V| psi_k",a,">")') &
           TRIM(q_names2(iq_prj_x)), TRIM(q_names(iq_prt_x)), TRIM(q_names2(iq_wfc_x))

      ENDIF
    ENDIF
    !
    ! Store restart info:
    CALL d3_check_restart('write')
    !
  ENDDO
  !
  CALL stop_clock('generate_dwfc')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE generate_dwfc2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter_d3q (irr, imode0, npe, iq_wfc, iq_prj, iq_prt, &
                             u_prt, xq_prt, unit_psidqvpsi, unit_dpsi,lmetq0,ldwfc)
  !-----------------------------------------------------------------------
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to the perturbation.
  !    It reads from a file the charge variation due to perturbation
  !    and calculates variation of the wavefunctions.
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat
  USE cell_base,      ONLY : tpiba2
  USE fft_base,       ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fft_interpolate
  USE io_global,      ONLY : stdout
  USE gvect,          ONLY : g, gstart
  USE ener,           ONLY : ef
  USE klist,          ONLY : xk, degauss, ngauss
  USE wvfct,          ONLY : nbnd, npwx, et
  USE uspp,           ONLY : nkb
  USE control_lr,     ONLY : nbnd_occ
  USE qpoint,         ONLY : nksq
  USE units_lr,       ONLY : iuwfc, lrwfc, lrdwf
  USE d3com,          ONLY : ethr_ph
  USE kplus3q,        ONLY : kplusq, q_sum_rule, nbnd_max
  USE mp_pools,       ONLY : inter_pool_comm, intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE d3_iofiles,     ONLY : lrpdqvp
  USE uspp_init,      ONLY : init_us_2
  !
  ! D3 subroutines called:
  USE dvdpsi_module
  USE incdrhoscf2_module
  USE d3_efermi_shift, ONLY : set_efsh
  USE dq_vscf_module
  USE d3_h_psi, ONLY : d3_ch_psi, d3_cg_psi, nproj

! debug: i.e. all the stuff that is used inside the phonon solver and that must be set correctly or it won't work
!   USE becmod,     ONLY : becp
!   USE lsda_mod,   ONLY : current_spin
!   USE gvecs,      ONLY : nls
!   USE fft_base,   ONLY : dffts
!   USE scf,        ONLY : vrs

  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: irr, npe, imode0
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes
  ! input: a switch
  INTEGER,INTENT(IN) :: iq_wfc, iq_prj, iq_prt
  ! input: the indexes of the q vectors to be used, they refer
  ! to the projector (prj), the perturbation (prt) and the unperturbet
  ! psi used (wfc). Their values are as following:
  ! 1,2,3 -> q_1, q_2, q_3
  ! 0 -> gamma
  ! -1,-2,-3 -> -q1, -q2, -q3
  !
  INTEGER,INTENT(IN)  :: unit_psidqvpsi, unit_dpsi
  REAL(DP),INTENT(IN) :: xq_prt (3)
  COMPLEX(DP),INTENT(IN) :: u_prt(3*nat,3*nat)
  LOGICAL,INTENT(IN) :: lmetq0 ! if true, fermi enegry shift is computed
  LOGICAL,INTENT(IN) :: ldwfc  ! if false dpsi is NOT computed, only the
                               ! <psi|dV|psi> terms are. 

  REAL (DP) :: thresh, wg1, wg2, wwg, deltae, anorm, averlt, &
               eprec1, aux_avg (2), tcpu,tcpu0, degaussm1=0._dp
  ! the convergence threshold
  ! weight for metals
  ! weight for metals
  ! weight for metals
  ! difference of energy
  ! the norm of the error
  ! average number of iterations
  ! cut-off for preconditioning
  ! auxiliary variable for avg. iter. coun

  REAL(DP),EXTERNAL :: w0gauss, wgauss, get_clock
  ! function computing the delta function
  ! function computing the theta function
  ! cpu time

  COMPLEX(DP) ::  ps (nbnd), psidvpsi
  ! the scalar products
  ! auxiliary dpsi dV matrix element between k+q  and  k wavefunctions
  COMPLEX(DP),EXTERNAL ::  ZDOTC

  REAL(DP),ALLOCATABLE :: h_diag (:,:)
  ! the diagonal part of the Hamiltonian
  COMPLEX(DP),ALLOCATABLE :: drhoscf (:,:), dvloc (:,:),  &
                             spsi (:), auxg (:), dpsiaux (:,:)
  ! the variation of the charge
  ! variation of local part of the potential
  ! the function spsi
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:), dpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: psidqvpsi(:,:)
  logical :: conv_root

  INTEGER :: ipert, ibnd, jbnd, lter, ltaver, lintercall, ik, ik_wfc, &
             ik_prj, ig, nrec, ios, mode
  !
  CHARACTER(len=16),PARAMETER :: sub='solve_linter_d3q'
  !CHARACTER(len=256) :: fmt_str
  !
  INTEGER         :: npw_wfc, npw_prj
  !
  LOGICAL :: lmetal, lhamiltonian
  !
  CALL start_clock ('solve_linter')
  tcpu0 = get_clock ('D3_toten')
  !
  ! Switches
  !
  lmetal = (degauss /= 0._dp)
  IF (lmetal) degaussm1 = 1._dp / degauss
  lhamiltonian = .true.
  !
  ! Stuff I don't know
  ltaver = 0
  lintercall = 0
  thresh = ethr_ph
  !
  ! FIXME!
  IF(allocated(vkb_prj)) DEALLOCATE(vkb_prj)
  !
  ALLOCATE( drhoscf(dfftp%nnr, npe) )
  ALLOCATE( dvloc(dfftp%nnr, npe) )
  ALLOCATE( spsi(npwx) )
  ALLOCATE( auxg(npwx) )
  !
  ALLOCATE(dpsi(npwx, nbnd))
  ALLOCATE(dvpsi(npwx, nbnd))
  !
  ALLOCATE(psi_wfc(npwx, nbnd), psi_prj(npwx, nbnd))
  ALLOCATE(vkb_wfc(npwx, nkb),  vkb_prj(npwx, nkb))
  ALLOCATE(igk_wfc(npwx),       igk_prj(npwx))
  !
  ALLOCATE(psidqvpsi(nbnd, nbnd))
  !
  IF (lmetal.and.ldwfc) ALLOCATE( dpsiaux(npwx, nbnd) )
  ALLOCATE( h_diag(npwx, nbnd) )
  !
  ! calculates the variation of the local part of the K-S potential
  !
  DO ipert = 1, npe
     mode = imode0 + ipert
     CALL dq_vscf(mode, dvloc(:,ipert), xq_prt, iq_prt, u_prt)
  ENDDO
  !
  ! fixme:
  drhoscf = (0._dp, 0._dp)
  !
  KPOINTS_LOOP : &
  DO ik = 1, nksq
    !
    ! Find the real index of the kpoints to be processed, they will be
    ! the same if (iq_wfc == iq_prj) 
    ! OR the two kpoint are equal by chance (<-- not yet implemented)
    !
    !
    ! read the number of plane waves and their ordering for the unperturbed wavefunction
    ! at k+q_wfc, then read the wavefunction itself
    !READ (kplusq(iq_wfc)%iunigkq, iostat = ios) npw_wfc, igk_wfc
    !print*, ik, allocated(kplusq(iq_wfc)%ngkq)
    ik_wfc  = kplusq(iq_wfc)%ikqs(ik)
    npw_wfc = kplusq(iq_wfc)%ngkq(ik)
    igk_wfc = kplusq(iq_wfc)%igkq(:,ik)
    !
    ik_prj  = kplusq(iq_prj)%ikqs(ik)
    npw_prj = kplusq(iq_prj)%ngkq(ik)
    igk_prj = kplusq(iq_prj)%igkq(:,ik)
    !
    CALL davcio (psi_wfc, lrwfc, iuwfc, ik_wfc, -1)
    !call get_buffer (evc, lrwfc, iuwfc, ikk)
    !
    ! calculate the variation of the non-local part of the K-S potential
    CALL init_us_2 (npw_wfc, igk_wfc, xk(1, ik_wfc), vkb_wfc)
    !
      ! If the two q vectors are equal just copy
    IF (ik_wfc == ik_prj) THEN
      psi_prj = psi_wfc
      vkb_prj = vkb_wfc
    ELSE
      ! Otherwise, we have to repeat for q_prj
      CALL davcio (psi_prj, lrwfc, iuwfc, ik_prj, -1)
      CALL init_us_2 (npw_prj, igk_prj, xk(1, ik_prj), vkb_prj)
      !
    ENDIF
    !
    !
    ! compute the kinetic energy used for H\psi and for preconditioning
    !
    FORALL(ig=1:npw_prj) &
      g2kin_prj(ig) = tpiba2 * SUM( (xk(:,ik_prj)+g(:,igk_prj(ig)))**2 )
    !
    PERTURBATION_LOOP : &
    DO ipert = 1, npe
      !
      ! calculates dvscf_q*psi_k in G_space, for all bands
      !
      mode = imode0 + ipert
      CALL dvdpsi( mode, u_prt, xq_prt, dvloc(:, ipert), &
                   npw_wfc, igk_wfc, npw_prj, igk_prj, vkb_wfc, vkb_prj, psi_wfc, dvpsi)
      !
      ! calculates matrix element of the perturbation to the hamiltonian:
      !   < psi_k+q_prj| d^q_prt H | psi_k+q_wfc >
      ! they will be written to file later on (in unit_psidqvpsi, if it is defined)
      !
      COMPUTE_PSI_dH_PSI : &
      IF (lhamiltonian) THEN
        !
        IF (lmetal.and.ldwfc) dpsiaux(:,:) = (0._dp, 0._dp)
        !
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
            !
            psidvpsi = ZDOTC(npw_prj, psi_prj(:, jbnd), 1, dvpsi(:, ibnd), 1)
            CALL mp_sum( psidvpsi, intra_pool_comm )
            psidqvpsi(jbnd, ibnd) = psidvpsi
            !
            IF (lmetal.and.ldwfc) THEN
                deltae = et(ibnd, ik_wfc) - et(jbnd, ik_prj)
                !            theta = 2.0d0*wgauss(deltae/degauss,0)
                IF (ABS(deltae) > 1.0d-5) THEN
                  wg1 = wgauss( (ef-et(ibnd, ik_wfc) ) * degaussm1, ngauss)
                  wg2 = wgauss( (ef-et(jbnd, ik_prj) ) * degaussm1, ngauss)
                  wwg = (wg1 - wg2) / deltae
                ELSE
                  wwg = - w0gauss( (ef-et(ibnd, ik_wfc))*degaussm1, ngauss)*degaussm1
                ENDIF
                !
                psidvpsi = 0.5_dp * wwg * psidvpsi
                CALL ZAXPY(npw_prj,psidvpsi,psi_prj(1,jbnd),1,dpsiaux(1,ibnd),1)
            ENDIF
            !
          ENDDO
        ENDDO
        !
      ENDIF &
      COMPUTE_PSI_dH_PSI
      !
      !
      ! writes psidqvpsi on file
      ! IMPORTANT: these terms are not summed inside pools!
      IF(unit_psidqvpsi>0) THEN
        nrec = imode0 + ipert + (ik - 1) * 3*nat
        CALL davcio (psidqvpsi, lrpdqvp, unit_psidqvpsi, nrec, +1)
!        write(stdout,'(a,i4)') "w psidvpsi:", unit_psidqvpsi
      ENDIF
      !
      ! if we only wanted to compute <psi|dV|psi>  we can skip the rest
      IF(.not. ldwfc) CYCLE PERTURBATION_LOOP
      !
      ! Ortogonalize dvpsi
      !
      CALL start_clock ('ortho')
      wwg = 1._dp
      DO ibnd = 1, nbnd_occ(ik_wfc)
          auxg = CMPLX(0._dp, 0._dp, kind=DP)
          DO jbnd = 1, nbnd_max
              ps (jbnd) = - wwg * ZDOTC(npw_prj, psi_prj(1,jbnd), 1, dvpsi(1,ibnd), 1)
          ENDDO
          !
          CALL mp_sum ( ps, intra_pool_comm )
          !
          DO jbnd = 1, nbnd_max !_occ(ik_prj)
              CALL ZAXPY(npw_prj, ps (jbnd), psi_prj(1, jbnd), 1, auxg, 1)
          ENDDO
          !
          CALL ZCOPY(npw_prj, auxg, 1, spsi, 1)
          CALL DAXPY(2 * npw_prj, 1._dp, spsi, 1, dvpsi(1, ibnd), 1)
      ENDDO
      CALL stop_clock ('ortho')
      !
      CALL DSCAL(2 * npwx * nbnd, -1._dp, dvpsi, 1) !dvpsi = - dvpsi
      !
      ! solution of the linear system (H-eS)*dpsi=dvpsi,
      ! dvpsi=-P_c^+ (dvscf)*psi
      !
      dpsi = (0._dp, 0._dp)
      DO ibnd = 1, nbnd_occ(ik_wfc)
          DO ig = 1, npw_prj
            auxg (ig) = g2kin_prj(ig) * psi_prj(ig, ibnd)
          ENDDO
          eprec1 = PRECONDITIONING_FACTOR &
                  *DBLE(ZDOTC(npw_prj, psi_prj(:, ibnd), 1, auxg, 1))
          !
          CALL mp_sum ( eprec1, intra_pool_comm )
          !
          IF(gstart==2) h_diag (1, ibnd) = 1._dp
          DO ig = gstart, npw_prj
            !  h_diag (ig, ibnd) = MIN(1._dp, eprec1/g2kin_prj(ig))
            h_diag(ig,ibnd)=1.d0/MAX(1.0d0,g2kin_prj(ig)/eprec1)
          ENDDO
      ENDDO
      !
      ! Solve the linear problem! Computes the first derivative of wfcs
      conv_root = .true.
      nproj = nbnd_max !_occ(ik_prj)
      CALL cgsolve_all (d3_ch_psi, d3_cg_psi, et(:, ik_wfc), dvpsi, dpsi, &
           h_diag, npwx, npw_prj, thresh, ik_prj, lter, conv_root, anorm, &
           nbnd_occ(ik_wfc), 1)
      !
      ltaver = ltaver + lter
      lintercall = lintercall + 1
      IF (.not.conv_root) THEN
          WRITE( *, '(5x,"iq_wfc/prt/prj",3i3," kpoint",i4," nbnd_occ",1i4, &
          & " linter: root not converged conv/thresh",2e10.3)') iq_wfc, iq_prt, iq_prj, ik_prj, &
          nbnd_occ(ik_prj), anorm, thresh
          CALL errore('solve_linter_d3q', 'root not converged', 1)
      ENDIF
      !
      ! writes dpsi 
      !
      nrec = (imode0 + ipert - 1) * nksq + ik
      CALL davcio (dpsi, lrdwf, unit_dpsi, nrec, + 1)
      !
      IF (lmetal) THEN
          DO ibnd = 1, nbnd !_occ(ik_wfc)
            wg1 = wgauss( (ef-et(ibnd,ik_wfc))*degaussm1, ngauss)
            CALL dscal(2*npw_prj, wg1, dpsi(1, ibnd), 1)
          ENDDO
          CALL daxpy(2*npw_prj * nbnd, 1._dp, dpsiaux, 1, dpsi, 1)
      ENDIF
      !
      ! This is used to calculate Fermi energy shift at q=0 in metals
      !
      IF (lmetq0) THEN
!           CALL incdrhoscf2 (drhoscf(1, ipert), npw_wfc, igk_wfc, psi_wfc, &
!                             npw_prj, igk_prj, dpsi, kplusq(iq_wfc)%wk(ik), &
!                             kplusq(iq_prj)%ikqs(ik), 1)
          CALL incdrhoscf2 (drhoscf(1, ipert), npw_prj, igk_prj, psi_prj, &
                            npw_prj, igk_prj, dpsi, kplusq(iq_wfc)%wk(ik), &
                            kplusq(iq_prj)%ikqs(ik), 1)
      ENDIF
      !
    ENDDO &
    PERTURBATION_LOOP
    !
  ENDDO &
  KPOINTS_LOOP
  !
  !
  IF (lmetq0) THEN
     DO ipert = 1, npe
        !CALL cinterpolate (drhoscf (1, ipert), drhoscf (1, ipert), 1)
        CALL fft_interpolate (dffts, drhoscf (:, ipert), dfftp,  drhoscf (:, ipert))
     ENDDO
  ENDIF
#ifdef __MPI
  CALL mp_sum( drhoscf, inter_pool_comm )
#endif
  !
  ! In the metal case, compute the shift of fermi energy
  IF (lmetq0) CALL set_efsh (drhoscf, imode0, irr, npe)
  !
  DEALLOCATE (h_diag)
  IF (lmetal.and.ldwfc) DEALLOCATE (dpsiaux)
  DEALLOCATE(psidqvpsi)
  DEALLOCATE(auxg)
  DEALLOCATE(spsi)
  DEALLOCATE(dvloc)
  DEALLOCATE(drhoscf)
  !
  DEALLOCATE(dpsi)
  DEALLOCATE(dvpsi)
  DEALLOCATE(psi_prj, psi_wfc)
  DEALLOCATE(igk_wfc, igk_prj)
  DEALLOCATE(vkb_wfc, vkb_prj)
  !
  CALL stop_clock ('solve_linter')
  !
  aux_avg (1) = DBLE (ltaver)
  aux_avg (2) = DBLE (lintercall)
  CALL mp_sum( aux_avg, inter_pool_comm )
  !
  averlt = 0
  IF(aux_avg(2)/=0._dp)THEN
    averlt = aux_avg(1) / aux_avg(2)
  ENDIF
  tcpu = get_clock ('D3_toten')

  IF (ldwfc) THEN
    !WRITE(fmt_str, '(a,"f",i1,".1",a,"f",i1,".1",a)') &
    !  '(9x,"pert=",i5,5x,"time taken/total="', &
    !  INT(log10(MAX((tcpu-tcpu0+.05_dp),1._dp))+3._dp), &
    !  '"/"', &
    !  INT(log10(MAX((tcpu+.05_dp),1._dp))+3._dp), &
    !  '" secs",5x,"av.it.=",f5.1)'
    WRITE( stdout, '(9x,"pert=",i5,5x,"time taken/total=",1f9.1,"/",1f9.1," secs",5x,"av.it.",1f6.1)') &
            irr, tcpu-tcpu0, tcpu, averlt
  ENDIF
  !
  FLUSH( stdout )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE solve_linter_d3q
!-----------------------------------------------------------------------

#undef PRECONDITIONING_FACTOR

END MODULE linter_d3q
