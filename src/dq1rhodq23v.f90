! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dq1rhodq23v_module
  !
!  USE kinds, ONLY : DP
!  COMPLEX(DP),ALLOCATABLE :: d3mrd(:,:,:)
!  LOGICAL :: first = .true.
  IMPLICIT NONE
  !
  PUBLIC  :: dq1rhodq23v
  PRIVATE :: dpsi_correction, dq23v_nonlocal, dq23v_local
  !
  LOGICAL,PARAMETER :: prof_dq1rhod2v = .false.
  !
  CONTAINS
!----------------------------------------------------------------------
SUBROUTINE dq1rhodq23v(iq_drho, iq_dva, iq_dvb, d3dyn)
  !-----------------------------------------------------------------------
  ! This subroutine reads the variation of charge density relative to q1 than passes the
  ! control to dq23v_(non)local that compute the second derivative of (non) local
  ! potential. Than it rotates the D3 matrix and ordinates it.
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE mp_pools,       ONLY : inter_pool_comm, my_pool_id
  USE mp,             ONLY : mp_sum, mp_bcast
  USE d3_iofiles,     ONLY : read_drho
  USE d3_basis,       ONLY : patq
  USE io_global,      ONLY : stdout
  USE kplus3q,        ONLY : q_names
  USE d3_debug,       ONLY : dbgwrite_d3dyn
  !
  IMPLICIT NONE
  COMPLEX(DP),VOLATILE,INTENT(inout) :: d3dyn( 3*nat, 3*nat, 3*nat) ! the D3 matrix
  INTEGER,INTENT(in) :: iq_drho ! q (index) for the perturbation of rho
  INTEGER,INTENT(in) :: iq_dva, iq_dvb ! index of qs w.r.t. we will derive the potential
  !
  COMPLEX(DP),ALLOCATABLE :: drhoscf (:)     ! delta rho (will be read from file)
  COMPLEX(DP),ALLOCATABLE :: d23_loc(:,:)    ! partial and not ordered D3 matrix, local part
  COMPLEX(DP),ALLOCATABLE :: d23_nlc(:,:)    ! partial and not ordered D3 matrix, non-local part

  COMPLEX(DP),ALLOCATABLE :: d3dyn_loc(:,:,:) ! contrib to D3, local
  COMPLEX(DP),ALLOCATABLE :: d3dyn_nlc(:,:,:) ! contrib to D3, non-local

  COMPLEX(DP) :: work_loc, work_nlc ! workspace
  INTEGER :: na, alpha, beta, na_alpha, na_beta ! atom, cartesian, cartesian+atom
  ! tricks for a correct and automatic ordering of D3
  ! (see dpsi1dv2dpsi3.f90 for a detailed discussion):
  INTEGER,VOLATILE,POINTER :: nu_dva=>null(), nu_dvb=>null(), nu_drho=>null()
  INTEGER,VOLATILE,TARGET  :: nu(3) = (/ 0, 0, 0 /)
  CHARACTER(len=11),PARAMETER :: sub='dq1rhodq23v'
  CHARACTER(len=20) :: filename="dq1rhodq23v.X.loc.d3"
  !
  CALL start_clock('dq1rhodq23v')
  !
  WRITE(stdout, '(7x,"<d_",a2," psi_k| d^2_(",a2,",",a2,") V|psi_k>")') &
      ADJUSTL(q_names(iq_drho)),ADJUSTL(q_names(iq_dva)), ADJUSTL(q_names(iq_dvb))
  !
  nu_drho => nu(ABS(iq_drho)); nu_dva => nu(ABS(iq_dva)); nu_dvb => nu(ABS(iq_dvb))
  ! check for consistency:
  nu = (/ -1,0,1 /)
  IF(nu_drho == nu_dva .or. nu_drho == nu_dvb .or. nu_dva == nu_dvb) &
    CALL errore(sub, "Invalid choice of iq's (repeated)", 2)
  !
  ALLOCATE(d23_loc(3*nat, 3*nat), d23_nlc(3*nat, 3*nat))
  ALLOCATE(d3dyn_loc(3*nat, 3*nat, 3*nat), d3dyn_nlc(3*nat, 3*nat, 3*nat))
  !
  d3dyn_loc = (0._dp, 0._dp)
  d3dyn_nlc = (0._dp, 0._dp)
  !
  RHO_PERTURBATION : &
  DO nu_drho = 1, 3 * nat
    ! prepare partial D3 matrix (only two indexes)
    d23_loc = (0._dp, 0._dp)
    d23_nlc = (0._dp, 0._dp)
#define DO_LOCAL_PART
#define DO_NONLOCAL_PART
    !
    ! FIRST: the local part, uses the G-space drho and is only computed 
    !        by the first pool (no sum of k-points).
#ifdef DO_LOCAL_PART
    IF ( my_pool_id == 0 ) THEN
      ! Read drho_q and takes it Fourier transform
      ALLOCATE(drhoscf(dfftp%nnr))
      CALL read_drho(drhoscf, iq_drho, nu_drho, with_core=.false., pool_only=.true.)
      CALL fwfft('Rho', drhoscf, dfftp)
      ! Compute the local part of \int dq1 rho dq2q3 v
      CALL dq23v_local(iq_drho, drhoscf, d23_loc)
      ! drho is no more necessary: we can free some ram
      DEALLOCATE(drhoscf)
      !
    ENDIF
    !
    ! Sum the partial D3 matrix among different pools
    !CALL mp_bcast( d23_loc, inter_pool_comm )
    !
#endif
    !
    ! SECOND: the non-local part is computed by all processors, it uses
    ! psi and dpsi, but not drho:
    !
#ifdef DO_NONLOCAL_PART
    CALL dq23v_nonlocal(nu_drho, iq_drho, d23_nlc)
    !
    ! Sum the partial D3 matrix among different pools
    CALL mp_sum( d23_nlc, inter_pool_comm )
#endif
    !
    ! Rotate the dynamical matrix on the basis of patterns.
    ! The index corresponding to iq_drho does not need to be rotated, as drho is already
    ! on the base of patterns.
    ! Furthermore, because dV/(dq1 dq2) = dV/(dq2 dq1) we have to copy each term twice with the indexes
    ! corresponding to d2V swapped.
    !
    DO nu_dva = 1, 3*nat
    DO nu_dvb = 1, 3*nat
      work_loc = (0._dp, 0._dp)
      work_nlc = (0._dp, 0._dp)
      DO na = 1, nat
        DO alpha = 1,3
          na_alpha = 3*(na-1) + alpha
          DO beta = 1,3
            na_beta = 3*(na-1) + beta
            !
            work_loc = work_loc &
                 +  patq(iq_dva)%u(na_alpha,nu_dva) &
                  * patq(iq_dvb)%u(na_beta, nu_dvb) &
                  * d23_loc(na_alpha,na_beta)
            work_nlc = work_nlc &
                 +  patq(iq_dva)%u(na_alpha,nu_dva) &
                  * patq(iq_dvb)%u(na_beta, nu_dvb) &
                  * d23_nlc(na_alpha,na_beta)
          ENDDO
        ENDDO
      ENDDO
      !
      d3dyn_loc(nu(1), nu(2), nu(3)) = d3dyn_loc(nu(1), nu(2), nu(3)) &
                                  + work_loc
      d3dyn_nlc(nu(1), nu(2), nu(3)) = d3dyn_nlc(nu(1), nu(2), nu(3)) &
                                  + work_nlc
    ENDDO
    ENDDO
    !
    !CALL print_clock('dq1rhodq23v')
    IF(prof_dq1rhod2v) THEN
      CALL print_clock('dq23v_loc')
      CALL print_clock('dq23v_nlc')
      CALL print_clock('dq23v_nlc:10')
      CALL print_clock('dq23v_nlc:20')
      CALL print_clock('dq23v_nlc:30')
      CALL print_clock('dq23v_nlc:40')
      CALL print_clock('dpsi_corr')
      CALL print_clock('dpsi_corr:10')
      CALL print_clock('dpsi_corr:20')
      CALL print_clock('dpsi_corr:30')
      CALL print_clock('dpsi_corr:40')
      CALL print_clock('dpsi_corr:50')
    ENDIF
    !
  ENDDO &
  RHO_PERTURBATION

  d3dyn = d3dyn + d3dyn_loc + d3dyn_nlc
  ! debug:
  write(filename,'("drd2v.loc.",3i1,".d3")') iq_drho, iq_dva, iq_dvb
  CALL dbgwrite_d3dyn (d3dyn_loc, filename, 1)  !
  write(filename,'("drd2v.nlc.",3i1,".d3")') iq_drho, iq_dva, iq_dvb
  CALL dbgwrite_d3dyn (d3dyn_nlc, filename, 1)  !
  !
  DEALLOCATE(d23_loc, d23_nlc, d3dyn_loc, d3dyn_nlc)
  !
  CALL stop_clock('dq1rhodq23v')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE dq1rhodq23v
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE dq23v_local (iq_drho, drhoscf_G, d3dyn_d23v)
  !-----------------------------------------------------------------------
  ! calculates the term containing the second variation of the potential
  ! and the first variation of the charge density with respect to a
  ! perturbation at a generic q
  !
  USE kinds,        ONLY : DP
  USE constants,    ONLY : tpi
  USE ions_base,    ONLY : nat, ityp, tau
  USE fft_base,     ONLY : dfftp
  USE gvect,        ONLY : ngm, g !, nl
  USE cell_base,    ONLY : tpiba2, omega
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE d3com,        ONLY : d3v
  USE kplus3q,      ONLY : kplusq

  IMPLICIT NONE
  ! I/O variables
  INTEGER,INTENT(in) :: iq_drho
  COMPLEX(DP),INTENT(in) :: drhoscf_G(dfftp%nnr)  ! fft of the variation of the charge density
  COMPLEX(DP),INTENT(inout) :: d3dyn_d23v( 3*nat, 3*nat)
  ! local variables
  COMPLEX(DP),ALLOCATABLE :: struct_factor(:)
  REAL(DP) :: gtau
  COMPLEX (DP), ALLOCATABLE :: d3dyn_wrk (:,:)
  REAL(DP) :: xq(3) ! 
  INTEGER :: na, icart, jcart, na_icart, na_jcart
  INTEGER :: ng
  !
  IF (prof_dq1rhod2v) CALL start_clock('dq23v_loc')
  !
  xq(:) = kplusq(iq_drho)%xq(:)
  !
  ALLOCATE  (d3dyn_wrk( 3 * nat, 3 * nat))
  d3dyn_wrk (:,:) = (0._dp, 0._dp)
  !
  ALLOCATE(struct_factor(ngm))
  !
  VLOC_ATOMS_LOOP : &
  DO na = 1, nat
    !
    ! Pre-calculate factor that depends only on atom:
    DO ng = 1, ngm
      gtau = tpi * SUM((xq+g(:,ng)) * tau(:,na))
      struct_factor(ng) = -0.5_dp * tpiba2 * omega *CMPLX(COS(gtau), SIN(gtau), kind=DP)
    ENDDO
    !
    DO icart = 1, 3
      na_icart = 3 * (na - 1) + icart
      DO jcart = 1, 3
        na_jcart = 3 * (na - 1) + jcart
        !
        DO ng = 1, ngm
          d3dyn_wrk(na_icart, na_jcart) = d3dyn_wrk(na_icart, na_jcart) &
                + (xq(icart) + g(icart, ng)) * (xq(jcart) + g(jcart, ng)) &
                *  d3v(iq_drho)%loc(ng, ityp(na)) * struct_factor(ng) * drhoscf_G(dfftp%nl(ng))
        ENDDO
        !
      ENDDO
    ENDDO
  ENDDO VLOC_ATOMS_LOOP
!  WRITE(10088, '(//)')
!  WRITE(10099, '(/)')
  !
  CALL mp_sum(  d3dyn_wrk, intra_pool_comm )
  d3dyn_d23v = d3dyn_d23v + d3dyn_wrk
  !
  DEALLOCATE(struct_factor)
  DEALLOCATE(d3dyn_wrk)
  IF (prof_dq1rhod2v) CALL stop_clock('dq23v_loc')
  !
  !-----------------------------------------------------------------------
END SUBROUTINE dq23v_local
!----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE dq23v_nonlocal(nu_drho, iq_drho, d3dyn_d23v)
  !-----------------------------------------------------------------------
  ! calculates the term containing the second variation of the potential
  ! and the first variation of the charge density with respect to a
  ! perturbation at a generic q
  !
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE cell_base,       ONLY : tpiba
  USE kinds,           ONLY : DP
  USE constants,       ONLY : tpi
  USE pwcom,           ONLY : lgauss, xk
  USE gvect,           ONLY : ngm, g
  USE wvfct,           ONLY : nbnd, npwx
  USE uspp,            ONLY : dvan, nkb
  USE uspp_param,      ONLY : nh
  USE units_lr,        ONLY : iuwfc, lrwfc, lrdwf
  USE control_lr,      ONLY : nbnd_occ
  USE mp_pools,        ONLY : intra_pool_comm
  USE mp,              ONLY : mp_sum
  USE kplus3q,         ONLY : nksq, kplusq, nbnd_max 
  USE d3_iofiles,      ONLY : iu_dwfc
  USE d3_efermi_shift, ONLY : read_efsh
  USE uspp_init,       ONLY : init_us_2
  !
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout) :: d3dyn_d23v( 3*nat, 3*nat)
  !
  INTEGER,INTENT(in) :: nu_drho ! index of the perturbation associated with drho
  INTEGER,INTENT(in) :: iq_drho ! -iq associated with the variation of rho
  !
  ! local variables
  INTEGER :: nt, na, icart, jcart, na_icart, na_jcart ! type, atom, cartesianx2, cartesian+atom
  INTEGER :: ng, ibnd, ih, jh, ikb, jkb, ijkb0       ! plane-waves, bands, projectors(x5)
  INTEGER :: ik, ikq_gamm, ikq_drho ! k-points
  INTEGER :: nrec, ios              ! auxiliary for i/o 

  COMPLEX (DP) :: ZDOTC
  !
  COMPLEX(DP),ALLOCATABLE :: d3dyn_wrk (:,:), alpha(:,:), beta(:,:)
  COMPLEX(DP),ALLOCATABLE :: psi(:,:), dpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: work1 (:), work2 (:), work3 (:), work4 (:), work5 (:), work6 (:)
  ! work space
  INTEGER :: npw_drho, npw_gamm
  INTEGER,ALLOCATABLE :: igk_drho(:), igk_gamm(:)
  COMPLEX(DP),ALLOCATABLE :: vkb_drho(:,:), vkb_gamm(:,:)
  ! for readability

  IF (prof_dq1rhod2v) CALL start_clock('dq23v_nlc')
  ALLOCATE(d3dyn_wrk(3*nat, 3*nat))
  ALLOCATE(vkb_drho (npwx, nkb), vkb_gamm(npwx, nkb))
  !
  ! FIXME: do we really need 6 copies of the wfcs???
  ALLOCATE(work1(npwx),work2(npwx),work3(npwx),&
           work4(npwx),work5(npwx),work6(npwx))
  !
  ALLOCATE(psi(npwx,nbnd),dpsi(npwx,nbnd))
  ALLOCATE(igk_gamm(npwx))
  ALLOCATE(igk_drho(npwx))

  d3dyn_wrk = (0.d0, 0.d0)
  !
  ! Here the contribution deriving from the local part of the potential;
  ! it is computed only by the first pool (no sum over k needed)
  !
  ! each pool contributes to next term
  !
  ! Here we compute the nonlocal (Kleinman-Bylander) contribution.
  !
!   REWIND (unit = kplusq(iq_drho)%iunigkq)
!   REWIND (unit = kplusq(0)%iunigkq)
  !
  IF(lgauss) CALL read_efsh()
  !
  LOOP_ON_KPOINTS : &
  DO ik = 1, nksq
    !
    IF (prof_dq1rhod2v) CALL start_clock('dq23v_nlc:10')
    ikq_drho = kplusq(iq_drho)%ikqs(ik)
    npw_drho = kplusq(iq_drho)%ngkq(ik)
    igk_drho = kplusq(iq_drho)%igkq(:,ik)
    ikq_gamm = kplusq(0)%ikqs(ik)
    npw_gamm = kplusq(0)%ngkq(ik)
    igk_gamm = kplusq(0)%igkq(:,ik)
    !
    ! Read the unperturbed wavefunctions at k (periodicity: k)
    CALL davcio (psi, lrwfc, iuwfc, ikq_gamm, -1)
    !
    CALL init_us_2(npw_drho, igk_drho, xk(1, ikq_drho), vkb_drho)
    CALL init_us_2(npw_gamm, igk_gamm, xk(1, ikq_gamm), vkb_gamm)
    !
    ! Reads the first variation of the wavefunction projected on conduction
    ! dpsi = d^q psi_k (periodicity: k+q)
    nrec = (nu_drho - 1) * nksq + ik
    CALL davcio (dpsi, lrdwf, iu_dwfc(0,iq_drho), nrec, -1)
    !
    ! In the metallic case corrects dpsi so as that the density matrix
    ! will be:   Sum_{k,nu} 2 * | dpsi > < psi |
    IF (prof_dq1rhod2v) CALL stop_clock('dq23v_nlc:10')

    IF (prof_dq1rhod2v) CALL start_clock('dq23v_nlc:20')
    IF (lgauss) CALL dpsi_correction(dpsi, ik, iq_drho, nu_drho, npw_gamm, npw_drho)
    IF (prof_dq1rhod2v) CALL stop_clock('dq23v_nlc:20')
    !
    DO icart = 1, 3
      DO jcart = 1, 3
          LOOP_ON_BANDS : &
          !DO ibnd = 1, MIN(nbnd_occ(ikq_gamm), nbnd_occ(ikq_drho))
          DO ibnd = 1, nbnd_max
            ! Note: in this loops, xk(*,ikq_drho) is actually k+q
            !
            IF (prof_dq1rhod2v) CALL start_clock('dq23v_nlc:30')
            ! w1 = 2pi/a \psi_k * (k+g)_icart
            FORALL(ng = 1:npw_gamm) work1(ng)= psi(ng,ibnd)*tpiba*(xk(icart,ikq_gamm)&
                                              +g(icart,igk_gamm(ng)))
            ! w2 = 2pi/a \psi_k * (k+g)_jcart
            FORALL(ng = 1:npw_gamm) work2(ng)= psi(ng,ibnd)*tpiba*(xk(jcart,ikq_gamm)&
                                              +g(jcart,igk_gamm(ng)))
            ! w5 = 2pi/a \psi_k * (k+g)_icart * (k+g)_jcart
            FORALL(ng = 1:npw_gamm) work5(ng)=    work1(ng)*tpiba*(xk(jcart,ikq_gamm)&
                                              +g(jcart,igk_gamm(ng)))
            !
            ! w3 = 2pi/a d^q \psi_k * (k+q+g)_icart
            FORALL(ng = 1:npw_drho) work3(ng)= dpsi(ng,ibnd)*tpiba*(xk(icart,ikq_drho)&
                                              +g(icart,igk_drho(ng)))
            ! w4 = 2pi/a d^q \psi_k * (k+q+g)_jcart
            FORALL(ng = 1:npw_drho) work4(ng)= dpsi(ng,ibnd)*tpiba*(xk(jcart,ikq_drho)&
                                              +g(jcart,igk_drho(ng)))
            ! w6 = 2pi/a d^q \psi_k * (k+q+g)_icart * (k+q+g)_jcart
            FORALL(ng = 1:npw_drho) work6(ng)=    work3(ng)*tpiba*(xk(jcart,ikq_drho)&
                                              +g(jcart,igk_drho(ng)))
            IF (prof_dq1rhod2v) CALL stop_clock('dq23v_nlc:30')
            !
            ijkb0 = 0
            !
            IF (prof_dq1rhod2v) CALL start_clock('dq23v_nlc:40')
            LOOP_ON_TYPES : &
            DO nt = 1, ntyp
              LOOP_ON_ATOMS : &
              DO na = 1, nat
              !
              CORRECT_TYPE : &
              IF (ityp (na) == nt) THEN
                  na_icart = 3 * (na - 1) + icart
                  na_jcart = 3 * (na - 1) + jcart
                  !
                  ALLOCATE(alpha(4,nh(nt)), beta(4,nh(nt)))
                  !
                  DO ih = 1, nh(nt)
                     ikb = ijkb0 + ih
                     alpha(1,ih) = ZDOTC(npw_gamm, work1,       1, vkb_gamm(:,ikb), 1)
                     alpha(2,ih) = ZDOTC(npw_gamm, work2,       1, vkb_gamm(:,ikb), 1)
                     alpha(3,ih) = ZDOTC(npw_gamm, work5,       1, vkb_gamm(:,ikb), 1)
                     alpha(4,ih) = ZDOTC(npw_gamm, psi(:,ibnd), 1, vkb_gamm(:,ikb), 1)
                     !
                  ENDDO
                     !
                  DO jh = 1, nh(nt)
                    jkb = ijkb0 + jh
                    beta(1,jh) = ZDOTC(npw_drho, vkb_drho(:,jkb), 1, work4,        1)
                    beta(2,jh) = ZDOTC(npw_drho, vkb_drho(:,jkb), 1, work3,        1)
                    beta(3,jh) = ZDOTC(npw_drho, vkb_drho(:,jkb), 1, dpsi(:,ibnd), 1)
                    beta(4,jh) = ZDOTC(npw_drho, vkb_drho(:,jkb), 1, work6,        1)
                  ENDDO
                  !
                  CALL mp_sum( alpha, intra_pool_comm )
                  CALL mp_sum( beta,  intra_pool_comm )
                  !
                  DO ih = 1, nh(nt)
                  DO jh = 1, nh(nt)
                    d3dyn_wrk (na_icart,na_jcart) = d3dyn_wrk (na_icart,na_jcart) &
                           +( alpha(1,ih) * beta(1,jh) + alpha(2,ih) * beta(2,jh)&
                             -alpha(3,ih) * beta(3,jh) - alpha(4,ih) * beta(4,jh) ) &
                            * dvan(ih, jh, nt) * kplusq(0)%wk(ik)
                  ENDDO
                  ENDDO
                  !
                  ijkb0 = ijkb0 + nh(nt)
                  !
                  DEALLOCATE(alpha,beta)
                  !
              ENDIF CORRECT_TYPE
              !
              ENDDO LOOP_ON_ATOMS
            END DO LOOP_ON_TYPES
            !            
            IF (prof_dq1rhod2v) CALL stop_clock('dq23v_nlc:40')
          END DO LOOP_ON_BANDS
      ENDDO ! jcart
    ENDDO ! icart
  ENDDO &
  LOOP_ON_KPOINTS
  !
  d3dyn_d23v = d3dyn_d23v + d3dyn_wrk
  !
  DEALLOCATE(work1,work2,work3,work4,work5,work6)
  DEALLOCATE(psi,dpsi)
  DEALLOCATE(vkb_drho, vkb_gamm)
  DEALLOCATE(d3dyn_wrk)
  IF (prof_dq1rhod2v) CALL stop_clock('dq23v_nlc')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE dq23v_nonlocal
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE dpsi_correction(dpsi, ik, iq, nu, npw, npwq)
  !-----------------------------------------------------------------------
  ! Used in the metallic case.
  ! If dpsi common variable contains the projection on the conduction
  ! states of the first variation of a wavefunction at a given k-point,
  ! this routine corrects dpsi in such a way that the density matrix
  ! is given by:   Sum_{k,nu} 2 * | dpsi > < psi |
  !
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat
  USE pwcom,           ONLY : npwx, nbnd, degauss, ngauss, et, ef
  USE units_lr,        ONLY : iuwfc, lrwfc
  USE control_lr,      ONLY : nbnd_occ
  USE d3_efermi_shift, ONLY : ef_sh
  USE io_global,       ONLY : stdout
  USE kplus3q,         ONLY : kplusq, nbnd_max
  USE d3_iofiles,      ONLY : iu_psi_dH_psi, lrpdqvp

  IMPLICIT NONE

  COMPLEX(DP),INTENT(inout) :: dpsi(npwx, nbnd)   ! d^q \psi_k to be corrected
  INTEGER,INTENT(in) :: ik, &        ! index of the k and k+q points
                        npw, npwq, & ! numer of plane waves at k and k+q
                        iq, &        ! index of q vector
                        nu           ! mode under consideration

  COMPLEX(DP),ALLOCATABLE :: psi_q(:,:)    ! k+q point wavefunction
  COMPLEX(DP),ALLOCATABLE :: psidvpsi(:,:) ! < psi_{k+q} | V(q) | psi_k >
  INTEGER  :: nrec, ik0, ikq
  INTEGER  :: ibnd, jbnd       ! counters on bands
  REAL(DP) :: wfshift, &      ! the shift coefficent for the wave function
              wgauss, &       ! function computing the theta function
              w0gauss, &      ! function computing the derivative of theta
              deltae, &       ! difference of energy
              wg1, wg2, wwg,& ! weights for metals
              degaussm1       ! 1/degauss
  COMPLEX(DP) :: work         ! workspace


  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr')
  ik0 = kplusq( 0)%ikqs(ik)
  ikq = kplusq(iq)%ikqs(ik)
  degaussm1 = 1._dp / degauss
  !
  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr:10')
  ! we need  the wave function at k+q point
  ALLOCATE(psi_q(npwx, nbnd))
  CALL davcio(psi_q, lrwfc, iuwfc, ikq, -1)
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr:10')

  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr:20')
  ! we also need < psi_{k+q} | V(q) | psi_k >
  ALLOCATE(psidvpsi(nbnd, nbnd))
  nrec = nu + (ik - 1) * 3 * nat
  CALL davcio(psidvpsi, lrpdqvp, iu_psi_dH_psi(0,iq), nrec, -1)
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr:20')
  !
  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr:30')
  ! Multiply dpsi by the theta function
  DO ibnd = 1, nbnd_max !_occ(ik0)
     wg1 = wgauss( (ef - et(ibnd, ik0) )*degaussm1, ngauss)
     CALL dscal(2 * npwq, wg1, dpsi(:, ibnd), 1)
  ENDDO
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr:30')
  !
  ! Adds to dpsi the term containing the valence wavefunctions
  !
  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr:40')
  DO ibnd = 1, nbnd_max !_occ(ik0)
     DO jbnd = 1, nbnd_max !_occ(ikq)
        deltae = et(ibnd, ik0) - et(jbnd, ikq)
        IF (ABS(deltae) > 1.0d-5) THEN
           wg1 = wgauss ( (ef - et(ibnd, ik0) )*degaussm1, ngauss)
           wg2 = wgauss ( (ef - et(jbnd, ikq) )*degaussm1, ngauss)
           wwg = (wg1 - wg2) / deltae
        ELSE
           wwg = - w0gauss( (ef - et(ibnd, ik0))*degaussm1, ngauss)*degaussm1
        ENDIF
        work = 0.5_dp * wwg * psidvpsi(jbnd, ibnd)
        CALL zaxpy (npwq, work, psi_q (1, jbnd), 1, dpsi (1, ibnd),1)
     ENDDO
  ENDDO
  !
  DEALLOCATE(psidvpsi)
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr:40')
  !
  ! If necessary corrects dpsi with a term depending on FermiEnergy shift
  !
  IF (prof_dq1rhod2v) CALL start_clock('dpsi_corr:50')
  IF (kplusq(iq)%lgamma) THEN
    !
    IF(ABS(ef_sh(nu))>1.d-12)THEN
      ! Next line is only used when memory-saving tricks are disabled in set_kplus3q
      IF(ik0/=ikq) THEN ! or, equivalently: .not.kplusq(iq)%lsame(0)
         IF(nu==1 .and. ik==1) WRITE(stdout, "(10x,a)") "-> reading |psi_k>"
        CALL davcio(psi_q, lrwfc, iuwfc, ik0, -1) ! psi_q are actually psi_Gamma from here on
      ENDIF
      !
      IF(nu==1 .and. ik==1) &
        WRITE(stdout, "(10x,a)") "-> including Efermi_shift correction for |d^G psi_k>"
      !
      DO ibnd = 1, nbnd_occ(ikq)
        wfshift = 0.5_dp * ef_sh(nu) * &
                  w0gauss((ef-et(ibnd, ikq))*degaussm1, ngauss)*degaussm1
        CALL daxpy(2*npw, wfshift, psi_q (1, ibnd), 1, dpsi (1, ibnd), 1) 
      ENDDO
    ELSE
      IF(nu==1 .and. ik==1) &
        WRITE(stdout, "(10x,a)") "-> skipping null Efermi_shift correction for |d^G psi_k>"
    ENDIF
  ENDIF
  !
  DEALLOCATE(psi_q)
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr:50')
  !
  IF (prof_dq1rhod2v) CALL stop_clock('dpsi_corr')
  RETURN
  !----------------------------------------------------------------------
END SUBROUTINE dpsi_correction
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
END MODULE dq1rhodq23v_module
!----------------------------------------------------------------------
