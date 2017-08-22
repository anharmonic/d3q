!
! Copyright(C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3_valence_module
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE d3com, ONLY : eps => eps_delta
  !REAL(DP),PARAMETER :: eps = 1.e-5_dp
  !
  !-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_valence_ijk(iq1, iq2, iq3, d3dyn, order)
  !-----------------------------------------------------------------------
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : nat
  use pwcom,        ONLY : degauss, ngauss, lgauss, nbnd, et, ef
  use qpoint,       ONLY : nksq
  USE io_global,    ONLY : stdout
  USE mp_pools,     ONLY : me_pool, inter_pool_comm
  USE mp,           ONLY : mp_sum
  USE mp_world,     ONLY : world_comm
  USE kplus3q,      ONLY : kplusq, q_sum_rule, nbnd_max
  USE d3_iofiles,   ONLY : iu_psi_dH_psi, lrpdqvp
  USE control_lr,   ONLY : nbnd_occ

  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout)   :: d3dyn( 3*nat, 3*nat, 3*nat)
  INTEGER,INTENT(in)          :: iq1,iq2,iq3
  LOGICAL,OPTIONAL,INTENT(in) :: order
  REAL(DP),PARAMETER :: twosixth = (2._dp / 6._dp)
  REAL(DP),PARAMETER :: onesixth = (2._dp / 6._dp)

  INTEGER :: ik, ik_i, ik_j, ik_k, &
             ibnd, jbnd, kbnd, &
             nrec
  REAL(DP) :: de_ij, de_ik, de_jk, &
              de_ji, de_ki, de_kj
  REAL(DP) :: wrk
  REAL(DP),ALLOCATABLE :: wg_i(:),  wg_j(:),  wg_k(:), &
                          w0g_i(:), w0g_j(:), w0g_k(:),&
                          w1g_i(:), w1g_j(:), w1g_k(:)
  REAL(DP) :: degaussm1
  REAL(DP),EXTERNAL ::  wgauss, w0gauss, w_1gauss
  COMPLEX(DP),ALLOCATABLE :: pdvp_i(:,:), pdvp_j(:,:), pdvp_k(:,:)
  COMPLEX(DP),VOLATILE,ALLOCATABLE :: d3dyn_aux(:,:,:)
  !
  ! perturbation indexes (see later)
  INTEGER,VOLATILE,TARGET  :: nu(3) = (/ 0,0,0 /)
  INTEGER,VOLATILE,POINTER :: nu_i, nu_j, nu_k
  LOGICAL :: order_internal
  !
  IF(.not.lgauss) RETURN
  !
  CALL start_clock('d3_smr_ijk')
  WRITE(stdout,'(5x,a,3i3)') "Computing 3-terms valence contribution ",iq1,iq2,iq3
  !
  IF(iq3 /= -q_sum_rule(iq1, iq2)) &
    CALL errore('d3_valence', "Invalid choice of iq's (sum coditions not respected)", 1)
  !
  ! nu(1) will contain the index of the mode associated with the perturbation
  ! at +/- q1, this may be the lft or rgt wavefunctions or the potential; in any case
  ! the perturbationpattern associated to q1 MUST correspond to the first index of D3.
  ! Accordingly, the pattern associated with q2 goes to the second index, q3 with the third.
  order_internal = .true.
  IF(present(order)) &
    order_internal = order
  IF(order_internal) THEN
    nu_i => nu(ABS(iq1))
    nu_j => nu(ABS(iq2))
    nu_k => nu(ABS(iq3))
  ELSE
    ! in this case, we do not use the correct ordering (useful for debugging)
    nu_i => nu(1)
    nu_j => nu(3)
    nu_k => nu(2)
  ENDIF
  ! check for consistency:
  nu = (/ -1,0,1 /)
  IF(nu_i == nu_j .or. nu_i == nu_k .or. nu_j == nu_k) &
    CALL errore('d3_valence', "Invalid choice of iq's (repeated)", 2)
  !
  ALLOCATE(d3dyn_aux(3*nat, 3*nat, 3*nat))
  d3dyn_aux = (0._dp, 0._dp)
  !
  ! Only one CPU per pool works (no sums over plane waves here)
  IF(me_pool==0)THEN
    !
    ALLOCATE(pdvp_i( nbnd, nbnd))
    ALLOCATE(pdvp_j( nbnd, nbnd))
    ALLOCATE(pdvp_k( nbnd, nbnd))
    !
    ALLOCATE(wg_i(nbnd),  wg_j(nbnd),  wg_k(nbnd))
    ALLOCATE(w0g_i(nbnd), w0g_j(nbnd), w0g_k(nbnd))
    ALLOCATE(w1g_i(nbnd), w1g_j(nbnd), w1g_k(nbnd))

    degaussm1 = 1._dp / degauss
    !
    K_POINTS_LOOP : &
    DO ik = 1, nksq
      !
      ik_i = kplusq(   0)%ikqs(ik)
      ik_j = kplusq( iq1)%ikqs(ik)
      ik_k = kplusq(-iq3)%ikqs(ik)
      !
      DO ibnd = 1, nbnd
        wg_i(ibnd) = wgauss((ef - et(ibnd, ik_i) ) * degaussm1, ngauss)
        wg_j(ibnd) = wgauss((ef - et(ibnd, ik_j) ) * degaussm1, ngauss)
        wg_k(ibnd) = wgauss((ef - et(ibnd, ik_k) ) * degaussm1, ngauss)

        w0g_i(ibnd) = w0gauss((ef - et(ibnd, ik_i) ) * degaussm1, ngauss) * degaussm1
        w0g_j(ibnd) = w0gauss((ef - et(ibnd, ik_j) ) * degaussm1, ngauss) * degaussm1
        w0g_k(ibnd) = w0gauss((ef - et(ibnd, ik_k) ) * degaussm1, ngauss) * degaussm1

        w1g_i(ibnd) = w_1gauss((ef - et(ibnd, ik_i) ) * degaussm1, ngauss) * degaussm1**2
        w1g_j(ibnd) = w_1gauss((ef - et(ibnd, ik_j) ) * degaussm1, ngauss) * degaussm1**2
        w1g_k(ibnd) = w_1gauss((ef - et(ibnd, ik_k) ) * degaussm1, ngauss) * degaussm1**2
      ENDDO

  !    write(stdout,'(a,i4)') "r psidvpsi:", iu_psi_dH_psi(0,iq1)
  !    write(stdout,'(a,i4)') "r psidvpsi:", iu_psi_dH_psi(iq1,iq2)
  !    write(stdout,'(a,i4)') "r psidvpsi:", iu_psi_dH_psi(-iq3,iq3)

      PERTURBATIONS_LOOPS : &
      DO nu_i = 1, 3 * nat
      nrec = nu_i + (ik-1) * 3*nat
      CALL davcio(pdvp_i, lrpdqvp, iu_psi_dH_psi(0,iq1), nrec, - 1)
      !
      DO nu_j = 1, 3 * nat
      nrec = nu_j + (ik-1) * 3*nat
      CALL davcio(pdvp_j, lrpdqvp, iu_psi_dH_psi(iq1,iq2), nrec, - 1)
      !
      DO nu_k = 1, 3 * nat
      nrec = nu_k + (ik-1) * 3*nat
      CALL davcio(pdvp_k, lrpdqvp, iu_psi_dH_psi(-iq3,iq3), nrec, - 1)
        !
        BANDS_LOOPS : &
        DO ibnd = 1,nbnd_max !_occ(ik_i)
        DO jbnd = 1,nbnd_max !_occ(ik_j)
          de_ij = et(ibnd, ik_i) - et(jbnd, ik_j)
          de_ji = -de_ij
          DO kbnd = 1,nbnd_max !_occ(ik_k)
            de_jk = et(jbnd, ik_j) - et(kbnd, ik_k)
            de_kj = -de_jk
            de_ik = et(ibnd, ik_i) - et(kbnd, ik_k)
            de_ki = -de_ik
            !
            wrk = 0._dp
            IF(       ABS(de_ij) < 2*eps &
                .and. ABS(de_jk) < 2*eps &
                .and. ABS(de_ki) < 2*eps ) &
            THEN
              !wrk = 0.5_dp * w1g_i(ibnd) !+w1g_j(jbnd)+w1g_k(kbnd))/3._dp
              wrk = onesixth*(w1g_i(ibnd)+w1g_j(jbnd)+w1g_k(kbnd))
            ELSEIF( ABS(de_ij) < eps) THEN
              wrk = ( (wg_i(ibnd)-wg_k(kbnd))/de_jk + w0g_i(ibnd) )/de_ki
            ELSEIF( ABS(de_jk) < eps) THEN
              wrk = ( (wg_j(jbnd)-wg_i(ibnd))/de_ki + w0g_j(jbnd) )/de_ij
            ELSEIF( ABS(de_ki) < eps) THEN
              wrk = ( (wg_k(kbnd)-wg_j(jbnd))/de_ij + w0g_k(kbnd) )/de_jk
            ELSE ! All three bigger than eps
              wrk = (wg_i(ibnd)*de_kj + wg_j(jbnd)*de_ik + wg_k(kbnd)*de_ji) &
                    / (de_jk * de_ki * de_ij)
            ENDIF
            !
            d3dyn_aux(nu(1),nu(2),nu(3)) = d3dyn_aux(nu(1),nu(2),nu(3)) &
                                  + twosixth *wrk *kplusq(0)%wk(ik) &
                                    *pdvp_i(jbnd, ibnd) &
                                    *pdvp_j(kbnd, jbnd) &
                                    *pdvp_k(ibnd, kbnd)
          ENDDO
        ENDDO
        ENDDO &
        BANDS_LOOPS
        !
      ENDDO
      ENDDO
      ENDDO &
      PERTURBATIONS_LOOPS
      !
    ENDDO &
    K_POINTS_LOOP
    !
    DEALLOCATE(pdvp_i,pdvp_j,pdvp_k)
    DEALLOCATE(wg_i,  wg_j,  wg_k  )
    DEALLOCATE(w0g_i, w0g_j, w0g_k )
    DEALLOCATE(w1g_i, w1g_j, w1g_k )
    !
  ENDIF
  !
  CALL mp_sum( d3dyn_aux, inter_pool_comm)
  d3dyn = d3dyn + d3dyn_aux
  !
  DEALLOCATE(d3dyn_aux)
  !
  CALL stop_clock('d3_smr_ijk')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_valence_ijk
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_valence_ij(iq_ef, iq_p, iq_m, d3dyn) !, order)
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat
  USE pwcom,           ONLY : degauss, ngauss, lgauss, nbnd, et, ef
  USE qpoint,          ONLY : nksq
  USE control_lr,      ONLY : nbnd_occ
  USE mp_pools,        ONLY : inter_pool_comm, me_pool
  USE io_global,       ONLY : stdout
  USE mp,              ONLY : mp_sum
  USE kplus3q,         ONLY : kplusq, q_sum_rule, nbnd_max
  USE d3_iofiles,      ONLY : iu_psi_dH_psi, iudpdvp, lrpdqvp, lrdpdvp
  USE d3_efermi_shift, ONLY : read_efsh, ef_sh
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: iq_ef, iq_p, iq_m
  COMPLEX(DP),VOLATILE,INTENT(inout)   :: d3dyn(3*nat, 3*nat, 3*nat)
  !
  COMPLEX(DP),ALLOCATABLE :: pdvp_pq(:,:),pdvp_mq(:,:), &
                             dpsidvpsi(:,:), d3dyn_aux(:,:,:)
  COMPLEX(DP) :: bsum, wrk
  REAL(DP) :: degaussm1, de
  !
  INTEGER :: nrec, &
             ik, ikG, ikpq, &
             ibnd, jbnd
  REAL(DP),EXTERNAL ::  w0gauss, w_1gauss
  !
  ! tricks for a correct and automatic ordering of D3
  ! (see dpsi1dv2dpsi3.f90 for a detailed discussion):
  INTEGER,POINTER :: nu_ef=>null(), nu_pq=>null(), nu_mq=>null()
  INTEGER,VOLATILE,TARGET  :: nu(3)
  CHARACTER(len=14),PARAMETER :: sub="d3_valence_ij"
  !
  ! This term only applies to metals
  IF (.not. lgauss) RETURN
  !
  ! This term is only non-zero when ef_sh is non-zero at this point:
  IF (.not.kplusq(iq_ef)%lgamma) RETURN
  !
  CALL start_clock('d3_smr_ij')
  !
  degaussm1 = 1._dp / degauss
  ALLOCATE(d3dyn_aux(3*nat, 3*nat, 3*nat))
  !
  !
  IF(    ANY( (/ iq_ef, iq_p, iq_m /) > 3) &
    .or. ANY( (/ iq_ef, iq_p, iq_m /) < 1) ) &
    CALL errore(sub, "Invalid choice of iq's (out of range)", 1)
  !
!   iq_m = -iq_p
!   iq_m = iq_m_
  nu_ef => nu(iq_ef)
  nu_pq => nu(iq_p)
  nu_mq => nu(iq_m)
  !
  nu = (/ -1, 0, 1 /)
  IF(nu_ef == nu_pq .or. nu_ef == nu_mq .or. nu_pq == nu_mq) &
    CALL errore(sub, "Invalid choice of iq's (repeated)", 2)
  !
  WRITE(stdout,'(5x,a,3i2)') "Doing d3_valence_ij: ", iq_ef, iq_p, iq_m
  !
  ! Read the Fermi energy shift (it is read into array ef_sh from module efermi_shift)
  CALL read_efsh()
  d3dyn_aux = (0._dp, 0._dp)
  !
  ! Return if all Fermi shifts are zero (i.e. in high symmetry materials)
  IF(ALL(ABS(ef_sh)<1.d-12)) RETURN
  !
  ! Only one CPU per pool works (no sums over plane waves here)
  IF(me_pool==0) THEN
    !
    ALLOCATE(pdvp_pq( nbnd, nbnd))
    ALLOCATE(pdvp_mq( nbnd, nbnd))
    ALLOCATE(dpsidvpsi( nbnd, nbnd))

    K_POINTS_LOOP : &
    DO ik = 1,nksq
      ikG  = kplusq(iq_ef)%ikqs(ik) ! k point
      ikpq = kplusq(iq_p)%ikqs(ik)   ! k+q point
      !
      PERT_MQ_LOOP : &
      DO nu_pq = 1,3*nat
        nrec = nu_pq + (ik-1) * 3*nat
        CALL davcio(pdvp_pq, lrpdqvp, iu_psi_dH_psi(0,iq_p),nrec, -1)
        !
        PERT_PQ_LOOP : &
        DO nu_mq = 1,3*nat
          nrec = nu_mq + (ik-1) * 3*nat
          CALL davcio(pdvp_mq, lrpdqvp, iu_psi_dH_psi(0,iq_p),nrec, -1)
          !
          nrec = nu_mq + (nu_pq-1)*3*nat + (ik-1)*(3*nat)**2   
          !         ^-- exchage pq and mq to get conjg at 0,q,-q
          CALL davcio(dpsidvpsi, lrdpdvp, iudpdvp(iq_p), nrec, -1)
          !
          PERT_GAMMA_LOOP : &
          DO nu_ef = 1,3*nat
            !
            bsum = (0._dp, 0._dp)
            DO ibnd = 1,nbnd_max !_occ(ikG)
            DO jbnd = 1,nbnd_max !_occ(ikpq)
              de = et(ibnd,ikG) - et(jbnd,ikpq)
              IF(ABS(de)>eps) THEN
                wrk = ( w0gauss((ef-et(ibnd,ikG)) *degaussm1,ngauss) &
                       -w0gauss((ef-et(jbnd,ikpq))*degaussm1,ngauss) )*degaussm1/de
              ELSE
                wrk = - w_1gauss((ef-et(ibnd,ikG))*degaussm1,ngauss)*degaussm1**2
              ENDIF
              bsum = bsum + kplusq(iq_ef)%wk(ik)*wrk*ef_sh(nu_ef) &
  !                          * pdvp_mq(ibnd,jbnd)*pdvp_pq(jbnd,ibnd)
                          * CONJG(pdvp_mq(jbnd,ibnd))*pdvp_pq(jbnd,ibnd)
            ENDDO
            ENDDO
            !
            d3dyn_aux(nu(1),nu(2),nu(3)) = d3dyn_aux(nu(1),nu(2),nu(3))+0.5_dp*bsum
            !
            bsum = (0._dp, 0._dp)
            DO ibnd = 1,nbnd_max !_occ(ikG)
              bsum = bsum + kplusq(iq_ef)%wk(ik)*ef_sh(nu_ef)*dpsidvpsi(ibnd,ibnd) &
                           *w0gauss((ef-et(ibnd,ikG))*degaussm1,ngauss)*degaussm1
            ENDDO
            !
            d3dyn_aux(nu(1),nu(2),nu(3))    = d3dyn_aux(nu(1),nu(2),nu(3)) +CONJG(bsum)
            !
          ENDDO &
          PERT_GAMMA_LOOP
        ENDDO &
        PERT_PQ_LOOP
      ENDDO &
      PERT_MQ_LOOP
    ENDDO &
    K_POINTS_LOOP
    !
    DEALLOCATE(pdvp_pq, pdvp_mq, dpsidvpsi)
    !
  ENDIF
  !
  CALL mp_sum(d3dyn_aux, inter_pool_comm)
  !
  d3dyn = d3dyn + d3dyn_aux
  DEALLOCATE(d3dyn_aux)
  !
  CALL stop_clock('d3_smr_ij')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_valence_ij
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_valence_gamma(d3dyn)
  !-----------------------------------------------------------------------
  ! The last two terms are only non-zero if q1=q2=q3=Gamma
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp,               ONLY : mp_sum
  USE kplus3q,          ONLY : kplusq, nbnd_max
  USE d3_efermi_shift,  ONLY : read_efsh, ef_sh
  USE io_global,        ONLY : stdout
  USE d3_iofiles,       ONLY : iu_psi_dH_psi, lrpdqvp
  USE pwcom,            ONLY : degauss, ngauss, lgauss, nbnd, et, ef
  USE qpoint,           ONLY : nksq
  USE d3_debug,         ONLY : dbgwrite_d3dyn
  USE control_lr,       ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(inout)   :: d3dyn(3*nat, 3*nat, 3*nat)
  !
  INTEGER :: nu_i, nu_j, nu_k, ibnd, ik, ikG, nrec
  REAL(DP) :: d_dos, degaussm1
  COMPLEX(DP),ALLOCATABLE :: aux(:), d3dyn_2efsh(:,:,:), d3dyn_3efsh(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: pdvp_i(:,:)
  !
  REAL(DP),EXTERNAL :: w_1gauss
  !
  IF (.not.kplusq(1)%lgamma .or. &
      .not.kplusq(2)%lgamma .or. &
      .not.kplusq(3)%lgamma ) &
    RETURN
  !
  IF (.not. lgauss) RETURN
  !
  CALL start_clock('d3_smr_gamma')
  !
  WRITE(stdout,'(7x,a)') "Computing Gamma-only terms"
  !
  ! read Fermi energy shift and store it in efermi_shift->ef_sh
  CALL read_efsh()
  !
  degaussm1 = 1._dp / degauss
  d_dos = 0.d0
  ALLOCATE(pdvp_i(nbnd, nbnd))
  ALLOCATE(aux(3*nat))
  aux(:) = (0.d0, 0.d0)
  !
  DO ik = 1, nksq
    ikG = kplusq(0)%ikqs(ik)
    !
    DO ibnd = 1, nbnd_max !_occ(ikG)
        d_dos = d_dos + kplusq(0)%wk(ik) * w_1gauss((ef-et(ibnd, ikG))*degaussm1, ngauss)*degaussm1**2
    ENDDO
    DO nu_i = 1, 3 * nat
        ! Note that the iu_psi_dH_psi(0,0) does not exist if the calculation does NOT include Gamma!!
        nrec = nu_i + (ik - 1) * 3 * nat
        CALL davcio (pdvp_i, lrpdqvp, iu_psi_dH_psi(0,1), nrec, -1)
        !
        DO ibnd = 1, nbnd_max !_occ(ikG)
          aux (nu_i) = aux(nu_i) + pdvp_i(ibnd, ibnd)*kplusq(0)%wk(ik) &
                                  *w_1gauss((ef-et(ibnd,ikG))*degaussm1, ngauss)*degaussm1**2
        ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE(pdvp_i)
  ALLOCATE(d3dyn_2efsh(3*nat,3*nat,3*nat))
  ALLOCATE(d3dyn_3efsh(3*nat,3*nat,3*nat))
  !
  d3dyn_2efsh = (0._dp, 0._dp)
  d3dyn_3efsh = (0._dp, 0._dp)
  !
  DO nu_i = 1, 3 * nat
  DO nu_j = 1, 3 * nat
  DO nu_k = 1, 3 * nat
      d3dyn_2efsh(nu_i, nu_j, nu_k) = d3dyn_2efsh(nu_i, nu_j, nu_k) + &
          ef_sh(nu_i) * ef_sh(nu_j) * aux(nu_k) + &
          ef_sh(nu_j) * ef_sh(nu_k) * aux(nu_i) + &
          ef_sh(nu_k) * ef_sh(nu_i) * aux(nu_j)
      d3dyn_3efsh(nu_i, nu_j, nu_k) = d3dyn_3efsh(nu_i, nu_j, nu_k) - &
          ef_sh(nu_i) * ef_sh(nu_j) * ef_sh(nu_k) * d_dos
  ENDDO
  ENDDO
  ENDDO
  !
  DEALLOCATE(aux)
  !
  CALL mp_sum( d3dyn_2efsh, inter_pool_comm )
  CALL mp_sum( d3dyn_3efsh, inter_pool_comm )
  !
  CALL dbgwrite_d3dyn (d3dyn_2efsh, 'd3_valence.3', 1)
  CALL dbgwrite_d3dyn (d3dyn_3efsh, 'd3_valence.4', 1)
  !
  d3dyn = d3dyn + d3dyn_2efsh + d3dyn_3efsh
  !
  DEALLOCATE(d3dyn_2efsh, d3dyn_3efsh)
  !
  CALL stop_clock('d3_smr_gamma')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_valence_gamma
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_valence_module
!-----------------------------------------------------------------------


