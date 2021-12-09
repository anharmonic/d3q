!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE d3_h_psi
! **********************************************************************
! *                   ~~~ IMPORTANT NOTE ~~~                           *
! * These subroutines are never CALLed explicitly in the D3 code!!!    *
! * They're instead USEd in solve_linter_d3q where they are passed     *
! * to the phonon subroutine PH/cgsolve_all.f90                        *
! * and used. Be extremely careful when changing global variables to   *
! * local ones as it may break something!!!                            *
! * In particular igkq=>igk_prj should NOT be touched lightheartedly.  *
! **********************************************************************

  USE kinds, ONLY : DP
  !
  !USE wvfct, ONLY : igk_prj => igk
  USE qpoint, ONLY : igkq
  USE wvfct,  ONLY : g2kin_prj => g2kin ! used in solve_linter_d3q
  USE uspp,   ONLY : vkb_prj => vkb
  !
  !PRIVATE
                         ! They compute:
  PUBLIC :: d3_ch_psi, & ! ( H - \epsilon S + alpha_pv P_v)
            d3_cg_psi    ! preconditioning

! note : h_psiq from phonon is used instead of the internal version
!   PUBLIC :: d3_h_psiq, & !  (H \psi) and (S \psi)
  !
  INTEGER :: nproj
  INTEGER,ALLOCATABLE,PUBLIC        :: igk_wfc(:)
  INTEGER,ALLOCATABLE,TARGET,PUBLIC :: igk_prj(:)
  COMPLEX(DP),ALLOCATABLE,PUBLIC :: vkb_wfc(:,:)!, vkb_prj(:,:)
  COMPLEX(DP),ALLOCATABLE,PUBLIC :: psi_prj(:,:), psi_wfc(:,:)
  !PUBLIC :: igk_prj, g2kin_prj, vkb_prj, nproj

CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_ch_psi(n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : npwx!, nbnd
! USE uspp,       ONLY : vkb
  USE becmod,     ONLY : becp, calbec
  USE control_lr, ONLY : alpha_pv
  USE mp_pools,   ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
!   USE linter_d3q, ONLY : igk_prj, wfc_prj, vkb_prj
  !
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  REAL(DP),INTENT(in) :: e (m)
  ! input: the eigenvalue

  COMPLEX(DP),INTENT(in)  :: h (npwx, m)
  ! input: the vector
  COMPLEX(DP),INTENT(out) :: ah (npwx, m)
  ! output: the operator applied to the vector
  !
  !   local variables
  !

  INTEGER :: ibnd,ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  COMPLEX (DP), ALLOCATABLE :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h


  CALL start_clock ('ch_psi')
  ALLOCATE (ps( nproj, m))
  ALLOCATE (hpsi( npwx, m))
  ALLOCATE (spsi( npwx, m))

  hpsi = (0.d0, 0.d0)
  spsi = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  igkq => igk_prj ! <-- igkq is hard-coded inside h_psiq
  CALL h_psiq (npwx, n, m, h, hpsi, spsi) !, igk_prj)
  NULLIFY(igkq)

  !CALL start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     ENDDO
  ENDDO
  !
  !   Here we compute the projector in the valence band
  !
  hpsi = (0.d0, 0.d0)
  !
  ps = (0.d0, 0.d0)

  CALL zgemm ('C', 'N', nproj, m, n, (1.d0, 0.d0) , psi_prj, npwx, spsi, &
       npwx, (0.d0, 0.d0) , ps, nproj)
  ps = ps * alpha_pv
#ifdef __MPI
  CALL mp_sum(  ps, intra_pool_comm )
#endif

  CALL zgemm ('N', 'N', n, m, nproj, (1.d0, 0.d0) , psi_prj, npwx, ps, &
       nproj, (1.d0, 0.d0) , hpsi, npwx)
  spsi = hpsi
  !
  !    And apply S again
  !
  CALL calbec (n, vkb_prj, hpsi, becp, m)
  CALL s_psi (npwx, n, m, hpsi, spsi)
  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     ENDDO
  ENDDO

  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)
  !CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
END SUBROUTINE d3_ch_psi

!-----------------------------------------------------------------
SUBROUTINE d3_cg_psi (lda, n, m, psi, h_diag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol !, noncolin
  IMPLICIT NONE

  INTEGER,INTENT(in) :: lda, n, m
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors

  COMPLEX(DP) :: psi (lda*npol, m)
  ! inp/out: the vector to be preconditioned

  REAL(DP) :: h_diag (lda*npol, m)
  ! input: the preconditioning vector

  INTEGER :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  DO k = 1, m
     DO i = 1, n
        psi (i, k) = psi (i, k) * h_diag (i, k)
     ENDDO
  ENDDO
!   IF (noncolin) THEN
!      DO k = 1, m
!         DO i = 1, n
!            psi (i+lda, k) = psi (i+lda, k) * h_diag (i+lda, k)
!         ENDDO
!      ENDDO
!   END IF
  RETURN
END SUBROUTINE d3_cg_psi

#ifdef DO_NOT_COMPILE_THIS_USE_H_PSIQ_FROM_PHONON
!-----------------------------------------------------------------------
SUBROUTINE d3_h_psiq (lda, n, m, psi, hpsi, spsi, igk_)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !

  USE kinds,                 ONLY : DP
  USE wavefunctions,         ONLY : psic, psic_nc
  USE becmod,                ONLY : bec_type, becp, calbec
  USE noncollin_module,      ONLY : npol,domag
  USE lsda_mod,              ONLY : current_spin
  USE gvecs,                 ONLY : nls
  USE fft_base,              ONLY : dffts
  !USE spin_orb,              ONLY : domag
  USE scf,                   ONLY : vrs
  USE uspp,                  ONLY : vkb
  USE wvfct,                 ONLY : npwx
  IMPLICIT NONE
  !
  !     Here the local variables
  !
  INTEGER :: ibnd
  ! counter on bands
  INTEGER,INTENT(in) :: igk_(npwx)
  ! input: plane-waves map
  INTEGER,INTENT(in) :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  INTEGER :: j
  ! DO loop index

  COMPLEX(DP),INTENT(in)  :: psi (lda*npol, m)
  COMPLEX(DP),INTENT(out) :: hpsi (lda*npol, m), spsi (lda*npol, m)
  COMPLEX(DP) :: sup, sdwn
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)

  CALL start_clock ('h_psiq')
  CALL start_clock ('init')

  CALL calbec ( n, vkb, psi, becp, m)
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  hpsi=(0.d0,0.d0)
  DO ibnd = 1, m
     DO j = 1, n
        hpsi (j, ibnd) = g2kin_prj(j) * psi(j, ibnd)
     ENDDO
  ENDDO
!   IF (noncolin) THEN
!      DO ibnd = 1, m
!         DO j = 1, n
!            hpsi (j+lda, ibnd) = g2kin (j) * psi (j+lda, ibnd)
!         ENDDO
!      ENDDO
!   ENDIF
  CALL stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  DO ibnd = 1, m
     CALL start_clock ('firstfft')
!      IF (noncolin) THEN
!         psic_nc = (0.d0, 0.d0)
!         DO j = 1, n
!            psic_nc(nls(igk_(j)),1) = psi (j, ibnd)
!            psic_nc(nls(igk_(j)),2) = psi (j+lda, ibnd)
!         ENDDO
!         CALL cft3s (psic_nc(1,1), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!         CALL cft3s (psic_nc(1,2), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!      ELSE
        psic(:) = (0.d0, 0.d0)
        DO j = 1, n
           psic (nls(igk_(j))) = psi (j, ibnd)
        ENDDO
        CALL invfft('Wave', psic, dffts)
        !CALL cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
!      END IF
!      CALL stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
!      if (noncolin) then
!         if (domag) then
!            DO j=1, nrxxs
!               sup = psic_nc(j,1) * (vrs(j,1)+vrs(j,4)) + &
!                     psic_nc(j,2) * (vrs(j,2)-(0.d0,1.d0)*vrs(j,3))
!               sdwn = psic_nc(j,2) * (vrs(j,1)-vrs(j,4)) + &
!                     psic_nc(j,1) * (vrs(j,2)+(0.d0,1.d0)*vrs(j,3))
!               psic_nc(j,1)=sup
!               psic_nc(j,2)=sdwn
!            end do
!         else
!            DO j=1, nrxxs
!               psic_nc(j,1)=psic_nc(j,1) * vrs(j,1)
!               psic_nc(j,2)=psic_nc(j,2) * vrs(j,1)
!            ENDDO
!         endif
!      else
        DO j = 1, nrxxs
           psic (j) = psic (j) * vrs (j, current_spin)
        ENDDO
!      endif
     !
     !   back to reciprocal space
     !
     CALL start_clock ('secondfft')
!      IF (noncolin) THEN
!         CALL cft3s(psic_nc(1,1),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
!         CALL cft3s(psic_nc(1,2),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
!      !
!      !   addition to the total product
!      !
!         DO j = 1, n
!            hpsi (j, ibnd) = hpsi (j, ibnd) + psic_nc (nls(igk_(j)), 1)
!            hpsi (j+lda, ibnd) = hpsi (j+lda, ibnd) + psic_nc (nls(igk_(j)), 2)
!         ENDDO
!      ELSE
        CALL fwfft('Wave', psic, dffts)
     !
     !   addition to the total product
     !
        DO j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igk_(j)))
        ENDDO
!      END IF
     CALL stop_clock ('secondfft')
  ENDDO
  !
  !  Here the product with the non local potential V_NL psi
  !

  CALL add_vuspsi (lda, n, m, hpsi)

  CALL s_psi (lda, n, m, psi, spsi)

  CALL stop_clock ('h_psiq')
  RETURN
END SUBROUTINe d3_h_psiq
#endif

END MODULE d3_h_psi
