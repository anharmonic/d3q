!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE d3_exc_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3_exc(d3dyn)
  !-----------------------------------------------------------------------
  !
  !    Calculates the contribution to the derivative of the dynamical
  !    matrix due to the third derivative of the exchange and correlation
  !    energy
  !
  USE ions_base,  ONLY : nat
  USE kinds,      ONLY : DP
  USE scf,        ONLY : rho, rho_core
  USE fft_base,   ONLY : dfftp
  USE cell_base,  ONLY : omega
  USE io_global,  ONLY : ionode_id
  USE mp_global,  ONLY : inter_pool_comm, my_pool_id, &
                         npool, intra_pool_comm
  USE mp,         ONLY : mp_bcast, mp_sum
  USE d3_iofiles, ONLY : read_drho
  !
  IMPLICIT NONE
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
  !
  INTEGER :: ir, ipert, jpert, kpert
  REAL(DP) :: d2mxc, rhotot
  COMPLEX(DP) :: aux
  REAL(DP)     :: pref
  REAL(DP),ALLOCATABLE :: d2muxc (:)
  COMPLEX (DP), ALLOCATABLE :: rho1 (:), rho2 (:), rho3 (:), &
                               d3dyn1 (:,:,:)
  !
  CALL start_clock('d3_exc')
  !
  ALLOCATE (d3dyn1(3*nat, 3*nat, 3*nat))
  !
  IF ( my_pool_id == 0 ) THEN
    !
    ALLOCATE (d2muxc(dfftp%nnr))
    ALLOCATE (rho1(dfftp%nnr))
    ALLOCATE (rho2(dfftp%nnr))
    ALLOCATE (rho3(dfftp%nnr))
    !
    pref = omega / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    !
    ! Calculates third derivative of Exc
    !
    d2muxc = 0._dp
    DO ir = 1, dfftp%nnr
      rhotot = rho%of_r(ir, 1) + rho_core (ir)
      IF (rhotot > 1.e-30_dp)   d2muxc (ir) =  d2mxc( rhotot)
      IF (rhotot < - 1.e-30_dp) d2muxc (ir) = -d2mxc(-rhotot)
    ENDDO
    !
    ! Calculates the contribution to d3dyn
    !
    d3dyn1 = (0._dp, 0._dp)
    DO ipert = 1, 3 * nat
      !
      CALL read_drho(rho1, 1, ipert, with_core=.true., pool_only=.true.)
      DO jpert = 1, 3 * nat
        !
        CALL read_drho(rho2, 2, jpert, with_core=.true., pool_only=.true.)
        DO kpert = 1, 3 * nat
          !
          CALL read_drho(rho3, 3, kpert, with_core=.true., pool_only=.true.)
          !
          ! Short version:
!          aux = SUM(d2muxc*rho1*rho2*rho3)
          !
          ! Exactly the same as above, but should be faster:
          DO ir = 1,dfftp%nnr
            rho3(ir) = rho3(ir)*rho2(ir)
          ENDDO
          DO ir = 1,dfftp%nnr
            rho3(ir) = rho3(ir)*rho1(ir)
          ENDDO
          DO ir = 1,dfftp%nnr
            rho3(ir) = rho3(ir)*d2muxc(ir)
          ENDDO
          aux = (0._dp, 0._dp)
          DO ir = 1,dfftp%nnr
            aux = aux + rho3(ir)
          ENDDO
!           !
          CALL mp_sum( aux, intra_pool_comm )
          !
          d3dyn1(ipert, jpert, kpert) = pref * aux
          !
        ENDDO
      ENDDO
    ENDDO
    !
    DEALLOCATE (d2muxc)
    DEALLOCATE (rho1, rho2, rho3)
    !
  END IF
  !
  IF ( npool > 1 ) CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
  !
  d3dyn = d3dyn  + d3dyn1
  DEALLOCATE (d3dyn1)
  !
  CALL stop_clock('d3_exc')
  !
  RETURN
  !
END SUBROUTINE d3_exc

END MODULE d3_exc_module
