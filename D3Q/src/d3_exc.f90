!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE d3_exc_module
    !
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
    USE mp_pools,   ONLY : inter_pool_comm, my_pool_id, &
                            npool, intra_pool_comm
    USE mp,         ONLY : mp_bcast, mp_sum
    USE d3_iofiles, ONLY : read_drho
    !
    IMPLICIT NONE
    COMPLEX(DP),INTENT(inout) :: d3dyn( 3 * nat, 3 * nat, 3 * nat)
    !
    INTEGER :: ir, ipert, jpert, kpert
    REAL(DP) :: rhotot
    COMPLEX(DP) :: aux
    REAL(DP)     :: domega
    REAL(DP),ALLOCATABLE :: d2muxc (:)
    COMPLEX (DP), ALLOCATABLE :: drho1d2muxc (:), dr1dr2d2muxc (:), drho3 (:), &
                                d3dyn1 (:,:,:)
    REAL(DP),EXTERNAL :: d2mxc
    !
    CALL start_clock('d3_exc')
    !
    ALLOCATE (d3dyn1(3*nat, 3*nat, 3*nat))
    !
    IF ( my_pool_id == 0 ) THEN
        !
        ALLOCATE (d2muxc(dfftp%nnr))
        ALLOCATE (drho1d2muxc(dfftp%nnr))
        ALLOCATE (dr1dr2d2muxc(dfftp%nnr))
        ALLOCATE (drho3(dfftp%nnr))
        !
        domega = omega / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
        !
        ! Calculates third derivative of Exc
        !
        d2muxc = 0._dp
        DO ir = 1, dfftp%nnr
          rhotot = rho%of_r(ir, 1) + rho_core (ir)
          IF (rhotot >   1.e-30_dp) d2muxc (ir) =  d2mxc( rhotot)
          IF (rhotot < - 1.e-30_dp) d2muxc (ir) = -d2mxc(-rhotot)
        ENDDO
        !d2muxc=1.e+10_dp
        !
        ! Calculates the contribution to d3dyn
        !
        d3dyn1 = (0._dp, 0._dp)
        DO ipert = 1, 3 * nat
        ! read drho1
        CALL read_drho(drho1d2muxc, 1, ipert, with_core=.true., pool_only=.true.)
        ! pre-multiply drho1*d2muxc
        FORALL (ir = 1:dfftp%nnr)
            drho1d2muxc(ir) = drho1d2muxc(ir)*d2muxc(ir)
        END FORALL
        !
        DO jpert = 1, 3 * nat
            ! read drho2
            CALL read_drho(dr1dr2d2muxc, 2, jpert, with_core=.true., pool_only=.true.)
            ! pre-multiply drho2*drho1*d2muxc
            FORALL (ir = 1:dfftp%nnr)
                dr1dr2d2muxc(ir) = dr1dr2d2muxc(ir)*drho1d2muxc(ir)
            END FORALL
            !
            DO kpert = 1, 3 * nat
            ! read drho3
            CALL read_drho(drho3, 3, kpert, with_core=.true., pool_only=.true.)
            ! integrate drho1*drho2*drho3*d2muxc
            aux = 0._dp
            DO ir = 1,dfftp%nnr
                aux = aux + drho3(ir)*dr1dr2d2muxc(ir) ! NOTE: dr1dr2d2muxc = rho2*rho1*d2muxc
            ENDDO
            !
            d3dyn1(ipert, jpert, kpert) = domega * aux
            !
            ENDDO
        ENDDO
        ENDDO
        !
        DEALLOCATE (d2muxc)
        DEALLOCATE (drho1d2muxc, dr1dr2d2muxc, drho3)
        !
        ! SUM inside pool, better to do it once here than (3*nat)**3 times above
        CALL mp_sum( d3dyn1,  intra_pool_comm )
        !
    ELSE
      d3dyn1 = 0._dp
    END IF
    !
    ! Let everybody get the correct matrix (this is not really needed as only cpu 1 writes)
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
    !
END MODULE d3_exc_module
