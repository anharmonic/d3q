!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
MODULE drho_add_phase_module
!---------------------------------------------------------------------
!
CONTAINS
!---------------------------------------------------------------------
SUBROUTINE drho_add_phase(drho, xq_old, xq_new)
  !---------------------------------------------------------------------
  ! In input drho is the derivative of the charge density
  ! wrt a perturbation having wavevector xq_old.
  ! In output drho is the charge density wrt to a different
  ! wavevector xq_new
  ! The routine should be used only when xq_new - xq_old = Gi
  !
  !
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE cell_base,        ONLY : at
  USE fft_base,         ONLY : dfftp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(inout) :: drho (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)! the derivative of the charge density
  REAL(DP),INTENT(in) :: xq_old(3), xq_new(3)        ! q-points: transform from old to new
!   COMPLEX(DP), ALLOCATABLE :: drho (:) ! auxiliary variable
  INTEGER :: i, j, k, ijk!, npp0           ! counters
  REAL(DP) :: gf(3), n(3)                 !  temp variables
  COMPLEX(DP) ::  term(3), phase          !  temp variables
  REAL(DP),PARAMETER :: gam(3) = (/ 0._dp, 0._dp, 0._dp /)
  LOGICAL,EXTERNAL   :: eqvect
  !
! #if defined (__MPI)
  call start_clock ('drho_add_phase')
!   ALLOCATE( drho(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
  !
  ! Gather all charge from all processors
!   CALL cgather_sym (drho, drho )
  !
  n(1) = tpi / DBLE (dfftp%nr1)
  n(2) = tpi / DBLE (dfftp%nr2)
  n(3) = tpi / DBLE (dfftp%nr3)
  !
  gf(:) = ( xq_old(1) - xq_new(1) ) * at(1,:) + &
          ( xq_old(2) - xq_new(2) ) * at(2,:) + &
          ( xq_old(3) - xq_new(3) ) * at(3,:)
  !
  IF( .not. eqvect(gf, gam, gam,1.d-5) ) &
    CALL errore('drho_add_phase', 'q_old and q_new are not equivalent minus a G vector!', 1)
  !
  gf(:) = gf(:)*n(:)
  !
  term(:) = CMPLX(cos(gf(:)), sin(gf(:)) ,kind=DP)
  !
!   print*, "terms:", term, gf, n, xq_new, xq_old, dfftp%nr1, dfftp%nr2, dfftp%nr3
  !
  phase = (1._dp, 0._dp)
  DO k = 1, dfftp%nr3
    DO j = 1, dfftp%nr2
      DO i = 1, dfftp%nr1
        ijk = i + (j-1)*dfftp%nr1x + (k-1)*dfftp%nr1x*dfftp%nr2x
        drho(ijk) = drho(ijk) * phase
        phase = phase * term (1)
      ENDDO
      phase = phase * term (2)
    ENDDO
    phase = phase * term (3)
  ENDDO
  !
  ! Find the first plane that belong to this processor
!   npp0 = 1
!   DO i = 1, me_pool
!      npp0 = npp0 + dfftp%npp(i) * dfftp%nnp
!   ENDDO
!   ! Each processor takes its own slice of charge
!   CALL zcopy (dfftp%npp(me_pool+1) * dfftp%nnp, drho(npp0), 1, drho, 1)
  !
!   DEALLOCATE (drho)
  !
  CALL stop_clock ('drho_add_phase')
! #else
! #error D3 serial version NOT implemented!
! #endif
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE drho_add_phase
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE drho_add_phase_module
!-----------------------------------------------------------------------
