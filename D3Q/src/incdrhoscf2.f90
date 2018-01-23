!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE incdrhoscf2_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE incdrhoscf2(drhoscf, npw, igk, psi, npwd, igkd, dpsi, weight, ikk, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
!   USE ions_base,  ONLY : nat
  USE kinds,          ONLY : DP
  USE pwcom,          ONLY : npwx, nbnd
  USE cell_base,      ONLY : omega
  USE control_lr,     ONLY : nbnd_occ
  !USE gvecs,          ONLY : nls
  USE fft_base,       ONLY : dffts
  USE fft_interfaces, ONLY : invfft
  USE uspp,           ONLY : okvan
  USE kplus3q,        ONLY : nbnd_max

  implicit none
  integer,intent(in) :: ikk
  ! input: the k point
  integer,intent(in) :: npw, npwd
  integer,intent(in) :: igkd(npwx), igk(npwx)
  ! orderign of plane waves
  real(DP),intent(in) :: weight
  complex (DP),intent(in) :: psi(npwx, nbnd), dpsi(npwx, nbnd)
  ! input: the weight of the k point
  complex (DP),intent(inout) :: drhoscf(dffts%nnr)
  ! output: the change of the charge densit
!   complex (DP),intent(in) :: dbecsum (nhm * (nhm + 1) / 2, nat)
  ! inp/out: the accumulated dbec
  integer,intent(in) :: mode !, flag
  ! flag =1 if dpsi is used (in solve_linte
  ! flag!=1 if dpsi is not used (in addusdd
  !
  !   here the local variable
  !

  real (DP) :: wgt
  ! the effective weight of the k point

  complex (DP), allocatable :: psic (:), dpsic (:)!, dpsif(:,:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  ALLOCATE(dpsic(dffts%nnr))
  ALLOCATE(psic(dffts%nnr))
!   ALLOCATE(dpsif(npwx, nbnd))

  wgt = 2.d0 * weight / omega
!  print*, "incdrhoscf", ik, weight, wgt
  !   if (lgamma) then
!      ikk = ikks(ik) !ik
!   else
!     ikk = kplusq(0)%ikqs(ik) !ikks(ik) !2 * ik - 1
!   endif
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_max !_occ(ikk)
     psic (:) = (0.d0, 0.d0)
     do ig = 1, npw
        psic (dffts%nl (igk (ig) ) ) = psi(ig, ibnd)
     enddo
     call invfft('Wave', psic, dffts)
     dpsic(:) =(0.d0, 0.d0)
     !
     !    here we add the term in the valence due to the change of the
     !    constraint. dpsif is used as work space, dpsi is unchanged
     !
!      if (flag == 1) then
!         dpsif(:, ibnd) = dpsi(:, ibnd)
!      else
!         dpsif(:, ibnd) = (0.d0, 0.d0)
!      endif
     !         call zgemm('N','N', npwd, nbnd, nbnd, (1.d0,0.d0),
     !     +              evq, npwx, prodval(1,1,mode),nbnd,
     !     +             (1.d0,0.d0),dpsif,npwx)
     if (okvan) then
        call errore ('incdrhoscf2', 'US not allowed', 1)
        !            do jbnd=1,nbnd
        !               call zaxpy(npwd,prodval(jbnd,ibnd,mode),
        !     +           evq(1,jbnd),1,dpsif(1,ibnd),1)
        !            enddo
     endif

     do ig = 1, npwd
        dpsic(dffts%nl(igkd(ig))) = dpsi(ig, ibnd)
     enddo
     call invfft('Wave', dpsic, dffts)

     do ir = 1, dffts%nnr
        drhoscf(ir) = drhoscf(ir) + wgt * CONJG(psic(ir)) * dpsic(ir)
        !            if (ir.lt.20) WRITE( stdout,*)   drhoscf(ir)
     enddo

  enddo

  ! USLTRASOFT ONLY:
!   call addusdbec (ik, wgt, dpsif, dbecsum)
  !      WRITE( stdout,*) '*********************'
  !      do ig=1,20
  !         WRITE( stdout,*) dbecsum(ig,1)
  !      enddo
  !      call stoallocate  (ph(.true.))    
  deallocate (psic)
  deallocate (dpsic)
!   deallocate (dpsif)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf2

end module incdrhoscf2_module
