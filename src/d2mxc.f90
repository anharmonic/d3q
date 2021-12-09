!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION d2mxc (rho)
  !-----------------------------------------------------------------------
  !
  !  second derivative of the xc potential with respect to the local density
  !  Analytical for Perdew and Zunger parameterization of the C.A. functional
  !  Numerical for any other functional
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : pi
!  USE dft_mod,     ONLY : get_dft_namea
!  USE dft_par_mod, ONLY: iexch,icorr
!  USE dft_mod,       ONLY : get_icorr, get_iexch
!  USE funct,       ONLY : init_lda_xc
!  USE qe_drivers_lda_lsda, ONLY :  xc_lda
  USE qe_drivers_lda_lsda, ONLY :  xc_lda
  USE xc_lib,       ONLY: xclib_get_id
  IMPLICIT NONE
  integer :: iexch, icorr
  real (DP) :: rho, d2mxc
  ! input: the charge density ( positive )
  ! output: the second derivative of the xc potent

  real (DP) :: b1, b2, gc, a, b, c, d, thofpi_3, fpioth_3, &
       thopi_3, tm1, tm2, tm3, tm4, tm5, tm6


  ! _ parameters defining the functionals
  ! /
  !    pi
  !    (3/4/pi)^0.333
  !    (4*pi/3)^0.333
  !      (3/pi)^0.333
  !    35*b1,
  !    76*b1*b1 + 64*b2
  !    35*b1*b1*b1 + 234*b1*b2
  !    140*b2*b1*b1 + 176*b2*b2
  !    175*b1*b2*b2
  !    64*b2*b2*b2

  parameter (b1 = 1.0529d0, b2 = 0.3334d0, gc = - 0.1423d0, a = &
       0.0311d0, b = - 0.0480d0, c = 0.0020d0, d = - 0.0116d0,  &
       fpioth_3 = 1.61199195401647d0, thofpi_3 = 0.620350490899400d0, &
       thopi_3 = 0.98474502184270d0, tm1 = 36.85150d0, tm2 =    &
       105.59107916d0, tm3 = 122.996139546115d0, tm4 =          &
       71.30831794516d0, tm5 = 20.4812455967d0, tm6 = 2.371792877056d0)

  real (DP) :: rs, x, den
  !
  real(DP) :: dr, varho(3), vx(3), vc(3), ex(3), ec(3)
  integer  :: i
  !
  ! Coefficients for finite differences second derivative, order of h**2
!   real(DP),parameter :: m_1over12 = 0.8333333333333333e-1_dp, &
!                         p_4over3  = 0.1333333333333333e+1_dp, &
!                         m_5over2  = -2.5_dp
  !real(DP),parameter ::  coeffs(-2:2) = (/ m_1over12, p_4over3, m_5over2, p_4over3, m_1over12 /)
  real(DP),parameter :: coeffs(3) = (/ 1._dp, -2._dp, 1._dp /)

  iexch = xclib_get_id('LDA','EXCH')
  icorr = xclib_get_id('LDA','CORR')
!  CALL init_lda_xc()
!print*, "mux", get_iexch() , get_icorr()
  if (iexch == 1 .and. icorr == 1) then
      ! First case: analitical derivative of PZ functional
      rs = thofpi_3 * (1._dp / rho) **0.3333333333333333d0
      if (rs > 1._dp) then
        x = sqrt (rs)
        den = 1._dp + x * b1 + b2 * x**2
        d2mxc = - gc * (tm1 * x + tm2 * x**2 + tm3 * x**3 + tm4 * x**4 &
              + tm5 * x**5 + tm6 * x**6) / ( (rho**2) * (den**4) * 216._dp)
      else
        d2mxc = (9._dp * a + (6._dp * c + 8._dp * d) * rs + 8._dp * c * rs &
              * log (rs) ) / (rho**2) / 27._dp

      endif

      rs = rs * fpioth_3
      d2mxc = d2mxc + (2._dp / 9._dp * thopi_3 * rs**5)

      d2mxc = 2._dp * d2mxc
  else
      !     second case: numerical derivatives
      dr = MIN(1.E-5_DP, 1.E-4_DP * ABS(rho))
      !do i = -1, 1
        varho = rho + (/-1._dp, 0._dp, 1._dp/) * dr
        !call xc (varho, ex, ec, vx(i), vc(i))
        !CALL xc( length, 1, 1, rhoaux, ex, ec, vx(:,1), vc(:,1) )

        CALL xc_lda( 3, varho, ex, ec, vx, vc )
      !enddo
      ! recompute dr, taking in account eventual loss of precision when rho >> 1.d-6
      dr = 0.5_dp * ((rho+dr)-(rho-dr))
      !
      d2mxc = SUM(coeffs * (vx+vc)) / (dr**2)
  endif
  return
END FUNCTION d2mxc
