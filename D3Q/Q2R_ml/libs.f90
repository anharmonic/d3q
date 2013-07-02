!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE parameters


  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: DP = kind(0.0d0)

END MODULE parameters
!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !     first the input variables
  !
  use parameters
  implicit none
  real(kind=DP) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
  ! input: first direct lattice vector
  ! input: second direct lattice vector
  ! input: third direct lattice vector
  ! output: first reciprocal lattice vector
  ! output: second reciprocal lattice vector
  ! output: third reciprocal lattice vector
  !
  !   then the local variables
  !
  real(kind=DP) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function rndm ()
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran1 of Num.Rec.
  !
  use parameters
  implicit none
  integer :: irand
  common/random_number/ irand
  real(kind=DP) :: rndm, shuffle (32)
  real(kind=DP), external :: rndx
  integer :: i
  logical :: first
  data first / .true. /
  save first, shuffle, i
  !
  if (first.or.irand.lt.0) then
     irand = - irand
     if (first) irand=1 ! starting seed, must be not be 0
     do i = 32 + 8, 1, - 1
        shuffle (min (i, 32) ) = rndx (irand)
     enddo
     i = 32 * shuffle (1) + 1
     first = .false.
  endif
  rndm = shuffle (i)
  shuffle (i) = rndx (irand)

  i = 32 * rndm + 1
  return

end function rndm

function rndx (irand)
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran0 of Num.Rec.
  !
  use parameters
  implicit none
  integer :: im, ia, iq, ir, irand, is, it
  real(kind=DP) :: rndx, obm
  logical :: first
  data first / .true. /
  save im, ia, iq, ir, obm, first
  !
  if (first) then
     ! this is 2**31-1 avoiding overflow
     im = 2 * (2**30 - 1) + 1
     obm = 1.0 / im
     ia = 7**5
     iq = im / ia
     ir = im - ia * iq
     ! starting seed, must be not be 0
!     irand = 1
     first = .false.
  endif
  is = irand / iq
  it = irand-is * iq
  irand = ia * it - is * ir
  if (irand.lt.0) irand = irand+im
  rndx = irand * obm
  return
end function rndx


subroutine set_rndm_seed(iseed)
!
! this subroutine initialize the random number with the given seed
!
  use parameters
  implicit none
  integer :: irand,iseed
  common/random_number/ irand
  real(kind=DP) :: dummy, rndm
  if (iseed.le.0) call errore('set_rndm_seed', &
                  'seed should be a positive integer',1)
! make sure rndm() has been called already once !
  dummy = rndm()
  irand = - iseed

  return
end subroutine set_rndm_seed

!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
subroutine latgen(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic p (sc)                8  orthorhombic p
  !       2  cubic f (fcc)               9  one face centered orthorhombic
  !       3  cubic i (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal p   11  body centered orthorhombic
  !       5  trigonal r                 12  monoclinic p
  !       6  tetragonal p (st)          13  one face centered monoclinic
  !       7  tetragonal i (bct)         14  triclinic p
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  use parameters, only:DP
  implicit none
  integer ibrav
  real(kind=DP) celldm(6), a1(3), a2(3), a3(3), omega
  !
  real(kind=DP), parameter:: sr2 = 1.414213562373d0, &
                             sr3 = 1.732050807569d0
  integer i,j,k,l,iperm,ir
  real(kind=DP) term, cbya, s, term1, term2, singam, sen
  !
  if(ibrav == 0) go to 100
  !
  !     user-supplied lattice
  !
  do ir=1,3
     a1(ir)=0.d0
     a2(ir)=0.d0
     a3(ir)=0.d0
  end do
  !
  if (ibrav == 1) then
     !
     !     simple cubic lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  else if (ibrav == 2) then
     !
     !     fcc lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  else if (ibrav == 3) then
     !
     !     bcc lattice
     !
     if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     term=celldm(1)/2.d0
     do ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     end do
     a2(1)=-term
     a3(1)=-term
     a3(2)=-term
     !
  else if (ibrav.eq.4) then
     !
     !     hexagonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  else if (ibrav == 5) then
     !
     !     trigonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (4) <= - 0.5d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     term1=sqrt(1.d0+2.d0*celldm(4))
     term2=sqrt(1.d0-celldm(4))
     a2(2)=sr2*celldm(1)*term2/sr3
     a2(3)=celldm(1)*term1/sr3
     a1(1)=celldm(1)*term2/sr2
     a1(2)=-a1(1)/sr3
     a1(3)= a2(3)
     a3(1)=-a1(1)
     a3(2)= a1(2)
     a3(3)= a2(3)
     !
  else if (ibrav == 6) then
     !
     !     tetragonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  else if (ibrav == 7) then
     !
     !     body centered tetragonal lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (3) <= 0.d0) &
          call errore ('latgen', 'wrong celldm', ibrav)
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  else if (ibrav == 8) then
     !
     !     Simple orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or. &
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  else if (ibrav == 9) then
     !
     !
     !     One face centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or. &
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1) = 0.5 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a3(3) = celldm(1) * celldm(3)
     !
  else if (ibrav == 10) then
     !
     !     All face centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.&
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a2(1) = 0.5 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 11) then
     !
     !     Body centered orthorhombic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.&
         celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm', ibrav)
     !
     a1(1) = 0.5 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 12) then
     !
     !     Simple monoclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.   &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  else if (ibrav == 13) then
     !
     !     One face centered monoclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.   &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = celldm(1) * celldm(4)
     a1(3) = celldm(1) * sen
     a2(1) = a1(1)
     a2(3) = - a1(3)
     a3(1) = celldm(1) * celldm(2)
     a3(2) = celldm(1) * celldm(3)
     !
  else if (ibrav == 14) then
     !
     !     Triclinic lattice
     !
     if (celldm (1) <= 0.d0 .or. celldm (2) <= 0.d0 .or.  &
         celldm (3) <= 0.d0 .or. abs (celldm (4) ) > 1.d0 .or. &
         abs (celldm (5) ) > 1.d0 .or. abs (celldm (6) ) > 1.d0) &
         call errore ('latgen', 'wrong celldm', ibrav)
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= sqrt((1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  else
     !
     call errore('latgen',' nonexistent bravais lattice',ibrav)
     !
  end if
  !
  !  calculate unit-cell volume omega
  !
100 omega=0.d0
  s=1.d0
  i=1
  j=2
  k=3
  !
101 do iperm=1,3
     omega=omega+s*a1(i)*a2(j)*a3(k)
     l=i
     i=j
     j=k
     k=l
  end do
!
  i=2
  j=1
  k=3
  s=-s
  if(s < 0.d0) go to 101
  omega=abs(omega)
  return
!
end subroutine latgen
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine volume (alat, a1, a2, a3, omega)
  !---------------------------------------------------------------------
  !
  !     Compute the volume of the unit cell
  !
  use parameters
  implicit none
  !
  !     First the I/O variables
  !
  real(kind=DP) :: alat, a1 (3), a2 (3), a3 (3), omega
  ! input:  lattice parameter (unit length)
  ! input: the first lattice vector
  ! input: the second lattice vector
  ! input: the third lattice vector
  ! input: the volume of the unit cell
  !
  !    Here the local variables required by the routine
  !

  real(kind=DP) :: s
  ! the sign of a permutation
  integer :: i, j, k, l, iperm
  !\
  ! \
  ! /   auxiliary indices
  !/
  ! counter on permutations
  !
  !   Compute the volume
  !
  omega = 0.d0
  s = 1.d0
  i = 1
  j = 2
  k = 3
101 do iperm = 1, 3
     omega = omega + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s

  if (s.lt.0.d0) goto 101

  omega = abs (omega) * alat**3
  return
end subroutine volume
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cryst_to_cart (nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates 
  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !
  use parameters
  implicit none
  !
  integer, intent(in) :: nvec, iflag
  ! nvec:  number of vectors (atomic positions or k-points)
  !        to be transformed from crystal to cartesian and vice versa
  ! iflag: gives the direction of the transformation
  real(kind=DP), intent(in) :: trmat (3, 3)
  ! trmat: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
  !          = bg ,  basis of the reciprocal-space lattice, for k-points
  ! if iflag=-1: the opposite
  real(kind=DP), intent(inout) :: vec (3, nvec)
  ! coordinates of the vector (atomic positions or k-points) to be
  ! transformed - overwritten on output
  !
  !    local variables
  !
  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  real(kind=DP) :: vau (3)
  ! workspace
  !
  !     Compute the cartesian coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
     if (iflag.eq.1) then
        do kpol = 1, 3
           vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
                * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
        enddo
     else
        do kpol = 1, 3
           vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
                * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
        enddo
     endif
     do kpol = 1, 3
        vec (kpol, nv) = vau (kpol)
     enddo
  enddo
  !
  return
end subroutine cryst_to_cart

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine errore (routin, messag, ierr)
  !----------------------------------------------------------------------
  !
  !    This is a simple routine which writes an error message to
  !    output. If ierr = 0 it does nothing, if ierr < 0 it
  !    writes the message but does not stop, if ierr > 0 stops.
  !
  !    **** Important note for parallel execution ***
  !    in parallel execution unit 6 is written only by the first node;
  !    all other nodes have unit 6 redirected to nothing (/dev/null).
  !    As a consequence an error not occurring on the first node
  !    will be invisible. For t3e and origin machines, this problem
  !    is solved by writing an error message to unit * instead of 6.
  !    For ibm sp machines, we write to the standard error, unit 0
  !    (this will appear in the error files produced by loadleveler).
  !
  use parameters
  implicit none
  character (len=*) :: routin, messag
  ! the name of the calling routine
  ! the output message

  integer :: ierr
  ! the error flag
  if (ierr.eq.0) return
  write (6, * ) ' '
  write (6, '(1x,78("%"))')
  write ( * , '(5x,"from ",a," : error #",i10)') routin, ierr
  write ( * , '(5x,a)') messag
  write (6, '(1x,78("%"))')

  if (ierr.gt.0) then
     write ( * , '("     stopping ...")')
     stop 2
  else
     write (6, * ) ' '
     return
  endif
end subroutine errore

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine rgd_blk (nax,nat,dyn,q,tau,epsil,zeu,bg,omega,sign)
  !-----------------------------------------------------------------------
  ! compute the rigid-ion (long-range) term for q 
  !
!!   #include "machine.h"
  implicit none
  real(kind=8), parameter :: e2=2.d0, pi=3.14159265358979d0, fpi=4.d0*pi
  integer ::  nax, nat         ! (maximum) number of atoms 
  complex(kind=8) :: dyn(3,3,nax,nax) ! dynamical matrix
  real(kind=8) &
       q(3),           &! q-vector
       tau(3,nax),     &! atomic positions
       epsil(3,3),     &! dielectric constant tensor
       zeu(3,3,nax),   &! effective charges tensor
       bg(3,3),        &! reciprocal lattice basis vectors
       omega,          &! unit cell volume
       sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term
  !
  ! local variables
  !
  real(kind=8) zag(3),zbg(3),zcg(3),&! eff. charges  times g-vector
       geg              !  <q+G| epsil | q+G>
  integer :: na,nb,nc, i,j
  !
  integer, parameter:: nrx1=16, nrx2=16, nrx3=16, &
       nmegax=(2*nrx1+1)*(2*nrx2+1)*(2*nrx3+1)
  integer im, m1, m2, m3
  complex(kind=8) facg, cmplx
  real(kind=8) gmega(3,nmegax), alph, fac,g1,g2,g3, fnat(3), facgd, arg
  save gmega
  !
  if (abs(sign).ne.1.0) &
       call errore ('rgd_blk',' wrong value for sign ',1)
  !
  fac = sign*e2*fpi/omega
  !
  ! DIAGONAL TERM FIRST (ONLY ONCE, INITIALIZE GMEGA)
  !
!  write (*,*) 'remember to fix Gmega'
!  write (*,*) 'remember to check diagonal-term symmetry'
  alph = 3.0
  im = 0
  do m1 = -nrx1,nrx1
     do m2 = -nrx2,nrx2
        do m3 = -nrx3,nrx3
           im = im + 1
           gmega(1,im) = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3)
           gmega(2,im) = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3)
           gmega(3,im) = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3)
        end do
     end do
  end do
  !
  do na = 1,nat
     do im =1,nmegax
        g1 = gmega(1,im)
        g2 = gmega(2,im)
        g3 = gmega(3,im)
        do i=1,3
           zag(i)=g1*zeu(1,i,na)+g2*zeu(2,i,na)+g3*zeu(3,i,na)
        end do
        !
        geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
               g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
               g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
        !
        if (geg.gt.0.0) then
           do j=1,3
              fnat(j) = 0.0
           end do
           do nc = 1,nat
              arg = 2*pi* (g1* (tau(1,na)-tau(1,nc))+  &
                           g2* (tau(2,na)-tau(2,nc))+  &
                           g3* (tau(3,na)-tau(3,nc)))
              do j=1,3
                 zcg(j) = g1*zeu(1,j,nc) + g2*zeu(2,j,nc) + g3*zeu(3,j,nc)
                 fnat(j) = fnat(j) + zcg(j)*cos(arg)
              end do
           end do
           facgd = fac*exp(-geg/alph/4.0)/geg
           do i = 1,3
              do j = 1,3
                 dyn(i,j,na,na) = dyn(i,j,na,na) - facgd * zag(i) * fnat(j) 
              end do
           end do
        end if
     end do
  end do
  !
  do na = 1,nat
     do nb = 1,nat
        !
        do im =1,nmegax
           !
           g1 = gmega(1,im) + q(1)
           g2 = gmega(2,im) + q(2)
           g3 = gmega(3,im) + q(3)
           !
           geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+   &
                  g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+   &
                  g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
           !
           if (geg.gt.0.0) then
              !
              do i=1,3
                 zag(i)=g1*zeu(1,i,na)+g2*zeu(2,i,na)+g3*zeu(3,i,na)
                 zbg(i)=g1*zeu(1,i,nb)+g2*zeu(2,i,nb)+g3*zeu(3,i,nb)
              end do
              !
              arg = 2*pi* (g1 * (tau(1,na)-tau(1,nb))+              &
                           g2 * (tau(2,na)-tau(2,nb))+              &
                           g3 * (tau(3,na)-tau(3,nb)))
              !
              facg = fac * exp(-geg/alph/4.0)/geg *                 &
                   cmplx(cos(arg),sin(arg))
              do i=1,3
                 do j=1,3 
                    dyn(i,j,na,nb) = dyn(i,j,na,nb) + facg * zag(i) * zbg(j) 
                 end do
              end do
           end if
           !
        end do
     end do
  end do
  return
  !
end subroutine rgd_blk
!
!-----------------------------------------------------------------------
subroutine nonanal(nax,nat,dyn,q,itau_blk,nax_blk,epsil,zeu,omega)
  !-----------------------------------------------------------------------
  !     add the term with macroscopic electric fields
  !
 implicit none
  real(kind=8), parameter :: e2=2.d0, pi=3.14159265358979d0, fpi=4.d0*pi
 integer     nax, nat,       &! (maximum) number of atoms 
      &      nax_blk,        &! maximum number of atoms in the bulk
      &      itau_blk(nat)    ! 
 !
 complex(kind=8) dyn(3,3,nax,nax) ! dynamical matrix
 real(kind=8) q(3),           &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nax_blk),&! effective charges tenso
      &       omega            ! unit cell volume
 !
 ! local variables
 !
 real(kind=8) zag(3),zbg(3),  &! eff. charges  times g-vector
      &       qeq              !  <q| epsil | q>
 integer na,nb,              &! counters on atoms
      &  na_blk,nb_blk,      &! counters on bulk atoms
      &  i,j                  ! counters on cartesian coordinates
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
 if(qeq.eq.0.0) return
 !
 do na = 1,nat
    na_blk = itau_blk(na)
    do nb = 1,nat
       nb_blk = itau_blk(nb)
       !
       do i=1,3
          !
          zag(i) = q(1)*zeu(1,i,na_blk) +  q(2)*zeu(2,i,na_blk) + &
                   q(3)*zeu(3,i,na_blk)
          zbg(i) = q(1)*zeu(1,i,nb_blk) +  q(2)*zeu(2,i,nb_blk) + &
                   q(3)*zeu(3,i,nb_blk)
       end do
       !
       do i = 1,3
          do j = 1,3
             dyn(i,j,na,nb) = dyn(i,j,na,nb)+ fpi*e2*zag(i)*zbg(j)/qeq/omega
          end do
       end do
    end do
 end do
 !
 return
end subroutine nonanal
!
!-----------------------------------------------------------------------
subroutine dyndiag (nax,nat,amass,ityp,dyn,w2,z)
  !-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On output: w2 = energies, z = displacements
  !
 implicit none
 ! input
 integer nax, nat, ityp(*)
 complex(kind=8) dyn(3,3,nax,nax)
 real(kind=8) amass(*)
 ! output
 real(kind=8) w2(3*nat)
 complex(kind=8) z(3*nax,3*nat)
 ! local
 integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
 complex(kind=8), allocatable :: dyn2(:,:)
 !
 !  fill the two-indices dynamical matrix
 !
 nat3 = 3*nat
 allocate(dyn2 (3*nax, nat3))
 !
 do na = 1,nat
    do nb = 1,nat
       do ipol = 1,3
          do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
          end do
       end do
    end do
 end do
 !
 !  impose hermiticity
 !
 do i = 1,nat3
    dyn2(i,i) = dcmplx(dreal(dyn2(i,i)),0.0)
    do j = 1,i - 1
       if (abs(dyn2(i,j)-conjg(dyn2(j,i))).gt.1.0e-6)              &
            &           write (6,9000) i,j,dyn2(i,j),dyn2(j,i)
       dyn2(j,i) = 0.5* (dyn2(i,j)+conjg(dyn2(j,i)))
       dyn2(i,j) = conjg(dyn2(j,i))
    end do
 end do
 !
 !  divide by the square root of masses
 !
 do na = 1,nat
    nta = ityp(na)
    do nb = 1,nat
       ntb = ityp(nb)
       do ipol = 1,3
          do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = &
                  dyn2((na-1)*3+ipol, (nb-1)*3+jpol) / &
                  sqrt(amass(nta)*amass(ntb))
          end do
       end do
    end do
 end do
 !
 !  diagonalisation
 !
 call cdiagh2(nat3,dyn2,3*nax,w2,z)
 !
 deallocate(dyn2)
 !
 !  displacements are eigenvectors divided by sqrt(amass)
 !
! do i = 1,nat3
!    do na = 1,nat
!       nta = ityp(na)
!       do ipol = 1,3
!          z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)/ sqrt(amass(nta))
!       end do
!    end do
! end do
 !
 !
9000 format ('  element ',2i3,' non symmetric ',2f12.6,/,35x,2f15.6)
 return
end subroutine dyndiag
!
!-----------------------------------------------------------------------
subroutine writemodes (nax,nat,q,w2,z,iout)
  !-----------------------------------------------------------------------
  !
  !   write modes on output file in a readable way
  !
 implicit none
 ! input
 integer nax, nat, iout
 real(kind=8) q(3), w2(3*nat)
 complex(kind=8) z(3*nax,3*nat)
 ! local
 integer nat3, na, nta, ipol, i, j
 real(kind=8):: freq(3*nat)
 real(kind=8):: rydthz,rydcm1,cm1thz,znorm
 !
 nat3=3*nat
 !
 !  conversion factors RYD=>THZ, RYD=>1/CM e 1/CM=>THZ
 !
 rydthz = 13.6058*241.796
 rydcm1 = 13.6058*8065.5
 cm1thz = 241.796/8065.5
 !
 !  write frequencies and normalised displacements
 !
 write(iout,'(5x,''diagonalizing the dynamical matrix ...''/)')
 write(iout,'(1x,''q = '',3f12.4)') q
 write(iout,'(1x,74(''*''))')
 do i = 1,nat3
    !
    freq(i)= sqrt(abs(w2(i)))*rydcm1
    if (w2(i).lt.0.0) freq(i) = -freq(i)
    write (iout,9010) i, freq(i)*cm1thz, freq(i)
    znorm = 0.0
    do j=1,nat3
       znorm=znorm+abs(z(j,i))**2
    end do
    znorm = sqrt(znorm)
    do na = 1,nat
       write (iout,9020) (z((na-1)*3+ipol,i)/znorm,ipol=1,3)
    end do
    !
 end do
 write(iout,'(1x,74(''*''))')
 !
 !      if (flvec.ne.' ') then
 !         open (unit=iout,file=flvec,status='unknown',form='unformatted')
 !         write(iout) nat, nat3, (ityp(i),i=1,nat), (q(i),i=1,3)
 !         write(iout) (freq(i),i=1,nat3), ((z(i,j),i=1,nat3),j=1,nat3)
 !         close(iout)
 !      end if
 !
 return
 !
9010 format(5x,'omega(',i2,') =',f15.6,' [THz] =',f15.6,' [cm-1]')
9020 format (1x,'(',3 (f10.6,1x,f10.6,3x),')')
 !
end subroutine writemodes
!
!-----------------------------------------------------------------------
subroutine writemolden(nax,nat,atm,a0,tau,ityp,w2,z,flmol)
  !-----------------------------------------------------------------------
  !
  !   write modes on output file in a molden-friendly way
  !
 implicit none
 ! input
 integer nax, nat, ityp(nat)
 real(kind=8) a0, tau(3,nat), w2(3*nat)
 complex(kind=8) z(3*nax,3*nat)
 character(len=50) flmol
 character(len=3) atm(*)
 ! local
 integer nat3, na, nta, ipol, i, j, iout
 real(kind=8) :: freq(3*nat)
 real(kind=8) :: rydcm1, znorm
 !
 if (flmol.eq.' ') then
    return
 else
    iout=4
    open (unit=iout,file=flmol,status='unknown',form='formatted')
 end if
 nat3=3*nat
 !
 rydcm1 = 13.6058*8065.5
 !
 !  write frequencies and normalised displacements
 !
 write(iout,'(''[Molden Format]'')')
 !
 write(iout,'(''[FREQ]'')')
 do i = 1,nat3
    freq(i)= sqrt(abs(w2(i)))*rydcm1
    if (w2(i).lt.0.0) freq(i) = 0.0
    write (iout,'(f8.2)') freq(i)
 end do
 !
 write(iout,'(''[FR-COORD]'')')
 do na = 1,nat
    write (iout,'(a6,1x,3f15.5)') atm(ityp(na)),  &
         a0*tau(1,na), a0*tau(2,na), a0*tau(3,na)
 end do
 !
 write(iout,'(''[FR-NORM-COORD]'')')
 do i = 1,nat3
    write(iout,'('' vibration'',i6)') i
    znorm = 0.0
    do j=1,nat3
       znorm=znorm+abs(z(j,i))**2
    end do
    znorm = sqrt(znorm)
    do na = 1,nat
       write (iout,'(3f10.5)') (real(z((na-1)*3+ipol,i))/znorm,ipol=1,3)
    end do
 end do
 !
 close(unit=iout)
 !
 return
 !
end subroutine writemolden
!
!-----------------------------------------------------------------------
subroutine cdiagh2 (n,h,ldh,e,v)
  !-----------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
 implicit none
 !
 ! on INPUT
 integer          n,       &! dimension of the matrix to be diagonalized
      &           ldh       ! leading dimension of h, as declared
 ! in the calling pgm unit
 complex(kind=8)  h(ldh,n)  ! matrix to be diagonalized
 !
 ! on OUTPUT
 real(kind=8)     e(n)      ! eigenvalues
 complex(kind=8)  v(ldh,n)  ! eigenvectors (column-wise)
 !
 ! LOCAL variables (LAPACK version)
 !
 integer          lwork,   &! aux. var.
      &           ILAENV,  &! function which gives block size
      &           nb,      &! block size
      &           info      ! flag saying if the exec. of libr. routines was ok
 !
 real(kind=8), allocatable::   rwork(:)
 complex(kind=8), allocatable:: work(:)
 !
 !     check for the block size
 !
 nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
 if (nb.lt.1) nb=max(1,n)
 if (nb.eq.1.or.nb.ge.n) then
    lwork=2*n-1
 else
    lwork = (nb+1)*n
 endif
 !
 ! allocate workspace
 !
 call ZCOPY(n*ldh,h,1,v,1)
 allocate(work (lwork))
 allocate(rwork (3*n-2))
 call ZHEEV('V','U',n,v,ldh,e,work,lwork,rwork,info)
 call errore ('cdiagh2','info =/= 0',abs(info))
 ! deallocate workspace
 deallocate(rwork)
 deallocate(work)
 !
 return
end subroutine cdiagh2
!
! Copyright (C) 2002-2003 PWSCF+CP group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
real(kind=8) function erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  implicit none  
  real(kind=8) :: x, x2, p1 (4), q1 (4)
  real(kind=8), external :: erfc  
  data p1 / 2.42667955230532d2, 2.19792616182942d1, &
       6.99638348861914d0, - 3.56098437018154d-2 /
  data q1 / 2.15058875869861d2, 9.11649054045149d1, &
       1.50827976304078d1, 1.00000000000000d0 /
  !
  if (abs (x) .gt.6.d0) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1 with 16-byte words
     !
     erf = sign (1.d0, x)  
  else  
     if (abs (x) .le.0.47d0) then  
        x2 = x**2  
        erf = x * (p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 ( &
             4) ) ) ) / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 ( &
             4) ) ) )
     else  
        erf = 1.d0 - erfc (x)  
     endif
  endif
  !
  return  
end function erf
!
!---------------------------------------------------------------------
real(kind=8) function erfc (x)  
  !---------------------------------------------------------------------
  !
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  implicit none  
  real(kind=8) :: x, ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
  real(kind=8), external :: erf  
  data p2 / 3.00459261020162d2, 4.51918953711873d2, &
       3.39320816734344d2, 1.52989285046940d2, 4.31622272220567d1, &
       7.21175825088309d0, 5.64195517478974d-1, - 1.36864857382717d-7 /
  data q2 / 3.00459260956983d2, 7.90950925327898d2, &
       9.31354094850610d2, 6.38980264465631d2, 2.77585444743988d2, &
       7.70001529352295d1, 1.27827273196294d1, 1.00000000000000d0 /
  data p3 / - 2.99610707703542d-3, - 4.94730910623251d-2, - &
       2.26956593539687d-1, - 2.78661308609648d-1, - 2.23192459734185d-2 &
       /
  data q3 / 1.06209230528468d-2, 1.91308926107830d-1, &
       1.05167510706793d0, 1.98733201817135d0, 1.00000000000000d0 /

  data pim1 / 0.564189583547756d0 /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax.gt.26.d0) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     erfc = 0.d0  
  elseif (ax.gt.4.d0) then  
     x2 = x**2  
     xm2 = (1.d0 / ax) **2  
     erfc = (1.d0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax.gt.0.47d0) then  
     x2 = x**2  
     erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     erfc = 1.d0 - erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x.lt.0.d0) erfc = 2.d0 - erfc  
  !
  return  
end function erfc
!
!---------------------------------------------------------------------
real(kind=8) function gauss_freq (x)
  !---------------------------------------------------------------------
  !
  !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
  !             - See comments in erf
  !
  real(kind=8) :: x
  real(kind=8), parameter :: c =  0.707106781186548d0
  !        ( c= sqrt(1/2) )
  real(kind=8), external :: erfc
  !
  gauss_freq = 0.5d0 * erfc ( - x * c)
  !
  return
end function gauss_freq

!-----------------------------------------------------------------------
subroutine dyndiag1 (nat,ntyp,amass,ityp,dyn,w2,z)
  !-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On output: w2 = energies, z = displacements
  !
  use constants, only: dp
  implicit none
  ! input
  integer nat, ntyp, ityp(nat)
  complex(DP) dyn(3,3,nat,nat)
  real(DP) amass(ntyp)
  ! output
  real(DP) w2(3*nat)
  complex(DP) z(3*nat,3*nat)
  ! local
  real(DP) diff, dif1, difrel
  integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
  complex(DP), allocatable :: dyn2(:,:)
  !
  !  fill the two-indices dynamical matrix
  !
  nat3 = 3*nat
  allocate(dyn2 (nat3, nat3))
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  !  impose hermiticity
  !
  diff = 0.d0
  difrel=0.d0
  do i = 1,nat3
     dyn2(i,i) = CMPLX( DBLE(dyn2(i,i)),0.d0,kind=DP)
     do j = 1,i - 1
        dif1 = abs(dyn2(i,j)-CONJG(dyn2(j,i)))
        if ( dif1 > diff .and. &
             max ( abs(dyn2(i,j)), abs(dyn2(j,i))) > 1.0d-6) then
           diff = dif1
           difrel=diff / min ( abs(dyn2(i,j)), abs(dyn2(j,i)))
        end if
        dyn2(i,j) = 0.5d0* (dyn2(i,j)+CONJG(dyn2(j,i)))
        dyn2(j,i) = CONJG(dyn2(i,j))
     end do
  end do
  if ( diff > 1.d-6 ) write (6,'(5x,"Max |d(i,j)-d*(j,i)| = ",f9.6,/,5x, &
       & "Max |d(i,j)-d*(j,i)|/|d(i,j)|: ",f8.4,"%")') diff, difrel*100
  !
  !  divide by the square root of masses
  !
  do na = 1,nat
     nta = ityp(na)
     do nb = 1,nat
        ntb = ityp(nb)
        do ipol = 1,3
           do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = &
                  dyn2((na-1)*3+ipol, (nb-1)*3+jpol) / &
                  sqrt(amass(nta)*amass(ntb))
          end do
       end do
    end do
 end do
 !
 !  diagonalisation
 !
 call cdiagh3(nat3,dyn2,nat3,w2,z)
 !
 deallocate(dyn2)
 !
 !  displacements are eigenvectors divided by sqrt(amass)
 !
 do i = 1,nat3
    do na = 1,nat
       nta = ityp(na)
       do ipol = 1,3
          z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)/ sqrt(amass(nta))
       end do
    end do
 end do
 !
 return
end subroutine dyndiag1
!-----------------------------------------------------------------------
subroutine cdiagh3 (n,h,ldh,e,v)
  !-----------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
  use constants, only: dp
  implicit none
  !
  ! on INPUT
  integer          n,       &! dimension of the matrix to be diagonalized
       &           ldh       ! leading dimension of h, as declared
  ! in the calling pgm unit
  complex(DP)  h(ldh,n)  ! matrix to be diagonalized
  !
  ! on OUTPUT
  real(DP)     e(n)      ! eigenvalues
  complex(DP)  v(ldh,n)  ! eigenvectors (column-wise)
  !
  ! LOCAL variables (LAPACK version)
  !
  integer          lwork,   &! aux. var.
       &           ILAENV,  &! function which gives block size
       &           nb,      &! block size
       &           info      ! flag saying if the exec. of libr. routines was ok
  !
  real(DP), allocatable::   rwork(:)
  complex(DP), allocatable:: work(:)
  !
  !     check for the block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  if (nb.lt.1) nb=max(1,n)
  if (nb.eq.1.or.nb.ge.n) then
     lwork=2*n-1
  else
     lwork = (nb+1)*n
  endif
  !
  ! allocate workspace
  !
  call zcopy(n*ldh,h,1,v,1)
  allocate(work (lwork))
  allocate(rwork (3*n-2))
  call ZHEEV('V','U',n,v,ldh,e,work,lwork,rwork,info)
  call errore ('cdiagh2','info =/= 0',abs(info))
  ! deallocate workspace
  deallocate(rwork)
  deallocate(work)
  !
  return
end subroutine cdiagh3
!-----------------------------------------------------------------------
subroutine wsinit(rws,nrwsx,nrws,atw)
!-----------------------------------------------------------------------
!

  implicit none
  integer i, ii, ir, jr, kr, nrws, nrwsx, nx
  real*8 :: eps, rws(0:3,nrwsx), atw(3,3)
  parameter (eps=1.0d-6,nx=2)
  ii = 1
  do ir=-nx,nx
     do jr=-nx,nx
        do kr=-nx,nx
           do i=1,3
              rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
           end do
           rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
                               rws(3,ii)*rws(3,ii)
           rws(0,ii)=0.5d0*rws(0,ii)
           if (rws(0,ii).gt.eps) ii = ii + 1
           if (ii.gt.nrwsx) call errore('wsinit', 'ii.gt.nrwsx',1)
        end do
     end do
  end do
  nrws = ii - 1
  return
end subroutine wsinit
!
!-----------------------------------------------------------------------
function wsweight(r,rws,nrws)
!-----------------------------------------------------------------------
!
  implicit none
  integer ir, nreq, nrws
  real*8 :: r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
  parameter (eps=1.0d-6)
!
  wsweight = 0.d0
  nreq = 1
  do ir =1,nrws
     rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
     ck = rrt-rws(0,ir)
     if ( ck .gt. eps ) return
     if ( abs(ck) .lt. eps ) nreq = nreq + 1
  end do
  wsweight = 1.d0/DBLE(nreq)
  return
end function wsweight
!
!
!------------------------------------------------------------------------
!subroutine addlist(Rout,nout,Rin,iout)
!-----------------------------------------------------------------------
!
!
!  implicit none
!  integer i, ii, ir, nout,iout
!  real*8 :: Rout(nrtot,nrtot) Rin(3,nrtot)
!  !parameter (eps=1.0d-6,nx=2)
!  ii = 1
!  do ir=-nx,nx
!     do jr=-nx,nx
!        do kr=-nx,nx
!           do i=1,3
!              rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
!           end do
!           rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
!                               rws(3,ii)*rws(3,ii)
!           rws(0,ii)=0.5d0*rws(0,ii)
!           if (rws(0,ii).gt.eps) ii = ii + 1
!           if (ii.gt.nrwsx) call errore('wsinit', 'ii.gt.nrwsx',1)
!        end do
!     end do
!  end do
!  nrws = ii - 1
!  return
!end subroutine addlist
!
