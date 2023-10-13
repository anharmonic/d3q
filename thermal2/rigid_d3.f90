!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE rigid_d3
  PUBLIC :: rgd_blk_d3, nonanal_d3
  PRIVATE
  CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE rgd_blk_d3(nr1, nr2, nr3, nat, dyn, q, tau, epsil, zeu, bg, omega, alat, loto_2d, sign) !, alpha)
  !-----------------------------------------------------------------------
  !! Compute the rigid-ion (long-range) term for q.  
  !! The long-range term used here, to be added to or subtracted from the
  !! dynamical matrices, is exactly the same of the formula introduced in:  
  !! X. Gonze et al, PRB 50. 13035 (1994).  
  !! Only the G-space term is implemented: the Ewald parameter alpha must
  !! be large enough to have negligible r-space contribution.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi, tpi, fpi, e2
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: loto_2d
  !! 2D LOTO correction
  INTEGER, INTENT(in) :: nr1, nr2, nr3
  !! FFT grid
  INTEGER, INTENT(in) :: nat
  !! Number of atoms  
  REAL(KIND = DP), INTENT(in) :: q(3)
  !! q-vector 
  REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
  !! dielectric constant tensor
  REAL(KIND = DP), INTENT(in) :: zeu(3, 3, nat)
  !! effective charges tensor
  REAL(KIND = DP), INTENT(in) :: sign
  !! sign =+/-1.0 ==> add/subtract rigid-ion term
!  REAL(DP), OPTIONAL :: alpha

  REAL(KIND = DP), INTENT(in) :: tau(3, nat)
  !! Atomic positions
  REAL(KIND = DP), INTENT(in) :: bg(3, 3)
  !! Reciprocal lattice basis vectors
  REAL(KIND = DP), INTENT(in) :: omega
  !! Unit cell volume 
  REAL(KIND = DP), INTENT(in) :: alat
  !! Cell dimension units
  COMPLEX(KIND = DP), INTENT(inout) :: dyn(3, 3, nat, nat)
  !! Dynamical matrix
  !
  ! Local variables
  INTEGER :: nr1x, nr2x, nr3x
  !! Max nr in direction 1, 2, 3
  INTEGER :: na
  !! Atom index 1 
  INTEGER :: nb
  !! Atom index 2
  INTEGER :: i
  !! Cartesian direction 1
  INTEGER :: j
  !! Cartesian direction 1
  INTEGER :: m1, m2, m3
  !! Loop over q-points
  REAL(KIND = DP):: geg
  !! <q+G| epsil | q+G>
  REAL(KIND = DP) :: alph
  !! Ewald parameter
  REAL(KIND = DP) :: fac
  !! Prefactor
  REAL(KIND = DP) :: g1, g2, g3
  !! G-vectors
  REAL(KIND = DP) :: facgd
  !! fac * EXP(-geg / (alph * 4.0d0)) / geg  
  REAL(KIND = DP) :: arg
  !! Argument of the function
  REAL(KIND = DP) :: gmax
  !! Maximum G
  REAL(KIND = DP) :: zag(3)
  !! Z * G
  REAL(KIND = DP) :: zbg(3)
  !! Z * G
  REAL(KIND = DP) :: zcg(3)
  !! Z * G
  REAL(KIND = DP) :: fnat(3)
  !! Z * G * cos(arg)
  REAL(KIND = DP) :: fmtx(3, 3)
  !! Z * Z * G * cos(arg)
  REAL(KIND = DP) :: reff(2, 2)
  !! Effective screening length for 2D materials  
  REAL(KIND = DP) :: grg
  !! G-vector * reff * G-vector for 2D loto
  COMPLEX(KIND = DP) :: facg
  !! Factor
  !
  ! alph is the Ewald parameter, geg is an estimate of G^2
  ! such that the G-space sum is convergent for that alph
  ! very rough estimate: geg/4/alph > gmax = 14
  ! (exp (-14) = 10^-6)
  !
  gmax = 20.d0
!  IF(present(alpha)) THEN
!    alph = alpha
!  ELSE
    alph = 1.d0
    !alph =  (0.5_dp * alat / tpi) ** 2
    IF(alph/=1.d0) print*, "alph", alph
!  ENDIF
  !
  geg  = gmax * alph * 4.0d0
  ! 
  ! Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
  ! Only for dimensions where periodicity is present, e.g. if nr1=1
  ! and nr2=1, then the G-vectors run along nr3 only.
  ! (useful if system is in vacuum, e.g. 1D or 2D)
  !
  IF (nr1 == 1) THEN
    nr1x = 0
  ELSE
    nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
  ENDIF
  IF (nr2 == 1) THEN
    nr2x = 0
  ELSE
    nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
  ENDIF
  IF (nr3 == 1) THEN
    nr3x=0
  ELSE
    nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
  ENDIF
  !
  IF (ABS(sign) /= 1.0_DP) CALL errore('rgd_blk',' wrong value for sign ',1)
  !
  IF (loto_2d) THEN 
    ! (e^2 * 2\pi) / Area     
    fac = (sign * e2 * tpi) / (omega * bg(3, 3) / alat)
    ! Effective screening length
    ! reff = (epsil - 1) * c/2 
    reff(:, :) = 0.0d0
    reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * tpi / bg(3, 3) ! (eps)*c/2 in 2pi/a units
    reff(1, 1) = reff(1, 1) - 0.5d0 * tpi / bg(3, 3) ! (-1)*c/2 in 2pi/a units
    reff(2, 2) = reff(2, 2) - 0.5d0 * tpi / bg(3, 3) ! (-1)*c/2 in 2pi/a units    
  ELSE
    ! (e^2 * 4\pi) / Volume
    fac = (sign * e2 * fpi) / omega          
  ENDIF
  DO m1 = -nr1x, nr1x
    DO m2 = -nr2x, nr2x
      DO m3 = -nr3x, nr3x
        !
        g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3)
        g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3)
        g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3)         
        !
        IF (loto_2d) THEN 
          geg = g1**2 + g2**2 + g3**2
          grg = 0.0d0
          IF (g1**2 + g2**2 > 1.0d-8) THEN
            grg = g1 * reff(1, 1) * g1 + g1 * reff(1, 2) * g2 + g2 * reff(2, 1) * g1 + g2 * reff(2, 2) * g2
            grg = grg / (g1**2 + g2**2)
          ENDIF
        ELSE
          geg = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))             
        ENDIF
        !
        IF (geg > 0.0d0 .AND. geg / (alph * 4) < gmax) THEN        
          !
          IF (loto_2d) THEN
            facgd = fac * (tpi / alat) * EXP(-geg / (alph * 4)) / (SQRT(geg) * (1.0 + grg * SQRT(geg)))       
          ELSE
            facgd = fac * EXP(-geg / (alph * 4)) / geg
          ENDIF
          !
!$OMP PARALLELDO DEFAULT(shared) PRIVATE(zcg,zag,arg,nb,na,i,j,fnat)
          DO na = 1, nat
            zag(:) = g1 * zeu(1, :, na) + g2 * zeu(2, :, na) + g3 * zeu(3, :, na)
            fnat(:) = 0.d0
            DO nb = 1, nat
              arg = 2.d0 * pi * (g1 * (tau(1, na) - tau(1, nb)) + &
                                 g2 * (tau(2, na) - tau(2, nb)) + &
                                 g3 * (tau(3, na) - tau(3, nb)))
              zcg(:)  = g1 * zeu(1, :, nb) + g2 * zeu(2, :, nb) + g3 * zeu(3, :, nb)
              fnat(:) = fnat(:) + zcg(:) * COS(arg)
            ENDDO
            ! Impose Hermiticity on long-range part of dynamical matrix
            ! Symmetrize long-range part of the on-site dynamical matrix by its conjugate transpose
            fmtx = MATMUL(RESHAPE(zag, (/3,1/)), RESHAPE(fnat, (/1,3/)))
            fmtx = (fmtx + TRANSPOSE(fmtx)) / 2.d0
            DO j = 1, 3
              DO i = 1, 3
                dyn(i, j, na, na) = dyn(i, j, na, na) - facgd * fmtx(i,j)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! nat
!$OMP END PARALLELDO
        ENDIF ! geg
        !
        g1 = g1 + q(1)
        g2 = g2 + q(2)
        g3 = g3 + q(3)
        !
        IF (loto_2d) THEN 
          geg = g1**2 + g2**2 + g3**2
          grg = 0.0d0
          IF (g1**2 + g2**2 > 1d-8) THEN
            grg = g1 * reff(1, 1) * g1 + g1 * reff(1, 2) * g2 + g2 * reff(2, 1) * g1 + g2 * reff(2, 2) * g2
            grg = grg / (g1**2 + g2**2)
          ENDIF
        ELSE
          geg = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))
        ENDIF
        !
        IF (geg > 0.0d0 .AND. geg / (alph * 4) < gmax) THEN        
          !
          IF (loto_2d) THEN
            facgd = fac * (tpi / alat) * EXP(-geg / (alph * 4)) / (SQRT(geg) * (1.0 + grg * SQRT(geg)))
          ELSE
            facgd = fac * EXP(-geg / (alph * 4)) / geg
          ENDIF
          !
!$OMP PARALLELDO DEFAULT(shared) PRIVATE(zbg,zag,arg,nb,na,i,j,facg)
          DO nb = 1, nat
            zbg(:) = g1 * zeu(1, :, nb) + g2 * zeu(2, :, nb) + g3 * zeu(3, :, nb)
            DO na = 1, nat
              zag(:) = g1 * zeu(1, :, na) + g2 * zeu(2, :, na) + g3 * zeu(3, :, na)
              arg = 2 * pi * (g1 * (tau(1, na) - tau(1 ,nb)) + &
                                 g2 * (tau(2, na) - tau(2, nb)) + &
                                 g3 * (tau(3, na) - tau(3, nb)) )
              facg = facgd * CMPLX(COS(arg), SIN(arg), KIND=DP)            
              !  
              DO j = 1, 3
                DO i = 1, 3
                  dyn(i, j, na, nb) = dyn(i, j, na, nb) + facg * zag(i) * zbg(j)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! na
          ENDDO ! nb
!$OMP END PARALLELDO
        ENDIF 
      ENDDO ! m3 
    ENDDO ! m2
  ENDDO ! m1
  !
  RETURN
!-----------------------------------------------------------------------
END SUBROUTINE rgd_blk_d3
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine nonanal_d3(nat, nat_blk, itau_blk, epsil, q, zeu, omega, dyn )
  !-----------------------------------------------------------------------
  !     add the nonanalytical term with macroscopic electric fields
  !
  use kinds, only: dp
  use constants, only: pi, fpi, e2
 implicit none
 integer, intent(in) :: nat, nat_blk, itau_blk(nat)
 !  nat: number of atoms in the cell (in the supercell in the case
 !       of a dyn.mat. constructed in the mass approximation)
 !  nat_blk: number of atoms in the original cell (the same as nat if
 !       we are not using the mass approximation to build a supercell)
 !  itau_blk(na): atom in the original cell corresponding to
 !                atom na in the supercell
 !
 complex(DP), intent(inout) :: dyn(3,3,nat,nat) ! dynamical matrix
 real(DP), intent(in) :: q(3),  &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nat_blk),   &! effective charges tensor
      &       omega            ! unit cell volume
 !
 ! local variables
 !
 real(DP) zag(3),zbg(3),  &! eff. charges  times g-vector
      &       qeq              !  <q| epsil | q>
 integer na,nb,              &! counters on atoms
      &  na_blk,nb_blk,      &! as above for the original cell
      &  i,j                  ! counters on cartesian coordinates
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
!print*, q(1), q(2), q(3)
 if (qeq < 1.d-8) then
    write(6,'(5x,"A direction for q was not specified:", &
      &          "TO-LO splitting will be absent")')
    return
 end if
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
!             print*, zag(i),zbg(j),qeq, fpi*e2*zag(i)*zbg(j)/qeq/omega
          end do
       end do
    end do
 end do
 !
 return
end subroutine nonanal_d3
END MODULE rigid_d3
