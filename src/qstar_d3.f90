!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qstar_d3_module
!
!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
 CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE star_3q(xq1, xq2, xq3, at, bg, nsym, s, invs, nst3q, sx3q, is3q, im3q )
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE kinds, only : DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nsym, s(3, 3, 48), invs(48)
  ! nsym matrices of symmetry operations of the crystal (lattic+atoms)
  ! invs: list of inverse operation indices: s(:,:,invs(j)) = s(:,:,j)^-1
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  ! xq*: triplet of q vectors
  REAL(DP),INTENT(in) :: at(3, 3), bg(3, 3)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  INTEGER,INTENT(out) :: nst3q
  ! n3q : degeneracy of the star of q triplets
  INTEGER,INTENT(out) :: is3q(48), im3q
  ! is3q : index of q in the star for a given sym
  ! im3q  : index of -q in the star (0 IF not present)
  !
  REAL(DP),INTENT(out) :: sx3q(3, 3, 48)
  ! list of vectors triples in the star of the 3 input q's
  !
  INTEGER :: nsq(48), isym, ism1, iq, istq, i
  ! number of symmetry ops. that generate a certain star (only used for sanity check)
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counter
  REAL(DP) :: x3q(3,3), sa3q (3,3, 48), a3q (3,3), ra3q (3,3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  REAL(DP),PARAMETER :: gam(3) = (/ 0._dp, 0._dp, 0._dp /)
  ! a gam vector: used in eqvect

  LOGICAL, EXTERNAL :: eqvect
  ! function used to compare two vectors
  !
  ! initialize
  x3q(:,1) = xq1(:)
  x3q(:,2) = xq2(:)
  x3q(:,3) = xq3(:)
  !
  ! go to  crystal coordinates
  DO iq = 1,3  ! <- index of the q vector
    DO i = 1,3 ! <- index of the crystal direction
      a3q(i,iq) =  x3q(1,iq) * at(1,i) &
                 + x3q(2,iq) * at(2,i) &
                 + x3q(3,iq) * at(3,i)
    ENDDO
  ENDDO
  !
  ! create the list of rotated q
  DO i = 1, 48
     nsq(i)  = 0
     is3q(i) = 0
  ENDDO
  !
  nst3q = 0
  DO isym = 1, nsym
      ism1 = invs(isym)
      !
      ! Compute the rotated triplet with this sym op
      DO iq = 1,3
        DO i = 1,3
          ra3q(i,iq) =  s(i,1,ism1) * a3q(1,iq) &
                      + s(i,2,ism1) * a3q(2,iq) &
                      + s(i,3,ism1) * a3q(3,iq)
        ENDDO
      ENDDO
      !
      ! check if this triplet has already been generated by another sym.op.
      DO istq = 1, nst3q ! <-- this loop is skipped when nst3q == 0
        IF ( eqvect(ra3q(:,1), sa3q(:,1,istq), gam, 1.d-5) .and. &
             eqvect(ra3q(:,2), sa3q(:,2,istq), gam, 1.d-5) .and. &
             eqvect(ra3q(:,3), sa3q(:,3,istq), gam, 1.d-5) ) THEN
            !
            is3q(isym)  = istq          ! <-- sym.op. isym produces the istq'th triplet
            nsq(istq)   = nsq(istq) + 1 ! <-- the istq'th triplet is generated by one more sym.op.
        ENDIF
      ENDDO
      !
      ! if no previous sym.op. produces this triplet, create a new set
      IF (is3q(isym) == 0) THEN
        !
        nst3q = nst3q + 1    ! <-- we have one new set
        nsq(nst3q) = 1       ! <-- generated by one sym.op. (the isym'th)
        is3q(isym) = nst3q   ! <-- isym'th sym.op. produces nst3q'th triplet
        ! save the new triplet in sa3q
        sa3q(:,:,nst3q) = ra3q(:,:)
        ! and take it back to cartesian axes in sx3q
        DO iq = 1, 3 
          DO i = 1, 3
              sx3q(i, iq, nst3q) =  bg(i, 1) * sa3q(1, iq, nst3q) &
                                  + bg(i, 2) * sa3q(2, iq, nst3q) &
                                  + bg(i, 3) * sa3q(3, iq, nst3q)
          ENDDO
        ENDDO
        !
      ENDIF
      !
  ENDDO
  !
  ! set im3q to the index of the star containing -q1,-q2,-q3 
  ! ...and check star degeneracy (just sanity check)
  !
  ra3q(:,:) = -a3q(:,:)
  im3q = 0
  DO istq = 1, nst3q
     IF (eqvect (ra3q(:,1), sa3q(1, 1, istq), gam, 1.d-5) .and. &
         eqvect (ra3q(:,2), sa3q(1, 2, istq), gam, 1.d-5) .and. &
         eqvect (ra3q(:,3), sa3q(1, 3, istq), gam, 1.d-5) ) im3q = istq
     IF(nsq(istq)*nst3q /= nsym) CALL errore ('star_q', 'wrong degeneracy', istq)
  ENDDO
  !
  ! writes star of q
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE star_3q
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE qstar_d3_module
!-----------------------------------------------------------------------

