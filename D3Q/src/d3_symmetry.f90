!
! Copyright (C) 2001-2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE d3_symmetry
  USE kinds, only: DP
  !
  PRIVATE
  !
  TYPE d3_symmetry_type
    INTEGER :: nirr    ! number of Irreducible Representations
    INTEGER :: nsymq   ! the order of the small group
    ! i.e. subset of nirr that preserves the symmetry of q
    INTEGER,ALLOCATABLE :: npert(:)! number of perturbations per IR
    INTEGER             :: npertx  ! MAXVAL(npert(:))
    !
    COMPLEX(DP),ALLOCATABLE :: t(:,:,:,:) ! ...
    ! the full set of symmetry operations, the first nsym
    ! are the ones that respect the symmetry of the crystal lattice,
    ! and atoms basis. The first nsymq (.le.nsym) also respect
    ! the symmetry of the q vector (Sq = q)
    LOGICAL :: sym(48)    ! .true. if the corresponding element of t() is a sym.op.
    LOGICAL :: minus_q    ! if true one symmetry send q -> -q+G
    COMPLEX(DP),ALLOCATABLE :: tmq(:,:,:)
    ! the symmetry q<->-q in the base of the pattern
    !
    INTEGER :: irotmq     ! the symmetry sending q -> -q+
    INTEGER :: irgq(48)   ! the small group of q
    REAL(DP):: gi(3, 48)  ! [S(irotq)*q - q]
    REAL(DP):: gimq(3)    ! [S(irotmq)*q + q]
  END TYPE d3_symmetry_type
  ! One of the above for each q vector, and one for Gamma.
  TYPE(d3_symmetry_type),TARGET  :: symq(1:3) ! symmetry of each q individually (the same as phonon at that q)
  ! symmetry of the crystal (at Gamma), will point to one of symq if appropriate (used in efermi_shift.f90)
  TYPE(d3_symmetry_type),POINTER :: sym_gamma=>null()
  !
  ! The public variables and subroutines
  PUBLIC :: symq, sym_gamma, d3_symmetry_type
  PUBLIC :: d3_set_sym_irr, d3_symmetrize, minus_3q
  PUBLIC :: allocate_d3_symmetry, deallocate_d3_symmetry
  !
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE allocate_d3_symmetry(nat, npertx, symq)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nat, npertx
  TYPE(d3_symmetry_type),INTENT(INOUT) :: symq
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  IF(npertx<=0) THEN
    ALLOCATE( symq%npert(3*nat) )
  ELSE
    ALLOCATE( symq%t(npertx, npertx, 48, 3*nat) )
    ALLOCATE( symq%tmq(npertx, npertx, 3*nat) )
  ENDIF
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE allocate_d3_symmetry
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE deallocate_d3_symmetry(symq)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(d3_symmetry_type),INTENT(INOUT) :: symq
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
!   IF(allocated(symq%npert)) DEALLOCATE( symq%npert )
!   IF(allocated(symq%t))     DEALLOCATE( symq%t )
!   IF(allocated(symq%tmq))   DEALLOCATE( symq%tmq )
  DEALLOCATE( symq%npert )
  DEALLOCATE( symq%t )
  DEALLOCATE( symq%tmq )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE deallocate_d3_symmetry
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
SUBROUTINE d3_set_sym_irr(nat, at, bg, xq, s, invs, nsym, rtau, irt, &
                          irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, &
                          npert, nirr, gi, gimq)
  !---------------------------------------------------------------------
  ! FIXME: make stuff only once in parallel, then broadcast
  !
  !     This subroutine computes a basis for all the irreducible
  !     representations of the small group of q, which are contained
  !     in the representation which has as basis the displacement vectors.
  !     This is achieved by building a random hermitean matrix,
  !     symmetrizing it and diagonalizing the result. The eigenvectors
  !     give a basis for the irreducible representations of the
  !     small group of q.
  !
  !     Furthermore it computes:
  !     1) the small group of q
  !     2) the possible G vectors associated to every symmetry operation
  !     3) the matrices which represent the small group of q on the
  !        pattern basis.
  !
  !     Original routine was from C. Bungaro.
  !     Revised Oct. 1995 by Andrea Dal Corso.
  !     April 1997: parallel stuff added (SdG)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
!  USE mp_global, ONLY : mpime, root
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  !   first the dummy variables
  !
  INTEGER,INTENT(in) ::  nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), nirr, npertx
  ! input: the number of atoms
  ! input: the number of symmetries
  ! input: the symmetry matrices
  ! input: the inverse of each matrix
  ! input: the rotated of each atom
  ! input: write control
  INTEGER,INTENT(in) :: npert (3 * nat)
  ! output: the dimension of each represe
  INTEGER,INTENT(out) :: irgq (48), nsymq, irotmq
  ! output: the small group of q
  ! output: the order of the small group
  ! output: the symmetry sending q -> -q+
  ! output: the number of irr. representa

  REAL(DP),INTENT(in)  :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3)
  ! input: the q point
  ! input: the R associated to each tau
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  REAL(DP),INTENT(out) :: gi (3, 48), gimq (3)
  ! output: [S(irotq)*q - q]
  ! output: [S(irotmq)*q + q]
  COMPLEX(DP),INTENT(in)  :: u (3 * nat, 3 * nat)
  ! input: the pattern vectors
  COMPLEX(DP),INTENT(out) :: t (npertx, npertx, 48, 3 * nat),   &
                             tmq (npertx, npertx, 3 * nat)
  ! output: the symmetry matrices
  ! output: the matrice sending q -> -q+G
  LOGICAL,INTENT(out) :: minus_q
  ! output: if true one symmetry send q -> -q+G
  !
  INTEGER :: na, imode, jmode, ipert, jpert, nsymtot, imode0, &
       irr, ipol, jpol, isymq, irot, sna
  ! counter on atoms
  ! counter on atoms
  ! counter on modes
  ! counter on modes
  ! counter on perturbations
  ! counter on perturbations
  ! total number of symmetries
  ! auxiliry variable for mode counting
  ! counter on irreducible representation
  ! counter on polarizations
  ! counter on polarizations
  ! counter on symmetries
  ! counter on rotations
  ! the rotated atom

  REAL(DP) :: arg
  ! the argument of the phase

  COMPLEX(DP) :: wrk_u (3, nat), wrk_ru (3, nat), fase
  ! the dynamical matrix
  ! the bi-dimensional dynamical ma
  ! one pattern
  ! the rotated of one pattern
  ! the phase factor

!  LOGICAL :: lgamma
  ! if true gamma point

!   IF ( mpime == root ) THEN
     !
     !   Allocate the necessary quantities
     !
!     lgamma = (xq(1).EQ.0._dp .AND. xq(2).EQ.0._dp .AND. xq(3).EQ.0._dp)
     !
     !   find the small group of q
     !
     CALL smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
     !
     !   And we compute the matrices which represent the symmetry transformations
     !   in the basis of the displacements
     !
     t(:,:,:,:) = (0._dp, 0._dp) 
     tmq(:,:,:) = (0._dp, 0._dp) 
     IF (minus_q) THEN
        nsymtot = nsymq + 1
     ELSE
        nsymtot = nsymq
     ENDIF
     !
     DO isymq = 1, nsymtot
        IF (isymq.LE.nsymq) THEN
           irot = irgq (isymq)
        ELSE
           irot = irotmq
        ENDIF
        imode0 = 0
        DO irr = 1, nirr
           DO ipert = 1, npert (irr)
              imode = imode0 + ipert
              DO na = 1, nat
                 DO ipol = 1, 3
                    jmode = 3 * (na - 1) + ipol
                    wrk_u (ipol, na) = u (jmode, imode)
                 ENDDO
              ENDDO
              !
              !     transform this pattern to crystal basis
              !
              DO na = 1, nat
                 CALL trnvecc (wrk_u (1, na), at, bg, - 1)
              ENDDO
              !
              !     the patterns are rotated with this symmetry
              !
              wrk_ru(:,:) = (0._dp, 0._dp)
              DO na = 1, nat
                 sna = irt (irot, na)
                 arg = 0._dp
                 DO ipol = 1, 3
                    arg = arg + xq (ipol) * rtau (ipol, irot, na)
                 ENDDO
                 arg = arg * tpi
                 IF (isymq.EQ.nsymtot.AND.minus_q) THEN
                    fase = CMPLX(COS (arg), SIN (arg) ,kind=DP)
                 ELSE
                    fase = CMPLX(COS (arg), - SIN (arg) ,kind=DP)
                 ENDIF
                 DO ipol = 1, 3
                    DO jpol = 1, 3
                       wrk_ru (ipol, sna) = wrk_ru (ipol, sna) + s (jpol, ipol, irot) &
                            * wrk_u (jpol, na) * fase
                    ENDDO
                 ENDDO
              ENDDO
              !
              !    Transform back the rotated pattern
              !
              DO na = 1, nat
                 CALL trnvecc (wrk_ru (1, na), at, bg, 1)
              ENDDO
              !
              !     Computes the symmetry matrices on the basis of the pattern
              !
              DO jpert = 1, npert (irr)
                 imode = imode0 + jpert
                 DO na = 1, nat
                    DO ipol = 1, 3
                       jmode = ipol + (na - 1) * 3
                       IF (isymq.EQ.nsymtot.AND.minus_q) THEN
                          tmq(jpert, ipert, irr) = tmq (jpert, ipert, irr) &
                                                  + CONJG(u(jmode, imode) * wrk_ru(ipol, na))
                       ELSE
!                           print*, "doing t", jpert,ipert,irot,irr,npert(irr),na,wrk_ru(ipol,na)
                          t(jpert, ipert, irot, irr) = t(jpert, ipert, irot, irr) &
                               + CONJG(u(jmode, imode)) * wrk_ru(ipol, na)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           imode0 = imode0 + npert (irr)
        ENDDO

     ENDDO
     !
!   END IF
  !
  ! parallel stuff: first node broadcasts everything to all nodes
  !
!   CALL mp_bcast (gi, root)
!   CALL mp_bcast (gimq, root)
!   CALL mp_bcast (t, root)
!   CALL mp_bcast (tmq, root)
!   CALL mp_bcast (u, root)
!   CALL mp_bcast (nsymq, root)
!   CALL mp_bcast (npert, root)
!   CALL mp_bcast (nirr, root)
!   CALL mp_bcast (irotmq, root)
!   CALL mp_bcast (irgq, root)
!   CALL mp_bcast (minus_q, root)
  !
  RETURN
END SUBROUTINE d3_set_sym_irr
!
!-----------------------------------------------------------------------
SUBROUTINE d3_symmetrize(d3dyn, xq1,xq2,xq3, s, invs, rtau, irt, irgq, &
                         at, bg, nsymq, nat, irotmq, minus_q )
  !-----------------------------------------------------------------------
  !
  !    This routine symmetrizes the dynamical matrix written in the basis
  !    of the modes
  !
  !
  USE kinds,       ONLY : DP
!  USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm
!  USE mp,          ONLY : mp_sum
  USE io_global,   ONLY : stdout
  USE d3_basis,    ONLY : d3_3idx_2_6idx, d3_6idx_2_3idx, &
                          d3_cart2crys, d3_crys2cart
  USE d3_shuffle,  ONLY : d3_shuffle_global
  USE kplus3q,     ONLY : kplusq
  IMPLICIT NONE
  INTEGER,INTENT(in) :: nat, s(3, 3, 48), irt(48, nat), irgq(48), invs(48), &
                        nsymq, irotmq !npert_i, npert_f,
  ! the number of atoms
  ! the symmetry matrices
  ! the rotated of each atom
  ! the small group of q
  ! the inverse of each matrix
  ! the order of the small group
  ! the symmetry q -> -q+G

  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3) ! the coordinates of q
  REAL(DP),INTENT(in) :: rtau(3, 48, nat), at(3, 3), bg(3, 3)
  ! the R associated at each r
  ! direct lattice vectors
  ! reciprocal lattice vectors

  LOGICAL,INTENT(in) :: minus_q ! if true symmetry sends q->-q

  COMPLEX(DP),INTENT(inout) :: d3dyn(3*nat, 3*nat, 3*nat) ! matrix to symmetrize
!   COMPLEX(DP),INTENT(in) ::u1(3*nat, 3*nat), u2(3*nat, 3*nat), u3(3*nat, 3*nat)! the patterns

!  INTEGER :: i, j, k, icart, jcart, kcart, na, nb, nc, nu1, nu2, nu3! counters

!  COMPLEX(DP) :: work, wrk(3, 3) ! auxiliary variables
  COMPLEX(DP),ALLOCATABLE :: p3tmp(:,:,:,:,:,:)! the dynamical matrix
  COMPLEX(DP),ALLOCATABLE :: d3tmp(:,:,:)      ! the dynamical matrix workspace
  LOGICAL :: any_gamma ! any of the q's is gamma
  !
  ! NOTE: imposing hemiticity seems to not be required as it is imposed 
  !       automatically by the shuffle mechanism in d3toten
  ! Start by imposing hermiticity 
  ! i.e. for (0,q,-q) matrix must be hermitean w.r.t exchange of index 2 and 3
  any_gamma = kplusq(1)%lgamma.or.kplusq(2)%lgamma.or.kplusq(3)%lgamma
  IF(any_gamma) ALLOCATE(d3tmp(3*nat,3*nat,3*nat))
  !
  IF (kplusq(1)%lgamma) THEN
    WRITE(stdout,'(7x,a)') "* imposing hermiticity on idx 2,3"
    d3tmp = d3dyn                                         !copy
    CALL d3_shuffle_global(nat,  1,2,3, 1,3,2, .true., d3tmp ) !shuffle and c.c.
    d3dyn = d3tmp+d3dyn                                   !add
    d3dyn = 0.5_dp * d3dyn                                !average
  ENDIF
  IF (kplusq(2)%lgamma) THEN
    WRITE(stdout,'(7x,a)') "* imposing hermiticity on idx 1,3"
    d3tmp = d3dyn
    CALL d3_shuffle_global(nat,  1,2,3, 3,2,1, .true., d3tmp )
    d3dyn = d3tmp+d3dyn
    d3dyn = 0.5_dp * d3dyn
  ENDIF
  IF (kplusq(3)%lgamma) THEN
    WRITE(stdout,'(7x,a)') "* imposing hermiticity on idx 1,2"
    d3tmp = d3dyn
    CALL d3_shuffle_global(nat,  1,2,3, 2,1,3, .true., d3tmp )
    d3dyn = d3tmp+d3dyn
    d3dyn = 0.5_dp * d3dyn
  ENDIF
  !
  IF(any_gamma) DEALLOCATE(d3tmp)
  !
  ALLOCATE(p3tmp( 3,3,3, nat,nat,nat))
  ! Repack the matrix with dimension (3,3,3, nat,nat,nat) ...this also make a copy of it
  CALL d3_3idx_2_6idx(nat, d3dyn, p3tmp)
  !
  ! Then we transform to the crystal axis
  CALL d3_cart2crys(nat, at, bg, p3tmp)
  !
  ! And symmetrize in this basis
  CALL d3_symdyn_core(xq1, xq2, xq3, p3tmp, s, invs, rtau, irt, irgq, nsymq, nat, &
                      irotmq, minus_q)
  !
  !  Back to cartesian coordinates
  CALL d3_crys2cart(nat, at, bg, p3tmp)
  !
  !  rewrite the dynamical matrix on the array dyn with dimension (3nat)^3
  CALL d3_6idx_2_3idx(nat, p3tmp, d3dyn)
  !
  DEALLOCATE (p3tmp)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_symmetrize
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_symdyn_core(xq1,xq2,xq3, phi, s, invs, rtau, irt, irgq, nsymq, &
                          nat, irotmq, minus_q)
  !-----------------------------------------------------------------------
  !
  !     This routine receives as input an unsymmetrized dynamical
  !     matrix expressed on the crystal axes and imposes the symmetry
  !     of the small group of q. Furthermore it imposes also the symmetry
  !     q -> -q+G if present.
  !
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER,INTENT(in) :: nat, s(3, 3, 48), irt(48, nat), irgq(48),&
                        invs(48), nsymq, irotmq
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each vector
  ! input: the small group of q
  ! input: the inverse of each matrix
  ! input: the order of the small gro
  ! input: the rotation sending q ->
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  ! input: the q point
  REAL(DP),INTENT(in) :: rtau(3, 48, nat)
  ! input: the R associated at each t

  LOGICAL,INTENT(in) :: minus_q
  ! input: true if a symmetry q->-q+G
  COMPLEX(DP),INTENT(inout) :: phi(3, 3, 3, nat, nat, nat)
  ! inp/out: the matrix to symmetrize
  !
  !   local variables
  !
  INTEGER :: isymq, irot, &
             na=0,  nb=0,  nc=0, &
             sna=0, snb=0, snc=0, &
             ipol, jpol, lpol, kpol, mpol, npol
  ! counters
  LOGICAL, ALLOCATABLE:: done_equiv(:,:,:)
  ! used to account for symmetrized elements
  REAL(DP) :: arg
  ! the argument of the phase
  COMPLEX(DP) :: faseq(48), fasemq, workmq
  COMPLEX(DP),ALLOCATABLE ::  work(:,:,:), phip(:,:,:, :,:,:)
  ! the phase factor
  ! the phases for each symmetry
  !
  !
  !
  !    If no symmetry is present we quit here
  IF((nsymq == 1) .and. (.not.minus_q) ) RETURN
  !
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
!   APPLY_3MINUSQ : &
!   IF(minus_q) THEN
!      ALLOCATE(phip(3, 3, 3, nat, nat, nat))
!      WRITE(stdout, "(7x, a, i3)") "* symmetrizing with q --> -q operation", irotmq
!      !
!      DO nc = 1, nat
!      DO na = 1, nat
!      DO nb = 1, nat
!         DO mpol = 1, 3
!         DO ipol = 1, 3
!         DO jpol = 1, 3
!           workmq =(0._dp, 0._dp)
!           arg = 0._dp
!           DO kpol = 1, 3
!               arg = arg - ( xq1(kpol) * rtau(kpol, irotmq, nc) + &
!                             xq2(kpol) * rtau(kpol, irotmq, na) + &
!                             xq3(kpol) * rtau(kpol, irotmq, nb) )
! 
!           ENDDO
!           arg = arg * tpi
!           fasemq = CMPLX(cos(arg), sin(arg) ,kind=DP)
!           !
!           snc = irt(irotmq, nc)
!           sna = irt(irotmq, na)
!           snb = irt(irotmq, nb)
!           DO npol = 1, 3
!           DO kpol = 1, 3
!           DO lpol = 1, 3
!               workmq = workmq + &
!                   fasemq * s(ipol, kpol, irotmq) * &
!                            s(jpol, lpol, irotmq) * &
!                            s(mpol, npol, irotmq) * &
!                           phi(npol, kpol, lpol, snc, sna, snb)* &
!                           fasemq
!           ENDDO
!           ENDDO
!           ENDDO
!           phip(mpol, ipol, jpol, nc, na, nb) &
!               = 0.5_dp*( phi(mpol, ipol, jpol, nc, na, nb) +CONJG(workmq) )
!         ENDDO
!         ENDDO
!         ENDDO
!     ENDDO
!     ENDDO
!     ENDDO
!     !
!     phi = phip
!     !
!     DEALLOCATE(phip)
!   ENDIF & ! (minus_q)
!   APPLY_3MINUSQ 
  !
  !
  !    Here we symmetrize with respect to the small group of q
  !
  IF(nsymq == 1) RETURN
  !
  WRITE(stdout,'(7x,a,i3,a)') "* symmetrizing with",nsymq," sym.ops."
  ALLOCATE(work(3, 3, 3))
  ALLOCATE(done_equiv( nat, nat, nat))
  done_equiv = .false.
  !
  DO nc = 1, nat
     DO na = 1, nat
        DO nb = 1, nat
          !  ===== nc ==== na ===== nb ====
          done_equiv_TODO : &
          IF( .not.done_equiv(nc, na, nb) )THEN
            work =(0._dp, 0._dp)
            !
            ISYMQ_LOOP_1 : &
            DO isymq = 1, nsymq
              irot = irgq(isymq)
              arg = 0._dp
              DO ipol = 1, 3
                arg = arg - ( xq1(ipol) * rtau(ipol, irot, nc) + &
                              xq2(ipol) * rtau(ipol, irot, na) + &
                              xq3(ipol) * rtau(ipol, irot, nb) )
              ENDDO
              arg = arg * tpi
              faseq(isymq) = CMPLX(cos(arg), sin(arg) ,kind=DP)
              !
              snc = irt(irot, nc)
              sna = irt(irot, na)
              snb = irt(irot, nb)
              DO mpol = 1, 3
              DO ipol = 1, 3
              DO jpol = 1, 3
                  DO npol = 1, 3
                  DO kpol = 1, 3
                  DO lpol = 1, 3
                    work(mpol, ipol, jpol) = work(mpol, ipol, jpol) + &
                          s(ipol, kpol, irot) * &
                          s(jpol, lpol, irot) * &
                          s(mpol, npol, irot) * &
                          phi(npol, kpol, lpol, snc, sna, snb) * &
                          faseq(isymq)

                  ENDDO
                  ENDDO
                  ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO &
            ISYMQ_LOOP_1
            !
            ISYMQ_LOOP_2 : &
            DO isymq = 1, nsymq
              irot = irgq(isymq)
              snc = irt(irot, nc)
              sna = irt(irot, na)
              snb = irt(irot, nb)
              DO mpol = 1, 3
              DO ipol = 1, 3
              DO jpol = 1, 3
                phi(mpol, ipol, jpol, snc, sna, snb) = (0._dp, 0._dp)
                DO npol = 1, 3
                DO kpol = 1, 3
                DO lpol = 1, 3
                    phi(mpol, ipol, jpol, snc, sna, snb) = &
                        phi(mpol, ipol, jpol, snc, sna, snb) +&
                        s(mpol, npol, invs(irot) ) * &
                        s(ipol, kpol, invs(irot) ) * &
                        s(jpol, lpol, invs(irot) ) * &
                        work(npol, kpol, lpol) * &
                        CONJG(faseq(isymq) )
                ENDDO
                ENDDO
                ENDDO
              ENDDO
              ENDDO
              ENDDO
              !
              done_equiv(snc, sna, snb) = .true.
              !
            ENDDO &
            ISYMQ_LOOP_2
            !
          ENDIF &
          done_equiv_TODO
          !  ===== nc ==== na ===== nb ====
        ENDDO
     ENDDO
  ENDDO
  !
  phi = phi / REAL(nsymq, kind=DP)
  !
  DEALLOCATE(done_equiv, work)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_symdyn_core
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE minus_3q(xq1, xq2, xq3, at, bg, nsym, s, minus_q, irotmq)
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE kinds, only : DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nsym, s(3, 3, 48)
  ! nsym matrices of symmetry operations of the crystal (lattic+atoms)
  REAL(DP),INTENT(in) :: xq1(3), xq2(3), xq3(3)
  ! xq*: triplet of q vectors
  REAL(DP),INTENT(in) :: at(3, 3), bg(3, 3)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  LOGICAL,INTENT(inout) :: minus_q
  INTEGER,INTENT(out)   :: irotmq
  ! sym op linking q1,q2,q3 to -q1,-q2,-q3
  !
  INTEGER :: isym, iq, i
  ! number of symmetry ops. that generate a certain star (only used for sanity check)
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counter
  REAL(DP) :: x3q(3,3), a3q(3,3), ra3q(3,3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  REAL(DP),PARAMETER :: gam(3) = (/ 0._dp, 0._dp, 0._dp /)
  ! a gam vector: used in eqvect
 
  LOGICAL, EXTERNAL :: eqvect
  ! function used to compare two vectors
  !
  irotmq = 0
  IF(.not.minus_q) RETURN ! do not even look for it in this case
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
  !
  DO isym = 1, nsym
      !
      ! Compute the rotated triplet with this sym op
      DO iq = 1,3
        DO i = 1,3
          ra3q(i,iq) =  s(i,1,isym) * a3q(1,iq) &
                      + s(i,2,isym) * a3q(2,iq) &
                      + s(i,3,isym) * a3q(3,iq)
        ENDDO
      ENDDO
      !
      ! check if this triplet is equivalent to -q1,-q2,-q3 
      IF ( eqvect(ra3q(:,1), -a3q(:,1), gam, 1.d-5) .and. &
           eqvect(ra3q(:,2), -a3q(:,2), gam, 1.d-5) .and. &
           eqvect(ra3q(:,3), -a3q(:,3), gam, 1.d-5) ) THEN
          !
          irotmq = isym
          RETURN
          !
      ENDIF
      !
      ! if no previous sym.op. produces this triplet, create a new set
      !
  ENDDO
  !
  ! We did not find a "minus_q" operation common to all three vectors..
  minus_q = .false.
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE minus_3q
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE
!-----------------------------------------------------------------------
