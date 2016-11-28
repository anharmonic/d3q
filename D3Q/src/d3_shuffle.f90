!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE d3_shuffle
!-----------------------------------------------------------------------
  !
  LOGICAL :: d3_shuffle_initialized = .false.
  !
  INTEGER,PARAMETER :: nperms = 6         ! total number of permutations
  INTEGER           :: dperms = 0         ! number of distinct permutations
  INTEGER           :: iperms(nperms) = 0 ! index of distinct permutations
  INTEGER,PARAMETER :: d3perms_order(3,nperms) = &
        RESHAPE( (/ &
          1,  2,  3, &  ! NOTE: it is important to keep the permutations in
          1,  3,  2, &  !       this order for the terms drhod2v and the second
          2,  1,  3, &  !       part of the valence contribution!!
          2,  3,  1, &  !       (Actually, other choices are possible, but
          3,  1,  2, &  !       d3toten must be changed accordingly)
          3,  2,  1  &
         /), (/3,6/) )
  INTEGER,PARAMETER :: d3perms_order2(3,nperms) = &
        RESHAPE( (/ &
          1,  2,  3, &  ! NOTE: it is important to keep the permutations in
          2,  3,  1, &  !       this order for the terms drhod2v and the second
          3,  1,  2, &  !       part of the valence contribution!!
          3,  2,  1, &  !       (Actually, other choices are possible, but
          2,  1,  3, &  !       d3toten must be changed accordingly)
          1,  3,  2  &
         /), (/3,6/) )
  TYPE perturbation_permutation
    INTEGER :: i,j,k ! the actual permutation of 1,2,3
    LOGICAL :: todo  ! .true. if this permutation has to be computed
    INTEGER :: shuffle_from  ! the index of the permutation we have to shuffle
                             ! to obtain this (only set if todo is false)
    LOGICAL :: shuffle_conjg ! if .true. shuffle and take the complex conjugate
    LOGICAL :: todo_first    ! if .true. it's the first permutation with this value of i
    INTEGER :: first_from    ! if todo_first is .false. this other permutation has the same value of i
    CHARACTER(len=3) :: name
  END TYPE
  TYPE(perturbation_permutation) :: d3perms(nperms)
  !
 CONTAINS
! 
! NOTE about the shuffle global and equiv:
! call d3_shuffle_global(nat, i1,i2,i3, j1,j2,j3,..) 
!   is equivalent to
! call d3_shuffle_equiv(nat, j1,j2,j3, i1,i2,i3,..) 
! And, contrary to what one may think, exchanging the is with the js 
! gives a different result.
!
!-----------------------------------------------------------------------
SUBROUTINE d3_reset_permutations()
  !-----------------------------------------------------------------------
  ! not really much to do here, as variables are initialized in d3_check_permutations
  IMPLICIT NONE
  dperms = 0
  iperms = 0
  d3_shuffle_initialized = .false.
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_reset_permutations
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_check_permutations()
  !-----------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE kplus3q,    ONLY : kplusq, q_names, q_special_cases_initialized
  USE constants,  ONLY : eps8
  !
  IMPLICIT NONE
  !
  INTEGER :: iperm, jperm
  INTEGER :: i1,i2,i3,  j1,j2,j3
  LOGICAL :: same_perm, conjg_perm
  CHARACTER(len=21),PARAMETER :: sub='d3_check_permutations'
  !
  IF(.not.q_special_cases_initialized) &
    CALL errore(sub,'you must call q_special_cases first', 1)
  !
  IF(d3_shuffle_initialized) CALL errore(sub, sub//' already called',1)
  d3_shuffle_initialized = .true.
  
  WRITE(stdout,'(5x,a)') "Looking for equivalent permutations of q vectors"
  !
  ! Set up special variables used for subroutines that only depend on the value of
  ! the first index (%i), namely d3_valence_ij, d3_nlcc_123 and the non-local part of dq1rhodq23v
  ! (at the moment, this tick is only used for d3_valence_ij, where it matters more)
  !
  ! NOTE: next loop is only valid for this specific ordering of the permutations
  DO iperm = 1,nperms,2
    d3perms(iperm)%todo_first = .true.
    d3perms(iperm)%first_from = iperm ! points to itself
    !
    d3perms(iperm+1)%todo_first = .false.
    d3perms(iperm+1)%first_from = iperm ! points to previous one
  ENDDO
  !
  ! Initialize with the most general case
  DO iperm = 1,nperms
    d3perms(iperm)%i = d3perms_order(1,iperm)
    d3perms(iperm)%j = d3perms_order(2,iperm)
    d3perms(iperm)%k = d3perms_order(3,iperm)
    !
    d3perms(iperm)%todo = .true.
    d3perms(iperm)%shuffle_from = iperm ! points to itself
    d3perms(iperm)%shuffle_conjg = .false.
    !
    WRITE(d3perms(iperm)%name, '(3i1)')  d3perms_order(:,iperm)
  ENDDO
  !
  dperms = 0
  !
  IPERM_LOOP : &
  DO iperm = 1,nperms
    !
    i1 = d3perms_order(1,iperm)
    i2 = d3perms_order(2,iperm)
    i3 = d3perms_order(3,iperm)
    !
    JPERM_LOOP : &
    DO jperm = 1,iperm
    !
    IF( iperm == jperm) THEN
      ! if no match was found we have one more permutation to compute explicitly
      dperms = dperms + 1     ! keep track of the number...
      iperms(dperms) = iperm  ! ...and index it
      !
      WRITE(stdout,'(7x,a,i2,a,3i1,a)') "Permutation", iperm, " (",i1,i2,i3,") has to be computed explicitly"
      CYCLE IPERM_LOOP ! no need to check in this case
    ENDIF
    !
    j1 = d3perms_order(1,jperm)
    j2 = d3perms_order(2,jperm)
    j3 = d3perms_order(3,jperm)
    !
! #define __DEBUG_D3_SHUFFLE
#ifdef __DEBUG_D3_SHUFFLE
    same_perm = .false.
    conjg_perm = .false.
#else
    same_perm  = kplusq(i1)%lsame(j1)  .and. kplusq(i2)%lsame(j2)  .and. kplusq(i3)%lsame(j3)
    conjg_perm = kplusq(i1)%lsame(-j1) .and. kplusq(i2)%lsame(-j2) .and. kplusq(i3)%lsame(-j3)
!    same_perm =  ( SUM(ABS(kplusq(i1)%xq-kplusq(j1)%xq)) < eps8 ) &
!            .and.( SUM(ABS(kplusq(i2)%xq-kplusq(j2)%xq)) < eps8 ) &
!            .and.( SUM(ABS(kplusq(i3)%xq-kplusq(j3)%xq)) < eps8 )
    !
!    conjg_perm = ( SUM(ABS(kplusq(i1)%xq+kplusq(j1)%xq)) < eps8 ) &
!            .and.( SUM(ABS(kplusq(i2)%xq+kplusq(j2)%xq)) < eps8 ) &
!            .and.( SUM(ABS(kplusq(i3)%xq+kplusq(j3)%xq)) < eps8 )
#endif
    !
    IF (same_perm) THEN
      d3perms(iperm)%todo = .false.
      d3perms(iperm)%shuffle_from = jperm
      d3perms(iperm)%shuffle_conjg = .false.
      !
      WRITE(stdout,'(7x,a,i2,a,3i1,a,i2,a,3i1,a)') "Permutation", iperm, " (",i1,i2,i3,") can be obtained from", &
                            jperm, " (",j1,j2,j3, ") "
      !
      CYCLE IPERM_LOOP ! one is enough
    ELSE IF (conjg_perm) THEN
      d3perms(iperm)%todo = .false.
      d3perms(iperm)%shuffle_from = jperm
      d3perms(iperm)%shuffle_conjg = .true.
      !
      WRITE(stdout,'(7x,a,i2,a,3i1,a,i2,a,3i1,a)') "Permutation", iperm, " (",i1,i2,i3,") can be obtained from", &
                            jperm, " (",j1,j2,j3, ") with c.c."
      !
      CYCLE IPERM_LOOP ! one is enough
    ENDIF
    !
    ENDDO JPERM_LOOP
    !
  ENDDO IPERM_LOOP
  !
  WRITE(stdout,'(5x,i1,a)') dperms, " inequivalent permutations were found"
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_check_permutations
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_shuffle_global(nat, i1,i2,i3, j1,j2,j3, conjugate, d3dyn, d3dyn_out)
  !-----------------------------------------------------------------------
  ! WARNING! what you want to do may be done by d3_shuffle_equiv! READ THE DESCRIPTION:
  !
  ! Nice and simple index shuffling routine.
  ! It assumes that d3dyn indexes are i1-th, i2-th and i3-th respectively
  ! and reshuffle to become the j1-th, j2-th and j3-th.
  ! Writes to d3dyn_out, if present, overwrites input otherwise.
  !
  ! This subroutine can be used to shuffle a *global* and *external* permutation
  ! of q1,q2,q3.
  !
  USE kinds,      ONLY : DP
  !USE ions_base,  ONLY : nat
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: nat, i1,i2,i3,  j1,j2,j3
  LOGICAL,INTENT(in) :: conjugate
  COMPLEX(DP),VOLATILE,INTENT(inout) :: d3dyn( 3*nat, 3*nat, 3*nat)
  COMPLEX(DP),OPTIONAL,INTENT(inout) :: d3dyn_out( 3*nat, 3*nat, 3*nat)
  !
  INTEGER,VOLATILE,TARGET  :: idx(3)
  INTEGER,VOLATILE,POINTER :: nu1, nu2, nu3
  INTEGER,VOLATILE,POINTER :: mu1, mu2, mu3
  COMPLEX(DP),VOLATILE,ALLOCATABLE :: d3aux( :,:,: )
  CHARACTER(len=17),PARAMETER :: sub='d3_shuffle_global'

  IF( ANY(ABS((/ i1,i2,i3, j1,j2,j3 /))>3) .or. &
      ANY(ABS((/ i1,i2,i3, j1,j2,j3 /))<1) ) &
    CALL errore(sub, 'Index out of range', 3)

  idx = (/ 0,-1,1 /)

  nu1 => idx(abs(i1));  nu2 => idx(abs(i2));  nu3 => idx(abs(i3))
  mu1 => idx(abs(j1));  mu2 => idx(abs(j2));  mu3 => idx(abs(j3))
!   IF( i1/=j1 .and. i2/=j2 .and. i3/=j3 ) THEN
!     WRITE(stdout, '(9x,a,3i1,a,3i1,a,l2)') "* global shuffling: ",i1,i2,i3," --> ",j1,j2,j3, &
!                                            " c.c.", conjugate
!   ENDIF

  IF(nu1==nu2 .or. nu1==nu3 .or. nu2==nu3) &
    CALL errore(sub, 'First indexes repeat', 1)
  IF(mu1==mu2 .or. mu1==mu3 .or. mu2==mu3) &
    CALL errore(sub, 'Second indexes repeat', 2)

  ALLOCATE(d3aux(3*nat, 3*nat, 3*nat))

  DO nu3 = 1,3*nat
    DO nu2 = 1,3*nat
      DO nu1 = 1,3*nat
        d3aux(mu1,mu2,mu3) = d3dyn(nu1,nu2,nu3)
      ENDDO
    ENDDO
  ENDDO
  !
  IF(conjugate) d3aux = CONJG(d3aux)
  !
  IF(present(d3dyn_out)) THEN
    d3dyn_out = d3aux
  ELSE
    d3dyn = d3aux
  ENDIF

  DEALLOCATE(d3aux)
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_shuffle_global
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_shuffle_equiv(nat, i1,i2,i3, j1,j2,j3, conjugate, d3dyn, d3dyn_out)
  !-----------------------------------------------------------------------
  ! Ugly and complicated reshuffling routine.
  !
  ! Assumes that d3dyn has been computed with a certain permutation of perturbation 
  ! q_i1,q_i2,q_i3. Nonetheless, first index of d3dyn contains the mode number of
  ! perturbation q1, its second perturbation q2 and its third perturbation q3. 
  !
  ! On output, d3dyn for the permutaion q_j1, q_j2 q_j3 equivalent to q_i1,q_i2,q_i3
  ! and again with q1, q2, q3 in respectively the first, second and third index of
  ! the matrix.
  !
  ! This subroutine shuffles and *internal* permutation of the perturbation and
  ! only makes sense for q_ix=+/-q_jx, x=1,2,3 (taking the complex conjugate if
  ! the sign is -).
  !
  USE kinds,      ONLY : DP
 ! USE ions_base,  ONLY : nat
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: nat, i1,i2,i3,  j1,j2,j3
  LOGICAL,INTENT(in) :: conjugate
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3*nat, 3*nat, 3*nat)
  COMPLEX(DP),OPTIONAL,INTENT(inout) :: d3dyn_out( 3*nat, 3*nat, 3*nat)
  !
  TYPE pointer_array
    INTEGER,POINTER :: x
  END TYPE
  TYPE(pointer_array),VOLATILE :: nu(3), mu(3)
  !
  INTEGER,VOLATILE,TARGET  :: idx1, idx2, idx3
  COMPLEX(DP),VOLATILE,ALLOCATABLE :: d3aux( :,:,: )
  CHARACTER(len=16),PARAMETER :: sub='d3_shuffle_equiv'

  IF( ANY((/ i1,i2,i3, j1,j2,j3 /)>3) .or. &
      ANY((/ i1,i2,i3, j1,j2,j3 /)<1) ) &
    CALL errore(sub, 'Index out of range', 3)

  idx1=0; idx2=1; idx3=-1

  nu(i1)%x => idx1; mu(j1)%x => idx1
  nu(i2)%x => idx2; mu(j2)%x => idx2
  nu(i3)%x => idx3; mu(j3)%x => idx3

  !WRITE(10099, '(6i3,l1)') i1,i2,i3,j1,j2,j3, conjugate
  !WRITE(10099,'(3(2f12.6,3x))') d3dyn 
  !WRITE(10099,'(/,a)')

!   WRITE(stdout, '(9x,a,3i1,a,3i1,a,l2)') "* shuffling: ",i1,i2,i3," --> ",j1,j2,j3, &
!                                           " - c.c.", conjugate

  IF(nu(1)%x==nu(2)%x .or. nu(1)%x==nu(3)%x .or. nu(2)%x==nu(3)%x) &
    CALL errore(sub, 'First indexes repeat', 1)
  IF(mu(1)%x==mu(2)%x .or. mu(1)%x==mu(3)%x .or. mu(2)%x==mu(3)%x) &
    CALL errore(sub, 'Second indexes repeat', 2)

  ALLOCATE(d3aux(3*nat, 3*nat, 3*nat))

  DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
      DO idx3 = 1,3*nat
        d3aux(mu(1)%x,mu(2)%x,mu(3)%x) = d3dyn(nu(1)%x,nu(2)%x,nu(3)%x)
        !WRITE(10099,'(3(3i3,3x))') idx1,idx2,idx3, nu(1)%x, nu(2)%x, nu(3)%x, mu(1)%x, mu(2)%x, mu(3)%x
      ENDDO
    ENDDO
  ENDDO
  !
  IF(conjugate) d3aux = CONJG(d3aux)
  !WRITE(10099,'(3(2f12.6,3x))') d3dyn_out
  !
  IF(present(d3dyn_out)) THEN
    d3dyn_out = d3aux
  ELSE
    d3dyn = d3aux
  ENDIF
  !
  DEALLOCATE(d3aux)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_shuffle_equiv
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
FUNCTION d3_shuffle_star(i,j,k, nqst, sxq, at, time_rev)
  !-----------------------------------------------------------------------
  ! return .true. if a certain permutation of a q triplet is already in the star,
  ! return .false. otherwise
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  LOGICAL :: d3_shuffle_star
  !
  INTEGER,INTENT(in)  :: i,j,k
  INTEGER,INTENT(in)  :: nqst
  REAL(DP),INTENT(in) :: sxq(3,3,48)
  REAL(DP),INTENT(in) :: at(3,3)
  LOGICAL,INTENT(in)  :: time_rev
  !
  REAL(DP) :: axq(3,3,nqst)
  INTEGER  :: istq, iq, ipol
  !
  LOGICAL,EXTERNAL :: eqvect
  !
  REAL(DP),PARAMETER :: gam(3) = (/ 0._dp, 0._dp, 0._dp /)
  INTEGER,PARAMETER  :: iorq = 1 ! Index of the original q's in the star
  CHARACTER(len=15),PARAMETER :: sub = 'd3_shuffle_star'
  !
  IF( ANY((/ i==j, i==k, j==k/)) .or. &
      ANY((/i,j,k/) > 3) .or. ANY((/i,j,k/) < 1) ) &
    CALL errore(sub, 'Invalid permutation.', 1)
  !
  d3_shuffle_star = .false.
  !
  ! We assume the original triplet to be the first in the star
  !
  DO istq = 1, nqst ! <- index of the triplet in the star
    DO iq = 1,3  ! <- index of the q vector
      DO ipol = 1,3 ! <- index of the crystal direction
        axq(ipol,iq,istq) =  sxq(1,iq,istq) * at(1,ipol) &
                           + sxq(2,iq,istq) * at(2,ipol) &
                           + sxq(3,iq,istq) * at(3,ipol)
      ENDDO
    ENDDO
  ENDDO
  !
  DO istq = 1, nqst
    ! Check if any q in the star is equivalent to this permutation
    d3_shuffle_star = &
            eqvect(axq(:,i,iorq), axq(:,1,istq), gam, 1.d-5) &
      .and. eqvect(axq(:,j,iorq), axq(:,2,istq), gam, 1.d-5) &
      .and. eqvect(axq(:,k,iorq), axq(:,3,istq), gam, 1.d-5)
    ! if there is one, we can return immediately
    write(*,'(L1,6(3f10.4,3x))') d3_shuffle_star, &
            axq(:,i,iorq), axq(:,j,iorq), axq(:,k,iorq), &
            axq(:,1,istq), axq(:,2,istq), axq(:,3,istq)
    IF(d3_shuffle_star) RETURN
    !
    ! If asked, also check if we can get a match by sending q_x->-q_x (time reversal used)
    IF(time_rev) THEN
      d3_shuffle_star = &
              eqvect(-axq(:,i,iorq), axq(:,1,istq), gam, 1.d-5) &
        .and. eqvect(-axq(:,j,iorq), axq(:,2,istq), gam, 1.d-5) &
        .and. eqvect(-axq(:,k,iorq), axq(:,3,istq), gam, 1.d-5)
      ! again, return immediately if possible
      IF(d3_shuffle_star) RETURN
    ENDIF
    !
  ENDDO
  !
  ! If not match is found, we return here (d3_shuffle_star is still .false.)
  RETURN
  !-----------------------------------------------------------------------
END FUNCTION d3_shuffle_star
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE d3_shuffle_stupid(nat, j1,j2,j3, conjugate, d3dyn)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !USE ions_base,  ONLY : nat
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: nat, j1,j2,j3
  LOGICAL,INTENT(in) :: conjugate
  COMPLEX(DP),INTENT(inout) :: d3dyn( 3*nat, 3*nat, 3*nat)
  !
  INTEGER :: idx1, idx2, idx3
  COMPLEX(DP),ALLOCATABLE :: d3aux( :,:,: )
  CHARACTER(len=16),PARAMETER :: sub='d3_shuffle_stupid'

  ALLOCATE(d3aux(3*nat,3*nat,3*nat))
  
  IF(j1==1.and.j2==2.and.j3==3)THEN
    d3aux = d3dyn
  ELSE &
  IF(j1==1.and.j2==3.and.j3==2)THEN
    DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
    DO idx3 = 1,3*nat
      !d3aux(idx1, idx2, idx3) = d3dyn(idx1,idx3,idx2)
      d3aux(idx1, idx3, idx2) = d3dyn(idx1,idx2,idx3)
    ENDDO
    ENDDO
    ENDDO
  ELSE &
  IF(j1==2.and.j2==1.and.j3==3)THEN
    DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
    DO idx3 = 1,3*nat
      !d3aux(idx1, idx2, idx3) = d3dyn(idx2,idx1,idx3)
      d3aux(idx2, idx1, idx3) = d3dyn(idx1,idx2,idx3)
    ENDDO
    ENDDO
    ENDDO
  ELSE &
  IF(j1==2.and.j2==3.and.j3==1)THEN
    DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
    DO idx3 = 1,3*nat
      !d3aux(idx1, idx2, idx3) = d3dyn(idx2,idx3,idx1)
      d3aux(idx2, idx3, idx1) = d3dyn(idx1,idx2,idx3)
    ENDDO
    ENDDO
    ENDDO
  ELSE &
  IF(j1==3.and.j2==1.and.j3==2)THEN
    DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
    DO idx3 = 1,3*nat
      !d3aux(idx1, idx2, idx3) = d3dyn(idx3,idx1,idx2)
      d3aux(idx3, idx1, idx2) = d3dyn(idx1,idx2,idx3)
    ENDDO
    ENDDO
    ENDDO
  ELSE &
  IF(j1==3.and.j2==2.and.j3==1)THEN
    DO idx1 = 1,3*nat
    DO idx2 = 1,3*nat
    DO idx3 = 1,3*nat
      !d3aux(idx1, idx2, idx3) = d3dyn(idx3,idx2,idx1)
      d3aux(idx3, idx2, idx1) = d3dyn(idx1,idx2,idx3)
    ENDDO
    ENDDO
    ENDDO
  ELSE
    CALL errore(sub, "unexpected permutation",1)
  ENDIF
  !
  IF(conjugate) d3aux = CONJG(d3aux)
  !
  d3dyn = d3aux
  
!   IF(present(d3dyn_out)) THEN
!     d3dyn_out = d3aux
!   ELSE
!     d3dyn = d3aux
!   ENDIF
  !
  DEALLOCATE(d3aux)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_shuffle_stupid
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
END MODULE d3_shuffle
!-----------------------------------------------------------------------
