
!
! Copyright (C) 2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE kplus3q
  !-----------------------------------------------------------------------
  ! This module takes the ikks and ikq2s mappings from the phonon code,
  ! furthermore it adds the additional mappings ikq1s and ikq3s corresponding
  ! to the grid of k+p and k-p-q points.
  !
  ! Inside the D3 code THIS modules should be used instead of taking the variables
  ! directly from qpoint, io_files, ...
  !
  USE kinds,          ONLY : DP
  USE qpoint,         ONLY : nksq
  USE control_lr,     ONLY : lgamma
  !
  INTEGER             :: orig_nks, &! number of points in the original grid (i.e. k)
                         tot_nks    ! total number of distinct points: k,k+q,...
  ! note that it can be less than 4*orig_nks if some points are equal to each other,
  ! or equal to Gamma
  !
  INTEGER :: nbnd_max
  !
  TYPE kplus3q_type
    REAL(DP) :: xq(3)      ! coordinates of the q-point from input (alat units)
    REAL(DP) :: xq_drho(3) ! coordinates of the q-point of the ph.x calculation (alat units)
    REAL(DP),ALLOCATABLE :: xk(:,:) ! k points xk+xq
    REAL(DP),ALLOCATABLE :: wk(:) ! weights of the kpoints
    LOGICAL :: lgamma   ! true if q==Gamma
    LOGICAL :: lsame(-3:3) ! true if this q is equal to q_lsame(#)
    LOGICAL :: lminus(-3:3) ! true if this q is equal to q_lsame(#)
    !
    LOGICAL  :: lstored  ! not lgamma, not lsame and not lopposite
    INTEGER  :: copy_of  ! if stored is false, the point/grid which this one is a copy of
    !
    LOGICAL  :: ldrho         ! this point has it's own drho file from phonon, if ldrho is false:
    INTEGER  :: drho_from     ! get drho from this other q-point
    LOGICAL  :: ldrho_cc      ! this point has it's own drho file from phonon, if ldrho is false:
    INTEGER  :: drho_cc_from  ! drho can obtained from this other point by complex-conjugate
    LOGICAL  :: ldrho_is_mine ! .ldrho. and drho_from points to itself
    !
    INTEGER  :: npwq    ! number of plane waves, will be overwritten (FIXME: is this necessary?)
    INTEGER  :: nksq    ! number of k+q points
    INTEGER  :: iunigkq ! unit where wavefunctions are written
    INTEGER  :: koffset !+1 is the index of the first kpoint belonging to this k+q grid (FIXME?)
    INTEGER,ALLOCATABLE :: ikqs(:) ! map/list of k+q points in the global list of kpoints
    INTEGER,ALLOCATABLE :: igkq(:,:) ! list of plane waves for each k+q point
    INTEGER,ALLOCATABLE :: ngkq(:)
  END TYPE
  !
  TYPE(kplus3q_type) :: kplusq(-3:3) ! -q3,-q2,-q1,Gamma,q1,q2,q3
  !
  LOGICAL :: q_special_cases_initialized = .false.
  LOGICAL :: kplus3q_grids_initialized = .false.
  !
  ! Units and file name suffixes (the full name will be _ph0$prefix.$suffix) for storing
  ! the plane-wave mapping at the various k and k+/-q_x points. Files are actually opened
  ! in set_kplus3q, where the data to be stored is computed.
  INTEGER,PARAMETER :: kplus3q_units(-3:3) = (/ 4027, 4028, 4029, 4020, 4021, 4022, 4023 /)
  CHARACTER(len=6),PARAMETER :: kplus3q_suffixes(-3:3) &
             = (/ 'igkmq3', 'igkmq2', 'igkmq1', 'igkgam', 'igkq1 ', 'igkq2 ', 'igkq3 ' /)
  !
  CHARACTER(len=3),PARAMETER :: q_names(-3:3) &
             = (/ '-q3', '-q2', '-q1', 'G  ', 'q1 ', 'q2 ', 'q3 ' /)
  CHARACTER(len=3),PARAMETER :: q_names2(-3:3) &
             = (/ '-q3', '-q2', '-q1', '   ', '+q1', '+q2', '+q3' /)
  !
  ! Note: we prefer to check the qpoint in this order: q1, q2, q3, Gamma, -q3, -q2, -q1
  !       (this used to be a problem but now it's purely estetical)
  INTEGER,PARAMETER :: q_sorting(1:7) = (/ 1, 2, 3, 0, -1, -2, -3 /)
  !
  !
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE reset_kplus3q(cleanup)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: cleanup
  INTEGER :: iq
  ! Allocate everything at the beginning
  DO iq = -3,3
    DEALLOCATE(kplusq(iq)%ikqs)
    DEALLOCATE(kplusq(iq)%xk)
    DEALLOCATE(kplusq(iq)%wk)
    DEALLOCATE(kplusq(iq)%igkq)
    DEALLOCATE(kplusq(iq)%ngkq)
  ENDDO
  !
  q_special_cases_initialized = .false.
  kplus3q_grids_initialized   = .false.
  !
  CALL kplus3q_close_iunigkq(cleanup)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE reset_kplus3q
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE kplus3q_close_iunigkq(cleanup)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: cleanup
  INTEGER :: iq
  CHARACTER(len=6) :: status

  status = 'KEEP'
  IF(cleanup) status= 'DELETE'''

  DO iq = -3,3
    IF( kplusq(iq)%lstored ) &
      CLOSE( kplusq(iq)%iunigkq, status=status)
  ENDDO
  !-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
FUNCTION q_sum_rule(iq_wfc, iq_prt, valid)
  !-----------------------------------------------------------------------
  ! when deriving the wavefunctions at k+q_(iq_wf), with respect to a perturbation
  ! at q_(iq_prt) the resulting dpsi will have periodicity q_prj, which is computed here.
  ! In principles all the combinations are possile, yet the values that do not make
  ! sense in the 3rd order theory are forbidden here.
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iq_wfc, iq_prt
  LOGICAL,OPTIONAL,INTENT(out) :: valid
  INTEGER :: iq_prj
  INTEGER :: q_sum_rule
  !
#define X 99
  INTEGER,PARAMETER :: q_sum_rule_table(-3:3, -3:3) = &
        RESHAPE( (/ &
          X,  1,  2, -3,  X,  X,  0, &
          1,  X,  3, -2,  X,  0,  X, &
          2,  3,  X, -1,  0,  X,  X, &
          X,  X,  X,  X,  X,  X,  X, &
          X,  X,  0,  1,  X, -3, -2, &
          X,  0,  X,  2, -3,  X, -1, &
          0,  X,  X,  3, -2, -1,  X  &
         /), (/7,7/) )

  iq_prj = q_sum_rule_table(iq_wfc, iq_prt)

  IF(present(valid)) THEN
    IF(iq_prj == X) THEN
      q_sum_rule = X
      valid=.false.
      RETURN
    ELSE
      valid = .true.
    ENDIF
  ELSE
    IF(iq_prj == X) &
      CALL errore('q_sum_rule', 'Violation of q-points sum rule', 1)
    !
    ! this test shoudn't not be necessary, if I compiled the q_sum_rule table correctly
    IF(SUM((kplusq(iq_wfc)%xq+kplusq(iq_prt)%xq-kplusq(iq_prj)%xq)**2) > 1.d-8) &
      CALL errore('q_sum_rule', 'Violation of q-points sum rule', 2)
  ENDIF
  !
  q_sum_rule = iq_prj
#undef X
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END FUNCTION
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
FUNCTION q_sum_check(iq_wfc, iq_prt)
  !-----------------------------------------------------------------------
  ! this fct returns .true. if iq_wfc & iq_prt are a valid combination
  INTEGER,INTENT(in) :: iq_wfc, iq_prt
  LOGICAL :: q_sum_check
  INTEGER :: dummy
  !
  dummy = q_sum_rule(iq_wfc, iq_prt, q_sum_check)
  !
  RETURN
  !-----------------------------------------------------------------------
END FUNCTION
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE q_special_cases(xq1,xq2,xq3, at)
  !-----------------------------------------------------------------------
  !     This routine sets the k and k+q points (with zero weight) used in
  !     the preparatory run for a linear response 3rd order calculation.
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the xk grid is multipleid in order to include xk, xk+/-q1, 
  !                xk+/-q2 and xk+/-q3. If some q points are equal or are zero
  !                they are NOT included twice. Only the original wk of the original xk
  !                are left unchanged, the others are all zero. 
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : eps => eps12
  USE io_global,  ONLY : stdout
#define INVALID -99
  IMPLICIT NONE
  !
  ! I/O variable. FIXME: comments!
  !
  REAL(DP),INTENT(in)    :: xq1 (3), xq2(3), xq3(3)
  REAL(DP),INTENT(in)    :: at(3,3)
  !
  INTEGER :: i, j, iq, jq, pos
  CHARACTER(len=15),PARAMETER :: sub='q_special_cases'
  !LOGICAL,EXTERNAL :: eqvect_cart
  CHARACTER(len=128) :: aux
  LOGICAL :: print_this
  !
  ! Sanity check
  IF(q_special_cases_initialized) &
    CALL errore(sub, 'Already initialized!', 1)
  q_special_cases_initialized = .true.
  !
  ! prepare the kplusq objects to hold the k+/-q_x lists:
  kplusq(-3)%xq = -xq3
  kplusq(-2)%xq = -xq2
  kplusq(-1)%xq = -xq1
  kplusq( 0)%xq =  0._dp
  kplusq( 1)%xq =  xq1
  kplusq( 2)%xq =  xq2
  kplusq( 3)%xq =  xq3
  !
  kplusq(:)%lstored = .true.
  !
  DO i = 1,7
    iq = q_sorting(i)
    kplusq(iq)%copy_of = iq
    IF (iq>0) THEN ! per default, q1,q2,q3 have a drho file from phonon
      kplusq(iq)%ldrho = .true.
      kplusq(iq)%drho_from = iq
      kplusq(iq)%ldrho_cc = .false.
      kplusq(iq)%drho_cc_from = INVALID
    ELSE IF (iq<0) THEN ! while -q1,-q2,-q3 have to do a complex-conjugate
      kplusq(iq)%ldrho = .false.
      kplusq(iq)%drho_from = INVALID
      kplusq(iq)%ldrho_cc = .true.
      kplusq(iq)%drho_cc_from = -iq
    ELSE ! Gamma does not have, nor need, a drho file
      kplusq(iq)%ldrho = .false.
      kplusq(iq)%drho_from = INVALID
      kplusq(iq)%ldrho_cc = .false.
      kplusq(iq)%drho_cc_from = INVALID
    ENDIF
    !
  ENDDO
  !________________________________________________________________
  !
!#define __DEBUG_KPLUSQ
#ifdef  __DEBUG_KPLUSQ
  ! ...do nothing...
  WRITE(stdout, '(5x,a)') "WARNING! q_special_cases compiled for debugging!"
#else
  ! look for equalities between q's, e.g. q_i = q_j, or q_i = -q_j or q_i = G
  ! only unique q's will later be used to build a grid of k+/-q points
  ! priority is given according to q_sorting: 1,2,3,0,-1,-2,-3
  IQ_LOOP : &
  DO i = 1,7
    iq = q_sorting(i)
    !
    ! Check if it is Gamma, or equal to some other point:
    ! note: lgamma is redundant with %lsame(0)
    kplusq(iq)%lgamma = (SUM(kplusq(iq)%xq**2) < eps)
    !
    JQ_LOOP : &
    DO j = 1,7 !i
      jq = q_sorting(j)
      !
      kplusq(iq)%lsame(jq)  = ( SUM((kplusq(iq)%xq - kplusq(jq)%xq)**2) < eps )
      kplusq(iq)%lminus(jq) = ( SUM((kplusq(iq)%xq + kplusq(jq)%xq)**2) < eps ) &
                               .and. .not. kplusq(iq)%lsame(jq)
!     this looked like a good idea, but does not work:
!       kplusq(iq)%lsame(jq)  = eqvect_cart(kplusq(iq)%xq,  kplusq(jq)%xq, kplusq( 0)%xq, at)
!       kplusq(iq)%lminus(jq) = eqvect_cart(kplusq(iq)%xq, -kplusq(jq)%xq, kplusq( 0)%xq, at) &
!                                .and. .not. kplusq(iq)%lsame(jq)
      !
      ! If I have several equal q points, only one of them have to be prepared
      ! and stored, the others can just point to it. Here we decide
      ! if the current point will be the "stored" one or just a copy
      IF(kplusq(iq)%lsame(jq) .and. kplusq(jq)%lstored &
         .and. kplusq(iq)%lstored .and. i > j) THEN
         !
         kplusq(iq)%lstored = .false.
         kplusq(iq)%copy_of = kplusq(jq)%copy_of
      ENDIF
      !
      ! If we are in a 0,q,-q configuration (and permutation) we can get the rho at -q
      ! from the rho at q with a simple complex conjugate. If q_iq = -q_jq and iq>0 and
      ! jq>0 we can take advantage of this.
      IF(kplusq(iq)%lsame(jq) .and. kplusq(jq)%ldrho &
         .and. i > j .and. iq /= 0) THEN
        !
        kplusq(iq)%ldrho = .true.
        kplusq(iq)%drho_from = kplusq(jq)%drho_from
        !
        kplusq(iq)%ldrho_cc = .false.
        kplusq(iq)%drho_cc_from = INVALID
      ENDIF
      !
!       IF(kplusq(iq)%lminus(jq) .and. kplusq(jq)%ldrho .and. &
!          kplusq(iq)%ldrho .and. i > j)THEN
      IF(kplusq(iq)%lminus(jq) .and. kplusq(jq)%ldrho &
         .and. i > j .and. iq /= 0) THEN
        ! 
        kplusq(iq)%ldrho = .false.
        kplusq(iq)%drho_from = INVALID
        !
        kplusq(iq)%ldrho_cc = .true.
        kplusq(iq)%drho_cc_from = kplusq(jq)%drho_from
      ENDIF
      !
    ENDDO JQ_LOOP
    !
    IF( kplusq(iq)%ldrho .and. kplusq(iq)%drho_from<0) &
      CALL errore(sub, 'This should never happen!', 2)
  ENDDO IQ_LOOP
  !
#endif
  ! Report:
  DO i = 1,7
    iq= q_sorting(i)
    kplusq(iq)%ldrho_is_mine = kplusq(iq)%ldrho .and. ( kplusq(iq)%drho_from == iq) 
  ENDDO

  WRITE(stdout, '(5x,a)') "Special cases for qs:"
  DO iq = 1,3
    print_this = .false.
    pos = 0
    aux=''
    !
    DO jq = 0,3
      IF(kplusq(iq)%lsame(jq).and. iq/=jq) THEN
        aux(5*pos+1:5*pos+5)=q_names(jq)//", "
        pos=pos+1
        print_this = .true.
      ENDIF
      IF(kplusq(iq)%lsame(-jq) .and. iq/=jq .and. jq/=0) THEN
        aux(5*pos+1:5*pos+5)=q_names(-jq)//", "
        pos=pos+1
        print_this = .true.
      ENDIF
    ENDDO
    !
    IF(print_this)THEN
      WRITE(stdout, '(7x,3a)') TRIM(q_names(iq)), " = ", TRIM(aux(1:5*pos-2))
    ELSE
      WRITE(stdout, '(7x,2a)') TRIM(q_names(iq)), " : no special case detected"
    ENDIF
    !
  ENDDO
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE q_special_cases
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE kplus3q_grids (xk, wk, nks, max_nks, nat, kunit)
  !-----------------------------------------------------------------------
  USE io_files, ONLY : seqopn
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in)     :: max_nks, nat
  INTEGER,INTENT(inout)  :: nks
  REAL(DP),INTENT(inout) :: xk (3, max_nks), wk (max_nks)
  INTEGER,INTENT(out)    :: kunit
  !
  INTEGER :: ik, i, iq, jq
  LOGICAL :: exst
  CHARACTER(len=13),PARAMETER :: sub='kplus3q_grids'
  REAL(DP),ALLOCATABLE :: xk0(:,:), wk0(:) ! the original kpoint and weight from the scf calculation
  !
  IF(.not.q_special_cases_initialized) &
    CALL errore(sub, 'q_special_cases must be called before!', 1)
  !
  IF(kplus3q_grids_initialized) &
    CALL errore(sub, 'Already initialized!', 1)
  kplus3q_grids_initialized = .true.

  ! Store a copy of the original kpoints and weights
  ALLOCATE(xk0(3,nks), wk0(nks))
  xk0(:,1:nks) = xk(:,1:nks)
  wk0(1:nks)   = wk(1:nks)
  !
  orig_nks  = nks
  nksq      = nks
  tot_nks = 0
  !
  ! Allocate everything at the beginning
  DO i = 1,7
    iq = q_sorting(i)
    ! Note: the number of kpoints per q point and relative weights may change in the future!
    !       using these variables in a consistent way should make the transition painless
    ALLOCATE(kplusq(iq)%ikqs(nks))
    ALLOCATE(kplusq(iq)%xk(3,nks))
    ALLOCATE(kplusq(iq)%wk(nks))
    kplusq(iq)%nksq = nks
    DO ik = 1,kplusq(iq)%nksq
      kplusq(iq)%xk(:,ik) = xk0(:,ik) + kplusq(iq)%xq(:)
      kplusq(iq)%wk(ik) = wk0(ik)
    ENDDO
  ENDDO
  !
  kunit = 0
  !
  IQ_LOOP : &
  DO i = 1, 7
    iq = q_sorting(i)
    !
    ! Some special cases:
    ! if q1 is gamma then q2 is -q3 and q3 is -q2,
    !    ... and all the permutations
    ! if q1 and q2 are gamma then also q3 is gamma,
    !    ... and so on
    !
    ! If this point is stored, i.e. it is unique or the first of a list of equal ones,
    ! i prepare the grid and update a few quantities:
    !
    NEW_KPLUSQ_GRID : &
    IF(kplusq(iq)%lstored) THEN
      ! set the offset for this grid:
      kplusq(iq)%koffset = tot_nks
      !
      kunit=kunit+1
      !
      ! than increase the total number of kpoints:
      tot_nks = tot_nks + kplusq(iq)%nksq
      IF(tot_nks>max_nks) CALL errore(sub, &
                'too many k points: increase npk in Modules/parameters.f90', tot_nks)
      !
      ! set up the new k points
      ! ********************************** NOTE **********************************
      ! The k points in xk and the relative weights in wk are global variables that must
      ! be set to perform the NSCF calculation of the unperturbed wavefunctions at k, k+/-q1, ...
      ! For this reason we must have wk /= 0 ONLY for the original set of kpoints used for the
      ! SCF calculation! However, using these weights is extremely confusing: only the weights
      ! stored in kplusq(iq)%wk should be used inside the D3q code!!
      ! ********************************** **** **********************************
      DO ik = 1,orig_nks
        xk(:,ik+kplusq(iq)%koffset) = xk0(:,ik) + kplusq(iq)%xq(:)
        !
#ifdef __DEBUG_KPLUSQ
        IF(iq==0) THEN
#else
        IF(kplusq(iq)%lgamma) THEN
#endif
          wk(ik+kplusq(iq)%koffset) = wk0(ik)
          wk0(ik) = 0._dp
        ELSE
          wk(ik+kplusq(iq)%koffset) = 0._dp
        ENDIF
        ! set up the index of k+q point 
        ! (it is just sequential, but could be more general)
        kplusq(iq)%ikqs(ik) = ik+kplusq(iq)%koffset
        !
      ENDDO
      !
      ! open the file to store the plane-wave mappings
      kplusq(iq)%iunigkq = kplus3q_units(iq)
      CALL seqopn(kplusq(iq)%iunigkq, TRIM(kplus3q_suffixes(iq))//"_", 'unformatted', exst)
      !
    ELSE NEW_KPLUSQ_GRID
      !
      jq = kplusq(iq)%copy_of
      ! copy another k+q grid
      kplusq(iq)%koffset = kplusq(jq)%koffset
      kplusq(iq)%ikqs    = kplusq(jq)%ikqs
      !
      ! set the unit where the wavefunctions are stored
      kplusq(iq)%iunigkq = kplus3q_units(jq)
      !
    ENDIF &
    NEW_KPLUSQ_GRID
    !
  ENDDO IQ_LOOP
  !
#ifdef __DEBUG_KPLUSQ
  WRITE(stdout, '("WARNING! not using q equivalence tricks!")')
#endif
  lgamma = ALL(kplusq(:)%lgamma)
  nks = tot_nks
  !
  DEALLOCATE(xk0, wk0)
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE kplus3q_grids
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_igkq_d3()
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : npwx
  USE gvecw,      ONLY : ecutwfc
  USE gvect,      ONLY : g, ngm
  USE cell_base,  ONLY : tpiba2
  USE klist,      ONLY : ngk, igk_k, xk
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  INTEGER  :: iq, ikq_global
  REAL(DP),ALLOCATABLE :: g2(:)
  REAL(DP) :: gcutwfc
  CHARACTER(len=13),PARAMETER :: sub = 'write_igkq_d3'
  !
  IF(.not.kplus3q_grids_initialized) &
    CALL errore(sub, 'set_kplus3q_grids must be called before!', 1)
  !
  gcutwfc = ecutwfc / tpiba2
  !
  ! FIXME!!!! lsda is currently not supported in D3, but when it will...
!   IF ( lsda ) current_spin = isk( ikk )
  !
  !
  ALLOCATE (g2(npwx))
  DO iq = -3,3
    !
    !IF(kplusq(iq)%lstored) THEN
    
      ALLOCATE(kplusq(iq)%igkq(npwx,nksq))
      ALLOCATE(kplusq(iq)%ngkq(nksq))
      ! DEBUG :
      kplusq(iq)%igkq = -10000000-iq
      kplusq(iq)%ngkq = -100000000-iq
      
      DO ik = 1, nksq
        !
        ikq_global = kplusq(iq)%ikqs(ik)
        ! This line would recompute:
!         CALL gk_sort( kplusq(iq)%xk(:,ik), ngm, g, gcutwfc, kplusq(iq)%ngkq(ik), &
!                       kplusq(iq)%igkq(:,ik), g2 )
        ! Instead we copy from the NSCF ariable:
        kplusq(iq)%ngkq(ik)   = ngk(ikq_global)
        kplusq(iq)%igkq(:,ik) = igk_k(:,ikq_global)
        !
        WRITE( kplusq(iq)%iunigkq ) kplusq(iq)%ngkq(ik), kplusq(iq)%igkq(:,ik)
  !       ! DEBUG!! :
!         WRITE( 90000+kplusq(iq)%iunigkq,* )  kplusq(iq)%ngkq(ik), kplusq(iq)%igkq(:,ik)
      ENDDO
      FLUSH(kplusq(iq)%iunigkq)
    !ENDIF
  ENDDO
  !
  DEALLOCATE(g2)
  !
  RETURN
END SUBROUTINE write_igkq_d3
!
END MODULE kplus3q
