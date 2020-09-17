!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE d3_grid
!-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  !
  LOGICAL :: grid_setup_done = .false.
  !
  INTEGER :: n_triplets = 0 ! number of triplets to compute
  TYPE d3_triplet
    REAL(DP) :: xq1(3)
    REAL(DP) :: xq2(3)
    REAL(DP) :: xq3(3)
  END TYPE d3_triplet
  TYPE(d3_triplet),ALLOCATABLE :: d3_triplets(:)
  !
  INTEGER :: i_triplet_first  = -1
  INTEGER :: i_triplet_last   = -1
  INTEGER :: i_triplet_step   = 1
  INTEGER :: i_triplet_offset = 0
  INTEGER :: n_triplets_todo

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE d3_grid_slice(first, last, step, offset)
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(in) :: first, last, step, offset
  INTEGER :: i
  CHARACTER(len=13),PARAMETER :: sub='d3_grid_slice'
  !
  IF(.not. grid_setup_done) &
    CALL errore(sub, 'grid must be setup first', 1)
  !
  i_triplet_first = 1
  IF(first>0) i_triplet_first=first
  !
  ! (do not stop if last > n_triplets, as it may not be known in advance)
  i_triplet_last = n_triplets
  IF( last>n_triplets) THEN
    WRITE(stdout, '(7x,a,i6)') "Remark: less triplets to compute then requested 'last' reset to", n_triplets
  ELSE IF (last>0) THEN
    i_triplet_last=last
  ENDIF
  !
  IF(offset >= step) &
    CALL errore(sub, 'offset >= step: offset must be smaller then step', 1)
  i_triplet_step   = step
  i_triplet_offset = offset
  !
  IF (i_triplet_offset>i_triplet_last-i_triplet_first)THEN
    WRITE(stdout, '(5x,a)') "Remark: no triplet will be computed in this run (offset is too high)"
  ENDIF
  !
  n_triplets_todo = 1+(i_triplet_last-i_triplet_first-i_triplet_offset)/i_triplet_step
  !
  IF(n_triplets_todo == n_triplets) RETURN
  !
  WRITE(stdout, '(5x,a,i6,a)') "Triplets to be computed in this run (", &
                               n_triplets_todo," total):"
  !
  IF(n_triplets<10)THEN
    WRITE(stdout, '(7x,10i2)') (i, i= i_triplet_first+i_triplet_offset, i_triplet_last, i_triplet_step)
  ELSE IF(n_triplets<1000) THEN
    WRITE(stdout, '(7x,12i4)') (i, i= i_triplet_first+i_triplet_offset, i_triplet_last, i_triplet_step)
  ELSE
    WRITE(stdout, '(7x,10i6)') (i, i= i_triplet_first+i_triplet_offset, i_triplet_last, i_triplet_step)
  ENDIF
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_grid_slice
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_single_point_init(xq1, xq2, xq3)
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  IMPLICIT NONE
  REAL(DP), INTENT(in) :: xq1(3), xq2(3), xq3(3)
  !
  IF(grid_setup_done) CALL errore('d3_single_point_init', 'grid already set up', 1)
  grid_setup_done = .true.
  !
  n_triplets = 1
  ALLOCATE(d3_triplets(1))
  d3_triplets(1)%xq1 = xq1
  d3_triplets(1)%xq2 = xq2
  d3_triplets(1)%xq3 = xq3
!  WRITE(stdout,"(5x,a)") "D3 calculation of triplet:"
!  WRITE(stdout,'(7x,3("(",3f8.4," )",3x))') xq1, xq2, xq3
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3_single_point_init
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE d3_grid_init(nq1,nq2,nq3, xq_base)
  !-----------------------------------------------------------------------
  USE kinds,           ONLY : DP
  USE symm_base,       ONLY : nsym, s, time_reversal, t_rev, invs
  USE cell_base,       ONLY : at, bg
  USE d3_shuffle,      ONLY : nperms, d3perms_order
  USE qstar_d3_module, ONLY : star_3q
  USE io_global,       ONLY : stdout
  USE d3matrix_io2,    ONLY : d3matrix_filename2
  !
  IMPLICIT NONE
  INTEGER,INTENT(in) :: nq1,nq2,nq3
  REAL(DP),INTENT(in),OPTIONAL :: xq_base(3)
  !
  INTEGER :: nqx, nxq1, nxq2, nqxx
  REAL(DP),ALLOCATABLE :: xq1(:,:), xq2(:,:), wk(:)
  !
  INTEGER :: n3q
  REAL(DP),ALLOCATABLE :: x3q(:,:,:), sx3q(:,:,:,:), a3q(:,:,:), sa3q(:,:,:,:)
  INTEGER,ALLOCATABLE  :: is3q(:,:), im3q(:), nst3q(:)
  REAL(DP) ::d(9)
  !
  INTEGER :: iq1, iq2, a,b,c, i3q, iperm, isq, iq, ipol
  !
  !LOGICAL,EXTERNAL :: eqvect_cart, eqvect
  REAL(DP),EXTERNAL :: get_clock
  REAL(DP),PARAMETER :: gamma(3) = (/ 0._dp, 0._dp, 0._dp /), accep=1.e-5_dp
  !
  IF(grid_setup_done) CALL errore('d3_grid_init', 'grid already set up', 1)
  grid_setup_done = .true.
  !
  nqx = nq1*nq2*nq3
  ALLOCATE(xq2(3,nqx), wk(nqx))
  !
  IF(present(xq_base)) THEN
    ! Do a partial grid calculation with q1 fixed
    ALLOCATE(xq1(3,1))
    xq1 = 0._dp
    xq1(:,1) = xq_base
    nxq1 = 1
  ELSE
    ! Do a full-grid calculation: generate the irreducible points,
    ! they will be the space of q1
    ALLOCATE(xq1(3,nqx))
    CALL kpoint_grid ( nsym, time_reversal, .false., s, t_rev, bg, nqx, &
                      0,0,0, nq1,nq2,nq3, nxq1, xq1, wk)
  ENDIF
  !
  ! Generate the full (no-symmetry) grid, this contains the space of q2
  CALL kpoint_grid ( 1, .false., .true., s, t_rev, bg, nqx, &
                     0,0,0, nq1,nq2,nq3, nxq2, xq2, wk)
  !
  nqxx = nxq1*nxq2
  ALLOCATE( x3q(3,3,nqxx), sx3q(3,3,48,nqxx))
  ALLOCATE( a3q(3,3,nqxx), sa3q(3,3,48,nqxx))
  ALLOCATE( is3q(48,nqxx), im3q(nqxx), nst3q(nqxx))
  !
  WRITE(stdout,"(5x,a,i8,a)") "Looking for irreducible triplets out of ", nqxx, " possibilities:"
  !
  n3q = 0
  Q1_LOOP : &
  DO iq1 = 1, nxq1
    Q2_LOOP : &
    DO iq2 = 1,nxq2
      !
      n3q = n3q+1
      x3q(:,1,n3q) = xq1(:,iq1)
      x3q(:,2,n3q) = xq2(:,iq2)
      x3q(:,3,n3q) = -xq1(:,iq1) -xq2(:,iq2)
      !
      DO iq = 1,3     ! <- index of the q vector
        DO ipol = 1,3 ! <- index of the direction
          a3q(ipol,iq,n3q) =  x3q(1,iq,n3q) * at(1,ipol) &
                            + x3q(2,iq,n3q) * at(2,ipol) &
                            + x3q(3,iq,n3q) * at(3,ipol)
        ENDDO
      ENDDO
      ! Check if this triplet is a permutation of the star of a previous triplet
      ! NOTE: we don't care about equivalent permutations, this wastes a bit of
      !       time but makes the code way easier to read
      DO i3q = 1, n3q-1
      DO isq = 1, nst3q(i3q)
      DO iperm = 1, nperms
          !
          a = d3perms_order(1,iperm)
          b = d3perms_order(2,iperm)
          c = d3perms_order(3,iperm)
          !
          d =  (/ &
               sa3q(1,a,isq,i3q)-a3q(1,1,n3q), &
               sa3q(2,a,isq,i3q)-a3q(2,1,n3q), &
               sa3q(3,a,isq,i3q)-a3q(3,1,n3q), &
               sa3q(1,b,isq,i3q)-a3q(1,2,n3q), &
               sa3q(2,b,isq,i3q)-a3q(2,2,n3q), &
               sa3q(3,b,isq,i3q)-a3q(3,2,n3q), &
               sa3q(1,c,isq,i3q)-a3q(1,3,n3q), &
               sa3q(2,c,isq,i3q)-a3q(2,3,n3q), &
               sa3q(3,c,isq,i3q)-a3q(3,3,n3q) /)
          IF (& !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
              ALL( ABS(d-NINT(d) )<accep)&
              ) &
          THEN  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            n3q = n3q - 1
            CYCLE Q2_LOOP
          ELSE  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            IF(im3q(i3q) == 0) THEN
              d =  (/ &
                  sa3q(1,a,isq,i3q)+a3q(1,1,n3q), &
                  sa3q(2,a,isq,i3q)+a3q(2,1,n3q), &
                  sa3q(3,a,isq,i3q)+a3q(3,1,n3q), &
                  sa3q(1,b,isq,i3q)+a3q(1,2,n3q), &
                  sa3q(2,b,isq,i3q)+a3q(2,2,n3q), &
                  sa3q(3,b,isq,i3q)+a3q(3,2,n3q), &
                  sa3q(1,c,isq,i3q)+a3q(1,3,n3q), &
                  sa3q(2,c,isq,i3q)+a3q(2,3,n3q), &
                  sa3q(3,c,isq,i3q)+a3q(3,3,n3q) /)
              IF (& !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                  ALL( ABS(d-NINT(d) )<accep)&
                  ) &
              THEN  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                n3q = n3q - 1
                CYCLE Q2_LOOP
              ENDIF

            ENDIF
          ENDIF !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          !
      ENDDO
      ENDDO
      ENDDO
      !
      ! NOTE: because of the CYCLE Q2_LOOP the code only arrives here
      !       when a new point is found.
      !
      ! Compute the star of the new triplet
      WRITE(stdout,'(7x,i5,":",3x,3("(",3f8.4," )",3x),a)') n3q, x3q(:,:,n3q), &
         TRIM(d3matrix_filename2(x3q(:,1,n3q), x3q(:,2,n3q), x3q(:,3,n3q), at, '-->'))
      !
      CALL star_3q(x3q(:,1,n3q), x3q(:,2,n3q), x3q(:,3,n3q), &
                   at, bg, nsym, s, invs, nst3q(n3q),        &
                   sx3q(:,:,:,n3q), is3q(:,n3q), im3q(n3q) )
      !
      DO isq = 1, nst3q(n3q) ! <- index of the triplet in the star
        DO iq = 1,3          ! <- index of the q vector
          DO ipol = 1,3      ! <- index of the direction
            sa3q(ipol,iq,isq,n3q) =  sx3q(1,iq,isq,n3q) * at(1,ipol) &
                                   + sx3q(2,iq,isq,n3q) * at(2,ipol) &
                                   + sx3q(3,iq,isq,n3q) * at(3,ipol)
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO &
    Q2_LOOP
  ENDDO &
  Q1_LOOP
  !
  WRITE(stdout,"(5x,a,i4)") "IRREDUCIBLE NUMBER OF TRIPLETS TO COMPUTE:", n3q
  !
  n_triplets = n3q
  ALLOCATE(d3_triplets(n3q))
  DO i3q = 1, n3q
    d3_triplets(i3q)%xq1 = x3q(:,1,i3q)
    d3_triplets(i3q)%xq2 = x3q(:,2,i3q)
    d3_triplets(i3q)%xq3 = x3q(:,3,i3q)
  ENDDO
  !
  DEALLOCATE(xq1, xq2, wk)
  DEALLOCATE(x3q, sx3q)
  DEALLOCATE(a3q, sa3q)
  DEALLOCATE(is3q, im3q, nst3q)
  !
  RETURN
  !
  !-----------------------------------------------------------------------
END SUBROUTINE d3_grid_init
!-----------------------------------------------------------------------
! !-----------------------------------------------------------------------
! SUBROUTINE kpoint_grids_difference(n1, xq1, n2, xq2)
!   !-----------------------------------------------------------------------
!   USE kinds, ONLY : DP
!   INTEGER,INTENT(in)     :: n1
!   REAL(DP),INTENT(in)    :: xq1(3,n1)
!   INTEGER,INTENT(inout)  :: n2
!   REAL(DP),INTENT(inout) :: xq2(3,n2)
!   !
!   REAL(DP),ALLOCATABLE :: xp(:,:)
!   INTEGER :: np
!   !
!   np = 
!   IF(n1 > n2) CALL errore
!   ALLOCATE(xp(3,n2-n1))
!   
!   !-----------------------------------------------------------------------
! END SUBROUTINE kpoint_grids_difference
! !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
END MODULE d3_grid
!-----------------------------------------------------------------------
