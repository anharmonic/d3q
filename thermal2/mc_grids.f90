!
! Written by Lorenzo Paulatto (2017-2018) IMPMC @ Universite' Sorbonne / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE mc_grids
  !
  USE q_grids
#include "mpi_thermal.h"

  REAL(DP),PRIVATE :: avg_npoints = 0._dp, avg_tested_npoints = 0._dp, saved_threshold
  INTEGER,PRIVATE  :: ngrids_optimized

  CONTAINS
  !
  SUBROUTINE setup_mcjdos_grid(input, S, fc, grid, xq0, nq_target, scatter)
    USE code_input,       ONLY : code_input_type
    USE input_fc,         ONLY : ph_system_info, forceconst2_grid
    USE ph_dos,           ONLY : joint_dos_q
    USE random_numbers,   ONLY : randy
    USE constants,        ONLY : RY_TO_CMM1
    USE mpi_thermal
    IMPLICIT NONE
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    !CHARACTER(len=*),INTENT(in)     :: grid_type
    !REAL(DP),INTENT(in)   :: bg(3,3) ! = System
    !INTEGER,INTENT(in) :: n1,n2,n3
    INTEGER,INTENT(in) :: nq_target
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),OPTIONAl,INTENT(in) :: xq0(3)
    LOGICAL,OPTIONAL,INTENT(in) :: scatter
    !
    REAL(DP) :: sigma_ry, avg_T
    !
    REAL(DP) :: xq_new(3), jdq_new, jdq_old
    INTEGER :: nq, iq
    LOGICAL :: accepted
    REAL(DP) :: acceptance, test
    TYPE(q_grid) :: grid0
    INTEGER :: warm_up

    !nq_target = 125
    warm_up = MIN(MAX(nq_target/10, 20), 100)
    nq = 1-warm_up
    jdq_old = -1._dp
    grid%type = 'mcjdos'

    sigma_ry = SUM(input%sigma)/RY_TO_CMM1/input%nconf
    avg_T = SUM(input%T)/input%nconf
    !print*, sigma_ry, avg_T

    ALLOCATE(grid%xq(3,nq_target))
    ALLOCATE(grid%w(nq_target))
    grid%scattered = .false.
    grid%nq    = nq_target
    grid%nqtot = nq_target

    CALL setup_grid("simple", S%bg, input%nk(1),input%nk(2),input%nk(3), &
                grid0, scatter=.true.)

    POPULATE_GRID : &
    DO
      xq_new = 0._dp
      IF(ionode)THEN
        IF(input%nk(1)>1) xq_new(1) = randy()/DBLE(input%nk(1))
        IF(input%nk(2)>1) xq_new(2) = randy()/DBLE(input%nk(2))
        IF(input%nk(3)>1) xq_new(3) = randy()/DBLE(input%nk(3))
        CALL cryst_to_cart(1,xq_new,S%bg, +1)
        test = randy()
      ENDIF
      CALL mpi_bcast_vec(3,xq_new)
      CALL mpi_bcast_scl(test)
      !
      jdq_new = 1._dp !joint_dos_q(grid0,sigma_ry,avg_T, S, fc, xq_new)
      !
      IF(jdq_old<=0._dp)THEN
        accepted = .true.
      ELSE
        acceptance = jdq_new/jdq_old
        accepted = (test <= acceptance)
      ENDIF
      !
      IF(accepted)THEN
        nq = nq+1
        jdq_old = jdq_new
        IF(nq>0)THEN
          grid%w(nq) = 1._dp
          grid%xq(:,nq) = xq_new
        ENDIF
        !ioWRITE(*,*) "accepted", nq, acceptance, jdq_old
      ELSE
        !ioWRITE(*,*) "discarded"
        IF(nq>0) grid%w(nq) = grid%w(nq)+1
      ENDIF
      !
      IF(nq==nq_target) EXIT POPULATE_GRID
      !
    ENDDO &
    POPULATE_GRID

    grid%w = 1._dp/grid%w
    grid%w = grid%w/SUM(grid%w)

    ioWRITE(7778,'(2x,"Setup a ",a," grid of",i9," q-points")') "mcjdos", grid%nqtot
    DO iq = 1,grid%nqtot
      ioWRITE(7778,'(3f12.6,f12.6)') grid%xq(:,iq), grid%w(iq)*grid%nqtot
    ENDDO
    IF(scatter) CALL grid%scatter()
    !
  END SUBROUTINE setup_mcjdos_grid

  ! As optimize grid, but tries to do everything in parallel to avoid memory bottleneck
  ! (currently sort is done in serial)
  SUBROUTINE setup_poptimized_grid(input, S, fc, grid, xq0, prec, scatter, fc3)
    USE code_input,       ONLY : code_input_type
    USE input_fc,         ONLY : ph_system_info, forceconst2_grid
    USE ph_dos,           ONLY : joint_dos_q
    USE random_numbers,   ONLY : randy
    USE constants,        ONLY : RY_TO_CMM1
    USE fc2_interpolate,  ONLY : freq_phq_safe, set_nu0, bose_phq
    USE fc3_interpolate,  ONLY : forceconst3
    USE functions,        ONLY : quicksort_idx
    USE linewidth,        ONLY : sum_linewidth_modes
    USE mpi_thermal
    USE timers
    IMPLICIT NONE
    TYPE(code_input_type),INTENT(in)  :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(inout) :: grid
    REAL(DP),INTENT(in) :: xq0(3)
    REAL(DP),INTENT(in) :: prec
    LOGICAL,INTENT(in) :: scatter
    CLASS(forceconst3),INTENT(in),OPTIONAL :: fc3
    !
    !
    TYPE(q_grid) :: grid0
    INTEGER :: i, iq, jq, nu0(3)
    INTEGER :: iband = 1
    REAL(DP) :: xq(3,3), totfklw(S%nat3), normfklw(S%nat3), partialfklw, targetfklw
    REAL(DP),ALLOCATABLE :: V3sq(:,:,:), contributions(:,:), contributions_sum(:), &
                            freq(:,:), bose(:,:), &
                            fklw(:), contributions_tot(:), xq_tot(:,:), w_tot(:), &
                             xq_sort(:,:), w_sort(:)
    COMPLEX(DP),ALLOCATABLE :: U(:,:,:), D3(:,:,:)
    INTEGER,ALLOCATABLE  :: idx(:)
    !
    LOGICAL,SAVE  :: first=.true.
    REAL(DP),SAVE :: saved_xq0(3) = 0._dp

    CALL t_optimize%start()

    grid%type = 'optimized'
    grid%scattered = .false.
    !
    ! Hack to always use the same random shift:
    IF(first)THEN
      first = .false.
      CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), &
                      grid0, xq0=input%xk0, scatter=.true., quiet=.true.)
      saved_xq0 = grid0%xq0
    ELSE
      CALL setup_grid("simple", S%bg, input%nk(1), input%nk(2), input%nk(3), &
                      grid0, xq0=saved_xq0, scatter=.true., quiet=.true.)
    ENDIF
    !ioWRITE(*,*) grid0%type, grid0%xq0
    !CALL setup_grid("simple", S%bg, input%nk(1),input%nk(2),input%nk(3), &
    !            grid0, scatter=.false.)

    ALLOCATE(contributions(S%nat3,grid0%nq))
    ALLOCATE(V3sq(S%nat3, S%nat3, S%nat3))
    ALLOCATE(U(S%nat3, S%nat3, 3))
    ALLOCATE(freq(S%nat3, 3))
    ALLOCATE(bose(S%nat3, 3))
    ALLOCATE(fklw(S%nat3))
    
    xq(:,1) = xq0
    
!     IF(present(fc3))THEN
!       ALLOCATE(D3(S%nat3, S%nat3, S%nat3))
!       xq(:,2) = -xq0/2
!       xq(:,3) = -xq0/2
!      CALL fc3%interpolate(xq(:,2), xq(:,3), S%nat3, D3)
!       V3sq = REAL( CONJG(D3)*D3 , kind=DP)
!       DEALLOCATE(D3)
!     ELSE
      V3sq = 1._dp
!     ENDIF
    
    nu0(1)  = set_nu0(xq(:,1), S%at)
    CALL freq_phq_safe(xq(:,1), S, fc, freq(:,1), U(:,:,1))
    CALL bose_phq(input%T(1),S%nat3, freq(:,1), bose(:,1))

    totfklw = 0._dp
    DO iq = 1, grid0%nq

      xq(:,2) = grid0%xq(:,iq)
      xq(:,3) = -(xq(:,2)+xq(:,1))
      DO jq = 2,3
        nu0(jq) = set_nu0(xq(:,jq), S%at)
        CALL freq_phq_safe(xq(:,jq), S, fc, freq(:,jq), U(:,:,jq))
        CALL bose_phq(MAXVAL(input%T),S%nat3, freq(:,jq), bose(:,jq))
      ENDDO

      fklw = sum_linewidth_modes(S, MAXVAL(input%sigma)/RY_TO_CMM1, freq, bose, V3sq, nu0)
      contributions(:,iq) = ABS(fklw)
      totfklw(:) =  totfklw(:) + contributions(:,iq)
    ENDDO
    DEALLOCATE(V3sq, U, freq, bose, fklw)
    !
    CALL mpi_bsum(S%nat3, totfklw)
    !
    ! We normalize the contribution of each band, otherwise the bands with a large joint dos
    ! dominate and can compromise the precision of the other bands
    ALLOCATE(contributions_sum(grid0%nq))
    WHERE(totfklw/=0._dp) normfklw = 1._dp/(totfklw)
    DO iq = 1,grid0%nq
      contributions(:,iq) = contributions(:,iq)*normfklw(:)
    ENDDO
    DO iq = 1, grid0%nq
      contributions_sum(iq) = contributions(iband,iq) !SUM(contributions(:,iq))
    ENDDO
    totfklw = totfklw*normfklw
    DEALLOCATE(contributions)
    !
    CALL gather_vec(grid0%nq, contributions_sum, contributions_tot)
    !
    DEALLOCATE(contributions_sum)
    !
    IF(ionode)THEN
      ALLOCATE(idx(grid0%nqtot))
      FORALL(iq=1:grid0%nqtot) idx(iq) = iq
      contributions_tot = -contributions_tot
      CALL quicksort_idx(contributions_tot, idx, 1, grid0%nqtot)
      contributions_tot = -contributions_tot
      !
      targetfklw = (1._dp-prec)*totfklw(iband)
      IF(targetfklw==0._dp)THEN
       grid%nqtot = 0
      ELSE
        partialfklw = 0._dp
        DO iq = 1,grid0%nqtot
            partialfklw = partialfklw + contributions_tot(iq)
            !print*, iq, contributions_tot(iq), partialfklw, targetfklw
            IF(partialfklw>=targetfklw) EXIT
        ENDDO
        IF(iq>grid0%nqtot) iq = grid0%nqtot !CALL errore("poptimize","wrong iq", 1)
        grid%nqtot = iq
      ENDIF
      
      !targetfklw = prec * totfklw
!       DO iq = 1,grid0%nqtot
!         partialfklw = partialfklw + contributions_tot(iq)
!         IF(partialfklw>prec) EXIT
!       ENDDO
!       grid%nqtot = grid0%nqtot-iq+1
    ENDIF
    !
    DEALLOCATE(contributions_tot)
    !
    CALL mpi_bcast_integer(grid%nqtot)
    !
    CALL gather_mat(3, grid0%nq, grid0%xq, xq_tot)
    CALL gather_vec(grid0%nq, grid0%w, w_tot)
    !
    IF(ionode)THEN
      ALLOCATE(xq_sort(3,grid%nqtot))
      ALLOCATE(w_sort(grid%nqtot))
      ! set point from less important to more, to reduce roundoff errors
      DO iq = 1,grid%nqtot
        xq_sort(:,iq) = xq_tot(:,idx(iq))
        w_sort(iq) = w_tot(idx(iq))
      ENDDO
    ELSE
      ALLOCATE(xq_sort(0,0))
      ALLOCATE(w_sort(0))
    ENDIF
    !
    DEALLOCATE(xq_tot, w_tot)
    !
    CALL scatter_vec(grid%nqtot, w_sort, grid%nq, grid%w, grid%iq0)
    CALL scatter_mat(3, grid%nqtot, xq_sort, grid%nq, grid%xq)
    grid%scattered = .true.
    !
    DEALLOCATE(xq_sort, w_sort)
    ! 
        avg_npoints = (avg_npoints*ngrids_optimized + grid%nqtot)&
                  /DBLE(ngrids_optimized+1)
    avg_tested_npoints = (avg_tested_npoints*ngrids_optimized + grid0%nqtot)&
                         /DBLE(ngrids_optimized+1)
    saved_threshold = input%optimize_grid_thr
    ngrids_optimized = ngrids_optimized+1
    !
    ioWRITE(*,'(a,f8.2,a,i8)') "Optimized grid:", &
         DBLE(grid%nqtot)/grid0%nqtot*100, "% ", grid%nqtot
    CALL grid0%destroy()

    CALL t_optimize%stop()
    !
    !
  END SUBROUTINE setup_poptimized_grid

  SUBROUTINE print_optimized_stats()
    USE nanoclock, ONLY : get_wall
    USE timers
    IMPLICIT NONE
    REAL(DP) :: wall_time, my_time, acceptance, estimated_time
    IF(avg_tested_npoints<1._dp) RETURN

    wall_time = get_wall()
    my_time = t_optimize%read()

    acceptance = avg_npoints/avg_tested_npoints
    ! estimated time required if we were not using optimize_grid:
    estimated_time = (wall_time-my_time)/acceptance

    ioWRITE(*,'(a)') "*** * Grid optimization statistics:"
    ioWRITE(*,'(2x," * ",a27," * ",f12.0," / ",f12.0," * ",a21," = ",f6.2," * ", a18, " = ", f6.2" *")') &
      "avg used/initial points:", avg_npoints, avg_tested_npoints, &
      "avg acceptance (%))", 100*acceptance, &
      "speedup (est.)", estimated_time/wall_time
  END SUBROUTINE

  
END MODULE mc_grids
