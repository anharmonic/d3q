!
! This module is rewritten from the tetra.f90 in PW/src
!
MODULE thtetra
   !
   ! Tetrahedron method, linear and optimized. opt is better for all purpose
   ! weights for delta integration are obtained by multiplying gi with Iik in
   ! https://iopscience.iop.org/article/10.1088/0022-3719/12/15/008
   !
   ! another useful link, for delta integration where the integrand is 1 (DOS)
   ! http://staff.ustc.edu.cn/~zqj/posts/LinearTetrahedronMethod/#fn:tet_weight
   !
   ! optimized tetrahedron method, used in QE for theta integration in opt_tetra_weights
   ! https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.094515
   ! they multiply ni with Jik in the first article, then they transform (fit) through wlsm matrices
   USE kinds, ONLY: DP
   USE mpi_thermal, ONLY: my_id, num_procs, mpi_bsum
   !
   IMPLICIT NONE
   !
   PRIVATE
   SAVE
   !
   INTEGER :: ntetra = 0
   !! number of tetrahedra
   INTEGER :: nntetra = 0
   !! k-points per tetrahedron used to compute weights.
   !! 4 for linear / 20 for optimized tetrahedron method
   INTEGER, ALLOCATABLE :: tetra(:,:)
   !! index of k-points in a given tetrahedron shape (nntetra,ntetra)
   REAL(DP), ALLOCATABLE :: wlsm(:,:)
   !! Weights for the optimized tetrahedron method
   INTEGER :: nqs = 0
   !! number of q-points
   INTEGER :: nbnd = 0
   !! number of bands
   INTEGER, ALLOCATABLE :: itetra(:,:,:)
   !! order index of vertices of each tetrahedron (4,nbnd,ntetra)
   REAL(DP), ALLOCATABLE :: ek_sort(:,:,:)
   !! sorted energies for each tetrahedron (4,nbnd,ntetra)
   ! INTEGER, allocatable :: which_tetra(:,:,:)
   !! inverse of tetra: given a q point, it gives all the tetrahedra that contain it

   REAL(DP), PARAMETER :: tet_cutoff = 1.0E-4_DP
   LOGICAL :: opt_flag
   !
   PUBLIC :: tetra, ntetra, nntetra, wlsm, tetra_weights_green
   PUBLIC :: tetra_init, deallocate_tetra, tetra_weights_delta
   PUBLIC :: tetra_weights_delta_vec, ek_sort
   ! PUBLIC :: tetra_weights_delta, tetra_weights_theta, tetra_delta

   EXTERNAL :: errore, hpsort

   !
CONTAINS
   !
   !--------------------------------------------------------------------------
   SUBROUTINE tetra_init(nq, bg, ek)
      !-----------------------------------------------------------------------------
      !! This rouotine sets the corners and additional points for each tetrahedron.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nq(3)
      !! number of q-points in each direction
      !
      REAL(DP), INTENT(IN) :: bg(3,3)
      !! Reciplocal lattice vectors [2 pi / a]
      REAL(DP), INTENT(IN) :: ek(:,:)
      !! energy in the form ek(ibnd, iq)
      ! LOGICAL, INTENT(IN) :: is_mpi
      !! if .true., the grid is scattered
      ! LOGICAL, INTENT(IN), OPTIONAL :: opt
      ! !! if .true., uses opt_tetra methods

      REAL(DP), PARAMETER :: eps = 1e-5_dp
      !
      INTEGER :: i1, i2, i3, itet, itettot, ii, ik,  &
         ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), ibnd
      ! integer :: tetra_ik(nq(1) * nq(2) * nq(3))
      !
      REAL(DP) :: l(4), bvec2(3,3), bvec3(3,4)

      IF(ntetra /= 0) CALL deallocate_tetra()
      !
      nbnd = SIZE(ek,1)
      nqs = nq(1) * nq(2) * nq(3)
      !
      IF(nqs /= SIZE(ek,2)) CALL errore("tetra_init", "n(1) * n(2) * n(3) /= SIZE(ek,2)", SIZE(ek,2))
      ntetra  = 6*nqs
      !
      ALLOCATE(ek_sort(4,nbnd,ntetra))
      ek_sort = 0.0_dp
      ALLOCATE(itetra (4,nbnd,ntetra))
      ! ALLOCATE(which_tetra (2,nqs,24))
      !
      opt_flag = .true.
      ! if(PRESENT(opt)) opt_flag = opt
      !
      ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
      !
      bvec2(1:3,1) = bg(1:3,1) / REAL(nq(1), dp)
      bvec2(1:3,2) = bg(1:3,2) / REAL(nq(2), dp)
      bvec2(1:3,3) = bg(1:3,3) / REAL(nq(3), dp)
      !
      bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
      bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
      bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
      !
      DO ii = 1, 4
         l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
      ENDDO
      !
      ii = MINLOC(l(1:4),1)
      !
      ivvec0(1:4) = (/ 0, 0, 0, 0 /)
      !
      divvec(1:4,1) = (/ 1, 0, 0, 0 /)
      divvec(1:4,2) = (/ 0, 1, 0, 0 /)
      divvec(1:4,3) = (/ 0, 0, 1, 0 /)
      divvec(1:4,4) = (/ 0, 0, 0, 1 /)
      !
      ivvec0(ii) = 1
      divvec(ii, ii) = - 1
      !
      ! Divide a subcell into 6 tetrahedra
      !
      itet = 0
      DO i1 = 1, 3
         DO i2 = 1, 3
            IF(i2 == i1) CYCLE
            DO i3 = 1, 3
               IF(i3 == i1 .OR. i3 == i2) CYCLE
               !
               itet = itet + 1
               !
               ivvec(1:3,1,itet) = ivvec0(1:3)
               ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
               ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
               ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
               !
            ENDDO
         ENDDO
      ENDDO
      !
      ! Additional points surrounding the tetrahedron
      !
      ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
      ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
      ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
      !
      ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
      ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
      ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
      !
      ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
      ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
      ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
      ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
      !
      ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
      ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
      ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
      ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
      !
      ! Set the weight for the each tetrahedron method
      !
      ! WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
      !
      IF (opt_flag) THEN
         !
         nntetra = 20
         IF (.NOT. ALLOCATED(tetra)) ALLOCATE( tetra(nntetra,ntetra) )
         IF (.NOT. ALLOCATED(wlsm))  ALLOCATE( wlsm(4,nntetra) )
         !
         wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
         wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
         wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
         wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
         !
         wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
         wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
         wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
         wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
         !
         wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
         wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
         wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
         wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
         !
         wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
         wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
         wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
         wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
         !
         wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
         wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
         wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
         wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
         !
         wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260.0_dp
         !
      ELSE
         !
         nntetra = 4
         IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
         IF(.NOT. ALLOCATED(wlsm))  ALLOCATE ( wlsm(4,nntetra) )
         wlsm(:,:) = 0.0_dp
         !
         wlsm(1,1) = 1.0_dp
         wlsm(2,2) = 1.0_dp
         wlsm(3,3) = 1.0_dp
         wlsm(4,4) = 1.0_dp
         !
      ENDIF

      !
      !  locate k-points of the uniform grid in the list of irreducible k-points
      !  that was previously calculated
      !
      !  bring irreducible k-points to crystal axis
      !
      ! Construct tetrahedra
      !
      itettot = 0
      ! tetra_ik = 0
      DO i1 = 1, nq(1)
         DO i2 = 1, nq(2)
            DO i3 = 1, nq(3)
               !
               DO itet = 1, 6
                  !
                  itettot = itettot + 1
                  DO ibnd = 1, nbnd
                     !
                     DO ii = 1, nntetra
                        !
                        ikv(1:3) = (/i1, i2, i3/) - 1
                        ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
                        ikv(1:3) = MODULO(ikv(1:3), (/nq(1), nq(2), nq(3)/))
                        !
                        ik = ikv(3) + nq(3) * (ikv(2) + nq(2) * ikv(1)) + 1
                        !
                        ! if(ii <= 4) then
                        !    tetra_ik(ik) = tetra_ik(ik) + 1
                        !    which_tetra(1, ik, tetra_ik(ik)) = itettot
                        !    which_tetra(2, ik, tetra_ik(ik)) = ii
                        ! endif
                        tetra(ii, itettot) = ik
                        !
                        ! IF(opt_flag) THEN
                           ek_sort(:,ibnd,itettot) = ek_sort(:,ibnd,itettot) + wlsm(:,ii) * ek(ibnd,ik)
                        ! ELSE
                        !    ek_sort(ii,ibnd,itettot) = ek(ibnd,ik)
                        ! ENDIF

                     END DO ! ii
                     !
                     !
                     itetra(1,ibnd,itettot) = 0 ! needed to initialize index inside hpsort
                     CALL hpsort( 4, ek_sort(:,ibnd,itettot), itetra(:,ibnd,itettot))
                  ENDDO ! ibnd
                  !
               ENDDO ! itet
               !
            ENDDO ! i3
         ENDDO ! i2
      ENDDO ! i1
      !
   END SUBROUTINE tetra_init
   !
   !--------------------------------------------------------------------------------------
   ! FUNCTION tetra_delta( nqs, nbnd, et, ef) RESULT(w)
   !   !-----------------------------------------------------------------------------------
   !   !! Calculate the area of the implicit surface for an integrals of the kind int(Ak delta(ef-ek))
   !   !! usage: create a grid of points around a point in which you know Ak, this function gives the area
   !   !! and it's the same as sum(tetra_weights_delta) over the grid
   !   !-----------------------------------------------------------------------------------
   !   INTEGER, INTENT(IN) :: nqs
   !   !! The total # of k in irr-BZ
   !   INTEGER, INTENT(IN) :: nbnd
   !   !! The # of bands
   !   REAL(DP), INTENT(IN) :: et(nbnd,nqs)
   !   !! Kohn Sham energy [Ry]
   !   REAL(DP) :: w(nbnd)
   !   !! Integration weight of each k
   !   REAL(DP), INTENT(IN) :: ef
   !   !! The Fermi energy
   !   !
   !   ! ... local variables
   !   !
   !   INTEGER :: ik, nt, ibnd, i, ii, itetra(4), my_id_, num_procs_
   !   REAL(DP) :: e(4), C, a(4,4)

   !   EXTERNAL hpsort
   !   !
   !   w = 0._dp
   !   !
   !   IF(is_mpi_flag) THEN
   !     my_id_ = my_id
   !     num_procs_ = num_procs
   !   ELSE
   !     my_id_ = 0
   !     num_procs_ = 1
   !   ENDIF
   !   DO nt = 1+my_id_, ntetra, num_procs_
   !     !
   !     DO ibnd = 1, nbnd
   !       !
   !       e(1:4) = 0.0_dp
   !       DO ii = 1, nntetra
   !         !
   !         ik = tetra(ii, nt)
   !         IF(opt_flag) THEN
   !           e(1:4) = e(1:4) + wlsm(1:4,ii) * et(ibnd,ik)
   !         ELSE
   !           e(ii) = et(ibnd,ik)
   !         ENDIF
   !         !
   !       ENDDO
   !       !
   !       IF(MAXVAL(e) < ef .or. MINVAL(e) > ef) CYCLE
   !       itetra(1) = 0
   !       CALL hpsort( 4, e, itetra )
   !       ! CALL quicksort_idx(e, itetra, 1, 4)
   !       DO ii = 1, 4
   !         DO i = 1, 4
   !           IF ( ABS(e(i)-e(ii)) < 1.d-12 ) THEN
   !             a(ii,i) = 0.0_dp
   !           ELSE
   !             a(ii,i) = ( ef - e(i) ) / (e(ii) - e(i) )
   !           END IF
   !         ENDDO
   !       ENDDO
   !       !
   !       IF( e(1) < ef .AND. ef < e(2) ) THEN
   !         !
   !         C = 3 * a(2,1) * a(3,1) * a(4,1) / (ef - e(1))
   !         !
   !       ELSEIF( e(2) <= ef .AND. ef < e(3)) THEN
   !         !
   !         C = 3 * ( a(2,3) * a(3,1) + a(3,2) * a(2,4) ) / (e(4) - e(1))
   !         !
   !       ELSEIF ( e(3) <= ef .AND. ef < e(4)) THEN
   !         !
   !         C = 3 * a(1,4) * a(2,4) * a(3,4) / (e(4) - ef)
   !         !
   !       ENDIF
   !       !
   !       w(ibnd) = w(ibnd) + C
   !       !
   !     ENDDO ! ibnd
   !     !
   !   ENDDO ! nt
   !   w = w / ntetra
   !   !
   !   ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
   ! END FUNCTION tetra_delta
   ! !--------------------------------------------------------------------------------------
   FUNCTION tetra_weights_delta_vec(nef, ef) !RESULT(wI)
      USE constants, ONLY : pi
      !-----------------------------------------------------------------------------------
      !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
      !! The resulting wg can be used as sum(Ak * wk)
      !-----------------------------------------------------------------------------------
      integer, intent(in) :: nef
      !! size of vector ef
      REAL(DP), INTENT(IN) :: ef(:)
      !! The Fermi energy
      REAL(DP) :: tetra_weights_delta_vec(nef, nbnd, nqs)
      REAL(DP),ALLOCATABLE :: wI(:,:,:)
      !! COMPLEX Integration weight of each k
      !
      !
      INTEGER :: ik, nt, ibnd, ii, ief
      REAL(DP) :: e(4), wI0(4)

      ! for real part calc
      ! REAL(DP) :: wR0(4), ef_e(4), log_ef_e(4), prod_a(4), sum_a(4), second_term(4)
      ! INTEGER :: i3, j3
      !
      ALLOCATE(wI(nef, nbnd, nqs))
      wI = 0._dp
      !
      DO nt = 1+my_id, ntetra, num_procs
         !
         DO ibnd = 1, nbnd
            DO ief = 1, SIZE(ef)
               !
               e = ek_sort(:,ibnd,nt)
               wI0 = delta_vertices(ef(ief), e)
               !
               DO ii = 1, nntetra
                  !
                  ik = tetra(ii, nt)
                  ! IF(opt_flag) THEN
                  wI(ief,ibnd,ik) = wI(ief,ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), wI0(:))
                  ! ELSE
                  !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
                  !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
                  !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
                  ! ENDIF
               ENDDO
               !
            ENDDO ! ief
         ENDDO ! ibnd
         !
      ENDDO ! nt
      ! wg = wg / REAL(ntetra, dp)
      wI = wI / (6.0_dp * nqs)
      !
      ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
      CALL mpi_bsum(nef, nbnd, nqs, wI)
      tetra_weights_delta_vec = wI
      DEALLOCATE(wI)
   END FUNCTION
   !
   FUNCTION tetra_weights_delta(ef) RESULT(wI)
      USE constants, ONLY : pi
      !-----------------------------------------------------------------------------------
      !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
      !! The resulting wg can be used as sum(Ak * wk)
      !-----------------------------------------------------------------------------------
      REAL(DP) :: wI(nbnd, nqs)
      !! COMPLEX Integration weight of each k
      REAL(DP), INTENT(IN) :: ef
      !! The Fermi energy
      INTEGER :: ik, nt, ibnd, ii
      REAL(DP) :: e(4), wI0(4)

      ! for real part calc
      ! REAL(DP) :: wR0(4), ef_e(4), log_ef_e(4), prod_a(4), sum_a(4), second_term(4)
      ! INTEGER :: i3, j3
      !
      wI = 0._dp
      !
      DO nt = 1+my_id, ntetra, num_procs
         !
         DO ibnd = 1, nbnd
            !
            e = ek_sort(:,ibnd,nt)
            wI0 = delta_vertices(ef, e)
            !
            DO ii = 1, nntetra
               !
               ik = tetra(ii, nt)
               ! IF(opt_flag) THEN
               wI(ibnd,ik) = wI(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), wI0(:))
               ! ELSE
               !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
               !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
               !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
               ! ENDIF
            ENDDO
            !
         ENDDO ! ibnd
         !
      ENDDO ! nt
      ! wg = wg / REAL(ntetra, dp)
      wI = wI / (6.0_dp * nqs)
      !
      ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
      CALL mpi_bsum(nbnd, nqs, wI)
   END FUNCTION

   subroutine rm_degen_vertices(hw, D)
      real(dp), INTENT(IN) :: hw
      real(dp), INTENT(INOUT) :: D(4)
      !
      real(dp) :: DAV, D_small_prev, D_large_prev
      !
      DAV = (D(2) + D(3))/2.0_dp
      if (abs((D(2) - D(3))/(DAV + hw)) < tet_cutoff) then
         D_small_prev = D(2); D_large_prev = D(3)
         D(3) = DAV + 0.5_dp*abs(DAV + hw)*tet_cutoff
         D(2) = DAV - 0.5_dp*abs(DAV + hw)*tet_cutoff
         if (D(1) > D(2)) D(1) = D(1) + (D(2) - D_small_prev)
         if (D(3) > D(4)) D(4) = D(4) + (D(3) - D_large_prev)
      endif
      DAV = (D(1) + D(2))/2.0_dp
      if (abs((D(1) - D(2))/(DAV + hw)) < tet_cutoff) then
         if (D(2) > 0) then
            D(1) = D(2)*(2.0_dp - tet_cutoff)/(2.0_dp + tet_cutoff)
         else
            D(1) = D(2)*(2.0_dp + tet_cutoff)/(2.0_dp - tet_cutoff)
         endif
      endif
      DAV = (D(3) + D(4))/2.0_dp
      if (abs((D(3) - D(4))/(DAV + hw)) < tet_cutoff) then
         if (D(3) > 0) then
            D(4) = D(3)*(2.0_dp + tet_cutoff)/(2.0_dp - tet_cutoff)
         else
            D(4) = D(3)*(2.0_dp - tet_cutoff)/(2.0_dp + tet_cutoff)
         endif
      endif
   end subroutine
   !
   PURE FUNCTION real_vertices(hw, D) result(wR0)
      !
      real(dp), INTENT(IN) :: hw
      real(dp), INTENT(IN) :: D(4)
      !
      real(dp) :: wR0(4)
      real(dp) :: dd(3), ll(3), ff, bb(4), cc(4, 3)
      integer :: a, b, c, i
      !
      ! intermediate variables for case 1 and 3
      !
      do i = 1, 3
         dd(i) = (D(4) - D(i))/(D(i) + hw)
         ll(i) = tetrahedron_log1p(dd(i))
      enddo
      !
      ff = 1.0_dp
      do i = 1, 3
         a = i
         b = mod(i, 3) + 1
         c = mod(i + 1, 3) + 1
         cc(a, a) = -(1.0_dp + dd(a))*(3.0_dp*dd(a)**2 - 2.0_dp*(dd(b) + dd(c))*dd(a) + dd(b)*dd(c)) &
            *((dd(b) - dd(c))*dd(b)*dd(c))**2
         cc(b, a) = -dd(a)*(1.0_dp + dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
         cc(c, a) = dd(a)*(1.0_dp + dd(c))*(dd(a) - dd(b))*((dd(b) - dd(c))*dd(b)*dd(c))**2
         cc(4, a) = -(dd(a) - dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
         bb(a) = cc(4, a)*dd(a)
         ff = ff*(1.0_dp + dd(a))/(dd(a)*(dd(a) - dd(b)))**2
      enddo
      bb(4) = -dd(1)*dd(2)*dd(3)*((dd(1) - dd(2))*(dd(2) - dd(3))*(dd(3) - dd(1)))**2
      ff = -ff/(D(4) + hw)
      !
      do i = 1, 4
         wR0(i) = (cc(i, 1)*ll(1) + cc(i, 2)*ll(2) + cc(i, 3)*ll(3) + bb(i))
      enddo
      !
      wR0 = wR0 * ff
   END FUNCTION

   FUNCTION delta_vertices(ef, e) result(wI0)
      !
      real(dp), INTENT(IN) :: ef
      real(dp), INTENT(IN) :: e(4)
      !
      real(dp) :: wI0(4)
      !
      real(dp) :: C, a(4,4)
      !
      integer :: i, ii
      !
      !
      !
      IF(ef > e(4) .or. ef < e(1)) THEN
         wI0 = 0.0_dp
         RETURN
      ENDIF
      !
      DO ii = 1, 4
         DO i = 1, 4
            IF ( ABS(e(i)-e(ii)) < 1.d-12 ) THEN
               a(ii,i) = 0.0_dp
            ELSE
               a(ii,i) = ( ef - e(i) ) / (e(ii) - e(i) )
            END IF
         ENDDO
      ENDDO
      !
      IF( e(1) < ef .AND. ef < e(2) ) THEN
         !
         C = a(2,1) * a(3,1) * a(4,1) / (ef - e(1))
         wI0(1) = a(1,2) + a(1,3) + a(1,4)
         wI0(2:4) = a(2:4,1)

         wI0 = wI0 * C
         !
      ELSEIF( e(2) <= ef .AND. ef < e(3)) THEN
         !
         C = a(2,3) * a(3,1) + a(3,2) * a(2,4)
         !
         wI0(1) = a(1,4) * C + a(1,3) * a(3,1) * a(2,3)
         wI0(2) = a(2,3) * C + a(2,4)**2 * a(3,2)
         wI0(3) = a(3,2) * C + a(3,1)**2 * a(2,3)
         wI0(4) = a(4,1) * C + a(4,2) * a(2,4) * a(3,2)

         wI0 = wI0 / (e(4) - e(1))
         !
      ELSEIF ( e(3) <= ef .AND. ef < e(4)) THEN
         !
         C = a(1,4) * a(2,4) * a(3,4) / (e(4) - ef)
         !
         wI0(1:3) = a(1:3,4)
         wI0(4) = a(4,1) + a(4,2) + a(4,3)
         !
         wI0 = wI0 * C
         !
      ENDIF
   END FUNCTION
   !
   !
   !
   FUNCTION tetra_weights_green(ef) RESULT(wg)
      USE constants, ONLY : pi
      !-----------------------------------------------------------------------------------
      !! Calculate weights for an integral of the kind int(Ak delta(ef-ek))
      !! The resulting wg can be used as sum(Ak * wk)
      !-----------------------------------------------------------------------------------
      COMPLEX(DP) :: wg(nbnd, nqs)
      !! COMPLEX Integration weight of each k
      REAL(DP), INTENT(IN) :: ef
      !! The Fermi energy
      !
      ! ... local variables
      !
      real(dp) :: D(4)
      !! end wannier90 tetra

      REAL(DP) :: wI(nbnd,nqs), wR(nbnd, nqs)

      INTEGER :: ik, nt, ibnd, ii
      REAL(DP) :: e(4), wI0(4), wR0(4)

      ! for real part calc
      ! REAL(DP) :: wR0(4), ef_e(4), log_ef_e(4), prod_a(4), sum_a(4), second_term(4)
      ! INTEGER :: i3, j3
      !
      wg = 0._dp
      wI = 0._dp
      wR = 0._dp
      !
      DO nt = 1+my_id, ntetra, num_procs
         !
         DO ibnd = 1, nbnd
            !
            e = ek_sort(:,ibnd,nt)
            wI0 = delta_vertices(ef, e)
            !
            D = -e
            CALL rm_degen_vertices(ef, D)
            wR0 = real_vertices(ef, D)
            !
            !
            DO ii = 1, nntetra
               !
               ik = tetra(ii, nt)
               ! IF(opt_flag) THEN
               wI(ibnd,ik) = wI(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), wI0(1:4))
               wR(ibnd,ik) = wR(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(:,ibnd,nt),ii), WR0(1:4))
               ! ELSE
               !   ik_s = tetra(itetra(ii,ibnd,nt), nt)
               !   wI(ibnd,ik_s) = wI(ibnd,ik_s) + wI0(ii)
               !   wR(ibnd,ik_s) = wR(ibnd,ik_s) + wR0(ii)
               ! ENDIF
            ENDDO
            !
         ENDDO ! ibnd
         !
      ENDDO ! nt
      ! wg = wg / REAL(ntetra, dp)
      wg = CMPLX(wR, -pi*wI, kind=DP) / (6.0_dp * nqs)
      !
      ! I LEFT OUT THE PART OF AVERAGING OF DEGENERACIES
      CALL mpi_bsum(nbnd, nqs, wg)
   END FUNCTION

   ! subroutine sort(list)
   !   !! Swap sort list of reals

   !   real(DP), intent(inout) :: list(:)
   !   real(DP) :: aux, tmp
   !   integer :: i, j, n

   !   n = size(list)

   !   do i = 1, n
   !     aux = list(i)
   !     do j = i + 1, n
   !       if (aux > list(j)) then
   !         tmp = list(j)
   !         list(j) = aux
   !         list(i) = tmp
   !         aux = tmp
   !       end if
   !     end do
   !   end do
   ! end subroutine sort
   !
   !--------------------------------------------------------------------
   SUBROUTINE deallocate_tetra( )
      !--------------------------------------------------------
      !! Deallocate tetra and wlsm
      !
      ntetra = 0
      nntetra = 0
      nqs = 0
      nbnd = 0
      IF (ALLOCATED(tetra  )) DEALLOCATE (tetra  )
      IF (ALLOCATED(wlsm   )) DEALLOCATE (wlsm   )
      IF (ALLOCATED(ek_sort)) DEALLOCATE (ek_sort)
      IF (ALLOCATED(itetra )) DEALLOCATE (itetra )
      !
   END SUBROUTINE deallocate_tetra
   !
   !
   ! FUNCTION tetra_weights_theta( nks, nbnd, et, ef) RESULT(wg)
   !   !--------------------------------------------------------------------
   !   !! Calculates weights with the tetrahedron method (P.E.Bloechl).
   !   !! Fermi energy has to be calculated in previous step.
   !   !! Generalization to noncollinear case courtesy of Iurii Timrov.
   !   !! @Note (P. Delugas 8/10/2019) Needs to be called only after initializations,
   !   !!       stops the program with an error call otherwise.
   !   !!
   !   !
   !   USE kinds
   !   !
   !   IMPLICIT NONE
   !   !
   !   INTEGER, INTENT(IN) :: nks
   !   !! Total # of k in irreducible BZ
   !   INTEGER, INTENT(IN) :: nbnd
   !   !! number of bands
   !   REAL(DP), INTENT(IN) :: et(nbnd,nks)
   !   !! eigenvalues of the hamiltonian
   !   REAL(DP) :: wg(nbnd,nks)
   !   !! the weight of each k point and band
   !   ! wg must be (inout) and not (out) because if is/=0 only terms for
   !   ! spin=is are initialized; the remaining terms should be kept, not lost.
   !   REAL(DP), INTENT(IN) :: ef
   !   !! Fermi energy

   !   ! ... local variables
   !   !
   !   REAL(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra(4), dosef
   !   INTEGER :: ibnd, nt, nk, i, kp1, kp2, kp3, kp4, itetra(4)
   !   !
   !   nk = 0
   !   wg = 0._dp
   !   !
   !   DO nt = 1, ntetra
   !     DO ibnd = 1, nbnd
   !       !
   !       ! etetra are the energies at the vertexes of the nt-th tetrahedron
   !       !
   !       DO i = 1, 4
   !         etetra(i) = et (ibnd, tetra(i,nt) + nk)
   !       ENDDO
   !       itetra (1) = 0
   !       CALL hpsort( 4, etetra, itetra )
   !       !
   !       ! ...sort in ascending order: e1 < e2 < e3 < e4
   !       !
   !       e1 = etetra(1)
   !       e2 = etetra(2)
   !       e3 = etetra(3)
   !       e4 = etetra(4)
   !       !
   !       ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
   !       !
   !       kp1 = tetra(itetra(1), nt) + nk
   !       kp2 = tetra(itetra(2), nt) + nk
   !       kp3 = tetra(itetra(3), nt) + nk
   !       kp4 = tetra(itetra(4), nt) + nk
   !       !
   !       ! calculate weights wg
   !       !
   !       IF (ef>=e4) THEN
   !         !
   !         wg(ibnd, kp1) = wg(ibnd, kp1) + 0.25d0 / ntetra
   !         wg(ibnd, kp2) = wg(ibnd, kp2) + 0.25d0 / ntetra
   !         wg(ibnd, kp3) = wg(ibnd, kp3) + 0.25d0 / ntetra
   !         wg(ibnd, kp4) = wg(ibnd, kp4) + 0.25d0 / ntetra
   !         !
   !       ELSEIF (ef<e4 .AND. ef>=e3) THEN
   !         !
   !         c4 = 0.25d0 / ntetra * (e4 - ef)**3 / (e4 - e1) / (e4 - e2) &
   !           / (e4 - e3)
   !         dosef = 3.d0 / ntetra * (e4 - ef)**2 / (e4 - e1) / (e4 - e2) &
   !           / (e4 - e3)
   !         wg(ibnd,kp1) = wg(ibnd,kp1) + 0.25d0 / ntetra - c4 * &
   !           (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
   !           (ibnd, kp1) ) / 40.d0
   !         wg(ibnd,kp2) = wg(ibnd,kp2) + 0.25d0 / ntetra - c4 * &
   !           (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
   !           (ibnd, kp2) ) / 40.d0
   !         wg(ibnd,kp3) = wg(ibnd,kp3) + 0.25d0 / ntetra - c4 * &
   !           (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
   !           (ibnd, kp3) ) / 40.d0
   !         wg(ibnd,kp4) = wg(ibnd,kp4) + 0.25d0 / ntetra - c4 * &
   !           (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
   !           + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
   !           et(ibnd,kp4) ) / 40.d0
   !         !
   !       ELSEIF (ef<e3 .AND. ef>=e2) THEN
   !         !
   !         c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
   !         c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
   !           / (e4 - e1) / (e3 - e2) / (e3 - e1)
   !         c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
   !           / (e3 - e2) / (e4 - e1)
   !         dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
   !           (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
   !           * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
   !         wg(ibnd, kp1) = wg(ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
   !           / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
   !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
   !         wg(ibnd, kp2) = wg(ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
   !           * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
   !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
   !         wg(ibnd, kp3) = wg(ibnd, kp3) + (c1 + c2) * (ef - e1) &
   !           / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
   !           (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
   !         wg(ibnd, kp4) = wg(ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
   !           / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
   !           e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
   !         !
   !       ELSEIF (ef<e2 .AND. ef>=e1) THEN
   !         !
   !         c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
   !           / (e4 - e1)
   !         dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
   !           / (e4 - e1)
   !         wg(ibnd, kp1) = wg(ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
   !           * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
   !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
   !         wg(ibnd, kp2) = wg(ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
   !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
   !         wg(ibnd, kp3) = wg(ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
   !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
   !         wg(ibnd, kp4) = wg(ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
   !           + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0

   !         ! c4 = (ef-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)/4/ntetra
   !         ! wg(ibnd,kp1) = wg(ibnd,kp1) + (1 + (ef-e2)/(e1-e2) + (ef-e3)/(e1-e3) + (ef-e4)/(e1-e4))*c4
   !         ! wg(ibnd,kp2) = wg(ibnd,kp2) + (ef-e1)/(e2-e1)*c4
   !         ! wg(ibnd,kp3) = wg(ibnd,kp3) + (ef-e1)/(e3-e1)*c4
   !         ! wg(ibnd,kp4) = wg(ibnd,kp4) + (ef-e1)/(e4-e1)*c4
   !       ENDIF
   !       !
   !     ENDDO
   !   ENDDO
   !   !
   ! END FUNCTION tetra_weights_theta

   PURE function tetrahedron_log1p(x)
      implicit none

      real(dp) :: tetrahedron_log1p
      real(dp), intent(in) :: x
      real(dp) :: y, z

      if (ABS(x) > 0.5_dp) then
         tetrahedron_log1p = LOG(ABS(1.0_dp + x))
      else
         y = 1.0_dp + x
         z = y - 1.0_dp
         if (z == 0) then
            tetrahedron_log1p = x
         else
            tetrahedron_log1p = x*LOG(y)/z
         endif
      endif

   end function tetrahedron_log1p

END MODULE thtetra
