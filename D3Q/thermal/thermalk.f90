! version 01\07\2013  
! contact info: giorgia.fugallo@gmail.com 
!================================================================================
Program r2q
  !================================================================================
  !
  USE constants
  USE common_variables, ONLY : &
    omega, at, nat, nat3, &
    asr2, F_0, velph, &
    freq1_, iq1_loc, nq1loc, nq1tot, nq1loc_max, deltaq1, &
    ntemp, tempm1, temp, tcond, const_cond, &
    sqrtm1, nRbig32, nRbig31t, nRbig2t, nRbig2_in, &
    iRbig2_in, mat2, mat2_in, mat3_in, mat3_rr, &
    alloc_dealloc_common
  USE mat_input
  USE mpi_base
!
  IMPLICIT NONE
  CALL start_parallel()
  !
  CALL read_input
  !
  !
  CALL read_mat2R
  CALL asr2_s(nRbig2_in, nat, mat2_in, iRbig2_in, asr2)
  CALL change_m2r_index
  CALL div_mass_mat2(nat3, nRbig2t, mat2, sqrtm1)
  !
  CALL read_mat3R
  CALL change_m3r_index_asr3(mat3_in)
  CALL change_m3q_index(nat3, nRbig31t, nRbig32, mat3_rr)
  CALL div_mass_mat3(nat3, nRbig31t, nRbig32, mat3_rr, sqrtm1)
  !
  ! deltaq1 and deltaq2 should be in crystal coords.
  ! if they are given in cartesian coords, use the following to transform
  CALL cryst_to_cart(1, deltaq1, at,-1)! from  cartesian to crystal
  !call cryst_to_cart(1,deltaq2,at,-1)
  !
  IF(nrank == 1) write(*,*) 'pre alloc'
  CALL alloc_dealloc_common(+1)
  !
  IF(nrank == 1) write(*,*) 'post alloc'
  IF(nrank == 1) CALL write_header
  CALL initializations
!
  CALL define_iq1_loc(nq1tot, iq1_loc, nq1loc, nq1loc_max, nrank, nsize)
!
  CALL calc_vel
!
!
  CALL loop_on_qpoints
!
  IF(nrank == 1) THEN
   CALL write_F0
   CALL calc_tcSMA
   CALL calc_conduct1(nat3, nq1tot, ntemp, const_cond, &
                       omega, at, velph, freq1_, tempm1, temp, F_0, tcond)
  ENDIF
!
!   CALL LANCZOS()
  
  CALL stop_parallel
!   stop 10

  CALL conjugate_gradient
!
  CALL alloc_dealloc_common(-1)
  !
  CALL stop_parallel
  !
  !================================================================================
END Program r2q
!================================================================================
!
!--------------------------------------------------------------------------------
SUBROUTINE conjugate_gradient
  !--------------------------------------------------------------------------------
  !
  USE constants
  USE common_variables, ONLY : &
    hh, gp, ggrd, ggt, matrix, hh_old, &
    Q_0rad, F_0, F_old, F_new, dim1, xcryst, nxcryst, &
    ntemp, temp, tempm1, &
    freq1_, velph, ltobedone, const_cond, dimrun, &
    tcond3, tcond, tcond2, &
    omega, nat3, at, &
    nq2tot, nq1tot
  USE mpi_base
  IMPLICIT NONE
  Integer :: iter, it, ix, im, ix_, irun, ii, jj
  Integer :: iq, i1
!
  Real(DP) :: F_oo(dim1, dimrun), F_nn(dim1, dimrun)
  Real(DP) :: rho(dimrun), rhold(dimrun), dgamma(dimrun), &
               dlambda(dimrun), aa, cc, pss(dimrun), tcond3old(3, 3, ntemp)
!
  F_old = 0.0d0
  F_new = 0.0d0
  ggrd = 0.0d0
  gp = 0.0d0
  ggt = 0.0d0
  hh = 0.0d0
  hh_old = 0.0d0
!
!
!  ltobedone(:)=.TRUE.
!
  IF(nrank == 1) THEN
     Write(*,*) '---------------------------------------------------------'
     Write(*,*) 'dim1,dimrun,nat3,nq1tot,ntemp'
     Write(*,*)  dim1, dimrun, nat3, nq1tot, ntemp
     Write(*,*)  'nxcryst, xcryst'
     Write(*,*)   nxcryst, xcryst
     Write(*,*) '---------------------------------------------------------'
  END IF
!
  DO it = 1, ntemp
     DO ix_ = 1, nxcryst
        ix = xcryst(ix_)
        i1 = ix + 3 *(it-1)
        DO iq = 1, nq1tot
           DO im = 1, nat3
              irun = im + nat3 *(iq-1)
              F_oo(irun, i1) = F_0(im, iq, ix, it) * Q_0rad(im, iq, ix, it)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!
! F_old la dai come variabile di input a llop_on_q_points ed a matrix_xxx
! Localmente dimensioni : F_oo(dim1,dimrun)
!
  IF(nrank == 1) write(*,*) 'Conjugate Gradient'
!
  DO iter = 1, 100
     IF(iter == 1) THEN
!         call loop_on_qpointsITcg(iter) ! questo mi restituisce A e Ax-b
        CALL loop_on_qpointsITcg(iter, F_oo, ggrd)! questo mi restituisce A e Ax-b
     END IF ! end loop 1
!
     gp(:, :) = ggrd(:, :)
!
     DO jj = 1, dimrun
        rho(jj) = 0.d0
        DO ii = 1, dim1
           rho(jj) = rho(jj) + gp(ii, jj) **2
        ENDDO
     ENDDO
!
!     IF(nrank == 1) write(*,*) 'rho', rho
!
  ! if(rho.lt.delta) nconv=nconv+1
     IF(iter == 1) THEN
        hh(:, :) = - ggrd(:, :)! compute h
     ELSE
!
        DO jj = 1, dimrun
           IF(ltobedone(jj)) THEN
              IF(rhold(jj)==0.d0)THEN
                print*, ">>>>>>>>> Zero rhold?", rho(jj), rhold(jj), jj
                dgamma(jj) = 0.d0
              ELSE
                dgamma(jj) = rho(jj) / rhold(jj)
              ENDIF
              hh(:, jj) = - ggrd(:, jj) + dgamma(jj) * hh_old(:, jj)! compute h
           ELSE
              hh(:, jj) = - ggrd(:, jj)
           END IF
        ENDDO
!
!        IF(nrank == 1) write(*,*) 'dgamma', dgamma
     END IF
!
     DO jj = 1, dimrun
        pss(jj) = 0.d0
        DO ii = 1, dim1
           pss(jj) = pss(jj) + ggrd(ii, jj) * hh_old(ii, jj)
        ENDDO
     ENDDO
!
!
!
!     IF(nrank == 1) write(*,*) 'grad dd h_old', pss
!
!
!        ggt(:)=0.0d0
     ggt(:, :) = 0.0d0
     matrix(:, :, :, :) = 0.0d0
!!        matrix(:,:)
     CALL matrixA_times_h(F_oo, ggt, hh)
!
     DO jj = 1, dimrun
        aa = 0.0d0
        DO ii = 1, dim1
           aa = aa + hh(ii, jj) * ggrd(ii, jj)
        ENDDO
        cc = 0.0d0
        DO ii = 1, dim1
           cc = cc + hh(ii, jj) * ggt(ii, jj)
        ENDDO
        dlambda(jj) = 0.d0
        IF(ltobedone(jj)) THEN
           IF(cc==0.d0) THEN
             print*, ">>>>>>>>>> Zero gradient?", aa, cc, jj
             dlambda(jj) =0.d0
           ELSE
             dlambda(jj) = - aa / cc ! compute lambda
           ENDIF
        ENDIF
     ENDDO
!
!
!     IF(nrank == 1) write(*,*) 'dlambda', dlambda
!
     DO jj = 1, dimrun
        F_nn(:, jj) = F_oo(:, jj) + dlambda(jj) * hh(:, jj)
        ggrd(:, jj) = gp(:, jj) + dlambda(jj) * ggt(:, jj)
     ENDDO
!
!
     rhold(:) = rho(:)
     hh_old(:, :) = hh(:, :)
!
!
     IF(nrank == 1) THEN
! we evaluated the conductivity of the previous step in order to not call  matrixA_times_h 2 times in a cycle
        IF(iter > 1) THEN
           Write(*,*) 'ITERATION', iter - 1
           CALL calc_conduct(nat3, nq1tot, ntemp, const_cond, omega, at, velph, freq1_, tempm1, temp, F_oo, tcond)
           CALL calc_conduct2(nat3, nq2tot, nq1tot, ntemp, const_cond, at, matrix, F_oo, tcond2, temp)
!
           CALL calc_conduct3(ntemp, temp, const_cond, tcond, tcond2, tcond3)
           DO it = 1, ntemp
!
              tcond3old(:, :, it) = tcond3(:, :, it)
!
           ENDDO
        END IF
     END IF
!
     F_oo(:, :) = F_nn(:, :)
!            F_old(:,:,:,:) = F_new(:,:,:,:)
     IF(nrank == 1) CALL write_FF(iter, F_nn)
!
  ENDDO !end iteration loop
!
  Return
    !--------------------------------------------------------------------------------
END SUBROUTINE conjugate_gradient
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE LANCZOS
  !--------------------------------------------------------------------------------
  !
  USE constants
  USE MPI_BASE
  USE common_variables, ONLY : F_0, ntemp, nat3, nq1tot, nxcryst
  USE mpi_base
  IMPLICIT NONE
!
  Integer,PARAMETER :: iter_max = 100
  Real(DP) :: alpha(3,ntemp,iter_max), beta(3,ntemp,iter_max), norm, norm2, value, inorm
  Integer  :: i,j, ix, it, ierr
  Real(DP),ALLOCATABLE :: vj(:,:,:,:), vjm(:,:,:,:), wj(:,:,:,:)
  Real(DP) :: mat(iter_max,iter_max), vec(iter_max,iter_max), eig(iter_max)

  Real,ALLOCATABLE :: oldvj(:,:,:,:,:)
  Logical,PARAMETER :: ortho = .true., random = .false., antisym = .false.
  Real(DP),EXTERNAL :: rndm

  ALLOCATE(vj(nat3, nq1tot, 3, ntemp))
  ALLOCATE(wj(nat3, nq1tot, 3, ntemp))
  ALLOCATE(vjm(nat3, nq1tot, 3, ntemp))

  IF (ortho) ALLOCATE(oldvj(nat3,nq1tot,3,ntemp,iter_max))
  
!   CALL antisymm_f(F_0)

  CALL A_times_f(F_0, wj)
  DO it = 1,ntemp
  DO ix = 1,nxcryst
    value = SUM( F_0(:,:,ix,it)*wj(:,:,ix,it) )
    norm2 =  SUM( F_0(:,:,ix,it)**2)
    IF(nrank==1) print*, "Approx 0 = ", value/norm2
  ENDDO
  ENDDO

  IF(random)THEN
    DO it = 1,ntemp
    DO ix = 1,nxcryst
    DO i = 1, nq1tot
    DO j = 1, nat3
      F_0(j,i,ix,it) = rndm()
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDIF


  DO it = 1,ntemp
  DO ix = 1,nxcryst
    norm = DSQRT( SUM(F_0(:,:,ix,it)**2) )
    inorm = 1/norm
    vj(:,:,ix,it) = inorm*F_0(:,:,ix,it)
  ENDDO
  ENDDO
  
  vjm   = 0.d0  !v_(j-1)
  wj    = 0.d0
  alpha = 0.d0
  beta  = 0.d0
  !lanczos start
  DO j = 1, iter_max
  
    IF(antisym) CALL antisymm_f(vj)
    IF(ortho) oldvj(:,:,:,:,j) = vj
    !
    CALL A_times_f(vj, wj)
    !
    DO it = 1,ntemp
    DO ix = 1,nxcryst
      alpha(ix,it,j) = SUM(wj(:,:,ix,it)*vj(:,:,ix,it))
    ENDDO
    ENDDO
    !
    WRITE(20000,*) j, alpha(1,1,j), beta(1,1,j)
    !
!     CALL MPI_ALLREDUCE(MPI_IN_PLACE, alpha(:,:,j), 3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
    
    DO it = 1,ntemp
    DO ix = 1,nxcryst
      wj(:,:,ix,it) = wj(:,:,ix,it) &
                     -alpha(ix,it,j)* vj(:,:,ix,it) &
                     -beta(ix,it,j) * vjm(:,:,ix,it)
    ENDDO
    ENDDO
    !
    IF(j<iter_max)THEN
      DO it = 1,ntemp
      DO ix = 1,nxcryst
        !
        beta(ix,it,j+1) = DSQRT( SUM(wj(:,:,ix,it)**2) )
        vjm(:,:,ix,it) = vj(:,:,ix,it)
        vj(:,:,ix,it)  = wj(:,:,ix,it) / beta(ix,it,j+1)
        !
        ! From here on, vj is v_(j+1)
        ! Graham-Schmidt orthogonalise to all previous vectors
        IF(ortho)THEN
          DO i = 1,j
            vj(:,:,ix,it) = vj(:,:,ix,it) - SUM(vj(:,:,ix,it)*oldvj(:,:,ix,it,i))*oldvj(:,:,ix,it,i)
          ENDDO
        ENDIF
        !
      ENDDO
      ENDDO
!       CALL MPI_ALLREDUCE(MPI_IN_PLACE, beta(:,:,j), 3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
!       CALL MPI_ALLREDUCE(MPI_IN_PLACE, vjm, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
!       CALL MPI_ALLREDUCE(MPI_IN_PLACE, vj, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
    ENDIF
    !     print*, "alphabeta", alpha(1,1,j), beta(1,1,j+1)
    !
    IF(nrank==1) THEN
    DO it = 1,ntemp
    DO ix = 1,nxcryst
      mat = 0.d0
      DO i=1,j
        mat(i,i) = alpha(ix,it,i)
        IF(i<j)THEN
          mat(i,i+1) = beta(ix,it,i+1)
          mat(i+1,i) = beta(ix,it,i+1)
        ENDIF
      ENDDO
      CALL rdiagh(j, mat, iter_max, eig, vec)
      
      write(10000+100*it+ix,'(1000e12.4)') eig(1:j)
      
    ENDDO
    ENDDO
    ENDIF
    !
  ENDDO
  
  write(20001,'(1e25.12)') eig

  !  
  DEALLOCATE(vj, wj, vjm)
  IF(ortho) DEALLOCATE(oldvj)

  Return
    !--------------------------------------------------------------------------------
END SUBROUTINE LANCZOS
!--------------------------------------------------------------------------------

SUBROUTINE antisymm_f(f)
  USE constants, ONLY: DP
  USE common_variables, ONLY : nq1tot, nq1loc, iq1_loc, nq1, &
                               ntemp, nat3, bg, &
                               ic_inn, ic_med, ic_out, deltaq1, deltaq2

  
  IMPLICIT NONE
  REAL(DP),INTENT(inout)  ::  f(nat3, nq1tot, 3, ntemp)
  Integer :: iq1_, mq1_, iq_init(3), dummy
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out
  Integer, Target :: iq1(3)
  Real (DP) :: q1d (3, nq1tot), aux(3), aux2(3)

  Logical,SAVE :: first = .true.
  Integer,ALLOCATABLE,SAVE :: idx_minus(:)
  
  Integer :: far = 2
  Integer :: ix, it
  Logical :: done(nq1tot)
  

  IF (first) THEN
    ALLOCATE(idx_minus(nq1tot))
    first = .false.
    iq1_out => iq1(ic_out)
    iq1_med => iq1(ic_med)
    iq1_inn => iq1(ic_inn)
  !
    IF(ANY(deltaq1/=0.d0).or.ANY(deltaq2/=0.d0)) CALL errore("antisymm_f","deltas detected",1)

    ! WARNING: I think q1d comes out in crystal coords, but you never know
    DO iq1_ = 1, nq1tot
      CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, -1)
      iq1(:) = iq_init(:)
      q1d(:,iq1_) = dfloat(iq1(:)) / dfloat(nq1(:))
    ENDDO

    !WARNING: I do not understand how to anti-symmetrize along a specific direction
    !
    idx_minus = -1
    PLUS : DO iq1_ = 1, nq1tot
    MINUS : DO mq1_ = 1, nq1tot
      aux(1) = q1d(1,iq1_)+q1d(1,mq1_)
      aux(2) = q1d(2,iq1_)+q1d(2,mq1_)
      aux(3) = q1d(3,iq1_)+q1d(3,mq1_)
      
      aux= aux-INT(aux)
      IF(  SUM(ABS(aux))<1.d-6  )THEN
        idx_minus(iq1_) = mq1_
        EXIT MINUS
      ENDIF
      
    ENDDO MINUS
    ENDDO PLUS
    
    IF(ANY(idx_minus==-1))THEN
      CALL errore("antisymm_f", "missing minus", 1)
    ENDIF
    
!     print*, idx_minus
  ENDIF
  !
  DO it = 1,ntemp
  DO ix = 1,3
    done = .false.
    DO iq1_ = 1, nq1tot
      mq1_ = idx_minus(iq1_)
      IF(.not.done(mq1_))THEN
        F(:,iq1_,ix,it) = 0.5d0*(F(:,iq1_,ix,it)-F(:,mq1_,ix,it))
        F(:,mq1_,ix,it) = -F(:,iq1_,ix,it)
        done(mq1_) = .true.
      ENDIF
    ENDDO
  ENDDO
  ENDDO


END SUBROUTINE antisymm_f

!--------------------------------------------------------------------------------
SUBROUTINE A_times_f(f,Af)
  !--------------------------------------------------------------------------------
  !
  USE constants, ONLY: DP, pi, sqrtpi, BOHR_TO_M, eps8
  USE common_variables, ONLY : &
     bos1, bos2, bos3, bos3b, &
     freq1, freq2, freq3, freq3b, &
     ntemp, tempm1, sigmam1, &
     nq1tot, nq1loc, iq1_loc, nq1, q1d, &
     iq1_, deltaq1, nq2tot, mq2d, nq2, q2d, deltaq2, q3d, q3db, &
     imq2_, ipq2_, ipq3_, ipq3b_, &
     Q_0, Q_0rad, Q_0radm1, &
     ic_inn, ic_med, ic_out, &
     nxcryst, xcryst, &
     nat, nat3, const_iso, &
     zz2_x, zz1_x, &
     d3mmx, d3mmx_1, ierr, &
     sum_IT, matrix
  USE mpi_base
  IMPLICIT NONE
!
  REAL(DP),INTENT(in)  ::  f(nat3, nq1tot, 3, ntemp)
  REAL(DP),INTENT(out) :: Af(nat3, nq1tot, 3, ntemp)

!
  Integer :: im3, ii, isw1, isw2, km, ix_
  Integer :: im, im2, ic, it, ix, iq_init(3)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out, iq2_out, iq2_med, iq2_inn
  Integer, Target :: iq1(3), iq2(3)
  Real(DP) :: matAm3, d3mm3(nat3), d3mm_13(nat3)
  Integer :: iqqmm1, iqqmm2
!
  Complex(DP) :: workc
  Real(DP) :: zz_12(nat3, nat3), norm_iso1
  Real(DP) :: bos_iso, dom_iso, ctm_iso, wimp_iso, prexp
!
  Integer :: ia
  Real(DP) :: freq1t(nat3, nq1tot)
  Real(DP) :: matrixloc(nat3, nq1tot, 3, ntemp)

!
  iq1_out => iq1(ic_out)
  iq1_med => iq1(ic_med)
  iq1_inn => iq1(ic_inn)
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
!
!
  isw1 = - 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  isw2 = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
!--------------------------------------------------------------------------------
!for isotopic scattering contribution
  norm_iso1 = pi * const_iso * 0.5d0 /(dfloat(nq2tot))! qui a differenza che nel loop dei qpoints la moltiplicazione per freq1 non e'
                                                      !in norm_iso ma dopo all'interno del loop
!     prexp= sigmam1(1) / sqrtpi
!--------------------------------------------------------------------------------
  matrixloc(:, :, :, :) = 0.0d0
!--------------------------------------------------------------------------------
!            EXTERNAL loop
!--------------------------------------------------------------------------------
  ! IF you change the order of these loop; you should change calc_vel accordingly
!
!
  LOOP_IQ1 : &
  DO ii = 1, nq1loc
     iq1_ = iq1_loc(ii)
!
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
!
     CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw1)
     iq1(:) = iq_init(:)
     q1d(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
!     iq_init(:) = iq1(:)
!--------------------------------------------------------------------------------
!
!
!      sum_IT3(:,:,:) = 0.0d0
     CALL setup(1)
!
!--------------------------------------------------------------------------------
!            INTERNAL loop
!--------------------------------------------------------------------------------
     LOOP_IQ2_OUT :  DO iq2_out = 0, nq2(ic_out) - 1
        CALL setup(2)
        LOOP_IQ2_MED : DO iq2_med = 0, nq2(ic_med) - 1
           CALL setup(3)
           LOOP_IQ2_INN : DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
              q2d(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
              iq_init(:) = iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw2)
!
              mq2d(:) = - q2d(:)
              iq_init(:) = - iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq2_, isw2)
!
!
              q3d(:) = - q1d(:) - q2d(:)
              iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) - 2 * Nint(deltaq1(:)*nq2(:)) - Nint(deltaq2(:)*nq2(:))

              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw2)
!
              q3db(:) = q2d(:) - q1d(:)
              iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint(deltaq1(:)*nq2(:))
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw2)
!--------------------------------------------------------------------------------
              CALL setup(4)
              !
              DO im = 1, nat3
                 freq1t(im, iq1_) = freq1(im)
              ENDDO
!
!
! write(*,*)'imq3',imq3_
! irun = 0     
              LOOP_TEMP : &
              DO it = 1, ntemp
                 prexp = sigmam1(it) / sqrtpi
                 DO im = 1, 3 * nat
                    bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(it))-1.d0)
                    bos2(im) = 1.d0 /(Exp(freq2(im)*tempm1(it))-1.d0)
                    bos3(im) = 1.d0 /(Exp(freq3(im)*tempm1(it))-1.d0)
                    bos3b(im) = 1.d0 /(Exp(freq3b(im)*tempm1(it))-1.d0)
                 ENDDO
                 !
                 LOOP_IX_ : &
                 DO ix_ = 1, nxcryst
                    ix = xcryst(ix_)
                    zz_12(:, :) = 0.0d0
                    !
                    LOOP_IM : &
                    DO im = 1, nat3
                       iqqmm1 = im + nat3 *(iq1_-1)
!
                       DO im2 = 1, 3 * nat !
                          DO ia = 1, nat
!
                             workc =(0.0d0, 0.0d0)
                             DO ic = 1, 3
!
                                km = ic + 3 *(ia-1)
                                workc = workc + CONJG(zz1_x(km, im, ix)) * zz2_x(km, im2, ix)
                             ENDDO
!
                             zz_12(im2, im) = zz_12(im2, im) + CONJG(workc) * workc
!
                          ENDDO
                       ENDDO
!
!
                       LOOP_IM2 : &
                       DO im2 = 1, nat3
                          iqqmm2 = im2 + nat3 *(ipq2_-1)
!
                          DO im3 = 1, nat3
                             d3mm3(im3) = d3mmx(im2, im3, im, ix)
                             d3mm_13(im3) = d3mmx_1(im2, im3, im, ix)
                          ENDDO
!
                          matAm3 = 0.0d0
!--------------------------------------------------------------------------------
! isotopic scattering contribution
                          bos_iso = bos1(im) * bos2(im2) + 0.5d0 *(bos1(im)+bos2(im2))
                          dom_iso =(freq1(im)-freq2(im2)) * sigmam1(it)
                          ctm_iso = prexp * Exp(-(dom_iso*dom_iso))
                          wimp_iso = norm_iso1 * freq1(im) * freq2(im2) * bos_iso * ctm_iso * zz_12(im2, im)
!--------------------------------------------------------------------------------
!
!
!                           CALL sum_modes3a2(nat3, nq2tot, sigmam1(it), freq1(im), freq2(im2), freq3, freq3b, &
!                                             bos1(im), bos2(im2), bos3, bos3b, d3mm3, d3mm_13, matAm3,+1)
! !
                          CALL sum_modes3a3(nat3, nq2tot, sigmam1(it), iq1_, ipq2_, im, im2, freq1(im), freq2(im2), freq3, &
                            & freq3b, bos1(im), bos2(im2), bos3, bos3b, d3mm3, d3mm_13, matAm3, Q_0(im2, ipq2_, ix, it), &
                            & wimp_iso)


                             IF((iq1_ == 1) .and.(im == 1)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 2)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 3)) matAm3 = 0
!
                             IF((ipq2_ == 1) .and.(im2 == 1)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 2)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 3)) matAm3 = 0


!                           matAm3 = matAm3 - wimp_iso
!
!                           IF(Q_0rad(im, iq1_, ix, it) /= 0.0d0 .and. Q_0rad(im2, ipq2_, ix, it) /= 0.0d0) THEN
!                             matAm3 = Q_0radm1(im, iq1_, ix, it) * matAm3 * Q_0radm1(im2, ipq2_, ix, it)
!                           ELSE
!                             matAm3 = 0.0d0
!                           END IF

                          matrixloc(im, iq1_, ix, it) = matrixloc(im, iq1_, ix, it) + matAm3 * f(im2, ipq2_, ix, it)
!
                       ENDDO LOOP_IM2 
                       
                    ENDDO LOOP_IM
                 ENDDO LOOP_IX_
                 !
              ENDDO LOOP_TEMP 
!
           ENDDO LOOP_IQ2_INN 
        ENDDO LOOP_IQ2_MED
     ENDDO LOOP_IQ2_OUT
     
  ENDDO LOOP_IQ1
!
!
!   Af = matrixloc
  CALL MPI_ALLREDUCE(matrixloc, Af, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)

  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE A_times_f
!--------------------------------------------------------------------------------
!






!
!--------------------------------------------------------------------------------
SUBROUTINE setupmat(nat, nrtot, q, phiq2, Rbig2, mat2)
  !--------------------------------------------------------------------------------
  USE constants
  IMPLICIT NONE
!
  Integer, Intent(In) :: nat, nrtot
  Integer :: iR
  Complex(DP) :: phiq2(3*nat, 3*nat), phase
  Complex(DP) :: mat2(3*nat, 3*nat, nrtot)
  Real(DP) :: q(3), arg, Rbig2(3, nrtot)
!
  phiq2 = cmplx(0.d0, 0.d0)
  DO iR = 1, nrtot
     arg = tpi *(q(1)*Rbig2(1, iR)+q(2)*Rbig2(2, iR)+q(3)*Rbig2(3, iR))
     phase = cmplx(Cos(arg),-Sin(arg), kind=DP)
     phiq2(:, :) = phiq2(:, :) + phase * mat2(:, :, iR)
  ENDDO
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE setupmat
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
SUBROUTINE setupmat3a(nat33, nRbig31t, nRbig32, qq, Rbig3, mat3_rr, mat3_r)
  !--------------------------------------------------------------------------------
  USE constants
  IMPLICIT NONE
!
  Integer, Intent(In) :: nat33, nRbig31t, nRbig32
  Real(DP) :: qq(3), Rbig3(3, nRbig32)
  Complex(DP) :: mat3_rr(nat33, nRbig31t, nRbig32)
  Complex(DP) :: mat3_r(nat33, nRbig31t)
!
  Integer :: iR
  Real(DP) :: arg
  Complex(DP) :: phase
!
  mat3_r = cmplx(0.d0, 0.d0)
  DO iR = 1, nRbig32
     arg = tpi *(qq(1)*Rbig3(1, iR)+qq(2)*Rbig3(2, iR)+qq(3)*Rbig3(3, iR))
     phase = cmplx(Cos(arg),-Sin(arg), kind=DP)
     mat3_r(:, :) = mat3_r(:, :) + phase * mat3_rr(:, :, iR)
  ENDDO
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE setupmat3a
!--------------------------------------------------------------------------------
SUBROUTINE setupmat3b(nat33, nRbig3, qq, Rbig3, mat3_r, mat3)
  !--------------------------------------------------------------------------------
  USE constants
  IMPLICIT NONE
!
  Integer, Intent(In) :: nat33, nRbig3
  Real(DP) :: qq(3), Rbig3(3, nRbig3)
  Complex(DP) :: mat3_r(nat33, nRbig3), mat3(nat33)
!
  Integer :: iR
  Real(DP) :: arg
  Complex(DP) :: phase
!
  mat3 = cmplx(0.d0, 0.d0)
  DO iR = 1, nRbig3
     arg = tpi *(qq(1)*Rbig3(1, iR)+qq(2)*Rbig3(2, iR)+qq(3)*Rbig3(3, iR))
     phase = cmplx(Cos(arg),-Sin(arg), kind=DP)
     mat3(:) = mat3(:) + phase * mat3_r(:, iR)
  ENDDO
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE setupmat3b
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_d3mm(d3mmx, mat3q, zz1, zz2, zz3, nat3, nat32, nat33)
  !--------------------------------------------------------------------------------
  !
  USE constants, ONLY: DP, cone, czero
  USE common_variables, ONLY: d3mm_aux, d3mm_aux1, d3mm_aux2
  IMPLICIT NONE
  Integer :: nat3, nat32, nat33
  Real(DP) :: d3mmx(nat3, nat3, nat3, 3)
  Complex(DP) :: mat3q(nat3, nat3, nat3), zz3(nat3, nat3, 3), zz2(nat3, nat3, 3), zz1(nat3, nat3, 3)
  Integer :: im1, im2, im3, ix
!
  DO ix = 1, 3
!
     DO im3 = 1, nat3
        DO im2 = 1, nat3
           CALL ZGEMV('T', nat3, nat3, cone, zz1(1, 1, ix), nat3, mat3q(im2, im3, 1), nat32, czero, d3mm_aux(im2, im3, 1, ix), &
          & nat32)
        ENDDO
     ENDDO
!
     DO im1 = 1, nat3
        DO im3 = 1, nat3
           CALL ZGEMV('T', nat3, nat3, cone, zz2(1, 1, ix), nat3, d3mm_aux(1, im3, im1, ix), 1, czero, d3mm_aux1(1, im3, im1, &
          & ix), 1)
        ENDDO
     ENDDO
!
     DO im1 = 1, nat3
        DO im2 = 1, nat3
           CALL ZGEMV('T', nat3, nat3, cone, zz3(1, 1, ix), nat3, d3mm_aux1(im2, 1, im1, ix), nat3, czero, d3mm_aux2(im2, 1, &
          & im1, ix), nat3)
        ENDDO
     ENDDO
!
!
     CALL square_mod_zvec(nat33, d3mm_aux2(1, 1, 1, ix), d3mmx(1, 1, 1, ix))
!
  ENDDO
!
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_d3mm
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE square_mod_zvec(nat33, d3mm_aux3, d3mm)
  !--------------------------------------------------------------------------------
  !
  USE constants, ONLY: DP
  IMPLICIT NONE
  Integer :: nat33
  Complex(DP) :: d3mm_aux3(nat33)
  Real(DP) :: d3mm(nat33)
!
  d3mm(:) = DBLE(d3mm_aux3(:)*CONJG(d3mm_aux3(:)))
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE square_mod_zvec
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes(nat3, nq2tot, sigmam1, freq1, freq2, freq3, bos1, bos2, bos3, d3mm, sum, isw)
  !--------------------------------------------------------------------------------
  USE constants
  IMPLICIT NONE
  Integer :: nat3, isw, nq2tot
  Real(DP) :: sum, sigmam1
  Real(DP) :: freq1, freq2(nat3), freq3(nat3), bos1, bos2(nat3), bos3(nat3), d3mm(nat3, nat3)
!
  Integer :: im1, im2
  Real(DP) :: norm, norm1, dom1, dom2, bos_a, bos_b, freqtot, delpi
  Real(DP) :: prexp, ctm1, ctm2
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))!hbar/(8*N0) where N0=N1*N2*N3
!
  DO im2 = 1, nat3
     DO im1 = 1, nat3
        bos_a = bos2(im1) + bos3(im2) + 1.d0
        bos_b = bos2(im1) - bos3(im2)
        dom1 =(freq1+freq2(im1)-freq3(im2)) * sigmam1
        dom2 =(freq1-freq2(im1)-freq3(im2)) * sigmam1
        ctm1 = 2.0d0 * bos_b * prexp * Exp(-(dom1*dom1))
        ctm2 = bos_a * prexp * Exp(-(dom2*dom2))
        delpi = ctm1 + ctm2
        freqtot = freq1 * freq2(im1) * freq3(im2)
        norm1 = norm / freqtot
        sum = sum + norm1 * delpi * d3mm(im1, im2)
     ENDDO
  ENDDO
!
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes2iso(nat3, nq2tot, const_iso, sigmam1, freq1, freq2, freq3, freq3b, &
                          bos1, bos2, bos3, bos3b, d3mm, d3mm1,  zz_12, sum, Q_0_, sum_iso)
  !--------------------------------------------------------------------------------
  !
  USE constants
  IMPLICIT NONE
  Integer :: nat3, nq2tot
  Real(DP) :: sum, sum_iso, sigmam1, sum_c, sum_d, const_iso
  Real(DP) :: freq1, freq2(nat3), freq3(nat3), freq3b(nat3), zz_12(nat3),&
               bos1, bos2(nat3), bos3(nat3), bos3b(nat3), &
               d3mm(nat3, nat3), d3mm1(nat3, nat3)
!
  Integer :: im2, im3
  Real(DP) :: norm, norm1, norm1b, dom1, dom2, bos_a, bos_b, freqtot, freqtotb, delpi
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d, Q_0_
!
  Real(DP) :: prexp, ctm1, ctm2
  Real(DP), Parameter :: freq_window = 0.0d0 !10.d0 / ry_to_cmm1
  Real(DP) :: norm_iso, bos_iso, dom_iso, ctm_iso, wimp_iso
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))!hbar/(8*N0) where N0=N1*N2*N3
  norm_iso = pi * const_iso * 0.5d0 * freq1 /(dfloat(nq2tot))! omega^2/(2*N0)
!
  wimp_iso = 0.0d0
!
! if( dabs(freq1) .GT. eps8  ) THEN
  DO im2 = 1, nat3
 !      if( dabs(freq3(im3)) .GT. eps8 ) THEN
     DO im3 = 1, nat3
  !           if( dabs(freq2(im2)) .GT. eps8 ) THEN
        bos_a = bos2(im2) + bos3(im3) + 1.d0
        bos_b = bos2(im2) - bos3(im3)
        dom1 =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        dom2 =(freq1-freq2(im2)-freq3(im3)) * sigmam1
        ctm1 = 2.0d0 * bos_b * prexp * Exp(-(dom1*dom1))
        ctm2 = bos_a * prexp * Exp(-(dom2*dom2))
        delpi = ctm1 + ctm2
        freqtot = freq1 * freq2(im2) * freq3(im3)
        norm1 = norm / freqtot
        sum = sum + norm1 * delpi * d3mm(im2, im3)
          ! end if
     ENDDO
     dom_iso =(freq1-freq2(im2)) * sigmam1
     ctm_iso = prexp * Exp(-(dom_iso**2))
     sum_iso = sum_iso + norm_iso * freq2(im2) * ctm_iso * zz_12(im2)! la moltiplicazione per freq1 e' in norm2
!
      !  end if
  ENDDO
  !end if
!
   !  if( dabs(freq1) .GT. eps8  ) THEN
  DO im2 = 1, nat3
!    if( dabs(freq3(im3)) .GT. eps8 ) THEN
     DO im3 = 1, nat3
 !    if( dabs(freq2(im2)) .GT. eps8 ) THEN
!
        bos_c = bos1 * bos2(im2) *(bos3(im3)+1)! da moltiplicare per le F
        dom_c =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
        bos_d =(bos1+1.0d0) * bos2(im2) * bos3b(im3)! da moltiplicare per le F
        dom_d =(-freq1+freq2(im2)+freq3b(im3)) * sigmam1
        ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
        freqtot = freq1 * freq2(im2) * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
        freqtotb = freq1 * freq2(im2) * freq3b(im3)! if freqtot lt 10cm-1 do not calculate anything!
     !  if(abs(freqtot).gt.eps8) THEN
        norm1 = norm / freqtot
        norm1b = norm / freqtotb
!
        sum_c = ctm_c * norm1 * d3mm(im2, im3)
        sum_d = ctm_d * norm1b * d3mm1(im2, im3)
!
        Q_0_ = Q_0_ + sum_c + 0.5d0 * sum_d
          !  end if
  !  end if
     ENDDO
!
     bos_iso = bos1 * bos2(im2) + 0.5d0 *(bos1+bos2(im2))
     dom_iso =(freq1-freq2(im2)) * sigmam1
     ctm_iso = prexp * Exp(-(dom_iso*dom_iso))
     wimp_iso = wimp_iso + norm_iso * freq2(im2) * bos_iso * ctm_iso * zz_12(im2)! la moltiplicazione per freq1 e' in norm2
!
 !end if
  ENDDO
  !  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  Q_0_ = Q_0_ + wimp_iso
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes2iso
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes2(nat3, nq2tot, sigmam1, freq1, freq2, freq3, bos1, bos2, bos3, d3mm, sum, Q_0_, isw)
  !--------------------------------------------------------------------------------
  !
  USE constants
  IMPLICIT NONE
  Integer :: nat3, isw, nq2tot
  Real(DP) :: sum, sigmam1, sum_c, sum_d
  Real(DP) :: freq1, freq2(nat3), freq3(nat3), bos1, bos2(nat3), bos3(nat3), d3mm(nat3, nat3)
!
  Integer :: im2, im3
  Real(DP) :: norm, norm1, dom1, dom2, bos_a, bos_b, freqtot, delpi
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d, Q_0_
!
  Real(DP) :: prexp, ctm1, ctm2
  Real(DP), Parameter :: freq_window = 0.0d0 !10.d0 / ry_to_cmm1
!
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))!hbar/(8*N0) where N0=N1*N2*N3
!
!  SELECT case(isw)
!  case(1)
! if( dabs(freq1) .GT. eps8  ) THEN
  DO im3 = 1, nat3
 !      if( dabs(freq3(im3)) .GT. eps8 ) THEN
     DO im2 = 1, nat3
  !           if( dabs(freq2(im2)) .GT. eps8 ) THEN
        bos_a = bos2(im2) + bos3(im3) + 1.d0
        bos_b = bos2(im2) - bos3(im3)
        dom1 =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        dom2 =(freq1-freq2(im2)-freq3(im3)) * sigmam1
        ctm1 = 2.0d0 * bos_b * prexp * Exp(-(dom1*dom1))
        ctm2 = bos_a * prexp * Exp(-(dom2*dom2))
        delpi = ctm1 + ctm2
        freqtot = freq1 * freq2(im2) * freq3(im3)
        norm1 = norm / freqtot
        sum = sum + norm1 * delpi * d3mm(im2, im3)
          ! end if
     ENDDO
      !  end if
  ENDDO
  !end if
!
!  case(2)
   !  if( dabs(freq1) .GT. eps8  ) THEN
  DO im3 = 1, nat3
!    if( dabs(freq3(im3)) .GT. eps8 ) THEN
     DO im2 = 1, nat3
 !    if( dabs(freq2(im2)) .GT. eps8 ) THEN
!
        bos_c = bos1 * bos2(im2) *(bos3(im3)+1)! da moltiplicare per le F
        dom_c =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
        bos_d = bos1 *(bos2(im2)+1) *(bos3(im3)+1)! da moltiplicare per le F
        dom_d =(freq1-freq2(im2)-freq3(im3)) * sigmam1
        ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
        freqtot = freq1 * freq2(im2) * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
     !  if(abs(freqtot).gt.eps8) THEN
        norm1 = norm / freqtot
!
        sum_c = ctm_c * norm1 * d3mm(im2, im3)
        sum_d = ctm_d * norm1 * d3mm(im2, im3)
!
!
        Q_0_ = Q_0_ + sum_c + 0.5d0 * sum_d
          !  end if
  !  end if
     ENDDO
!
 !end if
  ENDDO
  !  end if
! case(3) ! calcola il Q=(sum_c+1/2 sum_d)
!end SELECT
!
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes2
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes3(nat3, nq2tot, sigmam1, freq1, freq2, freq3, &
                       bos1, bos2, bos3, d3mm, fp2, fm2, fm3, sumPF, isw)
  !--------------------------------------------------------------------------------
  !
  USE constants
  IMPLICIT NONE
  Integer,Intent(IN) :: nat3, isw, nq2tot
  Real(DP),Intent(IN) :: fp2(nat3), fm2(nat3), fm3(nat3)
  Real(DP),Intent(INOUT) :: sumPF
  Real(DP) :: sum_c, sum_d
  Real(DP),Intent(IN) :: sigmam1
  Real(DP),Intent(IN) :: freq1, freq2(nat3), freq3(nat3), bos1, bos2(nat3), bos3(nat3), d3mm(nat3, nat3)
!
  Integer :: im2, im3
  Real(DP) :: norm, norm1, freqtot
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d
!
  Real(DP) :: prexp
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))
!
!  SELECT case(isw)
!  case(1)
!     sumPF(:) = 0.0d0
  DO im3 = 1, nat3
!
     DO im2 = 1, nat3
!
        bos_c = bos1 * bos2(im2) *(bos3(im3)+1)
        dom_c =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
        bos_d = bos1 *(bos2(im2)+1) *(bos3(im3)+1)
        dom_d =(freq1-freq2(im2)-freq3(im3)) * sigmam1
        ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
        freqtot = freq1 * freq2(im2) * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
        norm1 = norm / freqtot
!
!
     !  if(abs(freqtot).gt.eps8) THEN
!
        sum_c = ctm_c * norm1 * d3mm(im2, im3)
        sum_d = ctm_d * norm1 * d3mm(im2, im3)
!
!
        sumPF = sumPF + sum_c *(fm3(im3)-fp2(im2)) + 0.5d0 * sum_d *(fm2(im2)+fm3(im3))! ricontrolla chi e' associato a im1 e chi a im2
!
      ! end if
 !         END IF
     ENDDO
!
 !        END IF
  ENDDO
 !         END IF
!
!end SELECT
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes3
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes3a(nat3, nq2tot, sigmam1, freq1, freq2, freq3, freq3b, &
                        bos1, bos2, bos3, bos3b, d3mm, d3mm_1, fp2, fm2, fm3, &
                        sumPF, isw)
  !--------------------------------------------------------------------------------
  USE constants
  IMPLICIT NONE
  Integer :: nat3, isw, nq2tot
  Real(DP) :: fp2(nat3), fm2(nat3), fm3(nat3)
  Real(DP) :: sumPF, sigmam1, sum_c, sum_d, sumPF_a
  Real(DP) :: freq1, freq2(nat3), freq3(nat3), freq3b(nat3),&
               bos1, bos2(nat3), bos3(nat3), bos3b(nat3), &
               d3mm(nat3, nat3), d3mm_1(nat3, nat3)
!
  Integer :: im2, im3
  Real(DP) :: norm, norm1, norm1b, freqtot, freqtotb
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d
  Real(DP) :: bos_c1, ctm_c1, dom_c1, sum_c1
!
  Real(DP) :: prexp
!
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))

  DO im2 = 1, nat3
     sumPF_a = 0.0d0
!
     DO im3 = 1, nat3
!
        bos_c = bos1 * bos2(im2) *(bos3(im3)+1.0d0)
        dom_c =(freq1+freq2(im2)-freq3(im3)) * sigmam1
        ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
        bos_c1 = bos1 * bos3b(im3) *(bos2(im2)+1.0d0)
        dom_c1 =(freq1+freq3b(im3)-freq2(im2)) * sigmam1
        ctm_c1 = 2.0d0 * pi * bos_c1 * prexp * Exp(-(dom_c1*dom_c1))

        bos_d =(bos1+1.0d0) * bos2(im2) * bos3b(im3)
        dom_d =(freq1-freq2(im2)-freq3b(im3)) * sigmam1
!
        ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
!
        freqtot = freq1 * freq2(im2) * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
        norm1 = norm / freqtot
!
        freqtotb = freq1 * freq2(im2) * freq3b(im3)! if freqtot lt 10cm-1 do not calculate anything!
        norm1b = norm / freqtotb
!
        sum_c = ctm_c * norm1 * d3mm(im2, im3)
        sum_c1 = ctm_c1 * norm1b * d3mm_1(im2, im3)
        sum_d = ctm_d * norm1b * d3mm_1(im2, im3)
        sumPF_a = sumPF_a + sum_c1 + sum_d - sum_c ! A

 
     ENDDO
!
     sumPF = sumPF + sumPF_a * fp2(im2)
!
!
  ENDDO

  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes3a
!--------------------------------------------------------------------------------!-----------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_modes3a2(nat3, nq2tot, sigmam1, freq1, freq2, freq3, freq3b, &
                         bos1, bos2, bos3, bos3b, d3mm, d3mm_1, sumPF_a, isw)
  !--------------------------------------------------------------------------------
  !
  USE constants, ONLY : DP, pi, sqrtpi
  IMPLICIT NONE
  Integer :: nat3, isw, nq2tot
  Real(DP) :: sigmam1, sum_c, sum_d, sumPF_a
  Real(DP) :: freq1, freq2, freq3(nat3), freq3b(nat3), &
               bos1, bos2, bos3(nat3), bos3b(nat3), d3mm(nat3), d3mm_1(nat3)
!
  Integer :: im3
  Real(DP) :: norm, norm1, norm1b, freqtot, freqtotb
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d
  Real(DP) :: bos_c1, ctm_c1, dom_c1, sum_c1
!
  Real(DP) :: prexp
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))
!
!
  DO im3 = 1, nat3
!
     bos_c = bos1 * bos2 *(bos3(im3)+1.0d0)
     dom_c =(freq1+freq2-freq3(im3)) * sigmam1
     ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
     bos_c1 = bos1 * bos3b(im3) *(bos2+1.0d0)
     dom_c1 =(freq1+freq3b(im3)-freq2) * sigmam1
     ctm_c1 = 2.0d0 * pi * bos_c1 * prexp * Exp(-(dom_c1*dom_c1))
!
!
     bos_d =(bos1+1.0d0) * bos2 * bos3b(im3)
     dom_d =(freq1-freq2-freq3b(im3)) * sigmam1 ! questo ad onor del vero dovrebbe essere freq2+freq3b-freq1
!
     ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
     freqtot = freq1 * freq2 * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
     norm1 = norm / freqtot

     freqtotb = freq1 * freq2 * freq3b(im3)! if freqtot lt 10cm-1 do not calculate anything!
     norm1b = norm / freqtotb
!
     sum_c = ctm_c * norm1 * d3mm(im3)
     sum_c1 = ctm_c1 * norm1b * d3mm_1(im3)
     sum_d = ctm_d * norm1b * d3mm_1(im3)
!
     sumPF_a = sumPF_a + sum_c1 + sum_d - sum_c ! A
!
  ENDDO
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes3a2
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------!-----------------------------------------------------
SUBROUTINE sum_modes3a3(nat3, nq2tot, sigmam1, iq1_, ipq2_, pim, pim2, freq1, freq2, freq3, freq3b,&
                         bos1, bos2, bos3, bos3b, d3mm,  d3mm_1, sumPF_a, Q_0_, wimp)
  !--------------------------------------------------------------------------------
  !
  USE constants
  IMPLICIT NONE
  Integer :: nat3, nq2tot
  Real(DP) :: sigmam1, sum_c, sum_d, sumPF_a
  Real(DP) :: freq1, freq2, freq3(nat3), freq3b(nat3),&
               bos1, bos2, bos3(nat3), bos3b(nat3),&
               d3mm(nat3), d3mm_1(nat3)
!
  Integer :: iq1_, ipq2_, pim, pim2
  Integer :: im3
  Real(DP) :: norm, norm1, norm1b, freqtot, freqtotb
  Real(DP) :: bos_c, bos_d, dom_c, dom_d, ctm_c, ctm_d
  Real(DP) :: bos_c1, ctm_c1, dom_c1, sum_c1, Q_0_
!
  Real(DP) :: prexp, wimp
  Real(DP), Parameter :: freq_window = eps8 !10.d0 / ry_to_cmm1
!
!
  prexp = sigmam1 / sqrtpi
  norm = 1.d0 /(8.0d0*dfloat(nq2tot))
!
  DO im3 = 1, nat3
!
     bos_c = bos1 * bos2 *(bos3(im3)+1.0d0)
     dom_c =(freq1+freq2-freq3(im3)) * sigmam1
     ctm_c = 2.0d0 * pi * bos_c * prexp * Exp(-(dom_c*dom_c))
!
     bos_c1 = bos1 * bos3b(im3) *(bos2+1.0d0)
     dom_c1 =(freq1+freq3b(im3)-freq2) * sigmam1
     ctm_c1 = 2.0d0 * pi * bos_c1 * prexp * Exp(-(dom_c1*dom_c1))
!
!
     bos_d =(bos1+1.0d0) * bos2 * bos3b(im3)
     dom_d =(freq1-freq2-freq3b(im3)) * sigmam1
     ctm_d = 2.0d0 * pi * bos_d * prexp * Exp(-(dom_d*dom_d))
!
     IF(dabs(freq3(im3)) > eps8) THEN
        freqtot = freq1 * freq2 * freq3(im3)! if freqtot lt 10cm-1 do not calculate anything!
        norm1 = norm / freqtot
     ELSE
        norm1 = 0.0d0
     END IF
     IF(dabs(freq3b(im3)) > eps8) THEN
        freqtotb = freq1 * freq2 * freq3b(im3)! if freqtot lt 10cm-1 do not calculate anything!
        norm1b = norm / freqtotb
     ELSE
        norm1b = 0.0d0
     END IF
     !  if(abs(freqtot).gt.eps8) THEN
!
     sum_c  = ctm_c  * norm1  * d3mm(im3)
     sum_c1 = ctm_c1 * norm1b * d3mm_1(im3)
     sum_d  = ctm_d  * norm1b * d3mm_1(im3)
!
     sumPF_a = sumPF_a + sum_c1 + sum_d - sum_c ! A
!
!
  ENDDO
!
  sumPF_a = sumPF_a + wimp
! costruisco qui la matrice M come Q_{ij} * \delta_{ij} - A_{ij}
!
  IF((iq1_ == ipq2_) .and.(pim == pim2)) THEN
     sumPF_a =(Q_0_-sumPF_a)
  ELSE
     sumPF_a = - sumPF_a
  END IF
!
  !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_modes3a3
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE dyndiag1b(nat3, dyn, w2)
  !--------------------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On output: w2 = energies, dyn = displacements
  !
  USE constants, ONLY: DP
  IMPLICIT NONE
  ! input
  Integer nat3
  Complex(DP) dyn(nat3, nat3)
  Real(DP) w2(nat3)
  ! local
  Integer :: im
!  logical :: on
!  if(on) THEN
  CALL cdiagh2b(nat3, dyn, w2)
  DO im = 1, nat3
     w2(im) = Sqrt(Abs(w2(im)))
! IF(w2(im) < 0.0d0) w2(im) = -w2(im)
  ENDDO
!  end if
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE dyndiag1b
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE cdiagh2b(n, h, e)
  !--------------------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
  USE constants, ONLY: DP
  USE cdiagh_work
  IMPLICIT NONE
 !
  Integer :: n
  Complex(Kind=DP) :: h(n, n)! in input: matrix to be diagonalized
                          ! in output: eigenvectors
  Real(Kind=DP) :: e(n)! eigenvalues
!
  Logical, Save :: first = .True.
  Integer :: ILAENV, nb, info
!
  IF(first) THEN
     nb = ILAENV(1, 'ZHETRD', 'U', n,-1,-1,-1)
     IF(nb < 1) nb = Max(1, n)
     IF(nb == 1 .Or. nb .Ge. n) THEN
        lwork = 2 * n - 1
     ELSE
        lwork =(nb+1) * n
     END IF
!
! allocate workspace
     ALLOCATE(work(lwork))
     ALLOCATE(rwork(3*n-2))
     n_first = n
  END IF
  first = .False.
!
  IF(n_first /= n) CALL errore('cdiagh2b', ' n has chaged', 1)
  CALL ZHEEV('V', 'U', n, h, n, e, work, lwork, rwork, info)
  CALL errore('cdiagh2', 'info =/= 0', Abs(info))
!
!
 !
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE cdiagh2b
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE change_m2r_index
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP
  USE common_variables, ONLY : &
    mat2, mat2_in, &
    nRbig2_in, iRbig2_in, nRbig2, iRbig2, nRbig2t, Rbig2, &
    ic_out, ic_med, ic_inn, nat3
  IMPLICIT NONE
  Logical, Allocatable :: ldone(:)
  Integer :: iR, jr
  Integer imin(3), imax(3), ii, ic, i1, i2, i3
!
  DO ic = 1, 3
     imin(ic) = iRbig2_in(ic, 1)
     imax(ic) = iRbig2_in(ic, 1)
     DO ii = 1, nRbig2_in
        imin(ic) = Min(iRbig2_in(ic, ii), imin(ic))
        imax(ic) = Max(iRbig2_in(ic, ii), imax(ic))
     ENDDO
     nRbig2(ic) = imax(ic) - imin(ic) + 1
  ENDDO
  nRbig2t = nRbig2(1) * nRbig2(2) * nRbig2(3)
  ALLOCATE(iRbig2(3, nRbig2t))
  ALLOCATE(Rbig2(3, nRbig2t))
  ALLOCATE(mat2(nat3, nat3, nRbig2t))
  ii = 0
  DO i1 = imin(ic_out), imax(ic_out)
     DO i2 = imin(ic_med), imax(ic_med)
        DO i3 = imin(ic_inn), imax(ic_inn)
           ii = ii + 1
           iRbig2(ic_out, ii) = i1
           iRbig2(ic_med, ii) = i2
           iRbig2(ic_inn, ii) = i3
        ENDDO
     ENDDO
  ENDDO
!
  ALLOCATE(ldone(nRbig2_in))
  mat2 = 0.d0
  ldone = .False.
  DO iR = 1, nRbig2t
     DO jr = 1, nRbig2_in
        IF(iRbig2(1, iR) == iRbig2_in(1, jr) .and. iRbig2(2, iR) == iRbig2_in(2, jr) .and. iRbig2(3, iR) == iRbig2_in(3, jr)) &
       & THEN
           mat2(:, :, iR) = mat2_in(:, :, jr)
           IF(ldone(jr)) CALL errore('change_mr2_index', 'something wrong 1', 1)
           ldone(jr) = .True.
           Go To 120
        END IF
     ENDDO
     mat2(:, :, iR) = 0.d0
120      Continue
  ENDDO
  DO jr = 1, nRbig2_in
     IF( .Not. ldone(jr)) CALL errore('change_nr2_index', 'something wrong 2', 1)
  ENDDO
!
  Rbig2(:, :) = dfloat(iRbig2(:, :))
!
  DEALLOCATE(iRbig2_in)
  DEALLOCATE(mat2_in)
  DEALLOCATE(ldone)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE change_m2r_index
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE change_m3q_index(nat3, nr1, nr2, mat3)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP
  IMPLICIT NONE
  Integer :: nat3, nr1, nr2
  Complex(DP) :: mat3(nat3, nat3, nat3, nr1, nr2)
  Integer :: im1, im2, im3, ir1, ir2
  Complex(DP), Allocatable :: aux(:, :, :, :, :)
  ALLOCATE(aux(nat3, nat3, nat3, nr1, nr2))
!
  aux = mat3
  DO ir2 = 1, nr2
     DO ir1 = 1, nr1
        DO im3 = 1, nat3
           DO im2 = 1, nat3
              DO im1 = 1, nat3
                 mat3(im2, im3, im1, ir1, ir2) = aux(im1, im2, im3, ir1, ir2)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE change_m3q_index
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE write_header
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP
  USE common_variables, ONLY : ispeed2, ispeed3, nq1, deltaq1, nq2, deltaq2
  IMPLICIT NONE
!
!
  Write(*,*) ' ispeed2:', ispeed2
  Write(*,*) ' ispeed3:', ispeed3
  Write(*, '(a,3i5,3f12.6)') '  external  grid:', nq1, deltaq1
  Write(*, '(a,3i5,3f12.6)') '  internal  grid:', nq2, deltaq2
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE write_header
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE initializations
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, tpi
  USE common_variables, ONLY : &
   T0, tinit_0, tinit_1, ttot_0, &
   tphase, tsetup2, tsetup3, tdindiag, tbose, td3mm, tmodes, &
   deltaq2, nRbig31, Rbig31, ispeed3, &
   nq2, ic_out, ic_med, ic_inn, &
   phasem3_q2_, phasem3_q2_init
  
  IMPLICIT NONE
  Integer :: ic, iR, jr
  Real(DP) :: arg
  Real(DP), EXTERNAL :: nanosec
!
  IF(ispeed3 == 3) THEN
     DO ic = 1, 3
        DO jr = 1, nRbig31(ic)
           IF(ic == ic_out) iR = 1 +(jr-1) * nRbig31(ic_inn) * nRbig31(ic_med)
           IF(ic == ic_med) iR = 1 +(jr-1) * nRbig31(ic_inn)
           IF(ic == ic_inn) iR = jr
           arg = tpi * Rbig31(ic, iR) / dfloat(nq2(ic))
           phasem3_q2_(jr, ic) = cmplx(Cos(arg),-Sin(arg))
           arg = tpi * Rbig31(ic, iR) * deltaq2(ic)
           phasem3_q2_init(jr, ic) = cmplx(Cos(arg),-Sin(arg))
        ENDDO
        phasem3_q2_init(:, ic) = phasem3_q2_init(:, ic) * CONJG(phasem3_q2_(:, ic))
!
     ENDDO
  END IF
!
  tphase = 0.d0
  tsetup2 = 0.d0
  tsetup3 = 0.d0
  tdindiag = 0.d0
  tbose = 0.0d0
  td3mm = 0.0d0
  tmodes = 0.0d0
!
  tinit_1 = nanosec(T0) - tinit_0
  ttot_0 = nanosec(T0)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE initializations
!--------------------------------------------------------------------------------
SUBROUTINE setup(isw)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, tpi, cone, czero, eps8
  USE common_variables
  USE mpi_base
  IMPLICIT NONE 
 integer :: isw

  integer :: ir, itmp,i,j,im,ix,jm
  real(DP) :: arg

  INTEGER :: nequiv1, nequiv2, nequiv3, nequiv3b
!  LOGICAL :: lklappq3, lklappq3b
  REAL(DP) :: q1d_rf(3), q2d_rf(3), q3d_rf(3), q3db_rf(3),qsum(3),mod1,mod2
  INTEGER :: iklapp1,iklapp2

 on=.true.
 select case (isw)
 case(1)

    if (ispeed2.eq.0) then
        
        freq1(:)= freq1_(:,iq1_)
        do ix=1,3
           zz1_x(:,:,ix) = zz1n(:,:,iq1_,ix)
           do im=1,nat3 
              if(dabs(freq1(im)).lt.eps8) then 
                 zz1_x(:,im,ix)=0.0d0
              !   if((nrank.eq.1).and.(iq1_.ne.1))  write(*,*) 'freq1  ',iq1_, q1d ,im ,freq1(im)
              end if
           end do
        end do
     end if



!    if (ispeed2.eq.1) then
!       T1 = nanosec(T0)
!       !        call setupmat(nat,nRbig2t,q1d,zz1,Rbig2,mat2)
!       call setupmat(nat,nRbig2t,q1d,zz1_x(1,1,1),Rbig2,mat2)
!       tsetup2  = tsetup2  + nanosec(T0) - T1
!       T1 = nanosec(T0)
!      CALL dyndiag1b(nat3,zz1_x(1,1,1),freq1)
!      do im=1,nat3 
!         if(dabs(freq1(im)).lt.eps8) zz1_x(:,im,1)=0.0d0
!      end do
!      do ix=2,3
!         zz1_x(:,:,ix)=zz1_x(:,:,1)
!      end do
!      tdindiag = tdindiag + nanosec(T0) - T1
!   end if
    if (ispeed3.eq.1) then
       T1 = nanosec(T0)
       call setupmat3a (nat33,nRbig31t,nRbig32,q1d,Rbig32,mat3_rr,mat3_r)
       tsetup3  = tsetup3  + nanosec(T0) - T1
    endif
    
    
    
    if (ispeed3.eq.2 .or. ispeed3.eq.3) then
       T1 = nanosec(T0) 
       do ir = 1, nRbig32
          arg = tpi*(q1d(1)*Rbig32(1,ir)+q1d(2)*Rbig32(2,ir)+q1d(3)*Rbig32(3,ir))
          phasem3_q1t(ir) = cmplx(cos(arg),-sin(arg))
       enddo
       tphase  = tphase + nanosec(T0) - T1
    endif
!    if (ispeed3.eq.2 ) then
!       T1 = nanosec(T0)
!      if (.not.ltest) then
!      phasem3_q2t(:) = phasem3_q2t_init(:)  ! DA ELIMINARE
!      else
!      phasem3_q2(:,ic_out) = phasem3_q2_init(:,ic_out)
!      endif
!      tphase  = tphase + nanosec(T0) - T1
!      T1 = nanosec(T0)
!      itmp = nat33*nRbig31t
!      call ZGEMV('N',itmp,nRbig32,cone,mat3_rr,itmp,phasem3_q1t,1,czero,mat3_r,1)
!      tsetup3  = tsetup3  + nanosec(T0) - T1
!
!   endif
    if (ispeed3.eq.3) then
        T1 = nanosec(T0)
        phasem3_q2(:,ic_out) = phasem3_q2_init(:,ic_out)
        tphase  = tphase + nanosec(T0) - T1 
        T1 = nanosec(T0)
        itmp = nat33*nRbig31t
        call ZGEMV('N',itmp,nRbig32,cone,mat3_rr,itmp,phasem3_q1t,1,czero,mat3_r,1)
        
        tsetup3  = tsetup3  + nanosec(T0) - T1
     endif
     


 
  case(2)
        
!     if (ispeed3.eq.2 ) then
!        T1 = nanosec(T0)   
!        if (.not.ltest) then
!           phasem3_q2t(:) = phasem3_q2t(:) * phasem3_q2t_(:,ic_out) !!!!DA ELIM
!       else
!       phasem3_q2(:,ic_out)  = phasem3_q2(:,ic_out) * phasem3_q2_(:,ic_out)
!       phasem3_q2(:,ic_med) = phasem3_q2_init(:,ic_med)
!       endif
!       tphase  = tphase + nanosec(T0) - T1
!    endif
     if (ispeed3.eq.3) then
        T1 = nanosec(T0)
        phasem3_q2(:,ic_out)  = phasem3_q2(:,ic_out) * phasem3_q2_(:,ic_out)
        phasem3_q2c(:,ic_out) = conjg(phasem3_q2(:,ic_out))  !
        tphase  = tphase + nanosec(T0) - T1
        T1 = nanosec(T0)
        itmp = nat33*nRbig31(ic_inn)*nRbig31(ic_med)
        call ZGEMV('N',itmp,nRbig31(ic_out),cone,mat3_r,itmp,phasem3_q2(1,ic_out),1,czero,mat3_rzy,1)
        !GF
        call ZGEMV('N',itmp,nRbig31(ic_out),cone,mat3_r,itmp,phasem3_q2c(1,ic_out),1,czero,mat3_rzyB,1)
        tsetup3  = tsetup3  + nanosec(T0) - T1
        T1 = nanosec(T0)
        phasem3_q2(:,ic_med) = phasem3_q2_init(:,ic_med)
        tphase  = tphase + nanosec(T0) - T1 
        
     endif
  case(3)
  
!    if (ispeed3.eq.2) then
!       T1 = nanosec(T0)
!       if (.not.ltest) then
!       phasem3_q2t(:) = phasem3_q2t(:) * phasem3_q2t_(:,ic_med)  ! DA ELIMINARE
!       else
!       phasem3_q2(:,ic_med)  = phasem3_q2(:,ic_med) * phasem3_q2_(:,ic_med)
!       phasem3_q2(:,ic_inn) = phasem3_q2_init(:,ic_inn)
!       endif
!       tphase  = tphase + nanosec(T0) - T1
!    endif
     
     if (ispeed3.eq.3) then
        T1 = nanosec(T0)
        phasem3_q2(:,ic_med)  = phasem3_q2(:,ic_med) * phasem3_q2_(:,ic_med)
        phasem3_q2c(:,ic_med) = conjg(phasem3_q2(:,ic_med))  !
        tphase  = tphase + nanosec(T0) - T1
        
        T1 = nanosec(T0)
        itmp = nat33*nRbig31(ic_inn)
        call ZGEMV('N',itmp,nRbig31(ic_med),cone,mat3_rzy,itmp,phasem3_q2(1,ic_med),1,czero,mat3_rz,1)
        
        call ZGEMV('N',itmp,nRbig31(ic_med),cone,mat3_rzyB,itmp,phasem3_q2c(1,ic_med),1,czero,mat3_rzB,1)
        tsetup3  = tsetup3  + nanosec(T0) - T1
        
        
        T1 = nanosec(T0)
        phasem3_q2(:,ic_inn) = phasem3_q2_init(:,ic_inn)
        tphase  = tphase + nanosec(T0) - T1
     endif


  case(4)

     if (ispeed2.eq.0) then
        
!        freq2(:)= freq1_(:,ipq2_)
        freq2(:)= freq2_(:,ipq2_)
        do ix=1,3
!           zz2_x(:,:,ix) = zz1n(:,:,ipq2_,ix)
           zz2_x(:,:,ix) = zz2n(:,:,ipq2_,ix)
           do im=1,nat3 
              if(dabs(freq2(im)).lt.eps8) then
                 zz2_x(:,im,ix)=0.0d0
             !    if((nrank.eq.1).and.(ipq2_.ne.1)) write(*,*) ipq2_,im,freq2(im)
              end if
           end do
        end do

!        freq3(:)= freq1_(:,ipq3_)
        freq3(:) =  freq3_(:,ipq3_)
        do ix=1,3
           !           zz3_x(:,:,ix) = zz1n(:,:,ipq3_,ix)
           zz3_x(:,:,ix) = zz3n(:,:,ipq3_,ix)
           do im=1,nat3 
              if(dabs(freq3(im)).lt.eps8)then
                 zz3_x(:,im,ix)=0.0d0
              !   if((nrank.eq.1).and.(ipq3_.ne.1)) write(*,*)'freq3 ', ipq3_,q3d,im,freq3(im)
              end if
           end do
        end do
        
!       freq3b(:)= freq1_(:,ipq3b_)
       freq3b(:)= freq3b_(:,ipq3b_)
        do ix=1,3
!           zz3b_x(:,:,ix) = zz1n(:,:,ipq3b_,ix)
           zz3b_x(:,:,ix) = zz3bn(:,:,ipq3b_,ix)
           do im=1,nat3 
              if(dabs(freq3b(im)).lt.eps8) then
                 zz3b_x(:,im,ix)=0.0d0
          !       if((nrank.eq.1) .and.(ipq3b_.ne.1)) write(*,*)'freq3b', ipq3b_,q3db,im,freq3b(im)
              end if
            end do
        end do

     end if


!    if ( ispeed2.eq.1 ) then
!       T1 = nanosec(T0)
!       call setupmat(nat,nRbig2t,q2d,zz2_x(1,1,1),Rbig2,mat2)
!       tsetup2  = tsetup2  + nanosec(T0) - T1
!       T1 = nanosec(T0)
!       CALL dyndiag1b(nat3,zz2_x(1,1,1),freq2)
!       do im=1,nat3
!          if(dabs(freq2(im)).lt.eps8) zz2_x(:,im,1)=0.0d0
!       end do
!       do ix=2,3
!          zz2_x(:,:,ix) = zz2_x(:,:,1)
!       end do
!       tdindiag = tdindiag + nanosec(T0) - T1
!       T1 = nanosec(T0)
!       call setupmat(nat,nRbig2t,q3d,zz3_x(1,1,1),Rbig2,mat2)
!       call setupmat(nat,nRbig2t,q3db,zz3b_x(1,1,1),Rbig2,mat2)
!       tsetup2  = tsetup2  + nanosec(T0) - T1
!       T1 = nanosec(T0)
!       CALL dyndiag1b(nat3,zz3_x(1,1,1),freq3)
!       do im=1,nat3
!          if(dabs(freq3(im)).lt.eps8) zz3_x(:,im,1)=0.0d0
!       end do
!       do ix=2,3
!          zz3_x(:,:,ix) = zz3_x(:,:,1)
!       end do
!       CALL dyndiag1b(nat3,zz3b_x(1,1,1),freq3b)
!       
!       do im=1,nat3
!          if(dabs(freq3b(im)).lt.eps8) zz3b_x(:,im,1)=0.0d0
!       end do
!       do ix=2,3
!          zz3b_x(:,:,ix) = zz3b_x(:,:,1)
!       end do
!    end if

     if ( ispeed3.eq.1 ) then
        tdindiag = tdindiag + nanosec(T0) - T1
        T1 = nanosec(T0)
        call setupmat3b (nat33,nRbig31t,q2d,Rbig31,mat3_r,mat3q)
        call setupmat3b(nat33,nRbig31t, mq2d,Rbig31,mat3_r,mat3qb)  
        
        tsetup3  = tsetup3  + nanosec(T0) - T1
        
     endif
     
     
!    if ( ispeed3.eq.2 ) then
!       if (.not.ltest) then
!       phasem3_q2t(:) = phasem3_q2t(:) * phasem3_q2t_(:,ic_inn) ! DA ELIMINARE
!       else
!       phasem3_q2(:,ic_inn)  = phasem3_q2(:,ic_inn) * phasem3_q2_(:,ic_inn)
!       phasem3_q2t(:) = phasem3_q2(:,ic_inn) * phasem3_q2(:,ic_med) * phasem3_q2(:,ic_out)
!       endif
!       phasem3_q2tc(:) = conjg( phasem3_q2t(:))
!       
!       T1 = nanosec(T0)
!       call ZGEMV('N',nat33,nRbig31t,cone,mat3_r,nat33,phasem3_q2t,1,czero,mat3q,1)
!       
!       call ZGEMV('N',nat33,nRbig31t,cone,mat3_r,nat33,phasem3_q2tc,1,czero,mat3qb,1)
!       tsetup3  = tsetup3  + nanosec(T0) - T1
!    endif
     
     if (ispeed3.eq.3) then
        
        T1 = nanosec(T0)
        phasem3_q2(:,ic_inn)  = phasem3_q2(:,ic_inn) * phasem3_q2_(:,ic_inn)
        phasem3_q2c(:,ic_inn) = conjg(phasem3_q2(:,ic_inn))  !
        
        tphase  = tphase + nanosec(T0) - T1
        
        
        T1 = nanosec(T0)
        call ZGEMV('N',nat33,nRbig31(ic_inn),cone,mat3_rz,nat33,phasem3_q2(1,ic_inn),1,czero,mat3q,1)
        call ZGEMV('N',nat33,nRbig31(ic_inn),cone,mat3_rzB,nat33,phasem3_q2c(1,ic_inn),1,czero,mat3qb,1)
        !
        tsetup3  = tsetup3  + nanosec(T0) - T1
     endif


     T1=nanosec(T0)
     call calc_d3mm (d3mmx,mat3q,zz1_x,zz2_x,zz3_x,nat3,nat32,nat33)
!     d3mm(:,:,:) =d3mmx(:,:,:,1)
     
     zz2b_c_x(:,:,:) =CONJG(zz2_x(:,:,:))
     call calc_d3mm (d3mmx_1,mat3qb,zz1_x,zz2b_c_x,zz3b_x,nat3,nat32,nat33)
     
!     d3mm_1(:,:,:) =d3mmx_1(:,:,:,1)
     

     
     td3mm = td3mm + nanosec(T0) - T1
     
  end select
  

end subroutine setup
!
!--------------------------------------------------------------------------------
SUBROUTINE loop_on_qpoints
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, pi, ry_to_cm1, K_BOLTZMANN_RY, RY_TO_WATT, BOHR_TO_M
  USE common_variables, ONLY : &
     bos1, bos2, bos3, bos3b, &
     freq1, freq2, freq3, freq3b, &
     ntemp, tempm1,temp, sigmam1, alpha, &
     nq1tot, nq1loc, iq1_loc, nq1, q1d, &
     iq1_, deltaq1, nq2tot, mq2d, nq2, q2d, deltaq2, q3d, q3db, &
     imq2_, ipq2_, ipq3_, ipq3b_, &
     Q_0, Q_0rad, Q_0radm1, F_0, &
     ic_inn, ic_med, ic_out, &
     nxcryst, xcryst, &
     nat, nat3, Lcasm1, lbordo, const_iso, &
     T1, T0, ttot_0, tsetup2, tmodes, tinit_1, tdindiag, td3mm, tbose, ttot_1, &
     tsetup3, tphase, lambdiso, tcondc_, &
     zz2_x, zz1_x, velph, lambd, &
     d3mmx, d3mmx_1, ierr, bg
     
  USE mpi_base
  IMPLICIT NONE
!
  Integer :: itemp, im, ic, it, ix, iq_init(3),  isw1, isw2, iqq1, km, im2, ia, ix_
  Real(DP) :: lambdloc(nat3, ntemp, 3, nq1tot), lambdloc2(nat3, ntemp, 3, nq1tot)
  Real(DP) :: Q_0loc(nat3, nq1tot, 3, ntemp)!,F_0loc(nat3,nq1tot,3,ntemp)
  Real(DP) :: xq(3,nq1tot)
!
  Real(DP) :: zz_12(nat3, nat3)
!
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out, iq2_out, iq2_med, iq2_inn
  Integer, Target :: iq1(3), iq2(3)
  Integer :: ii
  Complex(DP) :: workc
  Real(DP),EXTERNAL :: nanosec
!
  iq1_out => iq1(ic_out)
  iq1_med => iq1(ic_med)
  iq1_inn => iq1(ic_inn)
!
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
  F_0 = 0.0d0
!  F_0loc = 0.0d0
  tcondc_(:, :, :) = 0.0d0
  isw1 = - 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  isw2 = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
!
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!            EXTERNAL loop
!--------------------------------------------------------------------------------
  ! IF you change the order of these loop; you should change calc_vel accordingly
!
  lambdloc = 0.0d0
  lambdloc2 = 0.0d0
  Q_0loc = 0.0d0
!
  DO ii = 1, nq1loc
!do iq1_=1,nq1(1)* nq1(2)* nq1(3)
     iq1_ = iq1_loc(ii)
!
!   write(*,*) iq1_,nq1loc
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
!
!
     CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw1)
     iq1(:) = iq_init(:)
     q1d(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
     xq(:,iq1_) = q1d
!      IF(nrank == 1) write(10,*) iq_init(:), iq1_
     CALL setup(1)
!--------------------------------------------------------------------------------
!            INTERNAL loop
!--------------------------------------------------------------------------------
     DO iq2_out = 0, nq2(ic_out) - 1
        CALL setup(2)
        DO iq2_med = 0, nq2(ic_med) - 1
           CALL setup(3)
           DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
              q2d(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
              iq_init(:) = iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw2)
!
              mq2d(:) = - q2d(:)
              iq_init(:) = - iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq2_, isw2)
!
              q3d(:) = - q1d(:) - q2d(:)
              iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) - 2 * Nint(deltaq1(:)*nq2(:)) - Nint(deltaq2(:)*nq2(:))

              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw2)
!
              q3db(:) = q2d(:) - q1d(:)
              iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint(deltaq1(:)*nq2(:))

              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw2)
!--------------------------------------------------------------------------------
              CALL setup(4)
!
              DO itemp = 1, ntemp
                 T1 = nanosec(T0)
                 DO im = 1, 3 * nat
                    bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(itemp))-1.d0)
                    bos2(im) = 1.d0 /(Exp(freq2(im)*tempm1(itemp))-1.d0)
                    bos3(im) = 1.d0 /(Exp(freq3(im)*tempm1(itemp))-1.d0)
                    bos3b(im) = 1.d0 /(Exp(freq3b(im)*tempm1(itemp))-1.d0)
                 ENDDO
                 tbose = tbose + nanosec(T0) - T1

                 T1 = nanosec(T0)
                 DO ix_ = 1, nxcryst
                    ix = xcryst(ix_)
                    zz_12(:, :) = 0.0d0
                    DO im = 1, nat3
                       DO im2 = 1, 3 * nat !
                          DO ia = 1, nat
!
                             workc =(0.0d0, 0.0d0)
                             DO ic = 1, 3
!
                                km = ic + 3 *(ia-1)
                                workc = workc + CONJG(zz1_x(km, im, ix)) * zz2_x(km, im2, ix)
                             ENDDO
!
                             zz_12(im2, im) = zz_12(im2, im) + CONJG(workc) * workc
!
                          ENDDO
                       ENDDO
!
                       CALL sum_modes2iso(nat3, nq2tot, const_iso, sigmam1(itemp), freq1(im), freq2, freq3, freq3b, bos1(im), &
                      & bos2, bos3, bos3b, d3mmx(1, 1, im, ix), d3mmx_1(1, 1, im, ix), zz_12(1, im), lambdloc(im, itemp, ix, &
                      & iq1_), Q_0loc(im, iq1_, ix, itemp), lambdloc2(im, itemp, ix, iq1_))
!
                    ENDDO
                 ENDDO
                 tmodes = tmodes + nanosec(T0) - T1
       !   enddo
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!
!
!
     IF(lbordo) THEN
        DO it = 1, ntemp
           DO im = 1, 3 * nat
              bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(it))-1.d0)
           ENDDO
           DO ix_ = 1, nxcryst
              ix = xcryst(ix_)
              DO im = 1, nat3
                 Q_0loc(im, iq1_, ix, it) = Q_0loc(im, iq1_, ix, it) + bos1(im) *(1.0d0+bos1(im)) &
                                          * Abs(velph(ix, im, iq1_)) * Lcasm1 * 2.0d0
 !              Q_0loc(im,iq1_,ix,it) =  Q_0loc(im,iq1_,ix,it) + bos1(im)*(1.0d0+bos1(im))* 1096892.13894d-8 * Lcasm1 *2.0d0
!
              ENDDO
           ENDDO
        ENDDO
     END IF
!
  ENDDO
  CALL MPI_ALLREDUCE(lambdloc, lambd, nat3*ntemp*3*nq1tot, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(lambdloc2, lambdiso, nat3*ntemp*3*nq1tot, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(Q_0loc, Q_0, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
!CALL MPI_ALLREDUCE(F_0loc,F_0,nat3*nq1tot*3*ntemp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!
  IF(.True.) THEN
!      IF(nrank == 1) THEN
! !!   write(*,*) x_iso,fracmass2,const_iso,lambdiso(1,1,1,1)
!         Open(Unit=222, File='LW.dat', Position='append')
!         Open(Unit=223, File='LWiso.dat', Position='append')
! !
! !
!         DO iq1_ = 1, nq1tot
!            DO itemp = 1, ntemp
!               Write(222,*) 'TEMPERATURE', temp(itemp), iq1_
!               Write(223,*) 'TEMPERATURE', temp(itemp), iq1_
!    !      do im=1,nat3
!               DO ix = 1, 3
!                  Write(222, '(6e14.6)')(pi*lambd(im, itemp, ix, iq1_)*ry_to_cm1, im=1, nat3)
!                  Write(223, '(6e14.6)')(lambdiso(im, itemp, ix, iq1_)*ry_to_cm1, im=1, nat3)
! !!            WRITE(23,'(6e14.6)')(pi*lambdiso(im,itemp,isigma,iq1_)*ry_to_cm1,isigma=1,nsigma)
!               ENDDO
!    !      enddo
!               Write(222,*) '----------------------------------------------------------------'
!               Write(223,*) '----------------------------------------------------------------'
!            ENDDO
!         ENDDO
!         Close(Unit=222)
!         Close(Unit=223)
! !
     CALL cryst_to_cart(nq1tot, xq, bg, +1)

     IF(nrank == 1) THEN
!!   write(*,*) x_iso,fracmass2,const_iso,lambdiso(1,1,1,1)
!         Open(Unit=222, File='LW.dat', Position='append')
!         Open(Unit=223, File='LWiso.dat', Position='append')
!
!
        DO iq1_ = 1, nq1tot
           DO itemp = 1, ntemp
!               Write(1000+itemp,*) 'TEMPERATURE', temp(itemp), iq1_
!               Write(2000+itemp,*) 'TEMPERATURE', temp(itemp), iq1_
   !      do im=1,nat3
!               DO ix = 1, 3
                 Write(1000+itemp, '(3f12.6,99e20.10)') xq(:,iq1_), pi*lambd(:, itemp, 1, iq1_)*ry_to_cm1
                 Write(2000+itemp, '(3f12.6,99e20.10)') xq(:,iq1_), lambdiso(:, itemp, 1, iq1_)*ry_to_cm1
!!            WRITE(23,'(6e14.6)')(pi*lambdiso(im,itemp,isigma,iq1_)*ry_to_cm1,isigma=1,nsigma)
!               ENDDO
   !      enddo
              Write(222,*) '----------------------------------------------------------------'
              Write(223,*) '----------------------------------------------------------------'
           ENDDO
        ENDDO
!         Close(Unit=222)
!         Close(Unit=223)
!
     END IF
  END IF
  IF(nrank == 1) write(*,*) 'fine loop qpoints'
!
!
!
!
  DO iqq1 = 1, nq1tot
      DO ix_ = 1, nxcryst
          ix = xcryst(ix_)
          DO it = 1, ntemp
            DO im = 1, nat3
                IF(Q_0(im, iqq1, ix, it) /= 0.0d0) THEN
                  Q_0rad(im, iqq1, ix, it) = Q_0(im, iqq1, ix, it) ** alpha
                  Q_0radm1(im, iqq1, ix, it) = 1.0d0 / Q_0rad(im, iqq1, ix, it)
                ELSE
                  Q_0rad(im, iqq1, ix, it) = 0.0d0
                  Q_0radm1(im, iqq1, ix, it) = 0.0d0 ! equal to zero ONLY because anyway I never consider the contribution when Q_0 =0.0d0
                END IF
            ENDDO
          ENDDO
      ENDDO
!
     CALL calc_f0(iqq1)
!
     CALL sum_conduct(iqq1)
  ENDDO
!--------------------------------------------------------------------------------
!
  ttot_1 = nanosec(T0) - ttot_0
  ttot_0 = tdindiag + tsetup3 + tsetup2 + tphase + tbose + td3mm + tmodes
!
  IF(.False.) THEN
     Write(*, '(a,e24.12)') '  inititaliz time:', tinit_1
     Write(*, '(a,e24.12)') '       total time:', ttot_1
     Write(*, '(a,e24.12)') '  sum of partials:', ttot_0
     Write(*, '(a,e24.12,f10.2,a)') '           tphase:', tphase,(tphase/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup2:', tsetup2,(tsetup2/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup3:', tsetup3,(tsetup3/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '          dindiag:', tdindiag,(tdindiag/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             bose:', tbose,(tbose/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             d3mm:', td3mm,(td3mm/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '            modes:', tmodes,(tmodes/ttot_0) * 100.d0, '  %'
  END IF
!
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE loop_on_qpoints
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE loop_on_qpointsIT(iter)
  !--------------------------------------------------------------------------------
  !
  USE constants, ONLY: DP, pi, sqrtpi, BOHR_TO_M, eps8
  USE common_variables, ONLY : &
     bos1, bos2, bos3, bos3b, &
     freq1, freq2, freq3, freq3b, &
     ntemp, tempm1, sigmam1, &
     nq1tot, nq1loc, iq1_loc, nq1, q1d, &
     iq1_, deltaq1, nq2tot, mq2d, nq2, q2d, deltaq2, q3d, q3db, &
     imq2_, ipq2_, ipq3_, ipq3b_, &
     Q_0, F_0, &
     ic_inn, ic_med, ic_out, &
     nxcryst, xcryst, &
     nat, nat3, const_iso, &
     T1, T0, ttot_0, tsetup2, tmodes, tinit_1, tdindiag, td3mm, tbose, ttot_1, &
     tsetup3, tphase, &
     zz2_x, zz1_x, &
     d3mmx, d3mmx_1, ierr, &
     sum_IT, matrix, F_new, F_old
  USE mpi_base
  IMPLICIT NONE
!
!
  Integer :: im3, ii, iter, isw1, isw2, km, ix_
  Integer :: im, im2, ic, it, ix, iq_init(3)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out, iq2_out, iq2_med, iq2_inn
  Integer, Target :: iq1(3), iq2(3)
  Real(DP) :: matAm3, d3mm3(nat3), d3mm_13(nat3)
  Integer :: iqqmm1, iqqmm2
!
  Complex(DP) :: workc
  Real(DP) :: zz_12(nat3, nat3), norm_iso1
  Real(DP) :: bos_iso, dom_iso, ctm_iso, wimp_iso, prexp
!
  Integer :: ia
  Real(DP) :: freq1t(nat3, nq1tot), F_newloc(nat3, nq1tot, 3, ntemp)
  Real(DP) :: matrixloc(nat3, nq1tot, 3, ntemp), sum_ITloc(nat3, nq1tot, 3, ntemp)
  Real(DP),EXTERNAL :: nanosec
!
!
  iq1_out => iq1(ic_out)
  iq1_med => iq1(ic_med)
  iq1_inn => iq1(ic_inn)
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
!
!
  isw1 = - 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  isw2 = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
!--------------------------------------------------------------------------------
!for isotopic scattering contribution
  norm_iso1 = pi * const_iso * 0.5d0 /(dfloat(nq2tot))! qui a differenza che nel loop dei qpoints la moltiplicazione per freq1 non e'
                                                      !in norm_iso ma dopo all'interno del loop
!     prexp= sigmam1(1) / sqrtpi
!--------------------------------------------------------------------------------
  sum_ITloc(:, :, :, :) = 0.0d0
  matrixloc(:, :, :, :) = 0.0d0
  F_newloc(:, :, :, :) = 0.0d0
!--------------------------------------------------------------------------------
!            EXTERNAL loop
!--------------------------------------------------------------------------------
  ! IF you change the order of these loop; you should change calc_vel accordingly
!
!
  DO ii = 1, nq1loc
     iq1_ = iq1_loc(ii)
!
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
!
     CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw1)
     iq1(:) = iq_init(:)
     q1d(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
!     iq_init(:) = iq1(:)
!--------------------------------------------------------------------------------
!
!
!      sum_IT3(:,:,:) = 0.0d0
     CALL setup(1)
!
!--------------------------------------------------------------------------------
!            INTERNAL loop
!--------------------------------------------------------------------------------
     DO iq2_out = 0, nq2(ic_out) - 1
        CALL setup(2)
        DO iq2_med = 0, nq2(ic_med) - 1
           CALL setup(3)
           DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
              q2d(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
              iq_init(:) = iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw2)
!
              mq2d(:) = - q2d(:)
              iq_init(:) = - iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq2_, isw2)
!
!
              q3d(:) = - q1d(:) - q2d(:)
              iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) - 2 * Nint(deltaq1(:)*nq2(:)) - Nint(deltaq2(:)*nq2(:))

              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw2)
!
              q3db(:) = q2d(:) - q1d(:)
              iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint(deltaq1(:)*nq2(:))
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw2)
!--------------------------------------------------------------------------------
              CALL setup(4)
!
!
              DO im = 1, nat3
                 freq1t(im, iq1_) = freq1(im)
              ENDDO
!
!
! write(*,*)'imq3',imq3_
! irun = 0
              DO it = 1, ntemp
                 T1 = nanosec(T0)
                 prexp = sigmam1(it) / sqrtpi
                 DO im = 1, 3 * nat
                    bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(it))-1.d0)
                    bos2(im) = 1.d0 /(Exp(freq2(im)*tempm1(it))-1.d0)
                    bos3(im) = 1.d0 /(Exp(freq3(im)*tempm1(it))-1.d0)
                    bos3b(im) = 1.d0 /(Exp(freq3b(im)*tempm1(it))-1.d0)
                 ENDDO
                 tbose = tbose + nanosec(T0) - T1
                 T1 = nanosec(T0)
                 DO ix_ = 1, nxcryst
                    ix = xcryst(ix_)
                    zz_12(:, :) = 0.0d0
                    DO im = 1, nat3
                       iqqmm1 = im + nat3 *(iq1_-1)
!
                       DO im2 = 1, 3 * nat !
                          DO ia = 1, nat
!
                             workc =(0.0d0, 0.0d0)
                             DO ic = 1, 3
!
                                km = ic + 3 *(ia-1)
                                workc = workc + CONJG(zz1_x(km, im, ix)) * zz2_x(km, im2, ix)
                             ENDDO
!
                             zz_12(im2, im) = zz_12(im2, im) + CONJG(workc) * workc
!
                          ENDDO
                       ENDDO
!
!
                       DO im2 = 1, nat3
                          iqqmm2 = im2 + nat3 *(ipq2_-1)
!
                          DO im3 = 1, nat3
                             d3mm3(im3) = d3mmx(im2, im3, im, ix)
                             d3mm_13(im3) = d3mmx_1(im2, im3, im, ix)
                          ENDDO
!
                          matAm3 = 0.0d0
!--------------------------------------------------------------------------------
! isotopic scattering contribution
                          bos_iso = bos1(im) * bos2(im2) + 0.5d0 *(bos1(im)+bos2(im2))
                          dom_iso =(freq1(im)-freq2(im2)) * sigmam1(it)
                          ctm_iso = prexp * Exp(-(dom_iso*dom_iso))
                          wimp_iso = norm_iso1 * freq1(im) * freq2(im2) * bos_iso * ctm_iso * zz_12(im2, im)
!--------------------------------------------------------------------------------
!
!
                          CALL sum_modes3a2(nat3, nq2tot, sigmam1(it), freq1(im), freq2(im2), freq3, freq3b, bos1(im), &
                         & bos2(im2), bos3, bos3b, d3mm3, d3mm_13, matAm3,+1)
!
                          matAm3 = matAm3 + wimp_iso
!
                          sum_ITloc(im, iq1_, ix, it) = sum_ITloc(im, iq1_, ix, it) + matAm3 * F_old(im2, ipq2_, ix, it)
!
!
                          IF((iq1_ == ipq2_) .and.(im == im2)) THEN
                             matAm3 =(Q_0(im2, ipq2_, ix, it)-matAm3)
                          ELSE
                             matAm3 = - matAm3
                          END IF
!
                          matrixloc(im, iq1_, ix, it) = matrixloc(im, iq1_, ix, it) + matAm3 * F_old(im2, ipq2_, ix, it)
!
                       ENDDO
                    ENDDO
                 ENDDO
                 tmodes = tmodes + nanosec(T0) - T1
              ENDDO
!
           ENDDO
        ENDDO
     ENDDO
!
!
     DO it = 1, ntemp
        DO ix_ = 1, nxcryst
           ix = xcryst(ix_)
           DO im = 1, nat3
              IF(Q_0(im, iq1_, ix, it) /= 0.d0) THEN
                 F_newloc(im, iq1_, ix, it) = F_0(im, iq1_, ix, it) + sum_ITloc(im, iq1_, ix, it) / Q_0(im, iq1_, ix, it)
              END IF
           ENDDO
        ENDDO
     ENDDO
  ENDDO ! iq1 loop
!
  CALL MPI_ALLREDUCE(matrixloc, matrix, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(sum_ITloc, sum_IT, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(F_newloc, F_new, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
!
!--------------------------------------------------------------------------------
!
  ttot_1 = nanosec(T0) - ttot_0
  ttot_0 = tdindiag + tsetup3 + tsetup2 + tphase + tbose + td3mm + tmodes
!
  IF(.False.) THEN
     Write(*, '(a,e24.12)') '  inititaliz time:', tinit_1
     Write(*, '(a,e24.12)') '       total time:', ttot_1
     Write(*, '(a,e24.12)') '  sum of partials:', ttot_0
     Write(*, '(a,e24.12,f10.2,a)') '           tphase:', tphase,(tphase/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup2:', tsetup2,(tsetup2/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup3:', tsetup3,(tsetup3/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '          dindiag:', tdindiag,(tdindiag/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             bose:', tbose,(tbose/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             d3mm:', td3mm,(td3mm/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '            modes:', tmodes,(tmodes/ttot_0) * 100.d0, '  %'
  END IF
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE loop_on_qpointsIT
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE loop_on_qpointsITcg(iter, F_old1, sum_ITq1)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, pi, ry_to_cm1, K_BOLTZMANN_RY, RY_TO_WATT, BOHR_TO_M, eps8, sqrtpi
  USE common_variables, ONLY : &
     bos1, bos2, bos3, bos3b, &
     freq1, freq2, freq3, freq3b, &
     ntemp, tempm1, sigmam1, &
     nq1tot, nq1loc, iq1_loc, nq1, q1d, &
     iq1_, deltaq1, nq2tot, mq2d, nq2, q2d, deltaq2, q3d, q3db, &
     imq2_, ipq2_, ipq3_, ipq3b_, &
     Q_0, F_0, &
     ic_inn, ic_med, ic_out, &
     nxcryst, xcryst, &
     nat, nat3, const_iso, &
     T1, T0, ttot_0, tsetup2, tmodes, tinit_1, tdindiag, td3mm, tbose, ttot_1, &
     tsetup3, tphase, &
     zz2_x, zz1_x, &
     d3mmx, d3mmx_1, ierr, &
     Q_0radm1, sum_ITq2, Q_0rad, ltobedone
   USE mpi_base
  IMPLICIT NONE
!
!
  Integer :: im3, ii, iq, iter, isw1, isw2, km, ix_
  Integer :: im, im2, ic, it, ix, irun, iq_init(3)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out, iq2_out, iq2_med, iq2_inn
  Integer, Target :: iq1(3), iq2(3)
  Real(DP) :: matAm3, d3mm3(nat3), d3mm_13(nat3)
 ! real(DP) :: matA(nat3*nq2tot, nat3*nq1tot),autA(nat3*nq2tot),eigA(nat3*nq2tot, nat3*nq1tot)
  Integer :: iqqmm1, iqqmm2
!
  Real(DP) :: zz_12(nat3, nat3), norm_iso1
  Real(DP) :: bos_iso, dom_iso, ctm_iso, wimp_iso, prexp
  Real(DP) :: F_old1(nat3, nq1tot, 3, ntemp), sum_ITq1(nat3, nq1tot, 3, ntemp)
!
  Integer :: ia
  Real(DP) :: freq1t(nat3, nq1tot), sum_ITq2loc(nat3, nq1tot, 3, ntemp)
  Complex(DP) :: workc
  Real(DP),EXTERNAL :: nanosec
!
  iq1_out => iq1(ic_out)
  iq1_med => iq1(ic_med)
  iq1_inn => iq1(ic_inn)
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
!
  isw1 = - 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  isw2 = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
!--------------------------------------------------------------------------------
!for isotopic scattering contribution
  norm_iso1 = pi * const_iso * 0.5d0 /(dfloat(nq2tot))! qui a differenza che nel loop dei qpoints la moltiplicazione per freq1 non e'
                                                      !in norm_iso ma dopo all'interno del loop
 !    prexp= sigmam1(1) / sqrtpi
!--------------------------------------------------------------------------------
  sum_ITq2loc(:, :, :, :) = 0.0d0
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!            EXTERNAL loop
!--------------------------------------------------------------------------------
  ! IF you change the order of these loop; you should change calc_vel accordingly
!
  DO ii = 1, nq1loc
     iq1_ = iq1_loc(ii)
!
! do iq1_=1,nq1(1)*nq1(2)*nq1(3)
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
!
     CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw1)
     iq1(:) = iq_init(:)
     q1d(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
!     iq_init(:) = iq1(:)
!--------------------------------------------------------------------------------
!
!
!      sum_IT3(:,:,:) = 0.0d0
     CALL setup(1)
!
!--------------------------------------------------------------------------------
!            INTERNAL loop
!--------------------------------------------------------------------------------
     DO iq2_out = 0, nq2(ic_out) - 1
        CALL setup(2)
        DO iq2_med = 0, nq2(ic_med) - 1
           CALL setup(3)
           DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
              q2d(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
              iq_init(:) = iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw2)
!
              mq2d(:) = - q2d(:)
              iq_init(:) = - iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq2_, isw2)
!
              q3d(:) = - q1d(:) - q2d(:)
              iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) - 2 * Nint(deltaq1(:)*nq2(:)) - Nint(deltaq2(:)*nq2(:))

              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw2)

!
              q3db(:) = q2d(:) - q1d(:)
              iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint(deltaq1(:)*nq2(:))
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw2)
!--------------------------------------------------------------------------------
              CALL setup(4)
!
!
              DO im = 1, nat3
                 freq1t(im, iq1_) = freq1(im)
              ENDDO

              DO it = 1, ntemp
                 T1 = nanosec(T0)
                 prexp = sigmam1(it) / sqrtpi
                 DO im = 1, 3 * nat
                    bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(it))-1.d0)
                    bos2(im) = 1.d0 /(Exp(freq2(im)*tempm1(it))-1.d0)
                    bos3(im) = 1.d0 /(Exp(freq3(im)*tempm1(it))-1.d0)
                    bos3b(im) = 1.d0 /(Exp(freq3b(im)*tempm1(it))-1.d0)
                 ENDDO
                 tbose = tbose + nanosec(T0) - T1
                 T1 = nanosec(T0)
! do isig = 1, nsigma
                 DO ix_ = 1, nxcryst
                    ix = xcryst(ix_)
! Ritorni al
!           DO ix = 1,3
                    irun = ix +(it-1) * 3
                    IF(ltobedone(irun)) THEN
!
                       zz_12(:, :) = 0.0d0
!
                       DO im = 1, nat3
                          iqqmm1 = im + nat3 *(iq1_-1)
                !     if(dabs(freq1(im)).gt.eps8) THEN
!
                          DO im2 = 1, 3 * nat !
                             DO ia = 1, nat
!
                                workc =(0.0d0, 0.0d0)
                                DO ic = 1, 3
!
                                   km = ic + 3 *(ia-1)
                                   workc = workc + CONJG(zz1_x(km, im, ix)) * zz2_x(km, im2, ix)
                                ENDDO
!
                                zz_12(im2, im) = zz_12(im2, im) + CONJG(workc) * workc
!
                             ENDDO
                          ENDDO

                          DO im2 = 1, nat3
                             iqqmm2 = im2 + nat3 *(ipq2_-1)
                   !    if(dabs(freq2(im2)).gt.eps8) THEN
!
                             DO im3 = 1, nat3
                                d3mm3(im3) = d3mmx(im2, im3, im, ix)
                                d3mm_13(im3) = d3mmx_1(im2, im3, im, ix)
                             ENDDO
!
                             matAm3 = 0.0d0
                   !--------------------------------------------------------------------------------
                   ! isotopic scattering contribution
                             bos_iso = bos1(im) * bos2(im2) + 0.5d0 *(bos1(im)+bos2(im2))
                             dom_iso =(freq1(im)-freq2(im2)) * sigmam1(it)
                             ctm_iso = prexp * Exp(-(dom_iso*dom_iso))
                             wimp_iso = norm_iso1 * freq1(im) * freq2(im2) * bos_iso * ctm_iso * zz_12(im2, im)
                   !--------------------------------------------------------------------------------
!
                             CALL sum_modes3a3(nat3, nq2tot, sigmam1(it), iq1_, ipq2_, im, im2, freq1(im), freq2(im2), freq3, &
                            & freq3b, bos1(im), bos2(im2), bos3, bos3b, d3mm3, d3mm_13, matAm3, Q_0(im2, ipq2_, ix, it), &
                            & wimp_iso)

                             IF((iq1_ == 1) .and.(im == 1)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 2)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 3)) matAm3 = 0
!
                             IF((ipq2_ == 1) .and.(im2 == 1)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 2)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 3)) matAm3 = 0
!
                             IF(Q_0rad(im, iq1_, ix, it) /= 0.0d0 .and. Q_0rad(im2, ipq2_, ix, it) /= 0.0d0) THEN
                                matAm3 = Q_0radm1(im, iq1_, ix, it) * matAm3 * Q_0radm1(im2, ipq2_, ix, it)
                             ELSE
                                matAm3 = 0.0d0
                             END IF
!

                             sum_ITq2loc(im, iq1_, ix, it) = sum_ITq2loc(im, iq1_, ix, it) &
                                                             + matAm3 * F_old1(im2, ipq2_, ix,  it)!Ax
                   ! end if ! if dabs(freq2)
                   !
                          ENDDO
                ! end if ! if dabs(freq1)
                       ENDDO
         !   write(15,*)sum_IT(:,it)
                    END IF
!
                 ENDDO
                 tmodes = tmodes + nanosec(T0) - T1
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!
  ENDDO !iq1_
!
  CALL MPI_ALLREDUCE(sum_ITq2loc, sum_ITq2, nat3*nq1tot*3*ntemp, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)

!if(nq1tot.lt.80)THEN
!   CALL MPI_ALLREDUCE(sum_IT2loc,sum_IT2,nat3*nat3*3*ntemp*nq2tot*nq1tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if(iter.eq.1) THEN
!      CALL MPI_ALLREDUCE(matAloc,matA,nat3*nq1tot*nat3*nq2tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   end if
!end if
!
!
!
  sum_ITq1 = 0.0d0
  DO it = 1, ntemp
     DO ix_ = 1, nxcryst
        ix = xcryst(ix_)
        DO iq = 1, nq1tot
           DO im = 1, nat3
              IF(Q_0(im, iq, ix, it) /= 0.0d0) THEN
                 sum_ITq1(im, iq, ix, it) = sum_ITq2(im, iq, ix, it) &
                                           - F_0(im, iq, ix, it) &
                                            * Q_0(im, iq, ix, it) * Q_0radm1(im, iq, ix, it)! Ax-b
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!
!
!--------------------------------------------------------------------------------
!
  ttot_1 = nanosec(T0) - ttot_0
  ttot_0 = tdindiag + tsetup3 + tsetup2 + tphase + tbose + td3mm + tmodes
!
  IF(.False.) THEN
     Write(*, '(a,e24.12)') '  inititaliz time:', tinit_1
     Write(*, '(a,e24.12)') '       total time:', ttot_1
     Write(*, '(a,e24.12)') '  sum of partials:', ttot_0
     Write(*, '(a,e24.12,f10.2,a)') '           tphase:', tphase,(tphase/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup2:', tsetup2,(tsetup2/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '           setup3:', tsetup3,(tsetup3/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '          dindiag:', tdindiag,(tdindiag/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             bose:', tbose,(tbose/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '             d3mm:', td3mm,(td3mm/ttot_0) * 100.d0, '  %'
     Write(*, '(a,e24.12,f10.2,a)') '            modes:', tmodes,(tmodes/ttot_0) * 100.d0, '  %'
  END IF
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE loop_on_qpointsITcg
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE matrixA_times_h(F_old1, ggt1, hh1)
  !--------------------------------------------------------------------------------
!
  USE constants, ONLY: DP, pi, ry_to_cm1, K_BOLTZMANN_RY, RY_TO_WATT, BOHR_TO_M, eps8, sqrtpi
  USE common_variables, ONLY : &
    ic_inn, ic_med, ic_out, &
    iq1_, ipq3b_, ipq3_, ipq2_, imq2_, &
    nq1, nq1tot, nq1loc, nq2tot, nq2, q1d, q2d, q3d, mq2d, q3db, deltaq1, deltaq2, &
    ntemp, tempm1, sigmam1, &
    ierr, freq1, freq2, freq3b, freq3, const_iso, ltobedone, &
    bos1, bos2, bos3b, bos3, &
    d3mmx, d3mmx_1, zz1_x, zz2_x, &
    nxcryst, xcryst, iq1_loc, &
    tmodes, tbose, T1, T0, &
    nat, nat3, matrix, dimG, &
    Q_0radm1, Q_0rad, Q_0
  USE mpi_base
  IMPLICIT NONE
!
!
  Integer :: im3, ii, isw1, isw2, km, ix_ !,dimG
  Integer :: im, im2, ic, it, ix, irun, iq_init(3)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out, iq2_out, iq2_med, iq2_inn
  Integer, Target :: iq1(3), iq2(3)
  Real(DP) :: matAm3, d3mm3(nat3), d3mm_13(nat3)
  real(DP), intent(out) :: ggt1(dimG)
  Real(DP), intent(in)  :: hh1(dimG)
 ! real(DP) :: matA(nat3*nq2tot, nat3*nq1tot),autA(nat3*nq2tot),eigA(nat3*nq2tot, nat3*nq1tot)
  Integer :: count, count1
!
  Real(DP) :: zz_12(nat3, nat3), norm_iso1
  Real(DP) :: bos_iso, dom_iso, ctm_iso, wimp_iso, prexp
  Complex(DP) :: workc
  Real(DP) :: F_old1(nat3, nq1tot, 3, ntemp)
!
  Real(DP) :: ggtloc(dimG), matrixloc(nat3, nq1tot, 3, ntemp)
!
  Integer :: ia
  Real(DP) :: freq1t(nat3, nq1tot)
  Real(DP),EXTERNAL :: nanosec
!
!
!
!--------------------------------------------------------------------------------
!for isotopic scattering contribution
  norm_iso1 = pi * const_iso * 0.5d0 /(dfloat(nq2tot))! qui a differenza che nel loop dei qpoints la moltiplicazione per freq1 non e'
                                                      !in norm_iso ma dopo all'interno del loop
!--------------------------------------------------------------------------------
!
  iq1_out => iq1(ic_out)
  iq1_med => iq1(ic_med)
  iq1_inn => iq1(ic_inn)
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
!
  isw1 = - 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  isw2 = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
 !--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!            EXTERNAL loop
!--------------------------------------------------------------------------------
  ! IF you change the order of these loop; you should change calc_vel accordingly
  ! IF you change the order of these loop; you should change calc_vel accordingly
  ggtloc(:) = 0.0d0
  matrixloc(:, :, :, :) = 0.0d0
  DO ii = 1, nq1loc
     iq1_ = iq1_loc(ii)
! do iq1_=1,nq1(1)*nq1(2)*nq1(3)
!
!  do iq1_out = 0, nq1(ic_out) - 1
!  do iq1_med = 0, nq1(ic_med) - 1
!  do iq1_inn = 0, nq1(ic_inn) - 1
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
!
     CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw1)
     iq1(:) = iq_init(:)
     q1d(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
!     iq_init(:) = iq1(:)
!--------------------------------------------------------------------------------
!
!
!      sum_IT3(:,:,:) = 0.0d0
     CALL setup(1)
!
!--------------------------------------------------------------------------------
!            INTERNAL loop
!--------------------------------------------------------------------------------
     DO iq2_out = 0, nq2(ic_out) - 1
        CALL setup(2)
        DO iq2_med = 0, nq2(ic_med) - 1
           CALL setup(3)
           DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
              q2d(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
              iq_init(:) = iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw2)
!
              mq2d(:) = - q2d(:)
              iq_init(:) = - iq2(:)
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq2_, isw2)
!
              q3d(:) = - q1d(:) - q2d(:)
              iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) - 2 * Nint(deltaq1(:)*nq2(:)) - Nint(deltaq2(:)*nq2(:))
!        iq_init(:) = -iq1(:)*nq2(:)/nq1(:) -iq2(:)  - nint((2*deltaq1(:)+deltaq2(:))*nq2(:))
!        iq_init(:) = -iq1(:)*nq2(:)/nq1(:) -iq2(:)  - 3*nint(deltaq2(:)*nq2(:))
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw2)
! controlla che non serva e poi eliminalo:
!        mq3d(:) =  q1d(:) + q2d(:)
!        iq_init(:) =  iq1(:) + iq2(:)
!        call calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, imq3_ ,isw2)
!
              q3db(:) = q2d(:) - q1d(:)
              iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint(deltaq1(:)*nq2(:))
              CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw2)
!--------------------------------------------------------------------------------
              CALL setup(4)
!
!
              DO im = 1, nat3
                 freq1t(im, iq1_) = freq1(im)
              ENDDO
!
!
! write(*,*)'imq3',imq3_
! irun = 0
              DO it = 1, ntemp
                 T1 = nanosec(T0)
                 prexp = sigmam1(it) / sqrtpi
                 DO im = 1, 3 * nat
                    bos1(im) = 1.d0 /(Exp(freq1(im)*tempm1(it))-1.d0)
                    bos2(im) = 1.d0 /(Exp(freq2(im)*tempm1(it))-1.d0)
                    bos3(im) = 1.d0 /(Exp(freq3(im)*tempm1(it))-1.d0)
                    bos3b(im) = 1.d0 /(Exp(freq3b(im)*tempm1(it))-1.d0)
                 ENDDO
                 tbose = tbose + nanosec(T0) - T1
                 T1 = nanosec(T0)
          ! do isig = 1, nsigma
                 DO ix_ = 1, nxcryst
          ! Ritorni al
                    ix = xcryst(ix_)
                    irun = ix +(it-1) * 3
                    IF(ltobedone(irun)) THEN
                !
!
                       zz_12(:, :) = 0.0d0
                       DO im = 1, nat3
!
                   !     if(dabs(freq1(im)).gt.eps8) THEN
                          count1 = im + nat3 *(iq1_-1) + nat3 * nq1tot *(ix-1) + nat3 * nq1tot * 3 *(it-1)
!
!
                          DO im2 = 1, 3 * nat !
                             DO ia = 1, nat
!
                                workc =(0.0d0, 0.0d0)
                                DO ic = 1, 3
!
                                   km = ic + 3 *(ia-1)
                                   workc = workc + CONJG(zz1_x(km, im, ix)) * zz2_x(km, im2, ix)
                                ENDDO
!
                                zz_12(im2, im) = zz_12(im2, im) + CONJG(workc) * workc
!
                             ENDDO
                          ENDDO
!
!
                          DO im2 = 1, nat3
                             count = im2 + nat3 *(ipq2_-1) + nat3 * nq2tot *(ix-1) + nat3 * nq2tot * 3 *(it-1)
!
                      !    if(dabs(freq2(im2)).gt.eps8) THEN
!
                             DO im3 = 1, nat3
                                d3mm3(im3) = d3mmx(im2, im3, im, ix)
                                d3mm_13(im3) = d3mmx_1(im2, im3, im, ix)
                             ENDDO
!
                             matAm3 = 0.0d0
                      !--------------------------------------------------------------------------------
                      ! isotopic scattering contribution
                             bos_iso = bos1(im) * bos2(im2) + 0.5d0 *(bos1(im)+bos2(im2))
                             dom_iso =(freq1(im)-freq2(im2)) * sigmam1(it)
                             ctm_iso = prexp * Exp(-(dom_iso*dom_iso))
                             wimp_iso = norm_iso1 * freq1(im) * freq2(im2) * bos_iso * ctm_iso * zz_12(im2, im)
                      !--------------------------------------------------------------------------------
!
                             CALL sum_modes3a3(nat3, nq2tot, sigmam1(it), iq1_, ipq2_, im, im2, freq1(im), freq2(im2), freq3, &
                            & freq3b, bos1(im), bos2(im2), bos3, bos3b, d3mm3, d3mm_13, matAm3, Q_0(im2, ipq2_, ix, it), &
                            & wimp_iso)
!
!
                             IF((iq1_ == 1) .and.(im == 1)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 2)) matAm3 = 0
                             IF((iq1_ == 1) .and.(im == 3)) matAm3 = 0
!
                             IF((ipq2_ == 1) .and.(im2 == 1)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 2)) matAm3 = 0
                             IF((ipq2_ == 1) .and.(im2 == 3)) matAm3 = 0
!
                             IF(Q_0rad(im, iq1_, ix, it) /= 0.0d0 .and. Q_0rad(im2, ipq2_, ix, it) /= 0.0d0) THEN
                                matAm3 = Q_0radm1(im, iq1_, ix, it) * matAm3 * Q_0radm1(im2, ipq2_, ix, it)
                             ELSE
                                matAm3 = 0.0d0
                             END IF
!
                             ggtloc(count1) = ggtloc(count1) + matAm3 * hh1(count)! Ah
!
                             matrixloc(im, iq1_, ix, it) = matrixloc(im, iq1_, ix, it) + matAm3 * F_old1(im2, ipq2_, ix, it)
!
!
!
                      ! end if ! if dabs(freq2)
                      !
                          ENDDO
                   ! end if ! if dabs(freq1)
                       ENDDO
         !   write(15,*)sum_IT(:,it)
                    END IF
                 ENDDO
!
                 tmodes = tmodes + nanosec(T0) - T1
              ENDDO
!
!
!
!
           ENDDO
        ENDDO
     ENDDO
!
!
!enddo
!enddo
!enddo
!
  ENDDO !iq1_
!
! nell'ottica poi di eliminare dimG sarebbe il caso qua sotto di sostituire dimG con dimrun*
  CALL MPI_ALLREDUCE(matrixloc, matrix, dimG, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(ggtloc, ggt1, dimG, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE matrixA_times_h
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_vel()
!--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_CMM1, tpi, eps8, eps12
  USE common_variables, ONLY : &
    q1ws, q2ws, q3ws, q3bws, &
    neqq1, neqq2, neqq3, neqq3b, &
    freq1_, freq2_, freq3_, freq3b_, &
    zz1n, zz2n, zz3n, zz3bn, &
    nq1, nq2, iq1_, ipq2_, ipq3_, ipq3b_,&
    ic_inn, ic_med, ic_out, nq1tot, nq2tot,&
    mat2, deltaq2, deltaq1, rbig2, nRbig2t, &
    nat, celldm, bg, at, nat3, &
    velph, nxcryst, xcryst
  USE mpi_base
  IMPLICIT NONE
!
  Logical :: allzz
!
  Integer :: ix, i, im, isw, ix_
  Real(DP) :: freqc(nat3), freqc3(nat3), freqc3b(nat3), freqc2(nat3)
  Real(DP) :: q1d_(3)
  Real(DP) :: velph_io(3, nat3, nq1tot)              
  Complex(DP) :: zz1(nat3, nat3), zz3(nat3, nat3), zz3b(nat3, nat3), zz2(nat3, nat3)
  Complex(DP) :: zz_inout(nat3, nat3)
  Integer :: iq_init(3)
  Real(DP) :: q1dv(3), q3dv(3), q3dbv(3), q2dv(3), velphSCR(nat3, 3, nq1tot)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out
  Integer, Pointer :: iq2_inn, iq2_med, iq2_out
  Integer, Target :: iq1(3)
  Integer, Target :: iq2(3)
!
!complex(DP) ::   mat2f(nat3, nat3,nRbig2t)
!
!
!
  iq1_inn => iq1(ic_inn)
  iq1_med => iq1(ic_med)
  iq1_out => iq1(ic_out)
!
  iq2_out => iq2(ic_out)
  iq2_med => iq2(ic_med)
  iq2_inn => iq2(ic_inn)
!
  isw = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
  IF(nq1tot == nq2tot) THEN
!
!--------------------------------------------------------------------------------
!  -allzz- has to be true when deltaq1.ne.deltaq2 ....and if it is true the code evaluates the
!   modes for q1,q3 and q3b applying the same rotation used for the velocities.
!--------------------------------------------------------------------------------
     allzz = .False. ! deve essere true quando la griglia e' shiftata
!
     IF(nrank == 1) write(*,*) 'allzz = ', allzz
     IF(nrank == 1) write(*,*) 'allzz has to be TRUE when dq1.ne.dq2'
!
     IF((deltaq1(1) /= deltaq2(1)) .Or.(deltaq1(2) /= deltaq2(2)) .Or.(deltaq1(2) /= deltaq2(2))) THEN
        IF( .Not. allzz) CALL errore('calc vel', 'deltaq1.ne.deltaq2 but same  zz modes', 2)
     END IF
!--------------------------------------------------------------------------------
!
     DO iq1_out = 0, nq1(ic_out) - 1
        DO iq1_med = 0, nq1(ic_med) - 1
           DO iq1_inn = 0, nq1(ic_inn) - 1
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
              q1dv(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
              iq_init(:) = iq1(:)
              CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw)
!--------------------------------------------------------------------------------
!umklapp test
              CALL refold_q_ws(q1dv, q1ws(1, iq1_), neqq1(iq1_), bg)
              IF(allzz) CALL errore('calc vel', ' devi calcolare q2ws etc...', 1)
!--------------------------------------------------------------------------------
!            D(q1)
!--------------------------------------------------------------------------------
!
!   write(*,*) deltaq1,q1d
              CALL setupmat(nat, nRbig2t, q1dv, zz1, Rbig2, mat2)
              CALL dyndiag1b(nat3, zz1, freq1_(1, iq1_))
              freqc(:) = freq1_(:, iq1_)
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!!!! q1
!--------------------------------------------------------------------------------
              DO ix = 1, 3
!!! velph_io(1,1,iq1_)
                 zz_inout(:, :) = zz1(:, :)
                 CALL rotate_vel_zz(q1dv, ix, zz1, zz_inout, freqc, velph_io)! togli il common  e controlla a chi e' associato velph_io...se a quella numerica
                                                             !o a quella analitica
!!zznew   call rotate_vel_zz( q1d(ix),zz1,zz_inout(1,1,ix),freqc)
                 zz1n(:, :, iq1_, ix) = zz_inout(:, :)
              ENDDO
!
! if you uncomment the following line you need to commet all the calls to cryst_to_cart2
              CALL cryst_to_cart(nat3, velph_io(1, 1, iq1_), at,+1)
              velph(:, :, iq1_) = velph_io(:, :, iq1_) * celldm(1) / tpi
!
!
!--------------------------------------------------------------------------------
!            q3 =-q1-q2   Index
!--------------------------------------------------------------------------------
!
              IF(allzz) THEN
!
                 q3dv(:) = - dfloat(iq1(:)) / dfloat(nq1(:)) - 2.0d0 * deltaq1(:)
!
                 iq_init(:) = - iq1(:) - 3 * Nint(deltaq1(:)*nq1(:))
                 CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw)
!--------------------------------------------------------------------------------
!            D(q3)
!--------------------------------------------------------------------------------
                 CALL setupmat(nat, nRbig2t, q3dv, zz3, Rbig2, mat2)
                 CALL dyndiag1b(nat3, zz3, freq3_(1, ipq3_))
                 freqc3(:) = freq3_(:, ipq3_)
!
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
!            q3b = q2-q1   Index
!--------------------------------------------------------------------------------
 ! questa parte e' associata nei vari loop on qpoints a iq2...visto che q3db=q2-q1
                 q3dbv(:) = dfloat(iq1(:)) / dfloat(nq1(:))
!
                 iq_init(:) = iq1(:) - Nint(deltaq1(:)*nq1(:))
                 CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw)
!--------------------------------------------------------------------------------
!            D(q3b)
!--------------------------------------------------------------------------------
                 CALL setupmat(nat, nRbig2t, q3dbv, zz3b, Rbig2, mat2)
                 CALL dyndiag1b(nat3, zz3b, freq3b_(1, ipq3b_))
!
                 freqc3b(:) = freq3b_(:, ipq3b_)
!
!--------------------------------------------------------------------------------
!!!! q3
!--------------------------------------------------------------------------------
                 DO ix = 1, 3
                    zz_inout(:, :) = zz3(:, :)
                    CALL rotate_vel_zz(q3dv, ix, zz3, zz_inout, freqc3, velph_io)
                    zz3n(:, :, ipq3_, ix) = zz_inout(:, :)
                 ENDDO
!--------------------------------------------------------------------------------
!!!! q3b
!--------------------------------------------------------------------------------
                 DO ix = 1, 3
                    zz_inout(:, :) = zz3b(:, :)
                    CALL rotate_vel_zz(q3dbv, ix, zz3b, zz_inout, freqc3b, velph_io)
                    zz3bn(:, :, ipq3b_, ix) = zz_inout(:, :)
                 ENDDO
!
              END IF
!
!
!
!--------------------------------------------------------------------------------
!
!
!
              IF(.True.) THEN
                 IF(nrank == 1) THEN
                    q1d_(:) = q1dv(:)
                    CALL cryst_to_cart(1, q1d_, bg,+1)
!
       !     open(unit=15,file='grvelA',position='append')
                    Open(Unit=17, File='grvel', Position='append')
                    Open(Unit=18, File='freq1', Position='append')
!
                    Write(17, '(10x,3f10.6)') q1d_(1), q1d_(2), q1d_(3)
!
                    DO ix_ = 1, nxcryst
                       ix = xcryst(ix_)
!
                       Write(17, '(6e24.12)')(velph(ix, im, iq1_), im=1, nat3)
                    ENDDO
  !
                    Write(18, '(10x,3f10.6)') q1d_(1), q1d_(2), q1d_(3)
                    Write(18, '(6e24.12)')(freq1_(im, iq1_)*RY_TO_CMM1, im=1, nat3)
!
                 END IF
              END IF
!
 !end do
!
           ENDDO
        ENDDO
     ENDDO
!
     IF( .Not. allzz) THEN
        q2ws(:, :) = q1ws(:, :)
        q3ws(:, :) = q1ws(:, :)
        q3bws(:, :) = q1ws(:, :)
!
        neqq2(:) = neqq1(:)
        neqq3(:) = neqq1(:)
        neqq3b(:) = neqq1(:)
!
        freq2_(:, :) = freq1_(:, :)
        freq3_(:, :) = freq1_(:, :)
        freq3b_(:, :) = freq1_(:, :)
        zz2n(:, :, :, :) = zz1n(:, :, :, :)
        zz3n(:, :, :, :) = zz1n(:, :, :, :)
        zz3bn(:, :, :, :) = zz1n(:, :, :, :)
     END IF
!
!
  ELSE IF((nq1tot == 1) .and.(nq2tot > nq1tot)) THEN ! if nq1tot .eq.1
!
     DO iq1_out = 0, nq1(ic_out) - 1
        DO iq1_med = 0, nq1(ic_med) - 1
           DO iq1_inn = 0, nq1(ic_inn) - 1
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
              q1dv(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
              iq_init(:) = iq1(:)
              CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw)
!
    !   write(*,*) deltaq1,q1d
              CALL setupmat(nat, nRbig2t, q1dv, zz1, Rbig2, mat2)
              CALL dyndiag1b(nat3, zz1, freq1_(1, iq1_))
              freqc(:) = freq1_(:, iq1_)
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
!!!! q1
    !--------------------------------------------------------------------------------
              DO ix = 1, 3
                 zz_inout(:, :) = zz1(:, :)
                 CALL rotate_vel_zz(q1dv, ix, zz1, zz_inout, freqc, velph_io)
!
       !!zznew   call rotate_vel_zz( q1d(ix),zz1,zz_inout(1,1,ix),freqc)
                 zz1n(:, :, iq1_, ix) = zz_inout(:, :)
              ENDDO
!
              CALL cryst_to_cart(nat3, velph_io(1, 1, iq1_), at,+1)
!
              velph(:, :, iq1_) = velph_io(:, :, iq1_) * celldm(1) / tpi
!--------------------------------------------------------------------------------
              q1d_(:) = q1dv(:)
              CALL cryst_to_cart(1, q1d_, bg,+1)
!
  !   write(*,*) 'celldm,tpi', celldm(1),tpi
              Open(Unit=15, File='grvelA', Position='append')
              Open(Unit=18, File='freq1', Position='append')
!
              Write(15, '(3f10.6)') q1d_(:)
              Write(18, '(3f10.6)') q1d_(:)
              DO ix_ = 1, nxcryst
                 ix = xcryst(ix_)
                 DO i = 1, nat3
                    velphSCR(i, ix, iq1_) = velph(ix, i, iq1_)
                 ENDDO
                 Write(15, '(6e14.6)')(velphSCR(i, ix, iq1_)*RY_TO_CMM1, i=1, nat3)
              ENDDO
              Write(18, '(6e14.6)')(freq1_(i, iq1_)*RY_TO_CMM1, i=1, nat3)
!
!
!
!--------------------------------------------------------------------------------
!
              DO iq2_out = 0, nq2(ic_out) - 1
                 DO iq2_med = 0, nq2(ic_med) - 1
                    DO iq2_inn = 0, nq2(ic_inn) - 1
!
!--------------------------------------------------------------------------------
!            Indeces
!--------------------------------------------------------------------------------
                       q2dv(:) = dfloat(iq2(:)) / dfloat(nq2(:)) + deltaq2(:)
                       iq_init(:) = iq2(:)
                       CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq2_, isw)
!
!--------------------------------------------------------------------------------
!            D(q2)
!--------------------------------------------------------------------------------
                       CALL setupmat(nat, nRbig2t, q2dv, zz2, Rbig2, mat2)
                       CALL dyndiag1b(nat3, zz2, freq2_(1, ipq2_))
                       freqc2(:) = freq2_(:, ipq2_)
!
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!!!! q2
!--------------------------------------------------------------------------------
                       DO ix = 1, 3
                          zz_inout(:, :) = zz2(:, :)
                          CALL rotate_vel_zz(q2dv, ix, zz2, zz_inout, freqc2, velph_io)
                          zz2n(:, :, ipq2_, ix) = zz_inout(:, :)
                       ENDDO
!--------------------------------------------------------------------------------
!mq2 is simply the cc of q2 and this is evaluated inside the setup subroutine.
    !
    !        mq2d(:) = -q2d(:)
    !        iq_init(:) = -iq2(:)
    !        call calc_index(iq_init, nq2, ic_inn,ic_med,ic_out, imq2_,isw )
    !
    !--------------------------------------------------------------------------------
    !            q3 =-q1-q2   Index
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!  q3dv(:) =  - dfloat(iq1(:))/dfloat(nq1(:)) - 2*(deltaq1(:))
!   iq_init(:) = -iq1(:) - 3*nint(deltaq1(:)*nq2(:))
!
                       q3dv(:) = - q1dv(:) - q2dv(:)
                       iq_init(:) = - iq1(:) * nq2(:) / nq1(:) - iq2(:) &
                                    - 2 * Nint(deltaq1(:)*nq2(:)) - Nint (deltaq2(:)*nq2(:))
                       CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3_, isw)
!--------------------------------------------------------------------------------
!            D(q3)
!--------------------------------------------------------------------------------
                       CALL setupmat(nat, nRbig2t, q3dv, zz3, Rbig2, mat2)
                       CALL dyndiag1b(nat3, zz3, freq3_(1, ipq3_))
                       freqc3(:) = freq3_(:, ipq3_)
!
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
!            q3b = q2-q1   Index
!--------------------------------------------------------------------------------
 ! questa parte e' associata nei vari loop on qpoints a iq2...visto che q3db=q2-q1
!  q3dbv(:) = dfloat(iq1(:))/dfloat(nq1(:))
!   iq_init(:) = iq1(:) - nint(deltaq1(:)*nq1(:))
!
                       q3dbv(:) = q2dv(:) - q1dv(:)
                       iq_init(:) = iq2(:) - iq1(:) * nq2(:) / nq1(:) + Nint(deltaq2(:)*nq2(:)) - 2 * Nint &
                      &(deltaq1(:)*nq2(:))
!
                       CALL calc_index(iq_init, nq2, ic_inn, ic_med, ic_out, ipq3b_, isw)
!--------------------------------------------------------------------------------
!            D(q3b)
!--------------------------------------------------------------------------------
                       CALL setupmat(nat, nRbig2t, q3dbv, zz3b, Rbig2, mat2)
                       CALL dyndiag1b(nat3, zz3b, freq3b_(1, ipq3b_))
!
                       freqc3b(:) = freq3b_(:, ipq3b_)
!
!--------------------------------------------------------------------------------
!!!! q3
!--------------------------------------------------------------------------------
                       DO ix = 1, 3
                          zz_inout(:, :) = zz3(:, :)
                          CALL rotate_vel_zz(q3dv, ix, zz3, zz_inout, freqc3, velph_io)
                          zz3n(:, :, ipq3_, ix) = zz_inout(:, :)
                       ENDDO
!--------------------------------------------------------------------------------
!!!! q3b
!--------------------------------------------------------------------------------
                       DO ix = 1, 3
                          zz_inout(:, :) = zz3b(:, :)
                          CALL rotate_vel_zz(q3dbv, ix, zz3b, zz_inout, freqc3b, velph_io)
                          zz3bn(:, :, ipq3b_, ix) = zz_inout(:, :)
                       ENDDO
!
                    ENDDO
                 ENDDO
              ENDDO
!
!
           ENDDO
        ENDDO
     ENDDO
!
!
  END IF
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_vel
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_conduct(iqq1)
!--------------------------------------------------------------------------------
!
  USE constants, ONLY: DP, pi, eps8, RY_TO_WATT, BOHR_TO_M, tpi
  USE common_variables, ONLY : tcondc_, nat3, lbordo, Lcasm1, &
                               ntemp, tempm1, lambdiso, freq1_, lambd, velph,nxcryst,xcryst
  USE mpi_base
  IMPLICIT NONE
!
!  integer :: iq1_
  Integer :: ic, jc, im, it, iqq1
  Real(DP) :: wrk, bos_(nat3), lambdaBE, lambtot
!
! qui l'indice isigma di lambd e' stato messo uguale a 1, se fai un test variando la mesh cambia anche questo!
!
  DO it = 1, ntemp
     DO im = 1, nat3
        bos_(im) = 1.d0 /(Exp(freq1_(im, iqq1)*tempm1(it))-1.d0)
     ENDDO
!
!              write(119,'(6e24.12)')(bos_(im),im=1,nat3)
!
     DO jc = 1, nxcryst
        DO ic = 1, nxcryst
           wrk = 0.0d0
           DO im = 1, nat3
!
              lambdaBE = 0.0d0
              IF(lbordo) lambdaBE = Abs(velph(xcryst(ic), im, iqq1)) * Lcasm1 * 2.0d0
!            lambdaBE=  1096892.13894E-8*Lcasm1*2.0d0
!
              lambtot = pi * lambd(im, it, xcryst(ic), iqq1) + lambdiso(im, it, xcryst(ic), iqq1) + lambdaBE
              IF((lambtot /= 0.d0) .and.(dabs(freq1_(im, iqq1)) > eps8)) THEN
                 wrk = wrk + velph(xcryst(ic), im, iqq1) *velph(xcryst(jc), im, iqq1) &
                            * freq1_(im, iqq1)**2 &
                            * bos_(im) *(bos_(im)+1.0d0) / lambtot
             END IF
!
           ENDDO
           tcondc_(xcryst(ic), xcryst(jc), it) = tcondc_(xcryst(ic), xcryst(jc), it) + wrk

        ENDDO
     ENDDO

  ENDDO
!
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_conduct
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_f0(iqq1)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP
  USE common_variables, ONLY : nat3, ntemp, nxcryst, xcryst, Q_0, &
                               tempm1, freq1_, velph, F_0
  IMPLICIT NONE
!
  Integer :: iqq1, ix_
  Integer :: it, im, ix
  Real(DP) :: bosT, bos_ 
!
  DO it = 1, ntemp
     DO im = 1, nat3
        bos_ = 1.d0 /(Exp(freq1_(im, iqq1)*tempm1(it))-1.d0)
!      bosT=bos_*(bos_ + 1.0d0)*tempm1(it)/temp(it) 
       !sposto la costante 1/( k_B T2) nella subroutine che calcola la conduttivita'
        bosT = bos_ *(bos_+1.0d0)
        DO ix_ = 1, nxcryst
           ix = xcryst(ix_)
           IF(Q_0(im, iqq1, ix, it) /= 0.d0) THEN 
             F_0(im, iqq1, ix, it) &
               = - velph(ix, im, iqq1) * freq1_(im, iqq1) * bosT /  Q_0(im, iqq1, ix, it)
           ENDIF
        ENDDO
     ENDDO
  ENDDO
!
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_f0
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE sum_conductFI
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, eps8
  USE common_variables, ONLY : nat3, ntemp, F_0, tcondFI_, iq1_, &
                               tempm1, freq1_, velph,nxcryst,xcryst
  IMPLICIT NONE
!
!  integer :: iq1_
  Integer :: ic, jc, im, it
  Real(DP) :: wrk, bos_(nat3)
!
!
  DO it = 1, ntemp
     DO im = 1, nat3
        bos_(im) = 1.d0 /(Exp(freq1_(im, iq1_)*tempm1(it))-1.d0)
     ENDDO
!
     DO jc = 1, nxcryst
        DO ic = 1, nxcryst
           wrk = 0.0d0
           DO im = 1, nat3
              IF(dabs(freq1_(im, iq1_)) > eps8) THEN
                 wrk = wrk + velph(xcryst(ic), im, iq1_) &
                             * freq1_(im, iq1_) &
                             * bos_(im) *(bos_(im)+1.0d0) * F_0(im, iq1_, xcryst(jc), it)
              END IF
           ENDDO
           tcondFI_(xcryst(ic), xcryst(jc), it) = tcondFI_(xcryst(ic), xcryst(jc), it) + wrk
        ENDDO
     ENDDO
  ENDDO
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE sum_conductFI
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_conduct(nat3, nq1tot, ntemp, const_cond, omega, at, velph, freq1_, &
                         tempm1, temp, FF, tcond)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_WATT, BOHR_TO_M, RY_TO_CMM1
  USE common_variables, ONLY: F_0, Q_0, Q_0radm1, sigma
  IMPLICIT NONE
  Integer :: nat3, nq1tot, ntemp
  Real(DP) :: omega, at(3, 3), velph(3, nat3, nq1tot), &
               freq1_(nat3, nq1tot), FF(nat3, nq1tot, 3, ntemp)
  Real(DP) :: tempm1(ntemp), tcond(3, 3, ntemp), temp(ntemp), const_cond(ntemp)
!
  Integer :: iq1_
  Integer :: ic, jc, im, it
  Real(DP) :: wrk
  Real(DP), Allocatable :: bos_(:)
!
  ALLOCATE(bos_(nat3))
!
  DO it = 1, ntemp
     Write(*,*) '1st (standard) definition  k = - \lambda b F ', temp(it), sigma(it)*RY_TO_CMM1
     tcond(:, :, it) = 0.d0
     DO iq1_ = 1, nq1tot
        DO im = 1, nat3
           bos_(im) = 1.d0 /(Exp(freq1_(im, iq1_)*tempm1(it))-1.d0)
        ENDDO
!
        DO jc = 1, 3
           DO ic = 1, 3
              wrk = 0.0d0
              DO im = 1, nat3
          !   wrk = wrk + velph(ic,im,iq1_)*freq1_(im,iq1_) * bos_(im)*(bos_(im)+1.0d0)* FF(im,iq1_,jc,it)
                 IF(Q_0(im, iq1_, ic, it) /= 0.0d0) wrk = wrk - F_0(im, iq1_, ic, it) * Q_0(im, iq1_, ic, it) * Q_0radm1 &
                &(im, iq1_, ic, it) * FF(im, iq1_, jc, it)
              ENDDO
              tcond(ic, jc, it) = tcond(ic, jc, it) + wrk
!
           ENDDO
        ENDDO
     ENDDO
   !  tcond(:,:) =   - tcond(:,:) /(dfloat(nq1tot)*omega)*RY_TO_WATT/BOHR_TO_M ! utilizziamo la costante definita in read_mat2R
!
!
     DO ic = 1, 3
        DO jc = 1, 3
           tcond(ic, jc, it) = - tcond(ic, jc, it) * const_cond(it) * RY_TO_WATT / BOHR_TO_M
        ENDDO
     ENDDO
!
!       call cryst_to_cart2(1,tcond(1,1,it),at,+1)
!
     DO jc = 1, 3
        Write(*, '(6e20.12)')(tcond(ic, jc, it), ic=1, 3)
     ENDDO
     Write(*,*) '------------------------------------------------------------'
  ENDDO
!
!
  DEALLOCATE(bos_)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_conduct
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_conduct1(nat3, nq1tot, ntemp, const_cond, omega, at, velph, freq1_, tempm1, temp, FF, tcond)
!--------------------------------------------------------------------------------
!
  USE constants, ONLY: DP, RY_TO_WATT, BOHR_TO_M, RY_TO_CMM1
  USE common_variables, ONLY: F_0, Q_0, sigma
  IMPLICIT NONE
  Integer :: nat3, nq1tot, ntemp
  Real(DP) :: omega, at(3, 3), velph(3, nat3, nq1tot), &
               freq1_(nat3, nq1tot), FF(nat3, nq1tot, 3, ntemp)
  Real(DP) :: tempm1(ntemp), tcond(3, 3, ntemp), temp(ntemp), const_cond(ntemp)
!
  Integer :: iq1_
  Integer :: ic, jc, im, it
  Real(DP) :: wrk
  Real(DP), Allocatable :: bos_(:)
!
  ALLOCATE(bos_(nat3))
!
  DO it = 1, ntemp
     Write(*,*) 'Zeroth order solution  k = - \lambda b F_0', temp(it), sigma(it)*RY_TO_CMM1
     tcond(:, :, it) = 0.d0
     DO iq1_ = 1, nq1tot
        DO im = 1, nat3
           bos_(im) = 1.d0 /(Exp(freq1_(im, iq1_)*tempm1(it))-1.d0)
        ENDDO
!
        DO jc = 1, 3
           DO ic = 1, 3
              wrk = 0.0d0
              DO im = 1, nat3
          !   wrk = wrk + velph(ic,im,iq1_)*freq1_(im,iq1_) * bos_(im)*(bos_(im)+1.0d0)* FF(im,iq1_,jc,it)
           ! if(Q_0(im,it,iq1_).ne.0.0d0) &
                 wrk = wrk - F_0(im, iq1_, ic, it) * Q_0(im, iq1_, ic, it) * FF(im, iq1_, jc, it)
              ENDDO
              tcond(ic, jc, it) = tcond(ic, jc, it) + wrk
!
           ENDDO
        ENDDO
     ENDDO
   !  tcond(:,:) =   - tcond(:,:) /(dfloat(nq1tot)*omega)*RY_TO_WATT/BOHR_TO_M ! utilizziamo la costante definita in read_mat2R
!
     DO ic = 1, 3
        DO jc = 1, 3
           tcond(ic, jc, it) = - tcond(ic, jc, it) * const_cond(it) * RY_TO_WATT / BOHR_TO_M
        ENDDO
     ENDDO
!
     DO jc = 1, 3
        Write(*, '(6e20.12)')(tcond(ic, jc, it), ic=1, 3)
     ENDDO
     Write(*,*) '------------------------------------------------------------'
  ENDDO
!
!
  DEALLOCATE(bos_)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_conduct1
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_conduct2(nat3, nq2tot, nq1tot, ntemp, const_cond, at, matrix, FF, tcond2, temp)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_WATT, BOHR_TO_M, RY_TO_CMM1
  USE common_variables, ONLY : sigma
  IMPLICIT NONE
  Integer :: nat3, nq2tot, nq1tot, ntemp
  Real(DP) :: at(3, 3), temp(ntemp), const_cond(ntemp), FF(nat3, nq1tot, 3, ntemp)
  Real(DP) :: matrix(nat3, nq1tot, 3, ntemp)
  Real(DP) :: tcond2(3, 3, ntemp)
!
  Integer :: iqq1_
  Integer :: ic, jc, im1, it
  Real(DP) :: wrk
!
!
  DO it = 1, ntemp
     Write(*,*) '2nd Definition  k = \lambda F A F ', temp(it), sigma(it)*RY_TO_CMM1
     DO jc = 1, 3
        DO ic = 1, 3
           tcond2(ic, jc, it) = 0.d0
        ENDDO
     ENDDO
!
     DO iqq1_ = 1, nq1tot
!
        DO jc = 1, 3
           DO ic = 1, 3
              wrk = 0.0d0
              DO im1 = 1, nat3
           !         do im2=1,nat3
                 wrk = wrk + FF(im1, iqq1_, ic, it) * matrix(im1, iqq1_, jc, it)
           !         end do
              ENDDO
              tcond2(ic, jc, it) = tcond2(ic, jc, it) + wrk
!
           ENDDO
        ENDDO
     ENDDO
!
     DO ic = 1, 3
        DO jc = 1, 3
           tcond2(ic, jc, it) = tcond2(ic, jc, it) * const_cond(it) * RY_TO_WATT / BOHR_TO_M
        ENDDO
     ENDDO
!
!       call cryst_to_cart2(1,tcond2(1,1,it),at,+1)
!
     DO jc = 1, 3
        Write(*, '(6e20.12)')(tcond2(ic, jc, it), ic=1, 3)
     ENDDO
     Write(*,*) '------------------------------------------------------------'
  ENDDO
!
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_conduct2
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_conduct3(ntemp, temp, const_cond, tcond, tcond2, tcond3)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_CMM1
  USE common_variables, ONLY : sigma
  IMPLICIT NONE
  Integer it, ic, jc, ntemp
  Real(DP) :: tcond(3, 3, ntemp), tcond2(3, 3, ntemp), &
               tcond3(3, 3, ntemp), temp(ntemp), const_cond(ntemp)
!
  DO it = 1, ntemp
     Write(*,*) 'Variational definition k = - 2 \lambda (1/2 FAF - bF)', temp(it), sigma(it)*RY_TO_CMM1
!
!
     DO jc = 1, 3
        DO ic = 1, 3
           tcond3(ic, jc, it) = - 2.0d0 *(0.5d0*tcond2(ic, jc, it)-tcond(ic, jc, it))
        ENDDO
     ENDDO
!          call cryst_to_cart2(1,tcond3(1,1,it),at,+1)
     DO jc = 1, 3
        Write(*, '(6e20.12)')(tcond3(ic, jc, it), ic=1, 3)
     ENDDO
     Write(*,*) '------------------------------------------------------------'
  ENDDO
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_conduct3
!--------------------------------------------------------------------------------
SUBROUTINE calc_tcSMA
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_WATT, BOHR_TO_M, RY_TO_CMM1
  USE common_variables, ONLY : tcondc_SMA, const_cond, tcondc_, temp, sigma, ntemp,nxcryst,xcryst
  IMPLICIT NONE
  Integer :: ic, jc, it
!
!
  DO it = 1, ntemp
     DO jc = 1,nxcryst
        DO ic = 1, nxcryst
           tcondc_SMA(xcryst(ic), xcryst(jc)) = const_cond(it) * tcondc_(xcryst(ic), xcryst(jc), it) * RY_TO_WATT / BOHR_TO_M
!
        ENDDO
     ENDDO
!     call cryst_to_cart2(1,tcondc_SMA,at,+1)
!
     Write(*,*) 'SMA', temp(it), sigma(it)*RY_TO_CMM1
     DO jc = 1, 3
        Write(*, '(6e20.12)')(tcondc_SMA(xcryst(ic), xcryst(jc)), ic=1, nxcryst)
     ENDDO
     Write(*,*) '------------------------------------------------------------'
  ENDDO
  !--------------------------------------------------------------------------------
END SUBROUTINE calc_tcSMA
!--------------------------------------------------------------------------------
FUNCTION mod_per(ii, nn)
  !--------------------------------------------------------------------------------
  IMPLICIT NONE
  Integer :: ii, nn, mod_per
  IF(ii .Ge. 0) THEN
     mod_per = Mod(ii, nn)
  ELSE
     mod_per = nn - 1 - Mod(Abs(ii+1), nn)
  END IF
  Return
  !--------------------------------------------------------------------------------
END FUNCTION mod_per
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE write_F0
  !--------------------------------------------------------------------------------
  USE constants, ONLY: RY_TO_CMM1
  USE common_variables, ONLY: F_0, ntemp, nq1tot, Q_0, freq1_,nxcryst,xcryst
  IMPLICIT NONE
  Integer :: it, ix_,ix, iqq1_
!
  Open(Unit=23, File='F_0', Position='append')
  Open(Unit=22, File='Q_0', Position='append')
  Open(Unit=21, File='freq', Position='append')
!
!
  DO it = 1, ntemp
   !   do im=1,nat3
     DO ix_ = 1,nxcryst
       ix=xcryst(ix_)
   !   irun=ix+3*(it-1)
        DO iqq1_ = 1, nq1tot
!
           Write(23, '(i5,99e14.6)') iqq1_, F_0(:, iqq1_, ix, it)
        ENDDO
  !        write(22,*) Q_0(im,it,iq1_)
  !   end do
     ENDDO
  ENDDO
!
  DO ix_ = 1, nxcryst
     ix=xcryst(ix_)
     DO iqq1_ = 1, nq1tot
        DO it = 1, ntemp
           Write(22, '(i5,99e14.6)') iqq1_, Q_0(:, iqq1_, ix, it)
        ENDDO
     ENDDO
  ENDDO
!
  DO iqq1_ = 1, nq1tot
     Write(21, '(i5,99e14.6)') iqq1_, freq1_(:, iqq1_) * RY_TO_CMM1
  ENDDO
!
  Close(21)
  Close(22)
  Close(23)
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE write_F0
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE write_FF(iter, FF)
  !--------------------------------------------------------------------------------
  USE constants, ONLY: DP
  USE common_variables, ONLY: ntemp, nq1tot, nat3, &
                              bg, nq1, ic_inn, ic_med, ic_out, deltaq1,nxcryst,xcryst
  IMPLICIT NONE
  Integer :: it, ix,ix_, iq1_, iter, iq_init(3), isw
  Real(DP) :: q1d_(3), FF(nat3, nq1tot, 3, ntemp)
  Integer, Pointer :: iq1_inn, iq1_med, iq1_out
  Integer, Target :: iq1(3)
!
  Character(Len=256) :: name(3), filename
  Character(Len=2) :: ctemp
!
  name(1) = 'FFx'
  name(2) = 'FFy'
  name(3) = 'FFz'
!
!
  DO it = 1, ntemp
     ctemp = '00'
     IF(it .Le. 9) THEN
        Write(ctemp(2:2), '(i1)') it
     ELSE
        Write(ctemp(1:2), '(i2)') it
     END IF
!
!
     DO ix_ = 1, nxcryst
        ix=xcryst(ix_)
!    filename = trim(name(ix)) // trim(chlab)
        filename = trim(name(ix)) // trim(ctemp)
        Open(11, File=filename)
!
        iq1_inn => iq1(1)
        iq1_med => iq1(2)
        iq1_out => iq1(3)
!
   !   do im=1,nat3
        isw = 1 ! -1 from the total iq1_ index to the 3 iq1(:); isw=1 from the three iq1(:) indeces to iq1_
        DO iq1_out = 0, nq1(ic_out) - 1
           DO iq1_med = 0, nq1(ic_med) - 1
              DO iq1_inn = 0, nq1(ic_inn) - 1
!--------------------------------------------------------------------------------
!            q1   Index
!--------------------------------------------------------------------------------
                 q1d_(:) = dfloat(iq1(:)) / dfloat(nq1(:)) + deltaq1(:)
                 iq_init(:) = iq1(:)
                 CALL calc_index(iq_init, nq1, ic_inn, ic_med, ic_out, iq1_, isw)
!--------------------------------------------------------------------------------
                 CALL cryst_to_cart(1, q1d_, bg,+1)
!
!                 if(abs(q1d_(3)).lt.eps8) THEN
                 Write(11, '(9e14.6)') q1d_(1), q1d_(2), q1d_(3), FF(:, iq1_, ix, it)
!                 end if
              ENDDO
           ENDDO
        ENDDO
!
 !   end do
     ENDDO
     Close(11)
  ENDDO
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE write_FF
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE rdiagh(n, h, ldh, e, v)
  !--------------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a real
  ! ... simmetric matrix H . On output, the matrix is unchanged
  !
  USE constants, ONLY: DP
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  Integer :: n, ldh
! dimension of the matrix to be diagonalized
! leading dimension of h, as declared in the calling pgm unit
  Real(DP) :: h(ldh, n)
! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  Real(DP) :: e(n)! eigenvalues
  Real(DP) :: v(ldh, n)! eigenvectors(column-wise)
  !
  !
  !
!  CALL rdiagh_aix()
!
  CALL rdiagh_lapack()
  !
  !
  Return
  !
 Contains
!
!................................................................................
SUBROUTINE rdiagh_lapack()
  !................................................................................
     IMPLICIT NONE
  !
  ! ... local variables(LAPACK version)
  !
     Integer :: lwork, nb, info
     Integer, External :: ILAENV ! ILAENV returns optimal block size "nb"
     Real(DP), Allocatable :: work(:)
  !
  ! ... check for the block size
     nb = ILAENV(1, 'DSYTRD', 'U', n,-1,-1,-1)
     IF(nb < 1 .Or. nb >= n) THEN
        lwork = 3 * n
     ELSE
        lwork =(nb+2) * n
     END IF
  ! ... allocate workspace
     v = h
  !
     ALLOCATE(work(lwork))
     CALL DSYEV('V', 'U', n, v, ldh, e, work, lwork, info)
     CALL errore('rdiagh', 'diagonalization(DSYEV) failed', Abs(info))
  ! ... DEALLOCATE workspace
  !
     DEALLOCATE(work)
  !
     Return
    !
    !................................................................................
  END SUBROUTINE rdiagh_lapack
  !................................................................................
!
END SUBROUTINE rdiagh
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE calc_index(iq_init, nq, ic_inn, ic_med, ic_out, iq_fin, isw)
  !--------------------------------------------------------------------------------
  IMPLICIT NONE
  Integer :: ic_inn, ic_med, ic_out, iq_init(3), iq_fin, nq(3)
  Integer :: iq_1, iq_2, iq_3, mod_per, isw
  Integer :: itmp
!
  IF(isw == 1) THEN
!
     iq_3 = mod_per(iq_init(3), nq(ic_out))
     iq_2 = mod_per(iq_init(2), nq(ic_med))
     iq_1 = mod_per(iq_init(1), nq(ic_inn))
     iq_fin = 1 + iq_1 + iq_2 * nq(ic_inn) + iq_3 * nq(ic_inn) * nq(ic_med)! index q2
  ELSE IF(isw ==-1) THEN
     itmp = iq_fin - 1
     iq_init(1) = mod_per(itmp, nq(ic_inn))
     itmp =(itmp-iq_init(1)+1) / nq(ic_inn)
     iq_init(2) = mod_per(itmp, nq(ic_med))
     itmp =(itmp-iq_init(2)+1) / nq(ic_med)
     iq_init(3) = mod_per(itmp, nq(ic_out))
!
  END IF
!
  Return
!
END SUBROUTINE calc_index
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE rotate_vel_zz(qqd, ix, zz1, zz_inout, freqc, velphA)!qui ho associato la velocita' analitica
!--------------------------------------------------------------------------------
  USE constants, ONLY: DP, RY_TO_CMM1, tpi, eps8, eps12
! togli il common
  USE common_variables
  IMPLICIT NONE
  Complex(DP) :: zz1(nat3, nat3)
  Real(DP) :: qqd(3)
  Complex(DP) :: zz_inout(nat3, nat3)
  Complex(DP) :: z_(nat3, nat3), zp(nat3, nat3), zm(nat3, nat3), diffM(nat3, nat3), diffM_(nat3, nat3)
  Complex(DP) :: zz1sub(nat3, nat3)
  Real(DP) :: d_q, pqi(3), mqi(3), freqcp(nat3), freqcm(nat3), freqc(nat3)

  Real(DP) :: velphA(3, nat3, nq1tot), velphN(3, nat3, nq1tot), &
               groupcA3(3, nat3, nat3), groupcA2(nat3, nat3), &
               groupcA2sub(nat3, nat3), groupcA4(3, nat3, nat3)
  Integer :: nf, ndeg(nat3), is, ndim_, ldim, list(nat3, nat3)
  Real(DP) :: freqdeg(nat3), autgd(nat3), eiggd(nat3, nat3)
  Integer :: ix, im, i, j, ip, jp, ii, jj, k
!
!
  d_q = 1.0E-8_DP
  groupcA3 = 0.0d0
  groupcA4 = 0.0d0
!
!
  pqi(:) = qqd(:)
  pqi(ix) = qqd(ix) + d_q
  mqi(:) = qqd(:)
  mqi(ix) = qqd(ix) - d_q
!
  CALL setupmat(nat, nRbig2t, pqi, zp, Rbig2, mat2)
  z_ = zp
!
  CALL dyndiag1b(nat3, z_, freqcp)
!
!
  CALL setupmat(nat, nRbig2t, mqi, zm, Rbig2, mat2)
  z_ = zm
  CALL dyndiag1b(nat3, z_, freqcm)
!
  DO i = 1, nat3
     velphN(ix, i, iq1_) =(freqcp(i)-freqcm(i)) /(2.0d0*d_q)
  ENDDO
!
  DO ip = 1, nat3
     DO jp = 1, nat3
        diffM_(ip, jp) =(zp(ip, jp)-zm(ip, jp)) /(2.0d0*d_q)
     ENDDO
  ENDDO
  DO ip = 1, nat3
     DO jp = 1, nat3
        diffM(ip, jp) = 0.5d0 *(diffM_(ip, jp)+CONJG(diffM_(jp, ip)))
     ENDDO
  ENDDO
!
  DO i = 1, nat3
     DO j = 1, nat3
        DO ip = 1, nat3
           DO jp = 1, nat3
              groupcA3(ix, i, j) = groupcA3(ix, i, j) + real((CONJG(zz1(ip, i)))*diffM(ip, jp)*zz1(jp, j)) / &
             &(2.0d0*freq1_(i, iq1_))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!
  DO i = 1, nat3
     DO j = 1, nat3
        groupcA2(j, i) = groupcA3(ix, j, i)
     ENDDO
  ENDDO
  ldim = nat3
 !
  ndeg(:) = 0
  nf = 1
  ndeg(nf) = 1
  freqdeg(nf) = freqc(1)
  list(nf, 1) = 1
!
  DO i = 2, nat3
     IF(Abs(freqc(i)-freqdeg(nf)) < eps12) THEN
        ndeg(nf) = ndeg(nf) + 1
!                 !  write(*,*) 'deg',i,freqc(i),freqnf
        list(nf, ndeg(nf)) = i
     ELSE
        nf = nf + 1
        freqdeg(nf) = freqc(i)
        ndeg(nf) = 1
        list(nf, ndeg(nf)) = i
     END IF
  ENDDO
!
  im = 0
!
!
  DO is = 1, nf
     IF(ndeg(is) == 1) THEN
        im = im + 1
     ELSE
        DO j = 1, ndeg(is)
           DO k = 1, ndeg(is)
              groupcA2sub(j, k) = groupcA2(list(is, j), list(is, k))
           ENDDO
        ENDDO
             !
        ndim_ = ndeg(is)
!
        CALL rdiagh(ndim_, groupcA2sub, ldim, autgd, eiggd)
             !
!
        DO j = 1, ndim_
           groupcA2(list(is, j), list(is, j)) = autgd(j)
        ENDDO
!
        zz1sub = 0.0d0
!
        DO ii = 1, ndim_
           DO jj = 1, ndim_
              zz1sub(:, list(is, ii)) = zz1sub(:, list(is, ii)) + zz1(:, list(is, jj)) * eiggd(jj, ii)
           ENDDO
        ENDDO
!
        DO ii = 1, ndim_
           zz_inout(:, list(is, ii)) = zz1sub(:, list(is, ii))
        ENDDO
!
     END IF
!
!
  ENDDO
  DO i = 1, nat3
     velphA(ix, i, iq1_) = groupcA2(i, i)!/(2.0d0*freqc(i))
  ENDDO
  Return
  !
  !--------------------------------------------------------------------------------
END SUBROUTINE rotate_vel_zz
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE define_iq1_loc(nq_tot, iqp, nq, nq_max, nrank, nsize)
  !--------------------------------------------------------------------------------
  IMPLICIT NONE
  Integer :: nq_tot, nq_max, iqp(*), nq, nrank, nsize
  Integer :: ii
!
  IF(nrank .Le. 0) CALL errore('define_ien', 'wrong nrank', 1)
  IF(nsize < nrank) CALL errore('define_ien', 'wrong nsize', 1)
  IF(nq_tot .Le. 0) CALL errore('define_ien', 'wrong nrank', 1)
!
  IF(nrank .Le. nq_tot-Int(nq_tot/nsize)*nsize) THEN
     nq = nq_tot / nsize + 1
  ELSE
     nq = nq_tot / nsize
  END IF
  IF(nq > nq_max) CALL errore('define_ien', 'nq too big', 1)
!
  iqp(1) = nrank
  DO ii = 2, nq
     iqp(ii) = iqp(ii-1) + nsize
  ENDDO
!
  Return
END SUBROUTINE define_iq1_loc
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
SUBROUTINE refold_q_ws(qq_in, qq_out, nequiv, bg)
  !--------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  Integer :: nequiv
  Real(Kind=8) :: qq_in(3), qq_out(3), bg(3, 3)
!
  Logical, Save :: first = .True.
  Integer, Parameter :: nGmax = 10000, nGX = 2
  Integer, Save :: nGvec
  Integer :: i1, i2, i3, iG, iws
  Real(Kind=8), Save :: Gvec(3, nGmax)
  Integer, Save :: iGvec(3, nGmax)
  Real(Kind=8) :: qq(3), Mod, mod_
  Real(Kind=8), Parameter :: thresh = 1.d-12
!
  IF(first) THEN
     nGvec = 1
     Gvec(:, nGvec) = 0.d0
     iGvec(:, nGvec) = 0
     DO i1 = - nGX, nGX
        DO i2 = - nGX, nGX
           DO i3 = - nGX, nGX
              IF(i1 /= 0 .Or. i2 /= 0 .Or. i3 /= 0) THEN
                 nGvec = nGvec + 1
                 IF(nGvec > nGmax) CALL errore('generate G', 'nGmax exceeded', 1)
                 Gvec(:, nGvec) = bg(:, 1) * dfloat(i1) + bg(:, 2) * dfloat(i2) + bg(:, 3) * dfloat(i3)
                 iGvec(1, nGvec) = i1
                 iGvec(2, nGvec) = i2
                 iGvec(3, nGvec) = i3
              END IF
           ENDDO
        ENDDO
     ENDDO
  END IF
  first = .False.
!
  qq(:) = qq_in(:)
  CALL cryst_to_cart(1, qq, bg,+1)! crystal --> cartestian
  Mod = qq(1) ** 2 + qq(2) ** 2 + qq(3) ** 2
  iws = 1 ! The first vector of the Gvec is zero
  nequiv = 1
!
  DO iG = 2, nGvec
     mod_ =(qq(1)+Gvec(1, iG)) ** 2 +(qq(2)+Gvec(2, iG)) ** 2 +(qq(3)+Gvec(3, iG)) ** 2
     IF(mod_ < Mod-thresh) THEN
        Mod = mod_
        iws = iG
        nequiv = 1
     ELSE IF(Abs(Mod-mod_) .Le. thresh) THEN
        nequiv = nequiv + 1
     END IF
  ENDDO
  qq_out(:) = qq_in(:) + dfloat(iGvec(:, iws))
!
  Return
  !--------------------------------------------------------------------------------
END SUBROUTINE refold_q_ws
!--------------------------------------------------------------------------------
