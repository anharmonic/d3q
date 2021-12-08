!============================================================================================
!    ACRS - Derivative-Free Adaptive Controlled Random Search algorithm for 
!           bound constrained global optimization problems 
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G. Liuzzi, S. Lucidi, F. Parasiliti, M. Villani. Multiobjective optimization techniques 
!    for the design of induction motors, IEEE Transactions on Magnetics, 39(3): 1261-1264 (2003)
!    DOI: 10.1109/TMAG.2003.810193
!
!    L. Cirio, S. Lucidi, F. Parasiliti, M. Villani. A global optimization approach for 
!    the synchronous motors design by finite element analysis, International Journal of 
!    Applied Electromagnetics and Mechanics, 16(1): 13-27 (2002)
!============================================================================================

SUBROUTINE ACRS0(N,XSOL,FOUT,FUNCT)
   IMPLICIT NONE
   INTEGER,      INTENT(IN)   :: N
   DOUBLE PRECISION, INTENT(INOUT) :: XSOL(N)
   DOUBLE PRECISION, INTENT(OUT) :: FOUT
   EXTERNAL                :: FUNCT

   INTEGER            :: IOUT, M
   INTEGER        :: MAXITER, PRNLEV
   REAL           :: CPU_LIMIT, MAXTGEN

   DOUBLE PRECISION           :: LB(N),UB(N), TOL
   INTEGER        :: IEXIT
   DOUBLE PRECISION,ALLOCATABLE    :: S(:,:), FVAL(:) !S(N,M), FVAL(M), 
   INTEGER         :: ITER, NFTOT
   DOUBLE PRECISION :: diff_initial, omega

   M = MAX(50, 25*N)
   ALLOCATE(S(N,M), FVAL(M))
   S = 0.d0
   FVAL = 0.d0
   UB=MAXVAL(ABS(XSOL))*10
   LB=-UB(1)
   ! max number of steps and function evaluations
   MAXITER =  500000
   MAXTGEN = 1000000
   CPU_LIMIT = 36000 ! 10 hours
   TOL = 1.d-12
   PRNLEV = 1
   diff_initial = 100.d0
   omega = 100.d0
   IOUT=6
   IEXIT = -1

   CALL ACRS(N,LB,UB,M,S,FVAL,FUNCT,cpu_limit,maxtgen,maxiter,prnlev,tol,IOUT,   &
         XSOL,FOUT,iter,nftot,iexit,diff_initial, omega)

   DEALLOCATE(S, FVAL)

END SUBROUTINE

SUBROUTINE ACRS(N,LB,UB,M,S,FVAL,FUNCT,cpu_limit,maxtgen,maxiter,prnlev,tol,IOUT,   &
      XSOL,FOUT,iter,nftot,iexit,diff_initial,omega            )

!---------------------------------------------------------------------------------
!
!  PRICE ALGORITHM
!  first implemented by G. Di Pillo, S. Lucidi and
!  successively modified by G. Liuzzi
!  November 14, 2005
!
!  Reference: P. Brachetti, M. De Felice Ciccoli, G. Di Pillo, S. Lucidi
!             "A new version of the Price's algorithm for global optimization"
!             Journal of Global Optimization, 10  pp.165-184, 1997.
!---------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,      INTENT(IN)   :: N, IOUT, M
   INTEGER,      INTENT(IN)   :: MAXITER, PRNLEV
   REAL,         INTENT(IN)   :: CPU_LIMIT, MAXTGEN

   DOUBLE PRECISION           :: LB(N),UB(N), TOL
   INTEGER,      INTENT(INOUT)   :: IEXIT
   DOUBLE PRECISION, INTENT(INOUT)  :: S(N,M), FVAL(M), XSOL(N)
   INTEGER,      INTENT(OUT)  :: ITER, NFTOT
   DOUBLE PRECISION, INTENT(OUT) :: FOUT
   EXTERNAL                :: FUNCT
!---------------------------------------------------------------------------------
!  N        : (ON ENTRY) The number of variables. Unchanged on exit
!             it must be >= 1
!
!  LB(N)    : (ON ENTRY) The lower bounds on the variables. Unchanged on exit
!             Every variable must have a proper lower bound
!
!  UB(N)    : (ON ENTRY) The upper bounds on the variables. Unchanged on exit
!             Every variable must have a proper upper bound
!
!  M        : (ON ENTRY) Dimension of the working set. (ON EXIT) Unchanged
!             It should be a bit more than N+1. The suggested value is MAX(50, 25*N).
!             If WKDIMEN <= N+1 then the routines terminates with an error
!
!  S(N,M)      : (ON ENTRY)
!             if IEXIT > -2 Unspecified.
!             if IEXIT = -2, -3 initial working set.
!             (ON EXIT) The last working set used by ACRS
!
!  FVAL(M)     : (ON ENTRY)
!             if IEXIT > -3 Unspecified.
!             if IEXIT = -3 then FVAL(i) = f(S(1:n,i)).
!             (ON EXIT) The last working set values used by ACRS
!
!  FUNCT    : (ON ENTRY) The actual name of the subroutine that computes the obj.
!             function value. It must be declared EXTERNAL in the calling driver.
!             The FUNCT subroutine must have the following interface
!
!             SUBROUTINE FUNCT (X, N, F)
!              INTEGER           :: N
!              DOUBLE PRECISION  :: X(N), F
!
!  CPU_LIMIT   : (ON ENTRY) The maximum allowed cpu time in seconds. Unchanged on exit
!
!  MAXTGEN     : (ON ENTRY) The maximum allowed cpu time for initial working
!                set generation. MAXTGEN should be <= CPU_LIMIT. If MAXTGEN > CPU_LIMIT then
!                CPU_LIMIT is used instead. (ON EXIT) Unchanged
!
!  MAXITER     : (ON ENTRY) The maximum allowed number of iterations
!
!  PRNLEV      : (ON ENTRY) The printing level. (ON EXIT) Unchanged
!             If PRNLEV < -1 then OUTLEV is automatically reset to -1
!             If PRNLEV >  2 then OUTLEV is automatically reset to  2
!             PRNLEV = -1 -->  NO OUTPUT
!                 PRNLEV =  0 -->  CONCISE OUTPUT
!                 PRNLEV = +1 -->  VERBOSE
!                 PRNLEV = +2 -->  VERY VERBOSE (for debug purpose only)
!
!  TOL         : (ON ENTRY) Tolerance in the stopping criterion. (ON EXIT) Unchanged
!
!  IOUT     : (ON ENTRY) File header onto which redirect the routine output. (ON EXIT) Unchanged
!
!  XSOL     : (ON ENTRY)
!             if IEXIT ==  0, unspecified
!             if IEXIT == -1, specifies the point to be added to the working set
!             (ON EXIT) The computed solution.
!
!  FOUT     : (ON ENTRY) Unspecified. (ON EXIT) The minimum computed function value.
!
!  ITER     : (ON ENTRY) Unspecified. (ON EXIT) The number of iterations to convergence
!
!  NFTOT    : (ON ENTRY) Unspecified. (ON EXIT) The number of objective function evaluations to convergence.
!
!  IEXIT    : (ON ENTRY)
!             0 --> Generate set S from scratch
!            -1 --> Generate set S then add point specified by XSOL to it
!            -2 --> Do not generate set S. Use the one provided
!            -3 --> Do not generate set S. Continue from previous set
!             (ON EXIT) An integer specifying the type of output.
!             0 --> Normal termination
!             Abnormal or early termination signals
!             1 --> N <= 0
!             2 --> maximum number of iterations reached
!             3 --> exceeded CPU time limit
!             4 --> unable to generate n+1 different points
!             5 --> exceeded number of constraints violations
!             6 --> exceeded max. number of f. evaluations
!             7 --> some var has an improper lower bound
!             8 --> some var has an improper upper bound
!             9 --> constraints are too tight to be easily satisfied
!            10 --> provided working set dimension too small
!---------------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Integer variables:
!-------------------------------------------------------------------------
   INTEGER  :: i, j, outlev, viol_count, collis_count, rej_count, ierr

!-------------------------------------------------------------------------
!  Single and Double prec. variables:
!-------------------------------------------------------------------------
   DOUBLE PRECISION :: diff, diff_initial, omega, ftemp
   REAL         :: timeset, timesol, timetot, time

!-------------------------------------------------------------------------
!  Logical variables
!-------------------------------------------------------------------------
   LOGICAL :: viol, collision, proj
!-------------------------------------------------------------------------

   CALL CPU_TIME(time)

   if( n <=  0) then
      iexit = 1
      goto 900
   endif

   if( M <= N+1 ) then
      iexit = 10
      goto 900
   endif

   outlev = prnlev
   if( prnlev < -1 ) outlev = -1
   if( prnlev >  2 ) outlev =  2

   DO I = 1,N
      IF(LB(I) <= -1.D+6) THEN
         iexit = 7
         exit
      endif
      if(UB(I) >= 1.D+16) THEN
         iexit = 8
         exit
      endif
   ENDDO

   if(iexit > 0) goto 1000

   !CALL RANDOM_SEED()

!========================================================
! Initializations of state variables
!========================================================

   proj       = .TRUE.
   if(iexit > -3) omega      = dble(n) !1.D+3

   iter       = 0
   nftot      = 0

   viol       = .FALSE.
   collision  = .FALSE.
   rej_count  = 0
   viol_count = 0

!========================================================
! Step 0 : Determine the initial set S = {x_1, ... ,x_m}
!          and evaluate f at each point
!========================================================
   IF ( OUTLEV >= 1 ) THEN
      WRITE(iout,*) '++++++++++++ Determine the initial set S +++++++++ '
   END IF

   if(iexit > -3) then
      CALL GENERATE_SET_S(MAXTGEN,N,M,NFTOT,OUTLEV,IOUT,LB,UB,TOL,S,FVAL,iexit,XSOL,FUNCT)
      IF(iexit > 0) GOTO 1000
   endif

   CALL CPU_TIME(timeset)
   timeset = timeset - time


!  open(11,file="cort.1",status="old")
!  do i = 1,m
!  do j = 1,n
!     read(11,*) S(j,i)
!  enddo
!  read(11,*) fval(i)
!  enddo
!  close(11)

   diff      = fval(m) - fval(1)
   if(iexit > -3) diff_initial = diff

   iexit = 0

   IF ( diff_initial > 1.D+10 ) THEN
      diff_initial = 1.D+10
      IF ( OUTLEV >= 1 ) THEN
         WRITE(iout,*) ' Modified diff_initial '
      END IF
   END IF


   do while (.true.)

      !---------------------
      ! do some printing
      !---------------------
      if(OUTLEV >= 0) then
         IF(MOD(ITER,30)==0) THEN
            WRITE(iout,*)
            WRITE(iout,3000)
            WRITE(iout,3010)
         ENDIF
         WRITE(iout,3020) ITER,NFTOT,FVAL(1),FVAL(M),diff
      endif

      IF ( OUTLEV >= 1 ) THEN
         WRITE(iout,*) ' '
         WRITE(iout,*) '****** Iteration ', iter , ' ****** '
      END IF

      !------------------------------------------
      ! Stopping criterion: stops either because
      !     tolerance    reached -- iexit = 0
      !     maxiter      reached -- iexit = 2
      !     max cpu time reached -- iexit = 3
      !------------------------------------------
      CALL CPU_TIME(timesol)
      timesol = timesol - time
      IF ( diff < tol      ) then
         iexit = 0
         exit
      END IF
      IF ( iter >= maxiter ) THEN
                        write(*,*) '********** ',iter,maxiter
         iexit = 2
         exit
      END IF
      IF ( timesol >= cpu_limit ) THEN
         iexit = 3
         exit
      END IF

      !-----------------------------------
      ! perform the basic Price iteration
      !-----------------------------------
      CALL PRICE_BASE_ITERATION(N,M,S,FVAL,IEXIT,NFTOT,IOUT,OUTLEV,LB,UB,OMEGA,DIFF,DIFF_INITIAL,  &
                 PROJ,COLLIS_COUNT,REJ_COUNT,VIOL_COUNT,FUNCT                    )

      !-----------------------------------
      ! if there were errors then STOP
      !-----------------------------------
      if(iexit > 0) exit

      iter = iter + 1

   enddo

!----------------------------------------------------------------
! Normal exit
!----------------------------------------------------------------

900   CONTINUE
        
   CALL CPU_TIME(timesol)
   timesol = timesol - time
   timetot = timeset + timesol

   if(OUTLEV >= -1) then

      SELECT CASE (iexit)
      CASE(0)

! Normal exit. Printing the final results

         WRITE(iout,2008)
         WRITE(iout,*) ' Modified Improved Price Algorithm '
         WRITE(iout,2010) n
         WRITE(iout,2009) m
         WRITE(iout,2006) tol

         WRITE(iout,2004) fval(1)

         write(iout,*) FVAL(m),FVAL(1)
         WRITE(iout,*) ' # iterations : ' , iter
         WRITE(iout,*) ' # function evaluations : ' , nftot
         WRITE(iout,2003) timetot
         WRITE(iout,2005)
         DO i=1,n
            WRITE(iout,2007) i , S(i,1)
         END DO
         WRITE(iout,2011)

! Error exit
      CASE(1)
         WRITE(iout,*) ' *** ERROR:  n  must be > 0        *** '
         WRITE(iout,*) ' *** Please, correct and resubmit. *** '

      CASE(2)
         WRITE(iout,*) ' *** WARNING: maximum number of iterations ***   (',iter,')'
         WRITE(iout,2003) time

      CASE(3)
         WRITE(iout,*) ' *** WARNING: exceeded CPU time limit ***'
         WRITE(iout,2003) time

      CASE(4)
         WRITE(iout,*) ' *** ERROR: unable to generate n+1 different points ***'

      CASE(5)
         WRITE(iout,*) ' *** ERROR: exceeded number of constraints violations ***'

      CASE(6)
         WRITE(iout,*) ' *** WARNING: exceeded max number f.evaluations ***'

      CASE(10)
         WRITE(iout,*) ' *** ERROR: working set dimension too small ***'

      END SELECT
    endif

XSOL(1:N) = S(1:N,1)
FOUT      = FVAL(1)

if(outlev >= 0) write(iout,2014) timeset,timesol,timetot

1000 CONTINUE

SELECT CASE (iexit)
CASE(7)
   if(OUTLEV >= 0) then
      WRITE(iout,*) ' *** ERROR: some var has an improper lower bound ***'
      WRITE(iout,*) ' *** Please, correct and resubmit.               ***'
   endif
CASE(8)
   if(OUTLEV >= 0) then
      WRITE(iout,*) ' *** ERROR: some var has an improper upper bound ***'
      WRITE(iout,*) ' *** Please, correct and resubmit.               ***'
   endif
CASE(9)
   if(OUTLEV >= 0) then
      WRITE(iout,*) ' *** ERROR: constraints are too tight to be easily satisfied ***'
      WRITE(iout,*) ' *** Please, correct and resubmit.                           ***'
   endif
END SELECT

!---------------------------------------------------------------
! Non executable statements
!---------------------------------------------------------------

 2000 FORMAT(I10)
 2001 FORMAT(D20.10)
 2002 FORMAT(2D20.10)
 2003 FORMAT(/, '  Total time = ', F10.2, ' seconds' ,/ )
 2004 FORMAT(/, '  Final function value = ', ES22.14 ,/ )
 2005 FORMAT(/, '  Minimum point : ',/ )
 2006 FORMAT('  tolerance used in the (global) stopping criterion = ', ES10.2 )
 2007 FORMAT('   X _',I2,'   = ' ,D22.14)
 2008 FORMAT(/, ' ********************* FINAL RESULTS ********************' ,/)
 2009 FORMAT('  number of initial points used (m) = ' , I6)
 2010 FORMAT(/,'  dimension of the problem (n) = ', I3)
 2011 FORMAT(/, ' ********************************************************' ,/)
 2014 FORMAT('Setup time = ', F10.2, ' seconds' ,/ &
             'Solve time = ', F10.2, ' seconds' ,/ &
             'Total time = ', F10.2, ' seconds')
 2020 FORMAT('  & ', F10.2)
 2021 FORMAT('  & ', ES10.3, ' \\ ')
 2030 FORMAT(A40,' & ',I2,' & ',I6,' & ',I6, $)
 2050 FORMAT(' & ',F13.4,$)
 2060 FORMAT(' & ',D9.3,$)
 2070 FORMAT(' & ',I4,$)
 2031 FORMAT(' & ',I2,' & ',I6,' & ',I6,' & ',F13.4,' & ',F12.2,' \\')
 2040 FORMAT(A40)
 3000 FORMAT('   Iter      NF        FMIN          FMAX        DIFF')
 3010 FORMAT('-----------------------------------------------------------')
!3020 FORMAT(   123456 | 123456 | -1.1234D+01 | -1.1234D+01 | -1.12D+01)
 3020 FORMAT( 2X, I6,' | ',I6,' | ', ES11.4,' | ', ES11.4,' | ',ES9.2)

!------------------------------------------------------------------------------------
END SUBROUTINE ACRS

!------------------------------------------------------------------------------------
! This subroutine implements the base price iteration that:
! - random selection of N+1 points of S
! - computing of weigthed centroid
! - weigthed reflection
! - updating of S, if necessary
!------------------------------------------------------------------------------------
SUBROUTINE PRICE_BASE_ITERATION(N,M,S,FVAL,IEXIT,NFTOT,IOUT,OUTLEV,LB,UB,OMEGA,DIFF,DIFF_INITIAL,  &
                                PROJ,COLLIS_COUNT,REJ_COUNT,VIOL_COUNT,FUNCT                      )
   IMPLICIT NONE

   INTEGER,       INTENT(IN)     :: N, M, IOUT, OUTLEV
   INTEGER,       INTENT(OUT)    :: IEXIT
   INTEGER,       INTENT(INOUT)  :: NFTOT, COLLIS_COUNT, REJ_COUNT, VIOL_COUNT
   DOUBLE PRECISION, INTENT(IN)     :: LB(N), UB(N)
   LOGICAL,       INTENT(IN)     :: PROJ
   DOUBLE PRECISION, INTENT(INOUT)  :: S(N,M), FVAL(M), OMEGA, DIFF, DIFF_INITIAL
   EXTERNAL                   :: FUNCT

   DOUBLE PRECISION, PARAMETER      :: ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0
   INTEGER,       PARAMETER      :: P   = 2
   REAL                       :: RR
   DOUBLE PRECISION              :: R, F_MIN, F_MAX, PHI, SUM, CENT(N), F_W, ALPHA
   DOUBLE PRECISION              :: X_TRIAL(N), F, ftemp
   INTEGER                       :: IND(N+1), I, J, MASS(1), MASS_I, ierr
   LOGICAL                       :: MASK(M), COLLISION, VIOL

   COLLISION = .FALSE.
   VIOL      = .FALSE.
   F_MIN     = FVAL(1)
   F_MAX     = FVAL(M)
   mask      = .false.

!==============================================================
! Step 2 : Choose at random  n+1  points  over  S and determine
!          the centroid
!==============================================================

! Random choice of n+1 points

   IF ( OUTLEV >= 2 ) THEN
      WRITE(iout,*) ' Random choice of ', n+1 , ' points among ', m
   END IF

!--------------------------------------------------------------------------
   CALL RANDOM_NUMBER(rr)

   r = DBLE(rr)
   r = (p**r-1)/(p-1)

   IND(1)       = INT(m*r) + 1
   MASK(IND(1)) =.TRUE.

   DO i=2,N+1

300      CALL RANDOM_NUMBER(rr)

      r = DBLE(rr)
      r = (p**r-1)/(p-1)

      IND(i)       = INT(m*r) + 1
      MASK(IND(i)) =.TRUE.
      DO j=1,i-1

         IF ( IND(i) == IND(j) ) THEN

            IF ( .NOT. collision ) collis_count = 0
            collision = .TRUE.
            collis_count = collis_count + 1

            IF ( collis_count > 1000 * n ) THEN
               iexit = 4
               return
            END IF

            GO TO 300
         END IF

      END DO

   END DO

!--------------------------------------------------------------------------

   collision = .FALSE.


!---------------------------------------------------------------------------
! Determine the trial point according to the IMPROVED+MODIFIED Price scheme:
!---------------------------------------------------------------------------
!
! - Determine the weighted centroid

   mass=MAXLOC(FVAL,MASK)
   mass_i=mass(1)

   IF ( rej_count >= 10 * m ) omega = omega / two

   phi = omega * ( ( diff * diff ) /  diff_initial )
   IF  (( diff_initial >= 1.D+4 ) .AND. ( phi < 1.D-6 ) ) THEN
      diff_initial = 1.D+1 * omega * diff
      IF ( OUTLEV >= 1 ) THEN
         WRITE(iout,*) ' Modified diff_initial determining the centroid '
      END IF
   END IF

   sum = zero
   DO j=1,N+1
!     IF(IND(J).NE.MASS_I) sum = sum + ( one / ( FVAL(IND(j)) - f_min + phi ) )
      sum = sum + ( one / ( FVAL(IND(j)) - f_min + phi ) )
   END DO

   DO i=1,n
      CENT(i) = zero
      DO j=1,N+1
 !       IF(IND(J).NE.MASS_I) CENT(i) = CENT(i) + S(i,IND(j)) / ( FVAL(IND(j)) - f_min + phi )
         CENT(i) = CENT(i) + S(i,IND(j)) / ( FVAL(IND(j)) - f_min + phi )
      END DO
      CENT(i) = CENT(i) / sum
   END DO

! - Determine the trial point by a weighted reflection

   f_w = zero
   DO j=1,N+1
!     IF(IND(J).NE.MASS_I) f_w = f_w + FVAL(IND(j)) / ( FVAL(IND(j)) - f_min + phi )
      f_w = f_w + FVAL(IND(j)) / ( FVAL(IND(j)) - f_min + phi )
   END DO
   f_w = f_w / sum

   alpha = one - ( ( FVAL(mass_i) - f_w ) / ( f_max - f_min + phi ) )

   DO i=1,n
      X_trial(i) = ( one + alpha ) * CENT(i) - alpha * S(i,mass_i)
   END DO

   DO j=1,N+1
      MASK(IND(j))=.false.
   END DO

   IF ( OUTLEV >= 2 ) THEN
      DO i=1,n
         WRITE(iout,*) ' X_trial(',i,') = ' , X_trial(i)
      END DO
   END IF

! Check the consistency of the new trial point with the box constraints
!  if (.false.) then
   IF ( proj ) THEN

      DO i=1,n

         IF ( ( X_trial(i) > UB(i) ) ) THEN

            X_trial(i) = UB(i)

            IF ( OUTLEV >= 2 ) THEN
               WRITE(iout,*) ' trial point does not satisfy UB constraints: projected '
            END IF

         END IF

         IF ( ( X_trial(i) < LB(i) ) ) THEN

            X_trial(i) = LB(i)

            IF ( OUTLEV >= 2 ) THEN
               WRITE(iout,*) ' trial point does not satisfy LB constraints: projected '
            END IF

         END IF

      END DO

   ELSE
      DO i=1,n
         IF ( ( X_trial(i) > UB(i) + 1.0d+4*epsilon(alpha)) .OR. ( X_trial(i) < LB(i) - 1.0d+4*epsilon(alpha)) ) THEN

            IF ( OUTLEV >= 2 ) THEN
               WRITE(iout,*) ' the trial point does not satisfy box constraints'
               WRITE(iout,*) lb(i),x_trial(i),ub(i)
            END IF

            IF ( .NOT. viol ) viol_count = 0
            viol       =.TRUE.
            viol_count = viol_count + 1

            IF ( viol_count > 1000 * n ) THEN
               iexit = 5
               return
            END IF

            return

         END IF
      END DO

   END IF
!  endif
!--------------------------------------------------------------------------


   viol = .FALSE.

   CALL FUNCT(x_trial,n,f)
   nftot = nftot + 1

   IF ( OUTLEV >= 2 ) THEN
      WRITE(iout,*) ' f_trial = ', f
   END IF

!======================================================================
! Step 3-4 : Comparison of f(X_trial) with f_max
!======================================================================

   IF ( f >= f_max ) THEN
      IF ( OUTLEV >= 1 ) THEN
         WRITE(iout,*) ' trial point rejected '
      END IF

      rej_count = rej_count + 1

   ELSE
      rej_count = 0

      IF ( OUTLEV >= 1 ) THEN
         WRITE(iout,*) ' trial point accepted '
         WRITE(iout,*) diff
      END IF

      DO i=1,n
         S(i,m) = X_trial(i)
      END DO

      FVAL(m) = f
      CALL INSERT (FVAL,S,m,m,n)

      diff = fval(m) - fval(1)

   END IF

   RETURN
!------------------------------------------------------------------------------------
END SUBROUTINE PRICE_BASE_ITERATION

SUBROUTINE INSERT(FVAL,S,i,m,n)
INTEGER i,m,n
DOUBLE PRECISION FVAL(m),S(n,m)
DOUBLE PRECISION ftemp
DOUBLE PRECISION, ALLOCATABLE :: XTEMP(:)
INTEGER j,l

ALLOCATE(XTEMP(N))

DO k=1,n
     XTEMP(k)=S(k,i)
END DO
ftemp=FVAL(i)
DO j=1,i-1
   IF (FVAL(j)>FVAL(i)) THEN
      DO l=i,j+1,-1
        DO k=1,n
             S(k,l)=S(k,l-1)
          END DO
          FVAL(l)=FVAL(l-1)
       END DO
       FVAL(j)=ftemp
      DO k=1,n
        S(k,j)=XTEMP(k)
      END DO
       GOTO 1
      END IF
END DO
1 CONTINUE
DEALLOCATE(XTEMP)
RETURN
END

SUBROUTINE GENERATE_SET_S(MAXTGEN,N,M,NFTOT,OUTLEV,IOUT,LB,UB,TOL,S,FVAL,ERRFLAG,XIN,FUNCT)
   IMPLICIT NONE

   INTEGER           :: N, M, NFTOT, OUTLEV, IOUT, ERRFLAG, ierr
   DOUBLE PRECISION  :: LB(n), UB(m), TOL, S(N,M), FVAL(M), W1(N), XIN(N)
   DOUBLE PRECISION  :: ST(N,M)
   INTEGER           :: I, J, I_MAX(1), VIOL_COUNT, VIOL_TOT, IT
   REAL           :: RR, timegen, time, MAXTGEN
   DOUBLE PRECISION  :: R, F
   DOUBLE PRECISION  :: F_MAX, F_MIN, diff
   LOGICAL           :: VIOL
   EXTERNAL       :: FUNCT

   VIOL_TOT = 0
   CALL CPU_TIME(time)
   !time = cputim(dum)
   ST = 0.d0
        FVAL = 1.d+30
   if((errflag==0).or.(errflag==-1)) S = 0.d0

   IT = 1
   DO I = 1, M
      VIOL_COUNT = 0
90    CONTINUE
      CALL CPU_TIME(timegen)
      timegen = timegen - time
      IF((timegen > MAXTGEN).OR.(VIOL_TOT > 100000000)) THEN
         ERRFLAG = 9
         RETURN
      ENDIF

      if((errflag == -2).and.(it <= M)) then
         !use the ith point provided in set S
         DO J = 1, N
            ST(J,I) = S(J,IT)
            W1(J) = S(J,IT)
            IF ( OUTLEV >= 2 ) THEN
               WRITE(IOUT,*) ' S(',J,I,')= ' , S(J,I)
            END IF
         ENDDO
         IT = IT + 1
         VIOL = .FALSE.
         !---------------------------------
         ! check whether the user point
         ! satisfies the box constraints
         !---------------------------------
         DO j=1,n
            IF ( ( W1(j) > UB(j) ) .OR. ( W1(j) < LB(j) ) ) THEN

               IF ( OUTLEV >= 2 ) THEN
                  WRITE(iout,*) ' the',i,'th initial point does not satisfy box constraints'
               END IF

               VIOL = .TRUE.

            END IF
         END DO
         if(viol) then
            DO J = 1, N
               CALL RANDOM_NUMBER(RR)
               R = DBLE(RR)
               S(J,I) = ( UB(J) - LB(J) ) * R + LB(J)
               W1(J)  = S(J,I)
               IF ( OUTLEV >= 2 ) THEN
                  WRITE(IOUT,*) ' S(',J,I,')= ' , S(J,I)
               END IF
            END DO
         endif
      else
         DO J = 1, N
            CALL RANDOM_NUMBER(RR)
            R = DBLE(RR)
            S(J,I) = ( UB(J) - LB(J) ) * R + LB(J)
            W1(J)  = S(J,I)
            IF ( OUTLEV >= 2 ) THEN
               WRITE(IOUT,*) ' S(',J,I,')= ' , S(J,I)
            END IF
         END DO
        endif

!--------------------------------------------------------------------------
      CALL FUNCT ( W1,N,F)
      NFTOT = NFTOT + 1

      IF(i>1) diff = fval(i-1) - fval(1)

      FVAL(I) = F


      if(errflag == -2) then
         CALL INSERT (FVAL,ST,I,M,N)
      else
         CALL INSERT (FVAL,S,I,M,N)
      endif
                
        if(outlev >= 0) then
            write(iout,3020) 0,nftot,fval(1),fval(i),0.d0
        endif
   END DO

!---------------------------------------------------------------------------
! If errflag < 0 then we have to add the user point to the working set
!---------------------------------------------------------------------------
   IF(ERRFLAG == -1) THEN
      VIOL = .FALSE.
      !---------------------------------
      ! check whether the user point
      ! satisfies the box constraints
      !---------------------------------
      DO i=1,n
         IF ( ( XIN(i) > UB(i) ) .OR. ( XIN(i) < LB(i) ) ) THEN

            IF ( OUTLEV >= 2 ) THEN
               WRITE(iout,*) ' the initial point does not satisfy box constraints'
            END IF

            VIOL = .TRUE.

         END IF
      END DO
      IF(.NOT.VIOL) THEN
         !---------------------------------
         ! if the user point is feasible
         ! then add it to S
         !---------------------------------
         CALL FUNCT ( XIN,N,F )
         NFTOT = NFTOT + 1

         I_MAX = MAXLOC(FVAL)
         DO J = 1, N
            S(J,I_MAX(1)) = XIN(J)
         END DO
         FVAL(I_MAX(1)) = F
         CALL INSERT (FVAL,S,I_MAX(1),M,N)
            if(outlev >= 0) then
                write(iout,3020) 0,nftot,fval(1),fval(m),0.d0
            endif
      ENDIF
   ENDIF

   F_MAX = MAXVAL(FVAL)
   F_MIN = MINVAL(FVAL)

   DO WHILE (F_MAX-F_MIN <= 1.0D+0*TOL)
      VIOL_COUNT = 0
100      CONTINUE
      CALL CPU_TIME(timegen)
      timegen = timegen - time
      IF((timegen > MAXTGEN).OR.(VIOL_TOT > 100000000)) THEN
         RETURN
      ENDIF
      DO J = 1, N
         CALL RANDOM_NUMBER(RR)
         R = DBLE(RR)
         W1(J) = ( UB(J) - LB(J) ) * R + LB(J)
      END DO

!--------------------------------------------------------------------------
      CALL FUNCT ( W1,N,F )
      NFTOT = NFTOT + 1

      IF(F<F_MAX) THEN
         I_MAX = MAXLOC(FVAL)

         DO J = 1, N
            S(J,I_MAX(1)) = W1(J)
         END DO
         FVAL(I_MAX(1)) = F
         if(errflag == -2) then
            CALL INSERT (FVAL,ST,I_MAX(1),M,N)
         else
            CALL INSERT (FVAL,S,I_MAX(1),M,N)
         endif
         F_MAX = MAXVAL(FVAL)
         F_MIN = MINVAL(FVAL)
            if(outlev >= 0) then
                write(iout,3020) 0,nftot,fval(1),fval(m),0.d0
            endif
      ENDIF
   END DO

   if(errflag == -2) then
      S = ST
   endif

   RETURN

3020 FORMAT( 2X, I6,' | ',I6,' | ', ES11.4,' | ', ES11.4,' | ',ES9.2)

END SUBROUTINE GENERATE_SET_S
