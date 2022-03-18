#include <stddef.h>

#include "minpack.h"

/*
 *     Subroutine hybrj
 *
 *     The purpose of hybrj is to find a zero of a system of 
 *     n nonlinear functions in n variables by a modification 
 *     of the Powell hybrid method.  The user must provide a 
 *     subroutine which calculates the functions and the jacobian. 
 *
 *     The subroutine statement is 
 *
 *       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag, 
 *                        mode,factor,nprint,info,nfev,njev,r,lr,qtf, 
 *                        wa1,wa2,wa3,wa4) 
 *
 *     Where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions and the jacobian.  fcn must 
 *         be declared in an external statement in the user 
 *         calling program, and should be written as follows. 
 *
 *         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag) 
 *         integer n,ldfjac,iflag 
 *         double precision x(n),fvec(n),fjac(ldfjac,n) 
 *         ---------- 
 *         If iflag = 1 calculate the functions at x and 
 *         return this vector in fvec.  Do not alter fjac. 
 *         If iflag = 2 calculate the jacobian at x and 
 *         return this matrix in fjac.  Do not alter fvec. 
 *         --------- 
 *         return 
 *         end 
 *
 *         The value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of hybrj. 
 *         In this case set iflag to a negative integer. 
 *
 *       n is a positive integer input variable set to the number 
 *         of functions and variables. 
 *
 *       x is an array of length n.  On input x must contain 
 *         an initial estimate of the solution vector.  On output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length n which contains 
 *         the functions evaluated at the output x. 
 *
 *       fjac is an output n by n array which contains the 
 *         orthogonal matrix q produced by the QR factorization 
 *         of the final approximate jacobian. 
 *
 *       ldfjac is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       xtol is a nonnegative input variable.  Termination 
 *         occurs when the relative error between two consecutive 
 *         iterates is at most xtol. 
 *
 *       maxfev is a positive integer input variable.  Termination 
 *         occurs when the number of calls to fcn with iflag = 1 
 *         has reached maxfev. 
 *
 *       diag is an array of length n.  If mode = 1 (see 
 *         below), diag is internally set.  If mode = 2, diag 
 *         must contain positive entries that serve as 
 *         multiplicative scale factors for the variables. 
 *
 *       mode is an integer input variable.  If mode = 1, the 
 *         variables will be scaled internally.  If mode = 2, 
 *         the scaling is specified by the input diag.  Other 
 *         values of mode are equivalent to mode = 1. 
 *
 *       factor is a positive input variable used in determining the 
 *         initial step bound.  This bound is set to the product of 
 *         factor and the euclidean norm of diag*x if nonzero, or else 
 *         to factor itself. in most cases factor should lie in the 
 *         interval (0.1, 100). 100 is a generally recommended value. 
 *
 *       nprint is an integer input variable that enables controlled 
 *         printing of iterates if it is positive.  In this case, 
 *         fcn is called with iflag = 0 at the beginning of the first 
 *         iteration and every nprint iterations thereafter and 
 *         immediately prior to return, with x and fvec available 
 *         for printing.  fvec and fjac should not be altered. 
 *         if nprint is not positive, no special calls of fcn 
 *         with iflag = 0 are made. 
 *
 *       info is an integer output variable.  If the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag. see description of fcn.  Otherwise, 
 *         info is set as follows. 
 *
 *         info = 0   improper input parameters. 
 *
 *         info = 1   relative error between two consecutive iterates 
 *                    is at most xtol. 
 *
 *         info = 2   number of calls to fcn with iflag = 1 has 
 *                    reached maxfev. 
 *
 *         info = 3   xtol is too small.  No further improvement in 
 *                    the approximate solution x is possible. 
 *
 *         info = 4   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    five jacobian evaluations. 
 *
 *         info = 5   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    ten iterations. 
 *
 *       nfev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 1. 
 *
 *       njev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 2. 
 *
 *       r is an output array of length lr which contains the 
 *         upper triangular matrix produced by the QR factorization 
 *         of the final approximate jacobian, stored rowwise. 
 *
 *       lr is a positive integer input variable not less than 
 *         (n*(n+1))/2. 
 *
 *       qtf is an output array of length n which contains 
 *         the vector (q transpose)*fvec. 
 *
 *       wa1, wa2, wa3, and wa4 are work arrays of length n. 
 * 
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
 */
void hybrj_c(minpack_func_nj fcn, const int *n, double *x, double *fvec,
	   double *fjac, const int *ldfjac, const double *xtol, const int *maxfev,
	   double *diag, const int *mode, const double *factor, const int *nprint,
	   int *info, int *nfev, int *njev, double *r, const int *lr, double *qtf, 
	   double *wa1, double *wa2, double *wa3, double *wa4)
{
	int mzero = 0;
	hybrbase(NULL, fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, &mzero, &mzero, NULL, diag, mode,
	      factor, nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, wa3, wa4);
}
