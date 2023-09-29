#include "minpack.h"

/*     Subroutine lmdif1
 *
 *     The purpose of lmdif1 is to minimize the sum of the squares of
 *     m nonlinear functions in n variables by a modification of the
 *     levenberg-marquardt algorithm.  This is done by using the more
 *     general least-squares solver lmdif.  The user must provide a
 *     subroutine which calculates the functions.  The jacobian is
 *     then calculated by a forward-difference approximation.
 *
 *     Rhe subroutine statement is
 *
 *       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
 *
 *     where
 *
 *       fcn is the name of the user-supplied subroutine which
 *         calculates the functions.  fcn must be declared
 *         in an external statement in the user calling
 *         program, and should be written as follows.
 *
 *         subroutine fcn(m,n,x,fvec,iflag)
 *         integer m,n,iflag
 *         double precision x(n),fvec(m)
 *         ----------
 *         calculate the functions at x and
 *         return this vector in fvec.
 *         ----------
 *         return
 *         end
 *
 *         The value of iflag should not be changed by fcn unless
 *         the user wants to terminate execution of lmdif1.
 *         In this case set iflag to a negative integer.
 *
 *       m is a positive integer input variable set to the number
 *         of functions.
 *
 *       n is a positive integer input variable set to the number
 *         of variables.  n must not exceed m.
 *
 *       x is an array of length n. on input x must contain
 *         an initial estimate of the solution vector.  On output x
 *         contains the final estimate of the solution vector.
 *
 *       fvec is an output array of length m which contains
 *         the functions evaluated at the output x.
 *qtf
 *       tol is a nonnegative input variable.  Termination occurs
 *         when the algorithm estimates either that the relative
 *         error in the sum of squares is at most tol or that
 *         the relative error between x and the solution is at
 *         most tol.
 *
 *       info is an integer output variable.  If the user has
 *         terminated execution, info is set to the (negative)
 *         value of iflag.  See description of fcn.  otherwise,
 *         info is set as follows.
 *
 *         info = 0  improper input parameters.
 *
 *         info = 1  algorithm estimates that the relative error
 *                   in the sum of squares is at most tol.
 *
 *         info = 2  algorithm estimates that the relative error
 *                   between x and the solution is at most tol.
 *
 *         info = 3  conditions for info = 1 and info = 2 both hold.
 *
 *         info = 4  fvec is orthogonal to the columns of the
 *                   jacobian to machine precision.
 *
 *         info = 5  number of calls to fcn has reached or
 *                   exceeded 200*(n+1).
 *
 *         info = 6  tol is too small. no further reduction in
 *                   the sum of squares is possible.
 *
 *         info = 7  tol is too small. no further improvement in
 *                   the approximate solution x is possible.
 *
 *       iwa is an integer work array of length n.
 *
 *       wa is a work array of length lwa.
 *
 *       lwa is a positive integer input variable not less than
 *         m*n+5*n+m.
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

void lmdif1_c(minpack_func_mn fcn, int const *m, int const *n, double *x, 
	double *fvec, double const *tol, int *info, int *iwa, 
	double *wa, int const *lwa)
{
	*info = 0;

	/* check the input parameters for errors. */
	if (*n <= 0 || *m < *n || *tol < 0 || *lwa < *m * *n + *n * 5 + *m)
		return;

	/* call lmdif. */
    double factor = 100.;
	int maxfev = (*n + 1) * 200;
	double ftol = *tol;
	double xtol = *tol;
	double gtol = 0;
	double epsfcn = 0;
	int mode = 1;
	int nprint = 1;
	int mp5n = *m + *n * 5;
    int nfev = 0;
	lmdif_c(fcn, m, n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn, 
        wa, &mode, &factor, &nprint, info, &nfev, &wa[mp5n], m, iwa, 
        &wa[*n], &wa[*n * 2], &wa[*n * 3], &wa[*n * 4], &wa[*n * 5]);
	if (*info == 8)
		*info = 4;
}
