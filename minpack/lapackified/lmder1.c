#include "minpack.h"

/*
 *     subroutine lmder1 
 *
 *     the purpose of lmder1 is to minimize the sum of the squares of 
 *     m nonlinear functions in n variables by a modification of the 
 *     levenberg-marquardt algorithm. this is done by using the more 
 *     general least-squares solver lmder. the user must provide a 
 *     subroutine which calculates the functions and the jacobian. 
 *
 *     the subroutine statement is 
 *
 *       subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info, 
 *                         ipvt,wa,lwa) 
 *
 *     where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions and the jacobian. fcn must 
 *         be declared in an external statement in the user 
 *         calling program, and should be written as follows. 
 *
 *         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag) 
 *         integer m,n,ldfjac,iflag 
 *         double precision x(n),fvec(m),fjac(ldfjac,n) 
 *         ---------- 
 *         if iflag = 1 calculate the functions at x and 
 *         return this vector in fvec. do not alter fjac. 
 *         if iflag = 2 calculate the jacobian at x and 
 *         return this matrix in fjac. do not alter fvec. 
 *         ---------- 
 *         return 
 *         end 
 *
 *         the value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of lmder1. 
 *         in this case set iflag to a negative integer. 
 *
 *       m is a positive integer input variable set to the number 
 *         of functions. 
 *
 *       n is a positive integer input variable set to the number 
 *         of variables. n must not exceed m. 
 *
 *       x is an array of length n. on input x must contain 
 *         an initial estimate of the solution vector. on output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length m which contains 
 *         the functions evaluated at the output x. 
 *
 *       fjac is an output m by n array. the upper n by n submatrix 
 *         of fjac contains an upper triangular matrix r with 
 *         diagonal elements of nonincreasing magnitude such that 
 *
 *                t     t           t 
 *               p *(jac *jac)*p = r *r, 
 *
 *         where p is a permutation matrix and jac is the final 
 *         calculated jacobian. column j of p is column ipvt(j) 
 *         (see below) of the identity matrix. the lower trapezoidal 
 *         part of fjac contains information generated during 
 *         the computation of r. 
 *
 *       ldfjac is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       tol is a nonnegative input variable. termination occurs 
 *         when the algorithm estimates either that the relative 
 *         error in the sum of squares is at most tol or that 
 *         the relative error between x and the solution is at 
 *         most tol. 
 *
 *       info is an integer output variable. if the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag. see description of fcn. otherwise, 
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
 *         info = 5  number of calls to fcn with iflag = 1 has 
 *                   reached 100*(n+1). 
 *
 *         info = 6  tol is too small. no further reduction in 
 *                   the sum of squares is possible. 
 *
 *         info = 7  tol is too small. no further improvement in 
 *                   the approximate solution x is possible. 
 *
 *       ipvt is an integer output array of length n. ipvt 
 *         defines a permutation matrix p such that jac*p = q*r, 
 *         where jac is the final calculated jacobian, q is 
 *         orthogonal (not stored), and r is upper triangular 
 *         with diagonal elements of nonincreasing magnitude. 
 *         column j of p is column ipvt(j) of the identity matrix. 
 *
 *       wa is a work array of length lwa. 
 *
 *       lwa is a positive integer input variable not less than 5*n+m. 
 *
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth e. Hillstrom, Jorge J. Moré 
 */
void lmder1_(minpack_func_mnj fcn, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *tol, 
	int *info, int *ipvt, double *wa, const int *lwa)
{
	*info = 0;

	/* check the input parameters for errors. */
	if (*n <= 0 || *m < *n || *ldfjac < *m || *tol < 0 || *lwa < *n * 5 + *m)
	return;

	/* Call lmder. */
	int maxfev = (*n + 1) * 100;
	double ftol = *tol;
	double xtol = *tol;
	double gtol = 0;
	int mode = 1;
	int nprint = 0;
	int nfev = 0;
	int njev = 0;
	double factor = 100.;
	lmder_c(fcn, m, n, x, fvec, fjac, ldfjac, &ftol, &xtol, &gtol, &maxfev, 
		wa, &mode, &factor, &nprint, info, &nfev, &njev, ipvt, &wa[*n], &wa[*n * 2], 
		&wa[*n * 3], &wa[*n * 4], &wa[*n * 5]);
	if (*info == 8)
	*info = 4;
}