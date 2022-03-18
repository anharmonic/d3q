#include "minpack.h"

/*
 *     Subroutine hybrj1 
 *
 *     The purpose of hybrj1 is to find a zero of a system of 
 *     n nonlinear functions in n variables by a modification 
 *     of the powell hybrid method.  This is done by using the 
 *     more general nonlinear equation solver hybrj.  The user 
 *     must provide a subroutine which calculates the functions 
 *     and the jacobian. 
 *
 *     The subroutine statement is 
 *
 *       subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa) 
 *
 *     where 
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
 *         the user wants to terminate execution of hybrj1. 
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
 *         orthogonal matrix Q produced by the QR factorization 
 *         of the final approximate jacobian. 
 *
 *       ldfjac is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       tol is a nonnegative input variable.  Termination occurs 
 *         when the algorithm estimates that the relative error 
 *         between x and the solution is at most tol. 
 *
 *       info is an integer output variable.  If the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag.  See description of fcn.  Otherwise, 
 *         info is set as follows. 
 *
 *         info = 0   improper input parameters. 
 *
 *         info = 1   algorithm estimates that the relative error 
 *                    between x and the solution is at most tol. 
 *
 *         info = 2   number of calls to fcn with iflag = 1 has 
 *                    reached 100*(n+1). 
 *
 *         info = 3   tol is too small.  No further improvement in 
 *                    the approximate solution x is possible. 
 *
 *         info = 4   iteration is not making good progress. 
 *
 *       wa is a work array of length lwa. 
 *
 *       lwa is a positive integer input variable not less than 
 *         (n*(n+13))/2. 
 *
 *     Argonne National Laboratory. MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
 */
void hybrj1_c(minpack_func_nj fcn, const int *n, double *x, double *fvec, double *fjac,
	     const int *ldfjac, const double *tol, int *info, double *wa, const int *lwa)
{
	*info = 0;

	/* check the input parameters for errors. */
	if (*n <= 0 || *ldfjac < *n || *tol < 0 || *lwa < *n * (*n + 13) / 2)
		return;

	/* Call hybrj */
	double factor = 100;
	int maxfev = (*n + 1) * 100;
	double xtol = *tol;
	int mode = 2;
	int nprint = 0;
	int lr = *n * (*n + 1) / 2;
	int nfev = 0;
	int njev = 0;
    double * diag = wa;
    for (long int j = 0; j < *n; ++j)
        diag[j] = 1;
    double *qtf = &wa[*n];
    double * wa1 = &wa[*n * 2];
    double * wa2 = &wa[*n * 3];
    double * wa3 = &wa[*n * 4];
    double * wa4 = &wa[*n * 5];
    double * r = &wa[*n * 6];
	
    hybrj_c(fcn, n, x, fvec, fjac, ldfjac, &xtol, &maxfev, wa, &mode, &factor, &nprint, info,
	       &nfev, &njev, r, &lr, qtf, wa1, wa2, wa3, wa4);
	
    if (*info == 5)
		*info = 4;
}
