#include <math.h>

#include "minpack.h"

/*
 *     subroutine fdjac1 
 *
 *     this subroutine computes a forward-difference approximation 
 *     to the n by n jacobian matrix associated with a specified 
 *     problem of n functions in n variables. if the jacobian has 
 *     a banded form, then function evaluations are saved by only 
 *     approximating the nonzero terms. 
 *
 *     the subroutine statement is 
 *
 *       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn, 
 *                         wa1,wa2) 
 *
 *     where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions. fcn must be declared 
 *         in an external statement in the user calling 
 *         program, and should be written as follows. 
 *
 *         subroutine fcn(n,x,fvec,iflag) 
 *         integer n,iflag 
 *         double precision x(n),fvec(n) 
 *         ---------- 
 *         calculate the functions at x and 
 *         return this vector in fvec. 
 *         ---------- 
 *         return 
 *         end 
 *
 *         the value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of fdjac1. 
 *         in this case set iflag to a negative integer. 
 *
 *       n is a positive integer input variable set to the number 
 *         of functions and variables. 
 *
 *       x is an input array of length n. 
 *
 *       fvec is an input array of length n which must contain the 
 *         functions evaluated at x. 
 *
 *       fjac is an output n by n array which contains the 
 *         approximation to the jacobian matrix evaluated at x. 
 *
 *       ldfjac is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       iflag is an integer variable which can be used to terminate 
 *         the execution of fdjac1. see description of fcn. 
 *
 *       ml is a nonnegative integer input variable which specifies 
 *         the number of subdiagonals within the band of the 
 *         jacobian matrix. if the jacobian is not banded, set 
 *         ml to at least n - 1. 
 *
 *       epsfcn is an input variable used in determining a suitable 
 *         step length for the forward-difference approximation. this 
 *         approximation assumes that the relative errors in the 
 *         functions are of the order of epsfcn. if epsfcn is less 
 *         than the machine precision, it is assumed that the relative 
 *         errors in the functions are of the order of the machine 
 *         precision. 
 *
 *       mu is a nonnegative integer input variable which specifies 
 *         the number of superdiagonals within the band of the 
 *         jacobian matrix. if the jacobian is not banded, set 
 *         mu to at least n - 1. 
 *
 *       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at 
 *         least n, then the jacobian is considered dense, and wa2 is 
 *         not referenced. 
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */
void fdjac1_(minpack_func_n fcn, const int *n, double *x, double *fvec,
	    double *fjac, const int *ldfjac, int *iflag, const int *ml,
	    const int *mu, const double *epsfcn, double *wa1, double *wa2)
{
	(void) ml;
    (void) mu;
    (void) wa2;
    int fjac_dim1 = *ldfjac;

	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;
	double eps = sqrt(fmax(*epsfcn, epsmch));
	
	/* computation of dense approximate jacobian. */
	for (long int j = 0; j < *n; ++j) {
		double temp = x[j];
		double h = eps * fabs(temp);
		if (h == 0)
			h = eps;
		x[j] = temp + h;
		(*fcn) (n, x, wa1, iflag);
		if (*iflag < 0)
			return;
		x[j] = temp;
		for (long int i = 0; i < *n; ++i)
			fjac[i + j * fjac_dim1] = (wa1[i] - fvec[i]) / h;
	}
	return;
}