#include <math.h>

#include "minpack.h"

 /*
  *     subroutine fdjac2 
  *
  *     this subroutine computes a forward-difference approximation 
  *     to the m by n jacobian matrix associated with a specified 
  *     problem of m functions in n variables. 
  *
  *     the subroutine statement is 
  *
  *       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa) 
  *
  *     where 
  *
  *       fcn is the name of the user-supplied subroutine which 
  *         calculates the functions. fcn must be declared 
  *         in an external statement in the user calling 
  *         program, and should be written as follows. 
  *
  *         subroutine fcn(m,n,x,fvec,iflag) 
  *         int m,n,iflag 
  *         double precision x(n),fvec(m) 
  *         ---------- 
  *         calculate the functions at x and 
  *         return this vector in fvec. 
  *         ---------- 
  *         return 
  *         end 
  *
  *         the value of iflag should not be changed by fcn unless 
  *         the user wants to terminate execution of fdjac2. 
  *         in this case set iflag to a negative int. 
  *
  *       m is a positive int input variable set to the number 
  *         of functions. 
  *
  *       n is a positive int input variable set to the number 
  *         of variables. n must not exceed m. 
  *
  *       x is an input array of length n. 
  *
  *       fvec is an input array of length m which must contain the 
  *         functions evaluated at x. 
  *
  *       fjac is an output m by n array which contains the 
  *         approximation to the jacobian matrix evaluated at x. 
  *
  *       ldfjac is a positive int input variable not less than m 
  *         which specifies the leading dimension of the array fjac. 
  *
  *       iflag is an int variable which can be used to terminate 
  *         the execution of fdjac2. see description of fcn. 
  *
  *       epsfcn is an input variable used in determining a suitable 
  *         step length for the forward-difference approximation. this 
  *         approximation assumes that the relative errors in the 
  *         functions are of the order of epsfcn. if epsfcn is less 
  *         than the machine precision, it is assumed that the relative 
  *         errors in the functions are of the order of the machine 
  *         precision. 
  *
  *       wa is a work array of length m. 
  *
  *     subprograms called 
  *
  *       user-supplied ...... fcn 
  *
  *       minpack-supplied ... dpmpar 
  *
  *       fortran-supplied ... dabs,dmax1,dsqrt 
  *
  *     argonne national laboratory. minpack project. march 1980. 
  *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
  */

void fdjac2_c(minpack_func_mn fcn, const int *m, const int *n, double *x,
	    double const *fvec, double *fjac, const int *ldfjac, int *iflag,
	    const double *epsfcn, double *wa)
{
	int fjac_dim1 = *ldfjac;

	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;

	double eps = sqrt(fmax(*epsfcn, epsmch));
	for (long int j = 0; j < *n; ++j) {
		double temp = x[j];
		double h = eps * fabs(x[j]);
		if (h == 0)
			h = eps;
		x[j] += h;
		(*fcn)(m, n, x, wa, iflag);
		if (*iflag < 0)
			return;
		x[j] = temp;              // restore x[j]
		for (long int i = 0; i < *m; ++i)
			fjac[i + j * fjac_dim1] = (wa[i] - fvec[i]) / h;
	}
}
