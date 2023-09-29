#include <math.h>

#include "minpack.h"

/*
 *     Subroutine chkder 
 *
 *     This subroutine checks the gradients of m nonlinear functions 
 *     in n variables, evaluated at a point x, for consistency with 
 *     the functions themselves.  The user must call chkder twice, 
 *     first with mode = 1 and then with mode = 2. 
 *
 *     mode = 1. On input, x must contain the point of evaluation. 
 *               On output, xp is set to a neighboring point. 
 *
 *     mode = 2. On input, fvec must contain the functions and the 
 *                         rows of fjac must contain the gradients 
 *                         of the respective functions each evaluated 
 *                         at x, and fvecp must contain the functions 
 *                         evaluated at xp. 
 *               On output, err contains measures of correctness of 
 *                          the respective gradients. 
 *
 *     The subroutine does not perform reliably if cancellation or 
 *     rounding errors cause a severe loss of significance in the 
 *     evaluation of a function.  Therefore, none of the components 
 *     of x should be unusually small (in particular, zero) or any 
 *     other value which may cause loss of significance. 
 *
 *     The subroutine statement is 
 *
 *       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err) 
 *
 *     where 
 *
 *       m is a positive int input variable set to the number 
 *         of functions. 
 *
 *       n is a positive int input variable set to the number 
 *         of variables. 
 *
 *       x is an input array of length n. 
 *
 *       fvec is an array of length m.  On input when mode = 2, 
 *         fvec must contain the functions evaluated at x. 
 *
 *       fjac is an m by n array.  On input when mode = 2, 
 *         the rows of fjac must contain the gradients of 
 *         the respective functions evaluated at x. 
 *
 *       ldfjac is a positive int input parameter not less than m 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       xp is an array of length n.  On output when mode = 1, 
 *         xp is set to a neighboring point of x. 
 *
 *       fvecp is an array of length m.  On input when mode = 2, 
 *         fvecp must contain the functions evaluated at xp. 
 *
 *       mode is an int input variable set to 1 on the first call 
 *         and 2 on the second.  Other values of mode are equivalent 
 *         to mode = 1. 
 *
 *       err is an array of length m.  On output when mode = 2, 
 *         err contains measures of correctness of the respective 
 *         gradients.  If there is no severe loss of significance, 
 *         then if err(i) is 1.0 the i-th gradient is correct, 
 *         while if err(i) is 0.0 the i-th gradient is incorrect. 
 *         For values of err between 0.0 and 1.0, the categorization 
 *         is less certain.  In general, a value of err(i) greater 
 *         than 0.5 indicates that the i-th gradient is probably 
 *         correct, while a value of err(i) less than 0.5 indicates 
 *         that the i-th gradient is probably incorrect. 
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */
void chkder_c(const int *m, const int *n, double *x, double *fvec, double *fjac, 
	const int *ldfjac, double *xp, double *fvecp, const int *mode, double *err)
{
	int fjac_dim1 = *ldfjac;

	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;
	double eps = sqrt(epsmch);
	double factor = 100;
	double epsf = factor * epsmch;
	double epslog = log10(eps);

	switch (*mode) {
	case 1:
		for (long int j = 0; j < *n; ++j) {
			double temp = eps * fabs(x[j]);
			if (temp == 0)
				temp = eps;
			xp[j] = x[j] + temp;
		}
		return;

	case 2:
		for (long int i = 0; i < *m; ++i)
			err[i] = 0;
		for (long int j = 0; j < *n; ++j) {
			double temp = fabs(x[j]);
			if (temp == 0)
				temp = 1;
			for (long int i = 0; i < *m; ++i)
				err[i] += temp * fjac[i + j * fjac_dim1];
		}
		for (long int i = 0; i < *m; ++i) {
			double temp = 1;
			double d2 = fvecp[i] - fvec[i];
			if (fvec[i] != 0 && fvecp[i] != 0 && fabs(d2) >= epsf * fabs(fvec[i])) {
				double d3 = d2 / eps - err[i];
				temp = eps * fabs(d3) / (fabs(fvec[i]) + fabs(fvecp[i]));
			}
			err[i] = 1;
			if (temp > epsmch && temp < eps)
				err[i] = (log10(temp) - epslog) / epslog;
			if (temp >= eps)
				err[i] = 0;
		}
	}
}
