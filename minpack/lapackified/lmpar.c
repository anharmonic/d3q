#include <math.h>

#include "minpack.h"

/*
 *     subroutine lmpar 
 *
 *     given an m by n matrix A, an n by n nonsingular diagonal 
 *     matrix D, an m-vector b, and a positive number delta, 
 *     the problem is to determine a value for the parameter 
 *     par such that if x solves the system 
 *
 *           A*x = b ,     sqrt(par)*D*x = 0 , 
 *
 *     in the least squares sense, and dxnorm is the euclidean 
 *     norm of D*x, then either par is zero and 
 *
 *           (dxnorm-delta) <= 0.1*delta, 
 *
 *     or par is positive and 
 *
 *           abs(dxnorm-delta) <= 0.1*delta. 
 *
 *     This subroutine completes the solution of the problem 
 *     if it is provided with the necessary information from the 
 *     QR factorization, with column pivoting, of A. That is, if 
 *     A*P = Q*R, where P is a permutation matrix, Q has orthogonal 
 *     columns, and R is an upper triangular matrix with diagonal 
 *     elements of nonincreasing magnitude, then lmpar expects 
 *     the full upper triangle of R, the permutation matrix P, 
 *     and the first n components of (Q transpose)*b. On output 
 *     lmpar also provides an upper triangular matrix S such that 
 *
 *            t   t                   t 
 *           P *(A *A + par*D*D)*P = S *S. 
 *
 *     S is employed within lmpar and may be of separate interest. 
 *
 *     only a few iterations are generally needed for convergence 
 *     of the algorithm. if, however, the limit of 10 iterations 
 *     is reached, then the output par will contain the best 
 *     value obtained so far. 
 *
 *     the subroutine statement is 
 *
 *       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, 
 *                        wa1,wa2) 
 *
 *     where 
 *
 *       n is a positive int input variable set to the order of r. 
 *
 *       r is an n by n array. on input the full upper triangle 
 *         must contain the full upper triangle of the matrix r. 
 *         on output the full upper triangle is unaltered, and the 
 *         strict lower triangle contains the strict upper triangle 
 *         (transposed) of the upper triangular matrix s. 
 *
 *       ldr is a positive int input variable not less than n 
 *         which specifies the leading dimension of the array r. 
 *
 *       ipvt is an int input array of length n which defines the 
 *         permutation matrix p such that a*p = q*r. column j of p 
 *         is column ipvt(j) of the identity matrix. 
 *
 *       diag is an input array of length n which must contain the 
 *         diagonal elements of the matrix d. 
 *
 *       qtb is an input array of length n which must contain the first 
 *         n elements of the vector (q transpose)*b. 
 *
 *       delta is a positive input variable which specifies an upper 
 *         bound on the euclidean norm of d*x. 
 *
 *       par is a nonnegative variable. on input par contains an 
 *         initial estimate of the levenberg-marquardt parameter. 
 *         on output par contains the final estimate. 
 *
 *       x is an output array of length n which contains the least 
 *         squares solution of the system a*x = b, sqrt(par)*d*x = 0, 
 *         for the output par. 
 *
 *       sdiag is an output array of length n which contains the 
 *         diagonal elements of the upper triangular matrix s. 
 *
 *       wa1 and wa2 are work arrays of length n. 
 *
 *     subprograms called 
 *
 *       minpack-supplied ... dpmpar,enorm,qrsolv 
 *
 *       fortran-supplied ... dabs,dmax1,dmin1,dsqrt 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */

void lmpar_(const int *n, double *r, const int *ldr, int *ipvt, double *diag, 
	double *qtb, double *delta, double *par, double *x, double *sdiag, double *wa1, double *wa2)
{
	/* dwarf is the smallest positive magnitude. */
	double dwarf = MINPACK_DWARF;
	int r_dim1 = *ldr;

	/* Compute and store in x the gauss-newton direction. If the
	   jacobian is rank-deficient, obtain a least squares solution. */
	int nsing = *n;
	for (long int j = 0; j < *n; ++j) {
		wa1[j] = qtb[j];
		if (r[j + j * r_dim1] == 0)
			nsing = j;
		if (nsing < *n)
			wa1[j] = 0;
	}

	int c1 = 1;
	dtrsv_("Up", "notrans", "nounit", &nsing, r, ldr, wa1, &c1);
	for (long int j = 0; j < *n; ++j) {
		int l = ipvt[j] - 1;
		x[l] = wa1[j];
	}

	/* initialize the iteration counter.
	   evaluate the function at the origin, and test
	   for acceptance of the gauss-newton direction. */
	int iter = 0;
	for (long int j = 0; j < *n; ++j)
		wa2[j] = diag[j] * x[j];
	double dxnorm = enorm_(n, wa2);
	double fp = dxnorm - *delta;
	if (fp <= 0.1 * *delta)
		goto fini;

	/* if the jacobian is not rank deficient, the newton
	   step provides a lower bound, parl, for the zero of
	   the function.  Otherwise set this bound to zero. */
	double parl = 0;
	if (nsing >= *n) {
		for (long int j = 0; j < *n; ++j) {
			int l = ipvt[j] - 1;
			wa1[j] = diag[l] * (wa2[l] / dxnorm);
		}
		dtrsv_("Up", "trans", "nounit", &nsing, r, ldr, wa1, &c1);
		double temp = enorm_(n, wa1);
		parl = fp / *delta / temp / temp;
	}

	/* Calculate an upper bound, paru, for the zero of the function. */
	for (long int j = 0; j < *n; ++j) {
		double sum = 0;
		for (long int i = 0; i < j; ++i)
			sum += r[i + j * r_dim1] * qtb[i];
		int l = ipvt[j] - 1;
		wa1[j] = sum / diag[l];
	}
	double gnorm = enorm_(n, wa1);
	double paru = gnorm / *delta;
	if (paru == 0)
		paru = dwarf / fmin(*delta, 0.1);

	/* if the input par lies outside of the interval (parl,paru),
	   set par to the closer endpoint. */
	*par = fmax(*par, parl);
	*par = fmin(*par, paru);
	if (*par == 0)
		*par = gnorm / dxnorm;

	/* beginning of an iteration. */
	for (;;) {
		++iter;
		/* evaluate the function at the current value of par. */
		if (*par == 0)
			*par = fmax(dwarf, 0.001 * paru);
		double temp = sqrt(*par);
		for (long int j = 0; j < *n; ++j)
			wa1[j] = temp * diag[j];
		qrsolv_(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
		for (long int j = 0; j < *n; ++j)
			wa2[j] = diag[j] * x[j];
		dxnorm = enorm_(n, wa2);
		temp = fp;
		fp = dxnorm - *delta;

		/* if the function is small enough, accept the current value
		   of par. also test for the exceptional cases where parl
		   is zero or the number of iterations has reached 10. */
		if (fabs(fp) <= 0.1 * *delta || (parl == 0 && fp <= temp && temp < 0) || iter == 100)
			break;

		/* compute the newton correction. */
		for (long int j = 0; j < *n; ++j) {
			int l = ipvt[j] - 1;
			wa1[j] = diag[l] * (wa2[l] / dxnorm);
		}
		/* cannot use dtrsv here because the diagonal is in sdiag, not in r */
		for (long int j = 0; j < *n; ++j) {
			wa1[j] /= sdiag[j];
			for (long int i = j + 1; i < *n; ++i)
				wa1[i] -= r[i + j * r_dim1] *  wa1[j];
		}
		temp = enorm_(n, wa1);
		double parc = fp / *delta / temp / temp;

		/* Depending on the sign of the function, update parl or paru. */
		if (fp > 0)
			parl = fmax(parl, *par);
		if (fp < 0)
			paru = fmin(paru, *par);

		/* Compute an improved estimate for par. */
		*par = fmax(parl, *par + parc);

	}
 fini:
	/* termination. */
	if (iter == 0)
		*par = 0;
}
