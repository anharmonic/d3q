#include <math.h>

#include "minpack.h"
/*
 *     subroutine dogleg 
 *
 *     Given an m by n matrix A, an n by n nonsingular diagonal 
 *     matrix D, an m-vector b, and a positive number delta, the 
 *     problem is to determine the convex combination x of the 
 *     gauss-newton and scaled gradient directions that minimizes 
 *     (A*x - b) in the least squares sense, subject to the 
 *     restriction that the euclidean norm of D*x be at most delta. 
 *
 *     This subroutine completes the solution of the problem 
 *     if it is provided with the necessary information from the 
 *     QR factorization of A. that is, if A = Q*R, where Q has 
 *     orthogonal columns and R is an upper triangular matrix, 
 *     then dogleg expects the full upper triangle of R and 
 *     the first n components of (Q transpose)*b. 
 *
 *     the subroutine statement is 
 *
 *       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2) 
 *
 *     where 
 *
 *       n is a positive integer input variable set to the order of r. 
 *
 *       r is an input array of length lr which must contain the upper 
 *         triangular matrix r stored by rows. 
 *
 *       lr is a positive integer input variable not less than 
 *         (n*(n+1))/2. 
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
 *       x is an output array of length n which contains the desired 
 *         convex combination of the gauss-newton direction and the 
 *         scaled gradient direction. 
 *
 *       wa1 and wa2 are work arrays of length n. 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void dogleg_(const int *n, double *r, const int *lr,
	     const double *diag, const double *qtb, const double *delta,
	     double *x, double *wa1, double *wa2)
{
	(void)lr;

	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;

	/* first, calculate the gauss-newton direction. */
	int jj = *n * (*n + 1) / 2;
	for (long int k = 1; k <= *n; ++k) {
		int j = *n - k;
		jj -= k;
		int l = jj + 1;
		double sum = 0;
		for (long int i = j+1; i < *n; ++i) {
			sum += r[l] * x[i];
			++l;
		}
		double temp = r[jj];
		if (temp == 0) {
			int l = j;
			for (long int i = 1; i < j; ++i) {
				temp = fmax(temp, fabs(r[l]));
				l += *n - i;
			}
			temp = epsmch * temp;
			if (temp == 0)
				temp = epsmch;
		}
		x[j] = (qtb[j] - sum) / temp;
	}

	/* test whether the gauss-newton direction is acceptable. */
	for (long int j = 0; j < *n; ++j) {
		wa1[j] = 0;
		wa2[j] = diag[j] * x[j];
	}
	double qnorm = enorm_(n, wa2);
	if (qnorm <= *delta)
		return;

	/* the gauss-newton direction is not acceptable. */
	/* next, calculate the scaled gradient direction. */
	int l = 0;
	for (long int j = 0; j < *n; ++j) {
		for (long int i = j; i < *n; ++i) {
			wa1[i] += r[l] * qtb[j];
			++l;
		}
		wa1[j] /= diag[j];
	}

	/* calculate the norm of the scaled gradient and test for */
	/* the special case in which the scaled gradient is zero. */
	double gnorm = enorm_(n, wa1);
	double sgnorm = 0;
	double alpha = *delta / qnorm;
	
	if (gnorm > 0) {
		/* calculate the point along the scaled gradient */
		/* at which the quadratic is minimized. */
		for (long int j = 0; j < *n; ++j)
			wa1[j] = wa1[j] / gnorm / diag[j];
		l = 0;
		for (long int j = 0; j < *n; ++j) {
			double sum = 0;
			for (long int i = j; i < *n; ++i) {
				sum += r[l] * wa1[i];
				++l;
			}
			wa2[j] = sum;
		}
		double temp = enorm_(n, wa2);
		sgnorm = gnorm / temp / temp;

		/* test whether the scaled gradient direction is acceptable. */
		alpha = 0;
		if (sgnorm < *delta) {
			/* the scaled gradient direction is not acceptable. */
			/* finally, calculate the point along the dogleg */
			/* at which the quadratic is minimized. */
			double bnorm = enorm_(n, qtb);
			double tmp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
			double d1 = sgnorm / *delta;
			double d2 = tmp - *delta / qnorm;
			double d3 = *delta / qnorm;
			double d4 = sgnorm / *delta;
			double tmp2 = tmp - *delta / qnorm * (d1 * d1) + sqrt(d2 * d2 + (1 - d3 * d3) * (1 - d4 * d4));
			double d5 = sgnorm / *delta;
			alpha = *delta / qnorm * (1 - d5 * d5) / tmp2;
		}
	}
	/* form appropriate convex combination of the gauss-newton */
	/* direction and the scaled gradient direction. */
	double tmp = (1 - alpha) * fmin(sgnorm, *delta);
	for (long int j = 0; j < *n; ++j)
		x[j] = tmp * wa1[j] + alpha * x[j];
}
