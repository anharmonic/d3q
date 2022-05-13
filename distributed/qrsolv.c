#include <math.h>
#include <cblas.h>
#include <stddef.h>

#include "pminpack.h"

/*
 *     subroutine qrsolv 
 *
 *     given an m by n matrix A, an n by n diagonal matrix D, 
 *     and an m-vector b, the problem is to determine an x which 
 *     solves the system 
 *
 *           A*x = b ,     D*x = 0, 
 *
 *     in the least squares sense. 
 *
 *     this subroutine completes the solution of the problem 
 *     if it is provided with the necessary information from the 
 *     QR factorization, with column pivoting, of A.  That is, if 
 *     A*P = Q*R, where P is a permutation matrix, Q has orthogonal 
 *     columns, and R is an upper triangular matrix with diagonal 
 *     elements of nonincreasing magnitude, then qrsolv expects 
 *     the full upper triangle of R, the permutation matrix P, 
 *     and the first n components of (Q transpose)*b.  The system 
 *     A*x = b, D*x = 0, is then equivalent to 
 *
 *                  t       t 
 *           R*z = Q *b ,  P *D*P*z = 0 , 
 *
 *     where x = P*z.  If this system does not have full rank, 
 *     then a least squares solution is obtained.  On output qrsolv 
 *     also provides an upper triangular matrix  such that 
 *
 *            t   t               t 
 *           P *(A *A + D*D)*P = S *S. 
 *
 *     S is computed within qrsolv and may be of separate interest. 
 *
 *     the subroutine statement is 
 *
 *       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa) 
 *
 *     where 
 *
 *       n is a positive integer input variable set to the order of R. 
 *
 *       R is an n by n array.  On input the full upper triangle 
 *         must contain the full upper triangle of the matrix R. 
 *         On output the full upper triangle is unaltered, and the 
 *         strict lower triangle contains the strict upper triangle 
 *         (transposed) of the upper triangular matrix S. 
 *
 *       ldr is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array R. 
 *
 *       ipvt is an integer input array of length n which defines the 
 *         permutation matrix P such that A*P = Q*T.  Column j of P 
 *         is column ipvt(j) of the identity matrix. 
 *
 *       diag is an input array of length n which must contain the 
 *         diagonal elements of the matrix D. 
 *
 *       qtb is an input array of length n which must contain the first 
 *         n elements of the vector (Q transpose)*b. 
 *
 *       x is an output array of length n which contains the least 
 *         squares solution of the system A*x = b, D*x = 0. 
 *
 *       sdiag is an output array of length n which contains the 
 *         diagonal elements of the upper triangular matrix s. 
 *
 *       wa is a work array of length n. 
 *
 *     subprograms called 
 *
 *       fortran-supplied ... dabs,dsqrt 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */

void qrsolv(int n, double *r, int ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa)
{
	ptrdiff_t r_dim1 = ldr;

	/* copy R and (Q transpose)*b to preserve input and initialize S. */
	/* in particular, save the diagonal elements of R in x. */
	for (int j = 0; j < n; ++j) {
 		for (int i = j; i < n; ++i)
 			r[i + j * r_dim1] = r[j + i * r_dim1];
		x[j] = r[j + j * r_dim1];
		wa[j] = qtb[j];
	}

	/* eliminate the diagonal matrix D using a givens rotation. */
	for (int j = 0; j < n; ++j) {

		/* prepare the row of D to be eliminated, locating the */
		/* diagonal element using P from the QR factorization. */
		int l = ipvt[j] - 1;   // LAPACK is one-based
		if (diag[l] == 0)
			continue;

		for (int k = j; k < n; ++k)
			sdiag[k] = 0;
		sdiag[j] = diag[l];

		/* the transformations to eliminate the row of d */
		/* modify only a single element of (q transpose)*b */
		/* beyond the first n, which is initially zero. */
		double qtbpj = 0;

		for (int k = j; k < n; ++k) {
			/* determine a givens rotation which eliminates the */
			/* appropriate element in the current row of d. */
			if (sdiag[k] == 0)
				continue;
			
			double cos, sin;
			double a = r[k + k * r_dim1];
			double b = sdiag[k];
			cblas_drotg(&a, &b, &cos, &sin);

			/* compute the modified diagonal element of r and */
			/* the modified element of ((q transpose)*b,0). */
			r[k + k * r_dim1] = cos * r[k + k * r_dim1] + sin * sdiag[k];
			double temp = cos * wa[k] + sin * qtbpj;
			qtbpj = -sin * wa[k] + cos * qtbpj;
			wa[k] = temp;

			/* accumulate the transformation in the row of s. */
			int len = n - (k+1);
			cblas_drot(len, &r[k+1 + k * r_dim1], 1, &sdiag[k+1], 1, cos, sin);
		}
	}

	/* solve the triangular system for z. if the system is
	   singular, then obtain a least squares solution. */
	int nsing = n;
	for (int j = 0; j < n; ++j) {
		if (r[j + j * r_dim1] == 0) {
			nsing = j;
			break;
		}
	}
	for (int j = nsing + 1; j < n; ++j)
		wa[j] = 0;
	cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, nsing, r, ldr, wa, 1);

	for (int j = 0; j < n; ++j) {
		/* store the diagonal element of s and restore
		   the corresponding diagonal element of R. */
		sdiag[j] = r[j + j * r_dim1];
		r[j + j * r_dim1] = x[j];
	}

	/* permute the components of z back to components of x. */
	for (int j = 0; j < n; ++j) {
		int l = ipvt[j] - 1;
		x[l] = wa[j];
	}
}
