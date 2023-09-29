#include <math.h>

#include "minpack.h"
/*
 *     Subroutine r1mpyq 
 *
 *     Given an m by n matrix A, this subroutine computes A*Q where 
 *     Q is the product of 2*(n - 1) transformations 
 *
 *           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) 
 *
 *     and gv(i), gw(i) are givens rotations in the (i,n) plane which 
 *     eliminate elements in the i-th and n-th planes, respectively. 
 *     Q itself is not given, rather the information to recover the 
 *     gv, gw rotations is supplied. 
 *
 *     The subroutine statement is 
 *
 *       subroutine r1mpyq(m,n,a,lda,v,w) 
 *
 *     where 
 *
 *       m is a positive integer input variable set to the number 
 *         of rows of a. 
 *
 *       n is a positive integer input variable set to the number 
 *         of columns of a. 
 *
 *       a is an m by n array.  On input a must contain the matrix 
 *         to be postmultiplied by the orthogonal matrix q 
 *         described above.  On output A*Q has replaced A. 
 *
 *       lda is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array A. 
 *
 *       v is an input array of length n.  v(i) must contain the 
 *         information necessary to recover the givens rotation gv(i) 
 *         described above. 
 *
 *       w is an input array of length n.  w(i) must contain the 
 *         information necessary to recover the givens rotation gw(i) 
 *         described above. 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void r1mpyq_(const int *m, const int *n, double *a, const int *lda, double *v, double *w)
{
	int a_dim1 = *lda;
	int c1 = 1;

	/* apply the first set of givens rotations to a. */	
	for (long int nmj = 1; nmj <= *n - 1; ++nmj) {
		int j = *n - 1 - nmj;
		double cos, sin;
		if (fabs(v[j]) > 1) {
			cos = 1 / v[j];
			sin = sqrt(1 - cos * cos);
		} else {
			sin = v[j];
			cos = sqrt(1 - sin * sin);
		}
		drot_(m, &a[(*n-1) * a_dim1], &c1, &a[j * a_dim1], &c1, &cos, &sin);
	}

	/* apply the second set of givens rotations to a. */
	for (long int j = 0; j < *n - 1; ++j) {
		double cos, sin;
		if (fabs(w[j]) > 1) {
			cos = 1 / w[j];
			sin = sqrt(1 - cos * cos);
		} else {
			sin = w[j];
			cos = sqrt(1 - sin * sin);
		}
		drot_(m, &a[j * a_dim1], &c1, &a[(*n-1) * a_dim1], &c1, &cos, &sin);
	}
}