#include "minpack.h"

/*
 *     function enorm
 *
 *     given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     the euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. the sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. non-destructive underflows are permitted. underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     the definitions of small, intermediate and large components
 *     depend on two constants, rdwarf and rgiant. the main
 *     restrictions on these constants are that rdwarf**2 not
 *     underflow and rgiant**2 not overflow. the constants
 *     given here are suitable for every known computer.
 *
 *     the function statement is
 *
 *       double precision function enorm(n,x)
 *
 *     where
 *
 *       n is a positive integer input variable.
 *
 *       x is an input array of length n.
 *
 *     subprograms called
 * 
 *       fortran-supplied ... dabs,dsqrt
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

double enorm_(const int *n, double const *x)
{
	int incx = 1;
	return dnrm2_(n, x, &incx);
}