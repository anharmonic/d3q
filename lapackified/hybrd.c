#include <stddef.h>

#include "minpack.h"

/*
 *     subroutine hybrd 
 *
 *     the purpose of hybrd is to find a zero of a system of 
 *     n nonlinear functions in n variables by a modification 
 *     of the powell hybrid method. the user must provide a 
 *     subroutine which calculates the functions. the jacobian is 
 *     then calculated by a forward-difference approximation. 
 *
 *     the subroutine statement is 
 *
 *       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn, 
 *                        diag,mode,factor,nprint,info,nfev,fjac, 
 *                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4) 
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
 *         --------- 
 *         return 
 *         end 
 *
 *         the value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of hybrd. 
 *         in this case set iflag to a negative integer. 
 *
 *       n is a positive integer input variable set to the number 
 *         of functions and variables. 
 *
 *       x is an array of length n. on input x must contain 
 *         an initial estimate of the solution vector. on output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length n which contains 
 *         the functions evaluated at the output x. 
 *
 *       xtol is a nonnegative input variable. termination 
 *         occurs when the relative error between two consecutive 
 *         iterates is at most xtol. 
 *
 *       maxfev is a positive integer input variable. termination 
 *         occurs when the number of calls to fcn is at least maxfev 
 *         by the end of an iteration. 
 *
 *       ml is a nonnegative integer input variable which specifies 
 *         the number of subdiagonals within the band of the 
 *         jacobian matrix. if the jacobian is not banded, set 
 *         ml to at least n - 1. 
 *
 *       mu is a nonnegative integer input variable which specifies 
 *         the number of superdiagonals within the band of the 
 *         jacobian matrix. if the jacobian is not banded, set 
 *         mu to at least n - 1. 
 *
 *       epsfcn is an input variable used in determining a suitable 
 *         step length for the forward-difference approximation. this 
 *         approximation assumes that the relative errors in the 
 *         functions are of the order of epsfcn. if epsfcn is less 
 *         than the machine precision, it is assumed that the relative 
 *         errors in the functions are of the order of the machine 
 *         precision. 
 *
 *       diag is an array of length n. if mode = 1 (see 
 *         below), diag is internally set. if mode = 2, diag 
 *         must contain positive entries that serve as 
 *         multiplicative scale factors for the variables. 
 *
 *       mode is an integer input variable. if mode = 1, the 
 *         variables will be scaled internally. if mode = 2, 
 *         the scaling is specified by the input diag. other 
 *         values of mode are equivalent to mode = 1. 
 *
 *       factor is a positive input variable used in determining the 
 *         initial step bound. this bound is set to the product of 
 *         factor and the euclidean norm of diag*x if nonzero, or else 
 *         to factor itself. in most cases factor should lie in the 
 *         interval (.1,100.). 100. is a generally recommended value. 
 *
 *       nprint is an integer input variable that enables controlled 
 *         printing of iterates if it is positive. in this case, 
 *         fcn is called with iflag = 0 at the beginning of the first 
 *         iteration and every nprint iterations thereafter and 
 *         immediately prior to return, with x and fvec available 
 *         for printing. if nprint is not positive, no special calls 
 *         of fcn with iflag = 0 are made. 
 *
 *       info is an integer output variable. if the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag. see description of fcn. otherwise, 
 *         info is set as follows. 
 *
 *         info = 0   improper input parameters. 
 *
 *         info = 1   relative error between two consecutive iterates 
 *                    is at most xtol. 
 *
 *         info = 2   number of calls to fcn has reached or exceeded 
 *                    maxfev. 
 *
 *         info = 3   xtol is too small. no further improvement in 
 *                    the approximate solution x is possible. 
 *
 *         info = 4   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    five jacobian evaluations. 
 *
 *         info = 5   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    ten iterations. 
 *
 *       nfev is an integer output variable set to the number of 
 *         calls to fcn. 
 *
 *       fjac is an output n by n array which contains the 
 *         orthogonal matrix q produced by the qr factorization 
 *         of the final approximate jacobian. 
 *
 *       ldfjac is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       r is an output array of length lr which contains the 
 *         upper triangular matrix produced by the qr factorization 
 *         of the final approximate jacobian, stored rowwise. 
 *
 *       lr is a positive integer input variable not less than 
 *         (n*(n+1))/2. 
 *
 *       qtf is an output array of length n which contains 
 *         the vector (q transpose)*fvec. 
 *
 *       wa1, wa2, wa3, and wa4 are work arrays of length n. 
 *
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void hybrd_c(minpack_func_n fcn, 
		  const int *n, double *x, double *fvec, const double *xtol, const int *maxfev,
		  const int *ml, const int *mu, const double *epsfcn, double *diag, const int *mode,
		  const double *factor, const int *nprint, int *info, int *nfev,
		  double *fjac, const int *ldfjac, double *r, const int *lr, double *qtf,
		  double *wa1, double *wa2, double *wa3, double *wa4)
{
	hybrbase(fcn, NULL, n, x, fvec, fjac, ldfjac, xtol, maxfev, ml, mu, epsfcn, diag, mode,
	      factor, nprint, info, nfev, NULL, r, lr, qtf, wa1, wa2, wa3, wa4);
}