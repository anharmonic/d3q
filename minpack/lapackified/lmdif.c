#include <stddef.h>

#include "minpack.h"

/*     subroutine lmdif
 *
 *     the purpose of lmdif is to minimize the sum of the squares of
 *     m nonlinear functions in n variables by a modification of
 *     the levenberg-marquardt algorithm. the user must provide a
 *     subroutine which calculates the functions. the jacobian is
 *     then calculated by a forward-difference approximation.
 *
 *     the subroutine statement is
 *
 *       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
 *                        diag,mode,factor,nprint,info,nfev,fjac,
 *                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
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
 *         the user wants to terminate execution of lmdif.
 *         in this case set iflag to a negative int.
 *
 *       m is a positive int input variable set to the number
 *         of functions.
 *
 *       n is a positive int input variable set to the number
 *         of variables. n must not exceed m.
 *
 *       x is an array of length n. on input x must contain
 *         an initial estimate of the solution vector. on output x
 *         contains the final estimate of the solution vector.
 *
 *       fvec is an output array of length m which contains
 *         the functions evaluated at the output x.
 *
 *       ftol is a nonnegative input variable. termination
 *         occurs when both the actual and predicted relative
 *         reductions in the sum of squares are at most ftol.
 *         therefore, ftol measures the relative error desired
 *         in the sum of squares.
 *
 *       xtol is a nonnegative input variable. termination
 *         occurs when the relative error between two consecutive
 *         iterates is at most xtol. therefore, xtol measures the
 *         relative error desired in the approximate solution.
 *
 *       gtol is a nonnegative input variable. termination
 *         occurs when the cosine of the angle between fvec and
 *         any column of the jacobian is at most gtol in absolute
 *         value. therefore, gtol measures the orthogonality
 *         desired between the function vector and the columns
 *         of the jacobian.
 *
 *       maxfev is a positive int input variable. termination
 *         occurs when the number of calls to fcn is at least
 *         maxfev by the end of an iteration.
 *
 *       epsfcn is an input variable used in determining a suitable
 *         step length for the forward-difference approximation.  This
 *         approximation assumes that the relative errors in the
 *         functions are of the order of epsfcn.  If epsfcn is less
 *         than the machine precision, it is assumed that the relative
 *         errors in the functions are of the order of the machine
 *         precision.
 *
 *       diag is an array of length n. if mode = 1 (see
 *         below), diag is internally set. if mode = 2, diag
 *         must contain positive entries that serve as
 *         multiplicative scale factors for the variables.
 *
 *       mode is an int input variable. if mode = 1, the
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
 *       nprint is an int input variable that enables controlled
 *         printing of iterates if it is positive. in this case,
 *         fcn is called with iflag = 0 at the beginning of the first
 *         iteration and every nprint iterations thereafter and
 *         immediately prior to return, with x and fvec available
 *         for printing. if nprint is not positive, no special calls
 *         of fcn with iflag = 0 are made.
 *
 *       info is an int output variable. if the user has
 *         terminated execution, info is set to the (negative)
 *         value of iflag. see description of fcn. otherwise,
 *         info is set as follows.
 *
 *         info = 0  improper input parameters.
 *
 *         info = 1  both actual and predicted relative reductions
 *                   in the sum of squares are at most ftol.
 *
 *         info = 2  relative error between two consecutive iterates
 *                   is at most xtol.
 *
 *         info = 3  conditions for info = 1 and info = 2 both hold.
 *
 *         info = 4  the cosine of the angle between fvec and any
 *                   column of the jacobian is at most gtol in
 *                   absolute value.
 *
 *         info = 5  number of calls to fcn has reached or
 *                   exceeded maxfev.
 *
 *         info = 6  ftol is too small. no further reduction in
 *                   the sum of squares is possible.
 *
 *         info = 7  xtol is too small. no further improvement in
 *                   the approximate solution x is possible.
 *
 *         info = 8  gtol is too small. fvec is orthogonal to the
 *                   columns of the jacobian to machine precision.
 *
 *       nfev is an int output variable set to the number of
 *         calls to fcn.
 *
 *       fjac is an output m by n array. The upper n by n submatrix
 *         of fjac contains an upper triangular matrix R with
 *         diagonal elements of nonincreasing magnitude such that
 *
 *                t     t           t
 *               P *(jac *jac)*P = R *R,
 *
 *         where P is a permutation matrix and jac is the final
 *         calculated jacobian. Column j of P is column ipvt(j)
 *         (see below) of the identity matrix. The lower trapezoidal
 *         part of fjac contains information generated during
 *         the computation of R.
 *
 *       ldfjac is a positive int input variable not less than m
 *         which specifies the leading dimension of the array fjac.
 *
 *       ipvt is an int output array of length n. ipvt
 *         defines a permutation matrix p such that jac*p = q*r,
 *         where jac is the final calculated jacobian, q is
 *         orthogonal (not stored), and r is upper triangular
 *         with diagonal elements of nonincreasing magnitude.
 *         column j of p is column ipvt(j) of the identity matrix.
 *
 *       qtf is an output array of length n which contains
 *         the first n elements of the vector (q transpose)*fvec.
 *
 *       wa1, wa2, and wa3 are work arrays of length n.
 *
 *       wa4 is a work array of length m.
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

void lmdif_c(minpack_func_mn fcn, const int * m, const int * n, double * x,
	   double * fvec, const double * ftol, const double * xtol, 
	   const double *gtol, const int * maxfev, const double * epsfcn, double * diag,
	   const int * mode, const double * factor, const int * nprint,
	   int * info, int * nfev, double * fjac, const int * ldfjac, 
	   int * ipvt, double * qtf, double * wa1, double * wa2, double * wa3, double * wa4)
{	
	lmbase(NULL, fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, epsfcn, 
    diag, mode, factor, nprint, info, nfev, NULL, ipvt, qtf, wa1, wa2, wa3, wa4);
}
