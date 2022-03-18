#include <stddef.h>

#include "minpack.h"

/*
 *     Subroutine lmder 
 *
 *     The purpose of lmder is to minimize the sum of the squares of 
 *     m nonlinear functions in n variables by a modification of 
 *     the Levenberg-Marquardt algorithm.  The user must provide a 
 *     subroutine which calculates the functions and the jacobian. 
 *
 *     The subroutine statement is 
 *
 *       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, 
 *                        maxfev,diag,mode,factor,nprint,info,nfev, 
 *                        njev,ipvt,qtf,wa1,wa2,wa3,wa4) 
 *
 *     where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions and the jacobian.  fcn must 
 *         be declared in an external statement in the user 
 *         calling program, and should be written as follows. 
 *
 *         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag) 
 *         integer m,n,ldfjac,iflag 
 *         double precision x(n),fvec(m),fjac(ldfjac,n) 
 *         ---------- 
 *         If iflag = 1 calculate the functions at x and 
 *         return this vector in fvec.  Do not alter fjac. 
 *         If iflag = 2 calculate the jacobian at x and 
 *         return this matrix in fjac.  Do not alter fvec. 
 *         ---------- 
 *         return 
 *         end 
 *
 *         The value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of lmder. 
 *         In this case set iflag to a negative integer. 
 *
 *       m is a positive integer input variable set to the number 
 *         of functions. 
 *
 *       n is a positive integer input variable set to the number 
 *         of variables. n must not exceed m. 
 *
 *       x is an array of length n.  On input x must contain 
 *         an initial estimate of the solution vector.  On output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length m which contains 
 *         the functions evaluated at the output x. 
 *
 *       fjac is an output m by n array.  The upper n by n submatrix 
 *         of fjac contains an upper triangular matrix R with 
 *         diagonal elements of nonincreasing magnitude such that 
 *
 *                t     t           t 
 *               P *(jac *jac)*P = R *R, 
 *
 *         where P is a permutation matrix and jac is the final 
 *         calculated jacobian.  Column j of P is column ipvt(j) 
 *         (see below) of the identity matrix.  The lower trapezoidal 
 *         part of fjac contains information generated during 
 *         the computation of R. 
 *
 *       ldfjac is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       ftol is a nonnegative input variable.  Termination 
 *         occurs when both the actual and predicted relative 
 *         reductions in the sum of squares are at most ftol. 
 *         Therefore, ftol measures the relative error desired 
 *         in the sum of squares. 
 *
 *       xtol is a nonnegative input variable. termination 
 *         occurs when the relative error between two consecutive 
 *         iterates is at most xtol.  Therefore, xtol measures the 
 *         relative error desired in the approximate solution. 
 *
 *       gtol is a nonnegative input variable.  Termination 
 *         occurs when the cosine of the angle between fvec and 
 *         any column of the jacobian is at most gtol in absolute 
 *         value.  Therefore, gtol measures the orthogonality 
 *         desired between the function vector and the columns 
 *         of the jacobian. 
 *
 *       maxfev is a positive integer input variable.  Termination 
 *         occurs when the number of calls to fcn with iflag = 1 
 *         has reached maxfev. 
 *
 *       diag is an array of length n.  If mode = 1 (see 
 *         below), diag is internally set.  If mode = 2, diag 
 *         must contain positive entries that serve as 
 *         multiplicative scale factors for the variables. 
 *
 *       mode is an integer input variable.  If mode = 1, the 
 *         variables will be scaled internally.  If mode = 2, 
 *         the scaling is specified by the input diag.  Other 
 *         values of mode are equivalent to mode = 1. 
 *
 *       factor is a positive input variable used in determining the 
 *         initial step bound. this bound is set to the product of 
 *         factor and the euclidean norm of diag*x if nonzero, or else 
 *         to factor itself. in most cases factor should lie in the 
 *         interval (.1,100.).100. is a generally recommended value. 
 *
 *       nprint is an integer input variable that enables controlled 
 *         printing of iterates if it is positive.  In this case, 
 *         fcn is called with iflag = 0 at the beginning of the first 
 *         iteration and every nprint iterations thereafter and 
 *         immediately prior to return, with x, fvec, and fjac 
 *         available for printing.  fvec and fjac should not be 
 *         altered.  If nprint is not positive, no special calls 
 *         of fcn with iflag = 0 are made. 
 *
 *       info is an integer output variable.  If the user has 
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
 *         info = 5  number of calls to fcn with iflag = 1 has 
 *                   reached maxfev. 
 *
 *         info = 6  ftol is too small.  No further reduction in 
 *                   the sum of squares is possible. 
 *
 *         info = 7  xtol is too small.  No further improvement in 
 *                   the approximate solution x is possible. 
 *
 *         info = 8  gtol is too small.  fvec is orthogonal to the 
 *                   columns of the jacobian to machine precision. 
 *
 *       nfev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 1. 
 *
 *       njev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 2. 
 *
 *       ipvt is an integer output array of length n.  ipvt 
 *         defines a permutation matrix p such that jac*P = Q*R, 
 *         where jac is the final calculated jacobian, Q is 
 *         orthogonal (not stored), and R is upper triangular 
 *         with diagonal elements of nonincreasing magnitude. 
 *         column j of P is column ipvt(j) of the identity matrix. 
 *
 *       qtf is an output array of length n which contains 
 *         the first n elements of the vector (Q transpose)*fvec. 
 *
 *       wa1, wa2, and wa3 are work arrays of length n. 
 *
 *       wa4 is a work array of length m. 
 *
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
*/
void lmder_c(minpack_func_mnj fcn, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *ftol,
	const double *xtol, const double *gtol, const int *maxfev, 
    double *diag, const int *mode, const double *factor, const int *nprint, 
    int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4)
{
    lmbase(fcn, NULL, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, NULL, 
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4);
}