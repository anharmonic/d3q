#include <float.h>

/************************* user-supplied functions *********************************/

/* 
 * for hybrd1 and hybrd:
 *    calculate the functions at x and return this vector in fvec.
 *   return a negative value to terminate hybrd1/hybrd.
 */
typedef void (*minpack_func_n)(const int * n, const double * x, double *fvec, int *iflag);

/*
 * for hybrj1 and hybrj
 *         if iflag = 1 calculate the functions at x and return this vector in fvec. do not alter fjac.
 *         if iflag = 2 calculate the jacobian  at x and return this matrix in fjac. do not alter fvec.
 * return a negative value to terminate hybrj1/hybrj
 */
typedef void (*minpack_func_nj)(const int *n, const double *x, double *fvec, double *fjac, const int *ldfjac, int *iflag);

/*
 * for lmdif1 and lmdif
 *         calculate the functions at x and return this vector in fvec.
 *         if iflag = 1 the result is used to compute the residuals.
 *         if iflag = 2 the result is used to compute the Jacobian by finite differences.
 *         One jacobian computation requires exactly n function calls with iflag = 2.
 * return a negative value to terminate lmdif1/lmdif
 */
typedef void (*minpack_func_mn)(const int * m, const int * n, double const *x, double *fvec, int *iflag);

/*
 * for lmder1 and lmder
 *         if iflag = 1 calculate the functions at x and return this vector in fvec. do not alter fjac.
 *         if iflag = 2 calculate the jacobian  at x and return this matrix in fjac. do not alter fvec.
 * return a negative value to terminate lmder1/lmder
 */
typedef void (*minpack_func_mnj)(const int *m, const int *n, const double *x, double *fvec, double *fjac, const int *ldfjac, int *iflag);


/************************** MINPACK routines **************************/

/* find a zero of a system of n nonlinear functions in n variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation) */
void hybrd1_c(minpack_func_n fcn, const int *n, double *x, double *fvec, 
    const double *tol, int *info, double *wa, const int *lwa);

/* find a zero of a system of n nonlinear functions in n variables by
   a modification of the Powell hybrid method (user-supplied Jacobian) */
void hybrj1_c(minpack_func_nj fcn, const int *n, double *x, double *fvec, 
    double *fjac, const int *ldfjac, const double *tol, int *info, 
    double *wa, const int *lwa);
     
/* minimize the sum of the squares of m nonlinear functions in n
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation) */
void lmdif1_c(minpack_func_mn fcn, int const *m, int const *n, double *x, 
	double *fvec, double const *tol, int *info, int *iwa, 
	double *wa, int const *lwa);

/* minimize the sum of the squares of m nonlinear functions in n
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian) */
void lmder1_(minpack_func_mnj fcn, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *tol, 
	int *info, int *ipvt, double *wa, const int *lwa);


/* find a zero of a system of n nonlinear functions in n variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */
void hybrd_c(minpack_func_n fcn, 
	      const int *n, double *x, double *fvec, const double *xtol, const int *maxfev,
	      const int *ml, const int *mu, const double *epsfcn, double *diag, const int *mode,
	      const double *factor, const int *nprint, int *info, int *nfev,
	      double *fjac, const int *ldfjac, double *r, const int *lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);
       
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */
void hybrj_c(minpack_func_nj fcn, const int *n, double *x, double *fvec, 
    double *fjac, const int *ldfjac, const double *xtol, const int *maxfev, 
    double *diag, const int *mode, const double *factor, const int *nprint, 
    int *info, int *nfev, int *njev, double *r, const int *lr, double *qtf, 
    double *wa1, double *wa2, double *wa3, double *wa4);

/* minimize the sum of the squares of m nonlinear functions in n
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */
void lmdif_c(minpack_func_mn fcn, const int * m, const int * n, double * x,
	   double * fvec, const double * ftol, const double * xtol, 
	   const double *gtol, const int * maxfev, const double * epsfcn, 
	   double * diag, const int * mode, const double * factor, 
	   const int * nprint, int * info, int * nfev, double * fjac, 
	   const int * ldfjac, int * ipvt, double * qtf, double * wa1, 
	   double * wa2, double * wa3, double * wa4);

/* minimize the sum of the squares of m nonlinear functions in n
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */
void lmder_c(minpack_func_mnj fcn, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *ftol,
	const double *xtol, const double *gtol, const int *maxfev, 
   double *diag, const int *mode, const double *factor, const int *nprint, 
   int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4);

void chkder_c(const int *m, const int *n, double *x,
	     double *fvec, double *fjac, const int *ldfjac, double *xp,
	     double *fvecp, const int *mode, double *err);


/***************************** internal MINPACK routines *****************************/


/* This replaces dpmpar */
#define MINPACK_EPSILON DBL_EPSILON
#define MINPACK_DWARF   DBL_MIN
#define MINPACK_GIANT   DBL_MAX

/* Original values in dpmpar */
// #define MINPACK_EPSILON  2.22044604926e-16 
// #define MINPACK_DWARF    2.22507385852e-308
// #define MINPACK_GIANT    1.79769313485e308

double enorm_(const int *n, double const *x);

void fdjac1_(minpack_func_n fcn, const int *n, double *x, double *fvec,
	    double *fjac, const int *ldfjac, int *iflag, const int *ml,
	    const int *mu, const double *epsfcn, double *wa1, double *wa2);

void fdjac2_c(minpack_func_mn fcn, const int *m, const int *n, double *x,
	    double const *fvec, double *fjac, const int *ldfjac, int *iflag,
	    const double *epsfcn, double *wa);

void lmpar_(const int *n, double *r, const int *ldr, int *ipvt, double *diag, 
	double *qtb, double *delta, double *par, double *x, double *sdiag, double *wa1, double *wa2);

void dogleg_(const int *n, double *r__, const int *lr, 
	const double *diag, const double *qtb, const double *delta, 
    double *x, double *wa1, double *wa2);

void qrsolv_(const int *n, double *r, const int *ldr, int *ipvt, 
	double *diag, double *qtb, double *x, double *sdiag, double *wa);

void r1updt_(const int *m, const int *n, double *s, const int *ls, 
	const double *u, double *v, double *w, int *sing);

void r1mpyq_(const int *m, const int *n, double *a, const int *lda, double *v, double *w);

void lmbase(minpack_func_mnj fcn_der, minpack_func_mn fcn_dif, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *ftol,
	const double *xtol, const double *gtol, const int *maxfev, const double * epsfcn, 
    double *diag, const int *mode, const double *factor, const int *nprint, 
    int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4);

void hybrbase(minpack_func_n fcn_dif, minpack_func_nj fcn_der,
	      const int *n, double *x, double *fvec, double *fjac, const int *ldfjac, 
	      const double *xtol, const int *maxfev,
	      const int *ml, const int *mu, const double *epsfcn, double *diag, const int *mode,
	      const double *factor, const int *nprint, int *info, int *nfev, int *njev,
	      double *r, const int *lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);


/**************************** BLAS ************************/

extern double dnrm2_(const int * n, double const * x, int * incx);

extern void drotg_(double *a, double *b, double *c, double *d);

extern void drot_(const int *n, double *dx, const int *incx, double *dy, 
	const int *incy, const double *c, const double *s);


extern void dtrsv_(const char *UPLO, const char *TRANS, const char *DIAG, const int *N, 
	const double *A, const int *LDA, double *X, const int *INCX);

extern void dtrmv_(const char *UPlo, const char * TRANS, const char *diag,
	const int *N, const double *A, const int *lda, const double *X, const int *incx);

extern void dgemv_(const char * TRANS, const int *M, const int *N, const double *alpha, 
	const double *A, const int *lda, const double *X, const int *incx, 
	const double *beta, double *Y, const int *incy);

/**************************** LAPACK ************************/

extern int dgeqp3_(const int *m, const int *n, double *a, const int *lda, 
	int *jpvt, double *tau, double *work, const int *lwork, int *info);

extern int dgeqrf_(const int *m, const int *n, double *a, const int *lda, 
	double *tau, double *work, const int *lwork, int *info);

extern int dormqr_(const char *side, const char *trans, 
	const int *m, const int *n, const int *k, 
	const double *a, const int *lda, const double *tau,
    double *c, const int *ldc, double *work, const int *lwork, int *info);

extern int dorgqr_(const int *m, const int *n, const int *k, 
	const double *a, const int *lda, const double *tau,
    double *work, const int *lwork, int *info);

