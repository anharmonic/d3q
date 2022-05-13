#include <float.h>

/*
 * User-supplied subroutine which calculates the functions.
 *
 * void fcn(void *farg, int m, int n, const double * x, double fvec)
 * {
 *         // calculate the functions at x and return this vector in fvec.
 * }
 *
 * Evaluating the functions is a collective operation.  All processes
 * participating in the minimization call fcn with the same arguments (except 
 * farg). All processes must receive the results in fvec.
 *
 * Argument list:
 *
 * - farg is an opaque user-provided pointer that can be used to provide 
 *   arguments to fcn.  Distinct values on all processes.
 *
 * - m is the number of functions to evaluate.  Same value on all processes.
 *
 * - n is the number of variables. Same value on all processes.
 *
 * - x is the vector of input variables.  Same value on all processes.
 *
 * - fvec is the output vector that receives the function evaluations.
 *   Must receive the same value on all processes.
 */
typedef void (*pminpack_func_mn)(void *arg, int m, int n, const double *x, double *fvec);


int plmdif1(pminpack_func_mn fcn, void *farg, int m, int n, double *x, 
            double *fvec, double tol, int ctx);

int plmdif(pminpack_func_mn fcn, void *farg, int m, int n, double *x, double *fvec, 
        double ftol, double xtol, double gtol, int maxfev, int *nfev, int ctx);

/* This replaces dpmpar */
#define MINPACK_EPSILON DBL_EPSILON
#define MINPACK_DWARF   DBL_MIN
#define MINPACK_GIANT   DBL_MAX

/* C equivalent of MINPACK routines */
double enorm(int n, const double * x);
double lmpar(int n, double *r, int ldr, int *ipvt, double *diag, double par,
        double *qtb, double delta, double *x, double *sdiag, double *wa1, double *wa2, int talk);
void qrsolv(int n, double *r, int ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa);

/* utils */
double pminpack_wtime();
void pminpack_human_format(char * target, long n);

/* ScaLAPACK wrappers */
int scalapack_numroc(int n, int nb, int rank, int srcrank, int nprocs);
int scalapack_descinit(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictx, int lld);
int scalapack_pdgeqpf(int M, int N, double * A, int IA, int JA, int * DESCA, int * IPIV, double * TAU, double * WORK, int LWORK);
int scalapack_pdormqr(char * SIDE, char * TRANS, int M, int N, int K, double * A, int IA, int JA, int * DESCA, double *TAU,
                      double * C, int IC, int JC, int * DESCC, double * WORK, int LWORK);

/* Extra "BLACS-like" routines */
void extrablacs_dgeld2d(const double *v, int vlen, double *A, int jA, const int *descA);
void extrablacs_dgedl2d(int m, int n, const double *A, int iA, int jA, const int *descA, double *B, int ldB);
void extrablacs_igedl2d(int m, int n, const int *A, int iA, int jA, const int *descA, int *B, int ldB);
void extrablacs_dtrdl2d(const char *uplo, const char *diag, int m, int n, 
                        const double *A, int iA, int jA, const int *descA, double *B, int ldB);

/* ScaLAPACK functions */
#define DTYPE_ 0
#define CTXT_  1 
#define M_     2
#define N_     3 
#define MB_    4 
#define NB_    5
#define RSRC_  6 
#define CSRC_  7 
#define LLD_   8

extern int numroc_(const int *n, const int *nb, const int *rank, const int *srcrank, const int *nprocs);
extern void descinit_(int * desc, const int * m, const int * n, const int * mb, const int * nb, 
                      const int * irsrc, const int * icsrc, const int * ictx, const int * lld, int * info);
extern void pdgeqpf_(int * M, int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, double * TAU, 
                     double * WORK, int * LWORK, int *INFO);
extern void pdormqr_(const char * SIDE, const char * TRANS, const int * M, const int *N, const int *K, const double * A, 
                     int * IA, int * JA, int * DESCA, double * TAU, double * C, int * IC, int * JC, int * DESCC, 
                     double * WORK, int * LWORK, int * INFO);
extern void Cpdgemr2d(int m, int n, const double *A, int IA, int JA, const int *descA, 
                      double *B, int IB, int JB, const int *descB, int gcontext);
extern void Cpigemr2d(int m, int n, const int *A, int IA, int JA, const int *descA, 
                      int *B, int IB, int JB, const int *descB, int gcontext);
extern void Cpdtrmr2d(const char *uplo, const char *diag, int m, int n, const double *A, int IA, int JA, const int *descA, 
                      double *B, int IB, int JB, const int *descB, int gcontext);

/* BLACS routines */
extern void Cblacs_get(int icontxt, int what, int *val);
extern void Cblacs_pinfo(int * rank, int * nprocs);
extern void Cblacs_gridinfo(int ictx, int * nprow, int * npcol, int * myrow, int * mycol);
extern void Cblacs_gridinit(int * ictx, const char * order, int nprow, int npcol);
extern void Cblacs_gridexit(int ictx);
extern void Cdgebs2d(int ctx, const char *scope, const char *top, int m, int n, const double *A, int ldA);
extern void Cdgebr2d(int ctx, const char *scope, const char *top, int m, int n, double *A, int ldA, int rsrc, int csrc);
extern void Cigebs2d(int ctx, const char *scope, const char *top, int m, int n, const int *A, int ldA);
extern void Cigebr2d(int ctx, const char *scope, const char *top, int m, int n, int *A, int ldA, int rsrc, int csrc);
extern void Cdtrbs2d(int ctx, const char *scope, const char *top, const char *uplo, const char *diag, int m, int n, const double *A, int ldA);
extern void Cdtrbr2d(int ctx, const char *scope, const char *top, const char *uplo, const char *diag, int m, int n, double *A, int ldA, int rsrc, int csrc);