#include <stdio.h>
#include <sys/time.h>

#include "pminpack.h"

/******************* utility functions ********************/

double pminpack_wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* represent n in <= 6 char  */
void pminpack_human_format(char * target, long n)
{
        if (n < 1000) {
                sprintf(target, "%ld", n);
                return;
        }
        if (n < 1000000) {
                sprintf(target, "%.1fK", n / 1e3);
                return;
        }
        if (n < 1000000000) {
                sprintf(target, "%.1fM", n / 1e6);
                return;
        }
        if (n < 1000000000000ll) {
                sprintf(target, "%.1fG", n / 1e9);
                return;
        }
        if (n < 1000000000000000ll) {
                sprintf(target, "%.1fT", n / 1e12);
                return;
        }
}

/****** ScaLAPACK wrappers */

int scalapack_numroc(int n, int nb, int rank, int srcrank, int nprocs)
{
        return numroc_(&n, &nb, &rank, &srcrank, &nprocs);
}

int scalapack_descinit(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictx, int lld)
{
        int info = 0;
        descinit_(desc, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictx, &lld, &info);
        return info;
}

int scalapack_pdgeqpf(int M, int N, double * A, int IA, int JA, int * DESCA, int * IPIV, double * TAU, double * WORK, int LWORK)
{
        int info = 0;
        pdgeqpf_(&M, &N, A, &IA, &JA, DESCA, IPIV, TAU, WORK, &LWORK, &info);       
        return info;
}

int scalapack_pdormqr(char * SIDE, char * TRANS, int M, int N, int K, double * A, int IA, int JA, int * DESCA, double *TAU,
                      double * C, int IC, int JC, int * DESCC, double * WORK, int LWORK)
{
        int info = 0;
        pdormqr_(SIDE, TRANS, &M, &N, &K, A, &IA, &JA, DESCA, TAU, C, &IC, &JC, DESCC, WORK, &LWORK, &info);       
        return info;
}