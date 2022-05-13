#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>

#include "pminpack.h"
#include "eq.h"

int rank, nprocs;	

void fcn(void *farg, int m, int n, const double * x, double *fvec)
{
	int * nprob = farg;
	(void) m;
	vecfcn(n, x, fvec, *nprob);
}

void do_test(int ctx, int nprob, int n, double factor)
{
	double tol = sqrt(MINPACK_EPSILON);
	printf("# %s with n=%d\n", problem_name[nprob - 1], n);

	double x[n];		// solution
	double * fvec = malloc(n * sizeof(*fvec));	    // residuals
	initpt(n, x, nprob, factor);	// set initial point

	double start = MPI_Wtime();
	int info = plmdif1(fcn, &nprob, n, n, x, fvec, tol, ctx);	// find solution
	double stop = MPI_Wtime();

	vecfcn(n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm(n, fvec);
	
	if (rank == 0) {
		printf("# LMDIF1: %.1fs\n", stop - start);
		printf("# Final norm of residual: %15.7e\n", fnorm2);
		printf("# info: %d\n", info);
		printf("\n");
	}
}


int main()
{
	/* MPI setup */
	MPI_Init(NULL, NULL);

	/* BLACS setup --- obtain system context */
	int nprow = 2;
	int npcol = 2;
	int sysctx;
	Cblacs_get(0, 0, &sysctx);

	/* obtain BLACS grid context */
	int ctx = sysctx;
	Cblacs_pinfo(&rank, &nprocs);
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);

	double start = MPI_Wtime();
	
	do_test(ctx, 13, 3000, 1);
	do_test(ctx, 11, 3200, 1);


	if (rank == 0)
		printf("# total time: %.1fs\n", MPI_Wtime() - start);

	/* MPI cleanup */
	MPI_Finalize();
}

