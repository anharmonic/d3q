/*
 * This program tests codes for the least-squares solution of
 * m nonlinear equations in n variables.
 *
 * Argonne National Laboratory. MINPACK project. march 1980.
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge j. Moré 
 */
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>

#include "pminpack.h"
#include "ls.h"

int lmdif1_known_failures[] = { 26, 27, 35, 38, 40, -1 };

/* global variables */
int nprob;
int nfev;
int njev;
int ma[60];
int na[60];
int nf[60];
int nj[60];
int np[60];
int nx[60];
double fnm[60];

/* This function is called by the solver */
void fcn(void * farg, int m, int n, const double *x, double *fvec)
{
	(void) farg;
	ssqfcn(m, n, x, fvec, nprob);
	nfev += 1;
}

void do_test(int ic, int rank, int ictx)
{
	double tol = sqrt(MINPACK_EPSILON);
	double ftol = 2e-5;

	// set global variables
	nprob = tests[ic].nprob;
	int n = tests[ic].n;
	int m = tests[ic].m;
	double factor = tests[ic].factor;

	double x[40];		// solution
	double fvec[100];	// residuals
	
	initpt(n, x, nprob, factor);	// set initial point
	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals

	nfev = 0;
	njev = 0;
	int info = plmdif1(fcn, NULL, m, n, x, fvec, tol, ictx);	// find solution

	if (rank > 0)
		return;

	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm(m, fvec);
	njev /= n;

	np[ic] = nprob;
	na[ic] = n;
	ma[ic] = m;
	nf[ic] = nfev;
	nj[ic] = njev;
	nx[ic] = info;
	fnm[ic] = fnorm2;

	// determine status and print it
	bool ok = (fnorm2 - tests[ic].fnorm2 < ftol);
	bool ok2 = (tests[ic].fnorm2_lastchance >= 0) && (fnorm2 - tests[ic].fnorm2_lastchance < ftol);

	if (ok || ok2) {
		printf("ok %d - %s (n = %d, m = %d, factor = %f)\n", ic + 1, problem_name[nprob - 1], n, m, factor);
	} else {
		bool TODO = false;
		for (int i = 0; lmdif1_known_failures[i] >= 0; i++)
			if (ic + 1 == lmdif1_known_failures[i])
				TODO = true;
		if (TODO)
			printf("not ok %d - %s (n = %d, m = %d, factor = %f) # TODO lmdif1 known failure\n", ic + 1, problem_name[nprob - 1], n, m, factor);
		else
			printf("not ok %d - %s (n = %d, m = %d, factor = %f)\n", ic + 1, problem_name[nprob - 1], n, m, factor);
	}
	
	commentator(ic, x, fvec, ftol, ftol, nfev, njev, info);
	
	for (int i = 0; i < n; i++)
		printf("\t# x[%2d] = %.20g\n", i, x[i]);
	for (int i = 0; i < m; i++)
		printf("\t# fvec[%2d] = %.20g\n", i, fvec[i]);
}


/* Run the usual collection of tests */
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
	int rank, nprocs;	
	Cblacs_pinfo(&rank, &nprocs);
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);

	/* TAP protocol */
	if (rank == 0)
		printf("1..53\n");

	for (int i = 0; i < 53; i++) {
		do_test(i, rank, ctx);
	}

	if (rank == 0) {
		printf("\n\n################################\n");
		printf("# Summary of 53 calls to lmdif1: \n\n");
		printf("#  test  nprob   n    m   nfev  njev  info  final l2 norm \n");
		for (int i = 0; i < 53; i++)
			printf("# %5d%5d%5d%5d%6d%6d%6d%16.7e\n", i + 1, np[i], na[i], ma[i], nf[i], nj[i], nx[i], fnm[i]);
	}

	/* MPI cleanup */
	MPI_Finalize();
	return EXIT_SUCCESS;
}
