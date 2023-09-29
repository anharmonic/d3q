#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>

#include "pminpack.h"

int main()
{
	/* MPI setup */
	MPI_Init(NULL, NULL);
	
	/* get BLACS system context */
	int sysctx;
	Cblacs_get(0, 0, &sysctx);

	/* BLACS setup */
	int rank, nprocs;	
	Cblacs_pinfo(&rank, &nprocs);
	
	int ctx = sysctx;
	int nprow = 2;
	int npcol = 2;
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);
	int myrow, mycol;
	Cblacs_gridinfo(ctx, &nprow, &npcol, &myrow, &mycol);

	/* setup distributed matrix */
	int n = 1000;
	int m = 1000;
	int nb = 31;
	int mb = 31;

	int A_nrow = scalapack_numroc(m, mb, myrow, 0, nprow);
	int A_ncol = scalapack_numroc(n, nb, mycol, 0, npcol);
	if (A_nrow == 0)
		A_nrow = 1;
	int ldA = A_nrow;
	printf("# process (%d, %d) owns local matrix of size %d x %d\n", 
		myrow, mycol, A_nrow, A_ncol);
	
	double *A = malloc(A_nrow * A_ncol * sizeof(*A));
	if (A == NULL) {
		printf("bail out! can't alloc A\n");
		exit(0);
	}
	
	int descA[9];
	int info = scalapack_descinit(descA, m, n, mb, nb, 0, 0, ctx, A_nrow);
	if (info < 0) {
		printf("bail out! descinit error\n");
		exit(0);
	}

	/* fill distributed matrix A with A[i,j] = (i-j) / (j + 1) */	
	for (int i = 0; i < m; i++) {
		int target_prow = (i / mb) % nprow;
		if (myrow != target_prow)
			continue;
		int locali = mb * (i / (nb * npcol)) + (i % mb); 
		
		for (int j = 0; j < n; j++) {
			int target_pcol = (j / nb) % npcol;
			if (mycol != target_pcol)
				continue;
			int localj = nb * (j / (nb * npcol)) + (j % nb); 
			A[locali + ldA * localj] = (i - j) / (j + 1);
		}
	}

	/* TAP protocol */
	if (rank == 0)
		printf("1..2\n");

	/********************** copy local vector to distributed matrix */
	double local[m];
	for (int i = 0; i < m; i++)
		local[i] = i;

	int modified_jA = 751;
	extrablacs_dgeld2d(local, m, A, modified_jA, descA);

	/* check copy */
	int ok = 1;
	for (int i = 0; i < m; i++) {
		int target_prow = (i / mb) % nprow;
		if (myrow != target_prow)
			continue;
		int locali = mb * (i / (mb * npcol)) + (i % mb); 
		
		for (int j = 0; j < n; j++) {
			int target_pcol = (j / nb) % npcol;
			if (mycol != target_pcol)
				continue;
			int localj = nb * (j / (nb * npcol)) + (j % nb);
			double Aij = A[locali + localj * ldA];
			if (j == modified_jA)
				ok &= (Aij == i);
			else
				ok &= (Aij == (i - j) / (j + 1));
		}
	}
	if (rank == 0)
		MPI_Reduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	else
		MPI_Reduce(&ok, NULL, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	if (rank == 0 && ok)
		printf("ok 1 - rvec2dmat\n");
	else if (rank == 0 && !ok)
		printf("nok 1 - rvec2dmat\n");

	/********************** copy distributed sub-matrix to local array */

	int nn = 20;
	int mm = 51;
	double L[nn * mm];

	int iA = 100;
	int jA = 740;
	extrablacs_dgedl2d(mm, nn, A, iA, jA, descA, L, mm);

	/* check copy */
	ok = 1;
	for (int ii = 0; ii < mm; ii++)
		for (int jj = 0; jj < nn; jj++) {
			int i = iA-1 + ii;
			int j = jA-1 + jj;
			double Aij = L[ii + jj * mm];
			if (j == modified_jA)
				ok &= (Aij == i);
			else
				ok &= (Aij == (i - j) / (j + 1));
		}
	if (rank == 0)
		MPI_Reduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	else
		MPI_Reduce(&ok, NULL, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	if (rank == 0 && ok)
		printf("ok 2 - dmat2rmat\n");
	else if (rank == 0 && !ok)
		printf("nok 2 - dmat2rmat\n");

	/* BLACS and MPI cleanup */
	MPI_Finalize();
	return EXIT_SUCCESS;
}
