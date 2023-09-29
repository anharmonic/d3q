#include <math.h>
#include <stdlib.h>

#include "minpack.h"

static int query_worksize(const int *n, const int *ldfjac)
{
        int c1 = 1;
        int lwork = -1;
        double work[1];
        int info;

        /* query LAPACK */
        dgeqrf_(n, n, NULL, ldfjac, NULL, work, &lwork, &info);
        if (info != 0)
                return -1;
        int needed_dgeqrf = work[0];
        if (needed_dgeqrf < *n)
                return -1;

        dormqr_("Left", "Transpose", n, &c1, n, NULL, ldfjac, NULL, NULL, n, work, &lwork, &info);
        if (info != 0)
                return -1;
        int needed_dormqr = work[0];
        if (needed_dormqr < *n)
                return -1;
        
	dorgqr_(n, n, n, NULL, ldfjac, NULL, work, &lwork, &info);
	if (info != 0)
                return -1;
        int needed_dorgqr = work[0];
        if (needed_dorgqr < *n)
                return -1;
        
        /* take max */
        int max = needed_dgeqrf;
        if (max < needed_dormqr)
        	max = needed_dormqr;
        if (max < needed_dorgqr)
        	max = needed_dorgqr;
        return max;
}

void hybrbase(minpack_func_n fcn_dif, minpack_func_nj fcn_der,
	      const int *n, double *x, double *fvec, double *fjac, const int *ldfjac, 
	      const double *xtol, const int *maxfev,
	      const int *ml, const int *mu, const double *epsfcn, double *diag, const int *mode,
	      const double *factor, const int *nprint, int *info, int *nfev, int *njev,
	      double *r, const int *lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4)
{
	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;
	int fjac_dim1 = *ldfjac;
	int c1 = 1;
	int msum = *ml + *mu + 1;
	if (msum > *n)
		msum = *n;
	int iflag = 0;
	*info = 0;
	*nfev = 0;
	double *work = NULL;

	/* check the input parameters for errors. */
	if (*n <= 0 || *ldfjac < *n || *xtol < 0 || *maxfev <= 0 || *factor <= 0 || *lr < *n * (*n + 1) / 2)
		goto fini;
	if ((fcn_dif == NULL) && (fcn_der == NULL))
		goto fini;
	if ((fcn_dif != NULL) && (fcn_der != NULL))
		goto fini;
	if (*mode == 2) {
		for (long int j = 0; j < *n; ++j)
			if (diag[j] <= 0)
				goto fini;
	}
	
	int lwork = query_worksize(n, ldfjac);
        if (lwork < 0)
                goto fini;
        work = malloc(lwork * sizeof(*work));
        if (work == NULL)
                goto fini;

	/* Evaluate the function at the starting point and calculate its norm. */
	iflag = 1;
	if (fcn_der)
		(*fcn_der) (n, x, fvec, fjac, ldfjac, &iflag);
	if (fcn_dif)
		(*fcn_dif) (n, x, fvec, &iflag);
	*nfev = 1;
	if (iflag < 0)
		goto fini;
	double fnorm = enorm_(n, fvec);

	/* initialize iteration counter and monitors. */
	int iter = 1;
	int ncsuc = 0;
	int ncfail = 0;
	int nslow1 = 0;
	int nslow2 = 0;

	/* beginning of the outer loop. */
	for (;;) {
		int jeval = 1;

		/* calculate the jacobian matrix. */
		iflag = 2;
		if (fcn_der) {
			(*fcn_der) (n, x, fvec, fjac, ldfjac, &iflag);
			*njev += 1;
		}
		if (fcn_dif) {
			fdjac1_(fcn_dif, n, x, fvec, fjac, ldfjac, &iflag, ml, mu, epsfcn, wa1, wa2);
			*nfev += msum;
		}
		if (iflag < 0)
			goto fini;

		/* on the first iteration and if mode is 1, scale according */
		/* to the norms of the columns of the initial jacobian. */
		double xnorm, delta;
		if (iter == 1) {
			if (*mode != 2)
				for (long int j = 0; j < *n; ++j) {
					diag[j] = wa2[j];
					if (wa2[j] == 0)
						diag[j] = 1;
				}

			/* on the first iteration, calculate the norm of the scaled x */
			/* and initialize the step bound delta. */
			for (long int j = 0; j < *n; ++j)
				wa3[j] = diag[j] * x[j];
			xnorm = enorm_(n, wa3);
			delta = *factor * xnorm;
			if (delta == 0)
				delta = *factor;
		}

		/* compute the QR factorization of the jacobian. */
		double *tau = wa1;
		int lapack_info;
		dgeqrf_(n, n, fjac, ldfjac, tau, work, &lwork, &lapack_info);	// LAPACK
		if (lapack_info != 0)
			goto fini;

		/* qtf <-- (Q transpose)*fvec */
		for (long int i = 0; i < *n; ++i)
			wa4[i] = fvec[i];
		dormqr_("Left", "Transpose", n, &c1, n, fjac, ldfjac, tau, wa4, n, work, &lwork, &lapack_info);
		if (lapack_info != 0)
			goto fini;
		for (long int j = 0; j < *n; ++j)
			qtf[j] = wa4[j];

		/* copy the triangular factor of the QR factorization into R. */
		int sing = 0;
		for (long int j = 0; j < *n; ++j) {
			int l = j;
			for (long int i = 0; i < j; ++i) {
				r[l] = fjac[i + j * fjac_dim1];
				l = l + *n - (i+1);
			}
			r[l] = fjac[j + j * fjac_dim1];
			if (r[l] == 0)
				sing = 1;
		}

		/* accumulate the orthogonal factor in fjac. */
		dorgqr_(n, n, n, fjac, ldfjac, tau, work, &lwork, &lapack_info);
		if (lapack_info != 0)
			goto fini;

		/* rescale if necessary. */
		if (*mode != 2)
			for (long int j = 0; j < *n; ++j)
				diag[j] = fmax(diag[j], wa2[j]);

		/* beginning of the inner loop. */
		for (;;) {

			/* if requested, call fcn to enable printing of iterates. */
			if (*nprint > 0) {
				iflag = 0;
				if ((iter - 1) % *nprint == 0) {
					if (fcn_der)
						(*fcn_der) (n, x, fvec, fjac, ldfjac, &iflag);
					if (fcn_dif)
						(*fcn_dif) (n, x, fvec, &iflag);
				}
				if (iflag < 0)
					goto fini;
			}

			/* determine the direction p. */
			dogleg_(n, r, lr, diag, qtf, &delta, wa1, wa2, wa3);

			/* store the direction p and x + p. calculate the norm of p. */
			for (long int j = 0; j < *n; ++j) {
				wa1[j] = -wa1[j];
				wa2[j] = x[j] + wa1[j];
				wa3[j] = diag[j] * wa1[j];
			}
			double pnorm = enorm_(n, wa3);

			/* on the first iteration, adjust the initial step bound. */
			if (iter == 1)
				delta = fmin(delta, pnorm);

			/* evaluate the function at x + p and calculate its norm. */
			iflag = 1;
			if (fcn_der)
				(*fcn_der) (n, wa2, wa4, fjac, ldfjac, &iflag);
			if (fcn_dif)
				(*fcn_dif) (n, wa2, wa4, &iflag);
			*nfev += 1;
			if (iflag < 0)
				goto fini;
			double fnorm1 = enorm_(n, wa4);

			/* compute the scaled actual reduction. */
			double actred = -1;
			if (fnorm1 < fnorm) {
				double d = fnorm1 / fnorm;
				actred = 1 - d * d;
			}

			/* compute the scaled predicted reduction. */
			int l = 0;
			for (long int i = 0; i < *n; ++i) {
				double sum = 0;
				for (long int j = i; j < *n; ++j) {
					sum += r[l] * wa1[j];
					++l;
				}
				wa3[i] = qtf[i] + sum;
			}
			double temp = enorm_(n, wa3);
			double prered = 0;
			if (temp < fnorm) {
				double d = temp / fnorm;
				prered = 1 - d * d;
			}

			/* compute the ratio of the actual to the predicted reduction. */
			double ratio = 0;
			if (prered > 0)
				ratio = actred / prered;

			/* update the step bound. */
			if (ratio < 0.1) {
				ncsuc = 0;
				++ncfail;
				delta = 0.5 * delta;
			} else {
				ncfail = 0;
				++ncsuc;
				if (ratio >= 0.5 || ncsuc > 1)
					delta = fmax(delta, pnorm / 0.5);
				if (fabs(ratio - 1) <= 0.1)
					delta = pnorm / 0.5;
			}

			/* test for successful iteration. */
			if (ratio >= 0.0001) {
				/* successful iteration. update x, fvec, and their norms. */
				for (long int j = 0; j < *n; ++j) {
					x[j] = wa2[j];
					wa2[j] = diag[j] * x[j];
					fvec[j] = wa4[j];
				}
				xnorm = enorm_(n, wa2);
				fnorm = fnorm1;
				++iter;
			}

			/* determine the progress of the iteration. */
			++nslow1;
			if (actred >= 0.001)
				nslow1 = 0;
			if (jeval)
				++nslow2;
			if (actred >= 0.1)
				nslow2 = 0;

			/* test for convergence. */
			if (delta <= *xtol * xnorm || fnorm == 0)
				*info = 1;
			if (*info != 0)
				goto fini;

			/* tests for termination and stringent tolerances. */
			if (*nfev >= *maxfev)
				*info = 2;
			if (0.1 * fmax(0.1 * delta, pnorm) <= epsmch * xnorm)
				*info = 3;
			if (nslow2 == 5)
				*info = 4;
			if (nslow1 == 10)
				*info = 5;
			if (*info != 0)
				goto fini;

			/* criterion for recalculating jacobian. */
			if (ncfail == 2)
				break;	/* exit the inner loop */

			/* calculate the rank one modification to the jacobian */
			/* and update qtf if necessary. */
			for (long int j = 0; j < *n; ++j) {
				double sum = 0;
				for (long int i = 0; i < *n; ++i)
					sum += fjac[i + j * fjac_dim1] * wa4[i];
				wa2[j] = (sum - wa3[j]) / pnorm;
				wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
				if (ratio >= 0.0001)
					qtf[j] = sum;
			}

			/* compute the QR factorization of the updated jacobian. */
			r1updt_(n, n, r, lr, wa1, wa2, wa3, &sing);
			r1mpyq_(n, n, fjac, ldfjac, wa2, wa3);
			r1mpyq_(&c1, n, qtf, &c1, wa2, wa3);
			jeval = 0;
		}		/* inner loop. */
	}			/* outer loop. */

 fini:
	/* termination, either normal or user imposed. */
 	free(work);
	if (iflag < 0)
		*info = iflag;
	iflag = 0;
	if (*nprint > 0) {
		if (fcn_der)
			(*fcn_der) (n, x, fvec, fjac, ldfjac, &iflag);
		if (fcn_dif)
			(*fcn_dif) (n, x, fvec, &iflag);
	}
}
