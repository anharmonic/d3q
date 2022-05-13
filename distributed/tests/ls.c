#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "pminpack.h"

#include "ls.h"

char *problem_name[] = {
	"Linear function --- full rank",
	"Linear function --- rank 1",
	"Linear function --- rank 1 with zero columns and rows",
	"Rosenbrock function",
	"Helical valley function",
	"Powell singular function",
	"Freudenstein and Roth function",
	"Bard function",
	"Kowalik and Osborne function",
	"Meyer function",
	"Watson function",
	"Box 3-dimensional function",
	"Jennrich and Sampson function",
	"Brown and Dennis function",
	"Chebyquad function",
	"Brown almost-linear function",
	"Osborne 1 function",
	"Osborne 2 function",
};

/* these are the test cases run by the MINPACK-1 test suite, along with the 
   values given by the original fortran-77 program. */
struct test_case tests[] = {
    // prob,    n,    m, factor,   nfev,  njev,  info,   fnorm2,   fnorm2_lastchance
    // Linear function --- full rank
    {     1,    5,   10,      1,      3,     2,     1,   2.236068, -1},   // 01
    {     1,    5,   50,      1,      3,     2,     1,   6.708204, -1},   // 02
    // Linear function --- rank 1
    {     2,    5,   10,      1,      3,     2,     1,   1.463850, -1},   // 03
    {     2,    5,   50,      1,      5,     2,     1,   3.482630, -1},   // 04
    // Linear function --- rank 1 with zero columns and rows
    {     3,    5,   10,      1,      3,     2,     1,   1.909727, -1},   // 05
    {     3,    5,   50,      1,      3,     2,     1,   3.691729, -1},   // 06
    // Rosenbrock function
    {     4,    2,    2,      1,     22,    16,     2,   0       , -1},   // 07
    {     4,    2,    2,     10,      8,     5,     2,   0       , -1},   // 08
    {     4,    2,    2,    100,      7,     5,     2,   0       , -1},   // 09
    // Helical valley function
    {     5,    3,    3,      1,     12,     9,     2,   0       , -1},   // 10
    {     5,    3,    3,     10,     20,    15,     2,   0       , -1},   // 11
    {     5,    3,    3,    100,     19,    16,     2,   0       , -1},   // 12
    // Powell singular function
    {     6,    4,    4,      1,    240,   190,     5,   0       , -1},   // 13
    {     6,    4,    4,     10,    236,   191,     5,   0       , -1},   // 14
    {     6,    4,    4,    100,    240,   190,     5,   0       , -1},   // 15
    // Freudenstein and Roth function
    {     7,    2,    2,      1,     14,     8,     1,   0       , 6.998875},   // 16
    {     7,    2,    2,     10,     19,    12,     1,   0       , 6.998875},   // 17
    {     7,    2,    2,    100,     24,    17,     1,   0       , 6.998875},   // 18
    // Bard function
    {     8,    3,   15,      1,      6,     5,     1,   0.090636, 4.174769},   // 19  
    {     8,    3,   15,     10,     38,    37,     1,   0.090636, 4.174769},   // 20
    {     8,    3,   15,    100,     15,    14,     1,   0.090636, 4.174769},   // 21
    // Kowalik and Osborne function
    {     9,    4,   11,      1,     18,    16,     1,   0.017536, 0.0320521},   // 22
    {     9,    4,   11,     10,     78,    69,     1,   0.017536, 0.0320521},   // 23
    {     9,    4,   11,    100,    255,   187,     5,   0.017536, 0.0320521},   // 24
    // Meyer function
    {    10,    3,   16,      1,    126,   116,     3,   9.377942, -1},   // 25
    {    10,    3,   16,     10,      4,     3,     4,   9.377942, -1},  // 26
    // Watson function
    {    11,    6,   31,      1,      8,     7,     1,   0.0478296, -1}, // 27
    {    11,    6,   31,     10,     14,    13,     1,   0.0478296, -1},   // 28
    {    11,    6,   31,    100,     15,    14,     1,   0.0478296, -1},   // 29

    {    11,    9,   31,      1,     10,     6,     1,   0.1183215e-02, -1},   // 30
    {    11,    9,   31,     10,     30,    16,     2,   0.1183115e-02, -1},   // 31
    {    11,    9,   31,    100,     30,    19,     2,   0.1183115e-02, -1},   // 32
    
    {    11,   12,   31,      1,     22,    11,     2,   2.1731015e-05, -1},   // 33
    {    11,   12,   31,     10,     34,    15,     2,   2.1731015e-05, -1},   // 34
    {    11,   12,   31,    100,     53,    31,     2,   2.1731015e-05, -1},   // 35
    // Box 3-dimensional function
    {    12,    3,   10,      1,      7,     6,     2,   0,           -1},   // 36
    // Jennrich and Sampson function
    {    13,    2,   10,      1,     21,    12,     1,   11.15177116, -1},   // 37
    // Brown and Dennis function
    {    14,    4,   20,      1,    216,   197,     5,   292.954262, -1}, // 38
    {    14,    4,   20,     10,     47,    36,     1,   292.954262, -1},   // 39
    {    14,    4,   20,    100,    212,   197,     5,   292.954262, -1}, // 40
    // Chebyquad function
    {    15,    1,    8,      1,      2,     1,     1,   1.886238     , -1},   // 41
    {    15,    1,    8,     10,     29,    28,     1,   1.884248     , -1},   // 42
    {    15,    1,    8,    100,     47,    46,     1,   1.884248     , -1},   // 43
    {    15,    8,    8,      1,     39,    20,     1,   0.05930320396, -1},   // 44
    {    15,    9,    9,      1,     12,     9,     2,   0            , -1},   // 45
    {    15,   10,   10,      1,     25,    12,     1,   0.08064707062, -1},   // 46
    // Brown almost-linear function
    {    16,   10,   10,      1,      7,     5,     1,   0            ,1},   // 47
    {    16,   10,   10,     10,     13,     8,     2,   0            ,1},   // 48
    {    16,   10,   10,    100,     34,    32,     2,   0            ,1},   // 49
    {    16,   30,   30,      1,     11,    10,     2,   0            ,1},   // 50
    {    16,   40,   40,      1,     11,    10,     2,   0            ,1},   // 51
    // Osborne 1 function
    {    17,    5,   33,      1,     18,    15,     1,   7.39248943e-3, -1}, // 52
    // Osborne 2 function
    {    18,   11,   65,      1,     16,    12,     1,   0.2003440E+00, -1}, // 53
};

void commentator(int ic, double *x, double *fvec, double ftol, double xtol, int nfev, int njev, int info)
{
	--fvec;
	--x;
	int n = tests[ic].n;
	int m = tests[ic].m;
	int nprob = tests[ic].nprob;

	double target = tests[ic].fnorm2;
	double target_alt = tests[ic].fnorm2_lastchance;

	double fnorm = enorm(m, &fvec[1]);
	double df = fnorm - target;
	double df2 = (target_alt >= 0) ? fnorm - target_alt : -1;
	double dx = -1;
	double dx2 = -1;
	double dx3 = -1;
	double tmp[n + 1];

	/* compute distance to known minimizers */
	switch (nprob) {
		/* linear function - full rank. */
	case 1:
		for (int j = 1; j <= n; j++)
			tmp[j] = x[j] + 1;
		dx = enorm(n, &tmp[1]);
		break;

		/* rosenbrock function. */
	case 4:
		tmp[1] = x[1] - 1;
		tmp[2] = x[2] - 1;
		dx = enorm(m, &tmp[1]);
		break;

		/* helical valley function. */
	case 5:
		tmp[1] = x[1] - 1;
		tmp[2] = x[2];
		tmp[3] = x[3];
		dx = enorm(n, &tmp[1]);
		break;

		/* powell singular function. */
	case 6:
		dx = enorm(n, &fvec[1]);
		break;

		/* freudenstein and roth function. */
	case 7:
		tmp[1] = x[1] - 5;
		tmp[2] = x[2] - 4;
		dx = enorm(n, &tmp[1]);
		tmp[1] = x[1] - 11.41;
		tmp[2] = x[2] + 0.8968;
		dx2 = enorm(n, &tmp[1]);
		break;

		/* bard function. */
	case 8:
		dx2 = fabs(x[1] - 0.8406);
		break;

		/* kowalik and osborne function. */
	case 9:
		dx2 = fabs(x[2] + 14.07);
		break;

		/* box 3-dimensional function. */
	case 12:
		tmp[1] = x[1] - 1;
		tmp[2] = x[2] - 10;
		tmp[3] = x[3] - 1;
		dx = enorm(n, &tmp[1]);
		tmp[1] = x[1] - 10;
		tmp[2] = x[2] - 1;
		tmp[3] = x[3] + 1;
		dx2 = enorm(n, &tmp[1]);
		tmp[1] = 0;
		tmp[2] = x[2] - x[1];
		tmp[3] = x[3];
		dx3 = enorm(n, &tmp[1]);
		break;

		/* jennrich and sampson function. */
	case 13:
		if (m == 10) {
			tmp[1] = x[1] - 0.2578;
			tmp[2] = x[2] - 0.2578;
			dx = enorm(n, &tmp[1]);
		}
		break;

		/* brown almost-linear function. */
	case 16:
		for (int j = 1; j <= m - 1; j++)
			tmp[j] = x[1];
		tmp[m] = pow(x[1], 1 - n);
		dx = enorm(n, &tmp[1]);
		for (int j = 1; j <= m - 1; j++)
			tmp[j] = x[j];
		tmp[m] = x[m] - (n + 1);
		dx2 = enorm(n, &tmp[1]);
		break;
	}

	/* print diagnostics */
	double rel_err = (fnorm - target) / target;
	double rel_err2 = (fnorm - target_alt) / target_alt;

	bool ok = df < ftol;
	bool ok2 = (target_alt > -1) && (df2 < ftol);

	if (info != tests[ic].info)
		printf("\t# NOTICE: mismatching status (got %d, expected %d)\n", info, tests[ic].info);

	if (info == 5)
		printf("\t# WARNING: unreliable result (exceeded maxfev; increasing it may work)\n");

	if (rel_err < -0.01 || df <= -ftol) {
		printf("\t# WONDERFUL: better final residual norm\n");
		printf("\t#        residual norm    : %15.7e\n", fnorm);
		printf("\t#        target minimum   : %15.7e\n", target);
		printf("\t#        relative distance: %15f\n", rel_err);
	}

	if (ok && dx > xtol)
		printf("\t# NOTICE: far from (known) minimizer: %15.7e\n", dx);

	if (!ok && target_alt < 0) {
		printf("\t# FATAL: missed (known) minimum\n");
		printf("\t#        residual norm    : %15.7e\n", fnorm);
		printf("\t#        target minimum   : %15.7e\n", target);
		printf("\t#        relative distance: %15f\n", rel_err);
	}

	if (!ok && ok2) {
		printf("\t# OK: missed (known) best minimum but hit other local minimum\n");
		printf("\t#     residual norm    : %15.7e\n", fnorm);
		printf("\t#     target minimum   : %15.7e\n", target);
		printf("\t#     relative distance: %f\n", rel_err);
		printf("\t#     other minimum    : %15.7e\n", target_alt);
		printf("\t#     relative distance: %15f\n", rel_err2);
		if (dx2 > xtol)
			printf("\t# NOTICE: far from other (known) minimizer: %.9e\n", dx2);
	}

	if (target_alt > -1 && !ok && !ok2) {
		printf("\t# FATAL: missed all (known) minima\n");
		printf("\t#        residual norm    : %15.7e\n", fnorm);
		printf("\t#        target minimum   : %15.7e\n", target);
		printf("\t#        relative distance: %15f\n", rel_err);
		printf("\t#        other minimum    : %15.7e\n", target_alt);
		printf("\t#        relative distance: %15f\n", rel_err2);
	}

	if (nprob == 12 && ok && dx > xtol && dx2 > xtol && dx3 > xtol)
		printf("\t# NOTICE: far from all (known) minimizer: %.9e\n", dx3);

	if (nfev < tests[ic].nfev)
		printf("\t# GOOD: less function evaluations (did %d, expected %d)\n", nfev, tests[ic].nfev);
	if (nfev > tests[ic].nfev)
		printf("\t# BAD: more function evaluations (did %d, expected %d)\n", nfev, tests[ic].nfev);
	if (njev < tests[ic].njev)
		printf("\t# GOOD: less jacobian evaluations (did %d, expected %d)\n", njev, tests[ic].njev);
	if (njev > tests[ic].njev)
		printf("\t# BAD: more jacobian evaluations (did %d, expected %d)\n", njev, tests[ic].njev);
}

/*
 * subroutine initpt
 * 
 * this subroutine specifies the standard starting points for the
 * functions defined by subroutine ssqfcn. the subroutine returns
 * in x a multiple (factor) of the standard starting point. for
 * the 11th function the standard starting point is zero, so in
 * this case, if factor is not unity, then the subroutine returns
 * the vector  x(j) = factor, j=1,...,n.
 *
 * the subroutine statement is
 *
 *     subroutine initpt(n,x,nprob,factor)
 *
 * where 
 *
 *     n is a positive integer input variable.
 *
 *     x is an output array of length n which contains the standard
 *     starting point for problem nprob multiplied by factor.
 *
 *     nprob is a positive integer input variable which defines the
 *     number of the problem. nprob must not exceed 18.
 *
 *     factor is an input variable which specifies the multiple of
 *     the standard starting point. if factor is unity, no
 *     multiplication is performed.
 *
 *  argonne national laboratory. minpack project. march 1980.
 *  burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void initpt(int n, double *x, int nprob, double factor)
{
	/* Parameter adjustments (fortran --> C) */
	--x;
	double h;
	
	switch (nprob) {
		/* linear function - full rank or rank 1. */
	case 1:
	case 2:
	case 3:
		for (int j = 1; j <= n; ++j)
			x[j] = 1.;
		break;

		/* rosenbrock function. */
	case 4:
		x[1] = -1.2;
		x[2] = 1.;
		break;

		/* helical valley function. */
	case 5:
		x[1] = -1.;
		x[2] = 0.;
		x[3] = 0.;
		break;

		/* powell singular function. */
	case 6:
		x[1] = 3.;
		x[2] = -1.;
		x[3] = 0.;
		x[4] = 1.;
		break;

		/* freudenstein and roth function. */
	case 7:
		x[1] = .5;
		x[2] = -2.;
		break;

		/* bard function. */
	case 8:
		x[1] = 1.;
		x[2] = 1.;
		x[3] = 1.;
		break;

		/* kowalik and osborne function. */
	case 9:
		x[1] = .25;
		x[2] = .39;
		x[3] = .415;
		x[4] = .39;
		break;

		/* meyer function. */
	case 10:
		x[1] = .02;
		x[2] = 4e3;
		x[3] = 250.;
		break;

		/* watson function. */
	case 11:
		for (int j = 1; j <= n; ++j)
			x[j] = 0.;
		break;

		/* box 3-dimensional function. */
	case 12:
		x[1] = 0.;
		x[2] = 10.;
		x[3] = 20.;
		break;

		/* jennrich and sampson function. */
	case 13:
		x[1] = .3;
		x[2] = .4;
		break;

		/* brown and dennis function. */
	case 14:
		x[1] = 25.;
		x[2] = 5.;
		x[3] = -5.;
		x[4] = -1.;
		break;

		/* chebyquad function. */
	case 15:
		h = 1. / (n + 1);
		for (int j = 1; j <= n; ++j)
			x[j] = j * h;
		break;

		/* brown almost-linear function. */
	case 16:
		for (int j = 1; j <= n; ++j)
			x[j] = .5;
		break;

		/* osborne 1 function. */
	case 17:
		x[1] = .5;
		x[2] = 1.5;
		x[3] = -1.;
		x[4] = .01;
		x[5] = .02;
		break;

		/* osborne 2 function. */
	case 18:
		x[1] = 1.3;
		x[2] = .65;
		x[3] = .65;
		x[4] = .7;
		x[5] = .6;
		x[6] = 3.;
		x[7] = 5.;
		x[8] = 7.;
		x[9] = 2.;
		x[10] = 4.5;
		x[11] = 5.5;
		break;
	}

	/* compute multiple of initial point. */
	if (factor != 1.) {
		if (nprob != 11) {
			for (int j = 1; j <= n; ++j)
				x[j] *= factor;
		} else {
			for (int j = 1; j <= n; ++j)
				x[j] = factor;
		}
	}
}

/*
 * subroutine ssqfcn
 *
 * This subroutine defines the functions of eighteen nonlinear
 * least squares problems.  The allowable values of (m,n) for
 * functions 1,2 and 3 are variable but with m >= n.
 * For functions 4,5,6,7,8,9 and 10 the values of (m,n) are
 * (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) and (16,3), respectively.
 * Function 11 (watson) has m = 31 with n usually 6 or 9.
 * However, any n, n = 2,...,31, is permitted.
 * Functions 12,13 and 14 have n = 3,2 and 4, respectively, but
 * allow any m >= n, with the usual choices being 10,10 and 20.
 * Function 15 (chebyquad) allows m and n variable with m >= n.
 * Function 16 (brown) allows n variable with m = n.
 * For functions 17 and 18, the values of (m,n) are
 * (33,5) and (65,11), respectively.
 *
 * The subroutine statement is
 *
 *       subroutine ssqfcn(m,n,x,fvec,nprob)
 *
 * where
 *
 *       m and n are positive integer input variables. n must not
 *         exceed m.
 *
 *       x is an input array of length n.
 *
 *       fvec is an output array of length m which contains the nprob
 *         function evaluated at x.
 *
 *       nprob is a positive integer input variable which defines the
 *         number of the problem. nprob must not exceed 18.
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. More
 */
void ssqfcn(int m, int n, const double *x, double *fvec, int nprob)
{
	/* Initialized data */
	const double v[11] = { 4., 2., 1., .5, .25, .167, .125, .1, .0833, .0714, .0625 };
	const double y1[15] = { .14, .18, .22, .25, .29, .32, .35, .39, .37, .58, .73,
		.96, 1.34, 2.1, 4.39
	};
	const double y2[11] = { .1957, .1947, .1735, .16, .0844, .0627, .0456,
		.0342, .0323, .0235, .0246
	};
	const double y3[16] = { 34780., 28610., 23650., 19630., 16370., 13720.,
		11540., 9744., 8261., 7030., 6005., 5147., 4427., 3820., 3307., 2872.
	};
	const double y4[33] = { .844, .908, .932, .936, .925, .908, .881, .85, .818,
		.784, .751, .718, .685, .658, .628, .603, .58, .558, .538, .522, .506, .49,
		.478, .467, .457, .448, .438, .431, .424, .42, .414, .411, .406
	};
	const double y5[65] = { 1.366, 1.191, 1.112, 1.013, .991, .885, .831, .847,
		.786, .725, .746, .679, .608, .655, .616, .606, .602, .626, .651, .724, .649,
		.649, .694, .644, .624, .661, .612, .558, .533, .495, .5, .423, .395, .375,
		.372, .391, .396, .405, .428, .429, .523, .562, .607, .653, .672, .708, .633,
		.668, .645, .632, .591, .559, .597, .625, .739, .71, .729, .72, .636, .581,
		.428, .292, .162, .098, .054
	};

	/* Local variables */
	double sum, tmp1, tmp2, tmp3, tmp4, temp;

	/* Parameter adjustments (fortran --> C) */
	--fvec;
	--x;
	double twopi;
	
	switch (nprob) {
		/* linear function - full rank. */
	case 1:
		sum = 0;
		for (int j = 1; j <= n; ++j)
			sum += x[j];
		temp = 2 * sum / m + 1;
		for (int i = 1; i <= m; ++i) {
			fvec[i] = -temp;
			if (i <= n)
				fvec[i] += x[i];
		}
		break;

		/* linear function - rank 1. */
	case 2:
		sum = 0;
		for (int j = 1; j <= n; ++j)
			sum += j * x[j];
		for (int i = 1; i <= m; ++i)
			fvec[i] = i * sum - 1;
		break;

		/* linear function - rank 1 with zero columns and rows. */
	case 3:
		sum = 0;
		for (int j = 2; j <= n - 1; ++j)
			sum += j * x[j];
		for (int i = 1; i <= m; ++i)
			fvec[i] = (i - 1) * sum - 1;
		fvec[m] = -1;
		break;

		/* rosenbrock function. */
	case 4:
		fvec[1] = 10 * (x[2] - x[1] * x[1]);
		fvec[2] = 1 - x[1];
		break;

		/* helical valley function. */
	case 5:
		twopi = 2 * M_PI;
		tmp1 = x[2] < 0 ? -0.25 : 0.25;
		if (x[1] > 0)
			tmp1 = atan(x[2] / x[1]) / twopi;
		if (x[1] < 0)
			tmp1 = atan(x[2] / x[1]) / twopi + .5;
		tmp2 = sqrt(x[1] * x[1] + x[2] * x[2]);
		fvec[1] = 10 * (x[3] - 10 * tmp1);
		fvec[2] = 10 * (tmp2 - 1);
		fvec[3] = x[3];
		break;

		/* powell singular function. */
	case 6:
		fvec[1] = x[1] + 10 * x[2];
		fvec[2] = sqrt(5) * (x[3] - x[4]);
		/* Computing 2nd power */
		tmp1 = x[2] - 2 * x[3];
		fvec[3] = tmp1 * tmp1;
		/* Computing 2nd power */
		tmp2 = x[1] - x[4];
		fvec[4] = sqrt(10) * (tmp2 * tmp2);
		break;

		/* freudenstein and roth function. */
	case 7:
		fvec[1] = -13 + x[1] + ((5 - x[2]) * x[2] - 2) * x[2];
		fvec[2] = -29 + x[1] + ((1 + x[2]) * x[2] - 14) * x[2];
		break;

		/* bard function. */
	case 8:
		for (int i = 1; i <= 15; ++i) {
			tmp1 = i;
			tmp2 = 16 - i;
			tmp3 = tmp1;
			if (i > 8)
				tmp3 = tmp2;
			fvec[i] = y1[i - 1] - (x[1] + tmp1 / (x[2] * tmp2 + x[3] * tmp3));
		}
		break;

		/* kowalik and osborne function. */
	case 9:
		for (int i = 1; i <= 11; ++i) {
			tmp1 = v[i - 1] * (v[i - 1] + x[2]);
			tmp2 = v[i - 1] * (v[i - 1] + x[3]) + x[4];
			fvec[i] = y2[i - 1] - x[1] * tmp1 / tmp2;
		}
		break;

		/* meyer function. */
	case 10:
		for (int i = 1; i <= 16; ++i) {
			temp = 5 * i + 45 + x[3];
			tmp1 = x[2] / temp;
			tmp2 = exp(tmp1);
			fvec[i] = x[1] * tmp2 - y3[i - 1];
		}
		break;

		/* watson function. */
	case 11:
		for (int i = 1; i <= 29; ++i) {
			double div = i / 29.;
			double s1 = 0;
			double dx = 1;
			for (int j = 2; j <= n; ++j) {
				s1 += (j - 1) * dx * x[j];
				dx = div * dx;
			}
			double s2 = 0;
			dx = 1;
			for (int j = 1; j <= n; ++j) {
				s2 += dx * x[j];
				dx = div * dx;
			}
			fvec[i] = s1 - s2 * s2 - 1;
		}
		fvec[30] = x[1];
		fvec[31] = x[2] - x[1] * x[1] - 1;
		break;

		/* box 3-dimensional function. */
	case 12:
		for (int i = 1; i <= m; ++i) {
			tmp1 = i / 10.;
			fvec[i] = exp(-tmp1 * x[1]) - exp(-tmp1 * x[2]) + (exp(-i) - exp(-tmp1)) * x[3];
		}
		break;

		/* jennrich and sampson function. */
	case 13:
		for (int i = 1; i <= m; ++i) {
			fvec[i] = 2 + 2 * i - exp(i * x[1]) - exp(i * x[2]);
		}
		break;

		/* brown and dennis function. */
	case 14:
		for (int i = 1; i <= m; ++i) {
			temp = i / 5.;
			tmp1 = x[1] + temp * x[2] - exp(temp);
			tmp2 = x[3] + sin(temp) * x[4] - cos(temp);
			fvec[i] = tmp1 * tmp1 + tmp2 * tmp2;
		}
		break;

		/* chebyquad function. */
	case 15:
		for (int i = 1; i <= m; ++i)
			fvec[i] = 0;
		for (int j = 1; j <= n; ++j) {
			tmp1 = 1;
			tmp2 = 2 * x[j] - 1;
			temp = 2 * tmp2;
			for (int i = 1; i <= m; ++i) {
				fvec[i] += tmp2;
				double ti = temp * tmp2 - tmp1;
				tmp1 = tmp2;
				tmp2 = ti;
			}
		}
		double dx = 1. / n;
		double iev = -1;
		for (int i = 1; i <= m; ++i) {
			fvec[i] = dx * fvec[i];
			if (iev > 0)
				fvec[i] += 1. / (i * i - 1);
			iev = -iev;
		}
		break;

		/* brown almost-linear function. */
	case 16:
		sum = -(n + 1);
		double prod = 1;
		for (int j = 1; j <= n; ++j) {
			sum += x[j];
			prod *= x[j];
		}
		for (int i = 1; i <= n; ++i)
			fvec[i] = x[i] + sum;
		fvec[n] = prod - 1;
		break;

		/* osborne 1 function. */
	case 17:
		for (int i = 1; i <= 33; ++i) {
			temp = 10. * (i - 1);
			tmp1 = exp(-x[4] * temp);
			tmp2 = exp(-x[5] * temp);
			fvec[i] = y4[i - 1] - (x[1] + x[2] * tmp1 + x[3] * tmp2);
		}
		break;

		/* osborne 2 function. */
	case 18:
		for (int i = 1; i <= 65; ++i) {
			temp = (i - 1) / 10.;
			tmp1 = exp(-x[5] * temp);
			/* Computing 2nd power */
			double e = temp - x[9];
			tmp2 = exp(-x[6] * (e * e));
			/* Computing 2nd power */
			double f = temp - x[10];
			tmp3 = exp(-x[7] * (f * f));
			/* Computing 2nd power */
			double g = temp - x[11];
			tmp4 = exp(-x[8] * (g * g));
			fvec[i] = y5[i - 1] - (x[1] * tmp1 + x[2] * tmp2 + x[3] * tmp3 + x[4] * tmp4);
		}
		break;
	}
}

/*
 *     Subroutine ssqjac
 *
 *     This subroutine defines the jacobian matrices of eighteen 
 *     nonlinear least squares problems.  The problem dimensions are 
 *     as described in the prologue comments of ssqfcn. 
 *
 *     The subroutine statement is 
 *
 *       subroutine ssqjac(m,n,x,fjac,ldfjac,nprob) 
 *
 *     where 
 *
 *       m and n are positive integer input variables.  n must not 
 *         exceed m. 
 *
 *       x is an input array of length n.
 *
 *       fjac is an m by n output array which contains the jacobian 
 *         matrix of the nprob function evaluated at x. 
 *
 *       ldfjac is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       nprob is a positive integer variable which defines the 
 *         number of the problem.  nprob must not exceed 18. 
 * 
 *     MINPACK. version of july 1978. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
 */

int ssqjac(int m, int n, const double *x, double *fjac, int ldfjac, int nprob)
{
	/* Initialized data */
	double v[11] = { 4., 2., 1., .5, .25, .167, .125, .1, .0833, .0714, .0625 };

	/* Parameter adjustments */
	--x;
	fjac -= 1 + ldfjac;
	double a, b, c, d, dx, prod;
	
	switch (nprob) {
	case 1:
		/* LINEAR FUNCTION - FULL RANK. */
		a = 2. / m;
		for (int j = 1; j <= n; ++j) {
			for (int i = 1; i <= m; ++i)
				fjac[i + j * ldfjac] = -a;
			fjac[j + j * ldfjac] += 1;
		}
		break;

	case 2:
		/* LINEAR FUNCTION - RANK 1. */
		for (int j = 1; j <= n; ++j)
			for (int i = 1; i <= m; ++i)
				fjac[i + j * ldfjac] = i * j;
		break;

	case 3:
		/* LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS. */
		for (int j = 1; j <= n; ++j)
			for (int i = 1; i <= m; ++i)
				fjac[i + j * ldfjac] = 0;
		for (int j = 2; j <= n - 1; ++j)
			for (int i = 2; i <= m - 1; ++i)
				fjac[i + j * ldfjac] = (i - 1) * j;
		break;

	case 4:
		/* ROSENBROCK FUNCTION. */
		fjac[ldfjac + 1] = -20 * x[1];
		fjac[2 * ldfjac + 1] = 10;
		fjac[ldfjac + 2] = -1;
		fjac[2 * ldfjac + 2] = 0;
		break;

	case 5:
		/* HELICAL VALLEY FUNCTION. */
		b = x[1] * x[1] + x[2] * x[2];
		c = 2 * M_PI * b;
		d = sqrt(b);
		fjac[ldfjac + 1] = 100 * x[2] / c;
		fjac[ldfjac * 2 + 1] = -100 * x[1] / c;
		fjac[ldfjac * 3 + 1] = 10;
		fjac[ldfjac + 2] = 10 * x[1] / d;
		fjac[ldfjac * 2 + 2] = 10 * x[2] / d;
		fjac[ldfjac * 3 + 2] = 0;
		fjac[ldfjac + 3] = 0;
		fjac[ldfjac * 2 + 3] = 0;
		fjac[ldfjac * 3 + 3] = 1;
		break;

	case 6:
		/* POWELL SINGULAR FUNCTION. */
		for (int j = 1; j <= 4; ++j)
			for (int i = 1; i <= 4; ++i)
				fjac[i + j * ldfjac] = 0;
		fjac[ldfjac + 1] = 1;
		fjac[ldfjac * 2 + 1] = 10;
		fjac[ldfjac * 3 + 2] = sqrt(5);
		fjac[ldfjac * 4 + 2] = -fjac[ldfjac * 3 + 2];
		fjac[ldfjac * 2 + 3] = 2 * (x[2] - 2 * x[3]);
		fjac[ldfjac * 3 + 3] = -2 * fjac[ldfjac * 2 + 3];
		fjac[ldfjac + 4] = 2 * sqrt(10) * (x[1] - x[4]);
		fjac[ldfjac * 4 + 4] = -fjac[ldfjac + 4];
		break;

	case 7:
		/* FREUDENSTEIN AND ROTH FUNCTION. */
		fjac[ldfjac + 1] = 1;
		fjac[ldfjac * 2 + 1] = x[2] * (10 - 3 * x[2]) - 2;
		fjac[ldfjac + 2] = 1;
		fjac[ldfjac * 2 + 2] = x[2] * (2 + 3 * x[2]) - 14;
		break;

	case 8:
		/* BARD FUNCTION. */
		for (int i = 1; i <= 15; ++i) {
			double tmp1 = i;
			double tmp2 = 16 - i;
			double tmp3 = tmp1;
			if (i > 8)
				tmp3 = tmp2;
			/* Computing 2nd power */
			double sq = x[2] * tmp2 + x[3] * tmp3;
			double tmp4 = sq * sq;
			fjac[i + ldfjac] = -1;
			fjac[i + ldfjac * 2] = tmp1 * tmp2 / tmp4;
			fjac[i + ldfjac * 3] = tmp1 * tmp3 / tmp4;
		}
		break;

	case 9:
		/* KOWALIK AND OSBORNE FUNCTION. */
		for (int i = 1; i <= 11; ++i) {
			double tmp1 = v[i - 1] * (v[i - 1] + x[2]);
			double tmp2 = v[i - 1] * (v[i - 1] + x[3]) + x[4];
			fjac[i + ldfjac] = -tmp1 / tmp2;
			fjac[i + ldfjac * 2] = -v[i - 1] * x[1] / tmp2;
			fjac[i + ldfjac * 3] = fjac[i + ldfjac] * fjac[i + ldfjac * 2];
			fjac[i + ldfjac * 4] = fjac[i + ldfjac * 3] / v[i - 1];
		}
		break;

	case 10:
		/* MEYER FUNCTION. */
		for (int i = 1; i <= 16; ++i) {
			double temp = 5 * i + 45 + x[3];
			double tmp1 = x[2] / temp;
			double tmp2 = exp(tmp1);
			fjac[i + ldfjac] = tmp2;
			fjac[i + ldfjac * 2] = x[1] * tmp2 / temp;
			fjac[i + ldfjac * 3] = -tmp1 * fjac[i + ldfjac * 2];
		}
		break;

	case 11:
		/* WATSON FUNCTION. */
		for (int i = 1; i <= 29; ++i) {
			double div = i / 29.;
			double s2 = 0;
			double dx = 1;
			for (int j = 1; j <= n; ++j) {
				s2 += dx * x[j];
				dx = div * dx;
			}
			double temp = 2 * div * s2;
			dx = 1 / div;
			for (int j = 1; j <= n; ++j) {
				fjac[i + j * ldfjac] = dx * (j - 1 - temp);
				dx = div * dx;
			}
		}
		for (int j = 1; j <= n; ++j)
			for (int i = 30; i <= 31; ++i)
				fjac[i + j * ldfjac] = 0;
		fjac[ldfjac + 30] = 1;
		fjac[ldfjac + 31] = -2 * x[1];
		fjac[ldfjac * 2 + 31] = 1;
		break;

	case 12:
		/* BOX 3-DIMENSIONAL FUNCTION. */
		for (int i = 1; i <= m; ++i) {
			double tmp1 = i / 10.;
			fjac[i + ldfjac] = -tmp1 * exp(-tmp1 * x[1]);
			fjac[i + ldfjac * 2] = tmp1 * exp(-tmp1 * x[2]);
			fjac[i + ldfjac * 3] = exp(-i) - exp(-tmp1);
		}
		break;

	case 13:
		/* JENNRICH AND SAMPSON FUNCTION. */
		for (int i = 1; i <= m; ++i) {
			fjac[i + ldfjac] = -i * exp(i * x[1]);
			fjac[i + ldfjac * 2] = -i * exp(i * x[2]);
		}
		break;

	case 14:
		/* BROWN AND DENNIS FUNCTION. */
		for (int i = 1; i <= m; ++i) {
			double temp = i / 5.;
			double ti = sin(temp);
			double tmp1 = x[1] + temp * x[2] - exp(temp);
			double tmp2 = x[3] + ti * x[4] - cos(temp);
			fjac[i + ldfjac] = 2 * tmp1;
			fjac[i + ldfjac * 2] = temp * fjac[i + ldfjac];
			fjac[i + ldfjac * 3] = 2 * tmp2;
			fjac[i + ldfjac * 4] = ti * fjac[i + ldfjac * 3];
		}
		break;

	case 15:
		/* CHEBYQUAD FUNCTION. */
		dx = 1. / n;
		for (int j = 1; j <= n; ++j) {
			double tmp1 = 1;
			double tmp2 = 2 * x[j] - 1;
			double temp = 2 * tmp2;
			double tmp3 = 0;
			double tmp4 = 2;
			for (int i = 1; i <= m; ++i) {
				fjac[i + j * ldfjac] = dx * tmp4;
				double foo = 4 * tmp2 + temp * tmp4 - tmp3;
				tmp3 = tmp4;
				tmp4 = foo;
				double bar = temp * tmp2 - tmp1;
				tmp1 = tmp2;
				tmp2 = bar;
			}
		}
		break;

	case 16:
		/* BROWN ALMOST-LINEAR FUNCTION. */
		prod = 1;
		for (int j = 1; j <= n; ++j) {
			prod = x[j] * prod;
			for (int i = 1; i <= n; ++i)
				fjac[i + j * ldfjac] = 1;
			fjac[j + j * ldfjac] = 2;
		}
		for (int j = 1; j <= n; ++j) {
			double temp = x[j];
			if (temp == 0) {
				temp = 1;
				prod = 2;
				for (int k = 1; k <= n; ++k)
					if (k != j)
						prod = x[k] * prod;
			}
			fjac[n + j * ldfjac] = prod / temp;
		}
		break;

	case 17:
		/* OSBORNE 1 FUNCTION. */
		for (int i = 1; i <= 33; ++i) {
			double temp = 10 * (i - 1);
			double tmp1 = exp(-x[4] * temp);
			double tmp2 = exp(-x[5] * temp);
			fjac[i + ldfjac] = -1;
			fjac[i + ldfjac * 2] = -tmp1;
			fjac[i + ldfjac * 3] = -tmp2;
			fjac[i + ldfjac * 4] = temp * x[2] * tmp1;
			fjac[i + ldfjac * 5] = temp * x[3] * tmp2;
		}
		break;

	case 18:
		/* OSBORNE 2 FUNCTION. */
		for (int i = 1; i <= 65; ++i) {
			double temp = (i - 1) / 10.;
			double tmp1 = exp(-x[5] * temp);
			double s1 = temp - x[9];
			double tmp2 = exp(-x[6] * (s1 * s1));
			double s2 = temp - x[10];
			double tmp3 = exp(-x[7] * (s2 * s2));
			double s3 = temp - x[11];
			double tmp4 = exp(-x[8] * (s3 * s3));
			fjac[i + ldfjac] = -tmp1;
			fjac[i + ldfjac * 2] = -tmp2;
			fjac[i + ldfjac * 3] = -tmp3;
			fjac[i + ldfjac * 4] = -tmp4;
			fjac[i + ldfjac * 5] = temp * x[1] * tmp1;
			double s4 = temp - x[9];
			fjac[i + ldfjac * 6] = x[2] * (s4 * s4) * tmp2;
			double s5 = temp - x[10];
			fjac[i + ldfjac * 7] = x[3] * (s5 * s5) * tmp3;
			double s6 = temp - x[11];
			fjac[i + ldfjac * 8] = x[4] * (s6 * s6) * tmp4;
			fjac[i + ldfjac * 9] = -2 * x[2] * x[6] * (temp - x[9]) * tmp2;
			fjac[i + ldfjac * 10] = -2 * x[3] * x[7] * (temp - x[10]) * tmp3;
			fjac[i + ldfjac * 11] = -2 * x[4] * x[8] * (temp - x[11]) * tmp4;
		}
		break;
	}
	return 0;
}
