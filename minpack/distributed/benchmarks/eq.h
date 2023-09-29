void initpt(int n, double * x, int nprob, double factor);
void vecfcn(int n, const double * x, double * fvec, int nprob);
void vecjac(int n, const double * x, double *fjac, int ldfjac, int nprob);
void errjac(int n, const double * x, double *fjac, int ldfjac, int nprob);

extern char * problem_name[];

struct test_case {
	int nprob;
	int n;
	double factor;
	int nfev;
	int info;
	double fnorm2;
};

extern struct test_case tests[];