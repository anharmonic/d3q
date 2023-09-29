struct test_case {
    int nprob;
    int n;
    int m;
    double factor;

    /* expected output */
    int nfev;
    int njev;
    int info;
    double fnorm2; 
    double fnorm2_lastchance; 
};

void initpt(int n, double *x, int nprob, double factor);
void ssqfcn(int m, int n, const double *x, double *fvec, int nprob);
int ssqjac(int m, int n, double *x, double *fjac, int ldfjac, int nprob);
void commentator(int ic, double *x, double *fvec, double ftol, double xtol, int nfev, int njev, int info);

extern struct test_case tests[];
extern char * problem_name[];