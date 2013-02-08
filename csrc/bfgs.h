void bfgs
(double x[], double * fopt, const int n, const int lm, const int m,
const double ftarget, const double gnormtol, const int maxit,
const int J, const double taux, const double taud,
const int echo, void(*testFunction)(double*, double*, double*, int), const char * datafilename, double info[]);
