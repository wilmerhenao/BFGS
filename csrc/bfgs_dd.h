void bfgs_dd
(dd_real x[], dd_real * fopt, const int n, const int lm, const int m,
const dd_real ftarget, const dd_real gnormtol, const int maxit,
const int J, const dd_real taux, const dd_real taud,
const int echo, void(*testFunction_dd)(dd_real*, dd_real*, dd_real*, int), const char * datafilename, double[]);
