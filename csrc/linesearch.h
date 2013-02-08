double linesearch_ww 
(double x[], double *f, double g[], double d[], double C1, double C2,
int n, void(*testFunction)(double*, double*, double*, int), int*nfeval, double ftarget,int*exitflag);
