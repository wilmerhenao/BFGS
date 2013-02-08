#include <qd/dd_real.h>

dd_real linesearch_ww 
(dd_real x[], dd_real *f, dd_real g[], dd_real d[], dd_real C1, dd_real C2,
int n, void(*testFunction)(dd_real*, dd_real*, dd_real*, int), int*nfeval, dd_real ftarget,int*exitflag);
