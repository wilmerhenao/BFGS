#include <cmath>

void absfun(double *f, double *g, double *x, int n)
{

    *f = fabs(*f);
    if (*f >= 0) {
	*g = 1;
    }

    else {
	*g = -1;
    }

}
