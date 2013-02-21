#include <cmath>

void testsquared(double *f, double *g, double *x, int n)
{

double sum = x[0]+x[1]+x[2];

*f = fabs(x[0]) + 10*pow(x[1],2) + (sin(sum))*(sin(sum));

if (x[0]>0) {
    g[0] = 1 + 2*cos(sum)*x[0];
}
else {
    g[0] = -1 + 2*cos(sum)*x[0];
}

g[1] = 20*x[1] + 2*cos(sum)*x[1];

g[2] = 2*cos(sum)*x[2];

}
