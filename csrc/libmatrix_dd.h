#include <qd/dd_real.h>

//void mxm(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k);
//void mtxm(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k);
//void mxmt(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k);
void mxv(dd_real y[], dd_real A[], dd_real x[], dd_real alpha, dd_real beta, int m, int n);
//int linsolve(double A[], double B[], int n, int m);
//int minv(double A[], const int n);
dd_real vecip(dd_real x[], dd_real y[], int n);
dd_real vecnorm(dd_real v[], int n);
void mat_set_eye(dd_real A[], int m, int n);
void vpv(dd_real y[], dd_real x[], dd_real a, int n);
void vcopy(dd_real y[], dd_real x[], int n);
void vcopyp(dd_real y[], dd_real x[], dd_real a, int n);
void vscal(dd_real x[], dd_real a, int n);
void mat_r1update(dd_real A[], dd_real x[], dd_real y[], dd_real alpha, int n);
//void changemajor(double A[], int m, int n);
