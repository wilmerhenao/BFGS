#ifndef _LIBMATRIX_TEMPLATE_HPP_
#define _LIBMATRIX_TEMPLATE_HPP_

#include <cstdlib>
#include <cmath>

template <class T> void mxv(T y[], T A[], T x[], T alpha, T beta, int m, int n);
template <class T> T vecip(T x[], T y[], int n);
template <class T> T vecnorm(T v[], int n);
template <class T> void mat_set_eye(T A[], int m, int n);
template <class T> void vpv(T y[], T x[], T a, int n);
template <class T> void vcopy(T y[], T x[], int n);
template <class T> void vcopyp(T y[], T x[], T a, int n);
template <class T> void vscal(T x[], T a, int n);
template <class T> void mat_r1update(T A[], T x[], T y[], T alpha, int n);

// Template implementations
template <class T> 
T vecip(T x[], T y[], int n) {
    T ip = 0;
    for(int i=0; i<n; i++) {
        ip += x[i] * y[i];
    }
    return ip;
}

template <class T> 
T vecnorm(T v[], int n) {
    T nv = 0;
    T nv_squared = 0;
    for(int i=0; i<n; i++) {
        nv_squared += v[i] * v[i];
    }
    nv = sqrt(nv_squared);
    return nv;
}


template <class T> 
void mxv(T y[], T A[], T x[], T alpha, T beta, int m, int n) {

    T Ax[m];
    for(int i=0; i<m; i++) {
	T row = 0;
        for(int j=0; j<n; j++) {
	    row += A[i*n+j]*x[j];
        }
	Ax[i] = row;
    }

    for(int i=0; i<m; i++) {
	y[i] = alpha*Ax[i] + beta*y[i];
    }

}

template <class T> 
void vpv(T y[], T x[], T a, int n) {
 /* add vectors (with a parameter):
  *  y = a*x + y
  */

    for(int i=0; i<n; i++) {
	y[i] = a*x[i] + y[i];
    }

}

template <class T> 
void vcopy(T y[], T x[], int n) {
    /* copies x into y */

    for(int i=0; i<n; i++) {
	y[i] = x[i];
    }

}

template <class T> 
void vcopyp(T y[], T x[], T a, int n) {
    /* copy with a parameter
     * y = a*x
     */

    for(int i=0; i<n; i++) {
	y[i] = a*x[i];
    }
   
}

template <class T> 
void vscal(T x[], T a, int n) {
    /* scale x by parameter a */

    for(int i=0; i<n; i++) {
	x[i] = a*x[i];
    }

}

template <class T>
void mat_r1update(T A[], T x[], T y[], T alpha, int n) {
 /*  performs the rank one update A = A + alpha*x*y' using BLAS_dger
     General cblas call:
  cblas_dger(CblasRowMajor, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
  *
  * In the case of square matrices nxn, this becomes:
  */

    int dim = n*n;
    T xy[dim];
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
	    xy[i*n+j] = x[i]*y[j];
        }
    }

    for(int i=0; i<dim; i++) {
	A[i] = A[i] + alpha*xy[i];
    }

}

template <class T>
void mat_set_eye(T A[], int m, int n) {
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j] = 0.0;
            if (i==j) {
                A[i*n+j] = 1; /* if on the diagonal, set = 1*/
            }
        }
    }
}

#endif // _LIBMATRIX_TEMPLATE_HPP_
