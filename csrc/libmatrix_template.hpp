#ifndef _LIBMATRIX_TEMPLATE_HPP_
#define _LIBMATRIX_TEMPLATE_HPP_

#include <cstdlib>
#include <cmath>

template<class T> void mxv(T y[], T [], T x[], T alpha, T beta, size_t m, size_t n);
template <class T> T vecip(T x[], T y[], size_t n);
template <class T> T vecnorm(T v[], size_t n);
template <class T> void mat_set_eye(T A[], size_t m, size_t n);
template <class T> void vpv(T y[], T x[], T a, size_t n);
template <class T> void vcopy(T y[], T x[], size_t n);
template <class T> void vcopyp(T y[], T x[], T a, size_t n);
template <class T> void vscal(T x[], T a, size_t n);
template <class T> void mat_r1update(T A[], T x[], T y[], T alpha, size_t n);

template <class T> 
T vecip(T x[], T y[], size_t n) {
  // Dot product between x and y vectors
  T ip = 0;
  for(size_t i=0; i < n; i++)
    ip += x[i] * y[i];
  return ip;
}

template <class T> 
T vecnorm(T v[], size_t n) {
  // calculates the euclidean norm of the vector
  T nv = 0.0;
  T nv_squared = 0.0;
  for(size_t i = 0; i < n; i++)
    nv_squared += v[i] * v[i];
  nv = sqrt(nv_squared);
  return nv;
}

template <class T> 
void mxv(T y[], T A[] , T x[], T alpha, T beta, size_t m, size_t n) {
  // First loop calculates p = -Hg or in this case Ax
  T * Ax = new T[m];
  for(size_t i = 0; i < m; i++) {
    T row = 0;
    for(size_t j = 0; j < n; j++) {
      row += A[i * n + j] * x[j];
    }
    Ax[i] = row;
  }
  // Second loop calculates \alpha * Hg + \beta * y
  // Beta seems to be always zero from what I've seen in the previous code
  for(size_t i = 0; i < m; i++)
    y[i] = alpha * Ax[i] + beta * y[i];
}

template <class T> 
void vpv(T y[], T x[], T a, size_t n) {
  /* add vectors (with a parameter):
     y = a*x + y
  */
  for(size_t i = 0; i < n; i++)
    y[i] = a * x[i] + y[i];
}

template <class T> 
void vcopy(T y[], T x[], size_t n) {
  /* copies x into y */
  for(size_t i = 0; i < n; i++)
    y[i] = x[i];
}

template <class T> 
void vcopyp(T y[], T x[], T a, size_t n) {
  /* copy with a parameter
     y = a*x
  */
  for(size_t i = 0; i < n; i++)
    y[i] = a * x[i];   
}

template <class T> 
void vscal(T x[], T a, size_t n) {
  /* scale x by parameter a */
  for(size_t i = 0; i < n; i++)
    x[i] = a * x[i];
}

template <class T>
void mat_r1update(T A[], T x[], T y[], T alpha, size_t n) {
  /*  performs the rank one update A = A + alpha*x*y' using BLAS_dger
      General cblas call:
      cblas_dger(CblasRowMajor, const size_t m, const size_t n,
      const double alpha, const double *X, const int incX,
      const double *Y, const int incY, double *A, const int lda);
      *
      * In the case of square matrices nxn, this becomes:
      */
  size_t dim = n * n;
  T * xy = new T[dim];
  for(size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < n; j++)
      xy[i * n + j] = x[i] * y[j];
  }
  for(size_t i = 0; i < dim; i++)
    A[i] = A[i] + alpha * xy[i];
}

template <class T>
void mat_set_eye(T A[], size_t m, size_t n) {
  // Set eye.  But it allows for non-square matrices (Please review)
  size_t i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A[i * n + j] = 0.0;
      if (i == j)
	A[i * n + j] = 1;
    }
  }
}

#endif // _LIBMATRIX_TEMPLATE_HPP_
