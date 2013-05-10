/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BFGS tools Copyright (C) 2013  Wilmer Henao
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _LIBMATRIX_TEMPLATE_HPP_
#define _LIBMATRIX_TEMPLATE_HPP_

#include <cstdlib>
#include <cmath>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
/*
template<class T> void mxv(T y[], T [], T x[], T alpha, T beta, int m, int n);
template <class T> T vecip(T x[], T y[], int n);
template <class T> double veciptd(double*&, double*&, int&);
template <class T> double veciptd(dd_real*&, double*&, int&);
template <class T> double veciptd(qd_real*&, double*&, int&);
template <class T> T vecnorm(T v[], int n);
template <class T> void mat_set_eye(T A[], int m, int n);
template <class T> void vpv(T y[], T x[], T a, int n);
template <class T> void vcopy(T y[], T x[], int n);
template <class T> void vcopyp(T y[], T x[], T a, int n);
template <class T> void vscal(T x[], T a, int n);
template <class T> void mat_r1update(T A[], T x[], T y[], T alpha, int n);
									      */
template <class T> 
double veciptd(double*& x, double*& y, int n){ 
// Dot product between x and y vectors
  double ip = 0;
  for(int i = 0; i < n; i++){
    ip += x[i] * y[i];
    std::cout << i << std::endl;
  }
  return ip;  
}

template <class T> 
double veciptd(dd_real*& x, double*& y, int n){
// Dot product between x and y vectors
  double ip = 0;
  for(int i=0; i < n; i++)
    ip += (x[i].x[0]) * y[i];
  return ip;  
}

template <class T> 
double veciptd(qd_real*& x, double*& y, int n){
// Dot product between x and y vectors
  double ip = 0;
  for(int i=0; i < n; i++)
    ip += (x[i][0]) * y[i];
  return ip;  
}

template <class T> 
T vecip(T x[], T y[], int n) {
  // Dot product between x and y vectors
  T ip = 0;
  for(int i=0; i < n; i++)
    ip += x[i] * y[i];
  return ip;
}

template <class T> 
T vecnorm(T v[], int n) {
  // calculates the euclidean norm of the vector
  T nv = 0.0;
  T nv_squared = 0.0;
  for(int i = 0; i < n; i++)
    nv_squared += v[i] * v[i];
  nv = sqrt(nv_squared);
  return nv;
}

template <class T> 
void mxv(T y[], T A[] , T x[], T alpha, T beta, int m, int n) {
  // First loop calculates p = -Hg or in this case Ax
  T * Ax = new T[m];
  for(int i = 0; i < m; i++) {
    T row = 0;
    for(int j = 0; j < n; j++) {
      row += A[i * n + j] * x[j];
    }
    Ax[i] = row;
  }
  // Second loop calculates \alpha * Hg + \beta * y
  // Beta seems to be always zero from what I've seen in the previous code
  for(int i = 0; i < m; i++)
    y[i] = alpha * Ax[i] + beta * y[i];
}

template <class T> 
void vpv(T y[], T x[], T a, int n) {
  /* add vectors (with a parameter):
     y = a*x + y
  */
  for(int i = 0; i < n; i++)
    y[i] = a * x[i] + y[i];
}

template <class T> 
void vcopy(T y[], T x[], int n) {
  /* copies x into y */
  for(int i = 0; i < n; i++)
    y[i] = x[i];
}

template <class T> 
void vcopyp(T y[], T x[], T a, int n) {
  /* copy with a parameter
     y = a*x
  */
  for(int i = 0; i < n; i++)
    y[i] = a * x[i];   
}

template <class T> 
void vscal(T x[], T a, int n) {
  /* scale x by parameter a */
  for(int i = 0; i < n; i++)
    x[i] = a * x[i];
}

template <typename T>
void mat_r1update(T A[], T x[], T y[], T alpha, int n) {
  /*  performs the rank one update A = A + alpha*x*y' using BLAS_dger
      General cblas call:
      cblas_dger(CblasRowMajor, const int m, const int n,
      const double alpha, const double *X, const int incX,
      const double *Y, const int incY, double *A, const int lda);
      *
      * In the case of square matrices nxn, this becomes:
      */
  
  //int dim = n * n;
  //FLAG();
  //T * xy = new T[dim];
  
  //FLAG();
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      A[i * n + j] = A[i * n + j] + alpha * x[i] * y[j];
  }
  
  //for(int i = 0; i < dim; i++)
  //  A[i] = A[i] + alpha * xy[i];
  
  //delete [] xy;
}

template <class T>
void mat_set_eye(T A[], int m, int n) {
  // Set eye.  But it allows for non-square matrices (Please review)
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A[i * n + j] = 0.0;
      if (i == j)
	A[i * n + j] = 1;
    }
  }
}

#endif // _LIBMATRIX_TEMPLATE_HPP_
