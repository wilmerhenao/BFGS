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

#ifndef _NUMMATRIX_HPP_
#define _NUMMATRIX_HPP_

#include<iostream>

// multiplication routine in lapack
extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);

/*
  The main idea of creatin this class is to use a library interacts with lapack and
  seamlessly can solve problems without having to worry for lapack's messy syntax.
 */

template<typename T>
class Matrix{
  // this class defines the typical matrix of dimensions m x n
protected:  
  int m; // vertical coordinate of the matrix
  int n; // horizontal coordinate of the matrix
  T * matrix;
public:
  Matrix();
  Matrix(T *, int , int );
  ~Matrix();
  void initializeToZero();
  bool isSquare();
  template<typename H> friend void matrixMultiply(Matrix<H> &, Matrix<H> &, 
						  Matrix<H> &, char transA = 'N',
						  char transB = 'N');
  
  T& operator()(int& );
  //operator=(T&);
};

template<typename T>
T& Matrix<T>::operator()(int& i){
  return matrix[i];
}
/*
template<typename T>
Matrix<T>::operator=(T&){
matrix
}
*/
template<typename T>
Matrix<T>::Matrix():m(1), n(1){
  matrix = new T[1];
  matrix[0] = 0.0;
}

template<typename T>
Matrix<T>::Matrix(T* A, int m0, int n0):m(m0), n(n0){
  if (0 >= m || 0 >= n)
    std::cerr << "Impossible to have a dimension zero or negative" << std::endl;
  matrix = new T[m * n];
  for(int i = 0; i < m * n; i++){
    matrix[i] = A[i];
  }
}

template<typename T>
Matrix<T>::~Matrix(){
  delete [] matrix;
}

template<typename T>
void Matrix<T>::initializeToZero(){
  for (int i = 0; i < m * n ; i++)
    matrix[i] = 0.0;
}

template<typename T>
bool Matrix<T>::isSquare(){
  bool retvalue;
  (m == n) ? retvalue = true : retvalue = false;
  return retvalue;
}

template<typename T>
void matrixMultiply(Matrix<T>& A, Matrix<T>& B, Matrix<T>& C, char transA = 'N',
		    char transB = 'N'){
  /*  Multiplies A * B = C */
  
  // Perform some basic checks (dependent on transA and transB)
  
  if('N' == transA && 'N' == transB){
    if(A.n != B.m)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.m != C.m)
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
    if(B.n != C.n)
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
  } else if('T' == transA && 'N' == transB){
    if(A.m != B.m)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.n != C.m)
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
    if(B.n != C.n)
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
  } else if('N' == transA && 'T' == transB){
    if(A.n != B.n)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.m != C.m)
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
    if(B.m != C.n)
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
  } else if('T' == transA && 'N' == transB){
    if(A.m != B.m)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.n != C.m)
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
    if(B.n != C.n)
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
  }

  // declare some local variables so that external variables do not get destroyed
  // Will make local copies of A and B.  This may sound like a waste of time.  But
  // avoids the problem of unintentionally
  // UPDATE: Not necessary for this function.  A and B will be unchanged on exit
  
  int m0;
  int n0;
  int k0;
  int LDA = A.m;
  int LDB = B.m;
  int LDC = A.n;
  double alpha = 1.0;
  double beta = 0.0;

  if('N' == transA || 'n' == transA){
    m0 = A.m;
    k0 = A.n;
  }

  if('T' == transA || 't' == transA){
    m0 = A.n;
    k0 = A.m;
  }

  if('N' == transB || 'n' == transB){
    n0 = B.n;
  } else{
    n0 = B.m;
  }
  
  dgemm_(&transA, &transB, &m0, &n0, &k0, &alpha, A.matrix, &LDA, B.matrix, 
	 &LDB, &beta, C.matrix, &LDC);

}

#endif // _NUMMATRIX_HPP_
