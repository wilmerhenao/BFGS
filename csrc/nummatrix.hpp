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

extern "C" int sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);

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
  Matrix(const Matrix<T>&); // copy constructor
  ~Matrix();
  void initializeToZero();
  bool isSquare();
  template<typename H> friend void matrixMultiply(Matrix<H> &, Matrix<H> &, 
						  Matrix<H> &, char transA = 'N',
						  char transB = 'N');
  template<typename H> friend void solver(Matrix<H>&, Matrix<H>&, Matrix<H>&);

  T& operator()(int& );
  Matrix<T>& operator=(const Matrix<T>& );
};

template<typename T>
T& Matrix<T>::operator()(int& i){
  return matrix[i];
}

//copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& other){
  m = other.m;
  n = other.n;
  matrix = new T[m * n];
  for(int i = 0; i < m * n; i++){
    matrix[i] = other(i);
  }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs){
  if(this != &rhs){
    int * newmatrix = new double[rhs.m * rhs * n];
    for (int i = 0; i < rhs.m * rhs.n; i++)
      newmatrix[i] = rhs(i);

    delete [] matrix;
    
    m = rhs.m;
    n = rhs.n;
    matrix = newmatrix;
  }
  return *this;
}

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
void solver(Matrix<T>& A, Matrix<T>& B, Matrix<T>& x){
  /*
    Solves the system Ax = B
  */
  
  // Perform some basic checks
  // A has to be nxn
  if((A.m != A.n) || A.m <= 0){
    std::cerr << "Dimensions of A are wrong.  Either not a square matrix or one of" << 
      "the dimensions of the matrix is wroong.  Aborting " << std::endl; 
  }
  if(A.m != B.m){
    std::cerr << "Size of matrix A and B are in disagreement" << std::endl;
  }
  if(A.n != x.m){
    std::cerr << "Size of matrix A and x are in disagreement" << std::endl;
  }
  
  // Variables to the equation
  int N = A.m;
  int NRHS = 1; //I always have to solve only one system
  double* A0 = new double[A.m * A.n];
  int* IPIV = new double[A.m]; // Not initialize it??? that is the question
  int info;
  
  // Use a copy of A0 that will be destroyed (do it in only one loop)
  for(int i = 0; i < a.m * A.n; i++)
    A0[i] = A(i);
  
  // Do not actually use matrix b since it will be destroyed.  Use "x" instead
  x = B;
    
  sgesv_(&N, &NRHS, A0, &N, IPIV, x.matrix, &N, &info);
  // A0 now contains the factors L and U from the LU factorization of A
  if(0 != info ){
    if(info < 0){
      std::cerr << "Argument " << -info << " had an illegal value" << std::endl;
    } else {
      std::cerr << "Position U(i,i) where i = " << info << " is exactly zero.  The " <<
	"factorization was completed but the U is exactly singular,  so the solution" <<
	" could not be computed" << std::endl;
    }
  }
  
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
