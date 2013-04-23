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

#include <iostream>
#include "../lib/qpspecial/lapackc.hpp"

// multiplication routine in lapack
/*extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);

extern "C" void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
*/
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
  Matrix(int, int); // Constructor without data
  Matrix(Matrix<T>&); // copy constructor
  ~Matrix();
  void initializeToZero();
  bool isSquare();
  template<typename H> friend void matrixMultiply(Matrix<H> &, Matrix<H> &, 
						  Matrix<H> &, char transA = 'N',
						  char transB = 'N');
  template<typename H> friend double squareForm(Matrix<H> &, Matrix<H> &, 
						Matrix<H>&);
  template<typename H> friend void GensquareForm(Matrix<H> &, Matrix<H> &, 
						 Matrix<H>&, Matrix<H> &);

  template<typename H> friend void bfgssolver(Matrix<H>&, Matrix<H>&, Matrix<H>&);

  T& operator()(int& );
  //T& operator()(int);
  Matrix<T>& operator=(Matrix<T>& );
  Matrix<T>& operator*=(double );
};

template<typename T>
T& Matrix<T>::operator()(int& i){
  return matrix[i];
}
/*
template<typename T>
T& Matrix<T>::operator()(int i){
  return matrix[i];
}
*/
//copy constructor
template<typename T>
Matrix<T>::Matrix(Matrix<T>& other){
  m = other.m;
  n = other.n;
  matrix = new T[m * n];
  for(int i = 0; i < m * n; i++){
    matrix[i] = other(i);
  }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& rhs){
  /*
    Please notice that this is not your average "equal" operator.  I am not creating any
    new memory and I am only using what I already had in the receiving Matrix<T>

    A check is performed... but you don't want this to crash while you run it so be
    careful
  */
  
  if(rhs.m != m)
    std::cerr << "This is an assignment to an element of different size" << std::endl;
  if(this != &rhs){
    m = rhs.m;
    n = rhs.n;   
    for (int i = 0; i < m * n; i++)
      matrix[i] = rhs(i);
  }
  return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(double rhs){
  for(int i = 0; i < n * m; i++){
    matrix[i] = rhs * matrix[i];
  }
  return (*this);
}

// constructors
template<typename T>
Matrix<T>::Matrix():m(1), n(1){
  matrix = new T[1];
  matrix[0] = 0.0;
}
// Initialize to zeroes
template<typename T>
Matrix<T>::Matrix(int m0, int n0):m(m0), n(n0){
  if (0 >= m || 0 >= n){
    std::cerr << "Impossible to have a dimension zero or negative" << std::endl;
    std::cerr << "m: " << m << " n:" << n << std::endl;
  }
  matrix = new T[m * n];
  for(int i = 0; i < m * n; i++){
    matrix[i] = 0;
  }
}


template<typename T>
Matrix<T>::Matrix(T* A, int m0, int n0):m(m0), n(n0){
  if (0 >= m || 0 >= n){
    std::cerr << "Impossible to have a dimension zero or negative" << std::endl;
    std::cerr << "m: " << m << " n:" << n << std::endl;
  }
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
void bfgssolver(Matrix<T>& A, Matrix<T>& B, Matrix<T>& x){
  /*
    Solves the system Ax = B
  */
  
  // Perform some basic checks
  // A has to be nxn
  if((A.m != A.n) || A.m <= 0){
    std::cerr << "Dimensions of A are wrong.  Either not a square matrix or one of" << 
      "the dimensions of the matrix is wrong.  Aborting " << std::endl;
    std::cerr << "Dimensions of A are: " << A.m << " by " << A.n << std::endl; 
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
  int* IPIV = new int[A.m]; // Not initialize it??? that is the question
  int info;
  
  // Use a copy of A0 that will be destroyed (do it in only one loop)
  for(int i = 0; i < A.m * A.n; i++)
    A0[i] = A(i);
  // Do not actually use matrix b since it will be destroyed.  Use "x" instead
  x = B;
    
  dgesv_(&N, &NRHS, A0, &N, IPIV, x.matrix, &N, &info);
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
    if(A.m != C.m){
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
      std::cerr << "A.m: " << A.m << " C.m: " << C.m << std::endl;
    }
    if(B.n != C.n){
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
      std::cerr << "B.n: " << B.n << " C.n: " << C.n << std::endl;
    }
  } else if('T' == transA && 'N' == transB){
    if(A.m != B.m)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.n != C.m){
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
      std::cerr << "A.n: " << A.n << " C.m: " << C.m << std::endl;
    }
    if(B.n != C.n){
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
      std::cerr << "B.n: " << B.n << " C.n: " << C.n << std::endl;
    }
  } else if('N' == transA && 'T' == transB){
    if(A.n != B.n)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.m != C.m){
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
      std::cerr << "A.m: " << A.m << " C.m: " << C.m << std::endl;
    }
    if(B.m != C.n){
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
      std::cerr << "B.m: " << B.m << " C.n: " << C.n << std::endl;
    }
  } else if('T' == transA && 'N' == transB){
    if(A.m != B.m)
      std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
	std::endl;
    if(A.n != C.m){
      std::cerr << "Size m of the result matrix is wrong" << std::endl;
      std::cerr << "A.n: " << A.n << " B.m: " << B.m << std::endl;
    }
    if(B.n != C.n){
      std::cerr << "Size n of the result matrix is wrong" << std::endl;
      std::cerr << "B.n: " << B.n << " C.n: " << C.n << std::endl;
    }
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

template<typename T>
double squareForm(Matrix<T> & A, Matrix<T> & B, Matrix<T>& C){
  /* 
     This one solves problems of the type x^TBz where x, and z are COLUMN vectors and
     B is a square matrix.  No warnings are checked here since they will show up in 
     the individual multiplications anyway
  */

  // Create a temporal instance for first multiplication result storage
  Matrix<T> temp(1, B.n);
  Matrix<T> finalresult(1, 1);
  int posit = 0;
  T squareFormResult;

  matrixMultiply(A, B, temp, 'T', 'N');
  matrixMultiply(temp, C, finalresult, 'N', 'N');
  squareFormResult = finalresult(posit);

  return squareFormResult;
}

template<typename H> 
void GensquareForm(Matrix<H> & A, Matrix<H> & B, Matrix<H>& C, Matrix<H>& res){
  /* 
     This one solves problems of the type x^TBz where x, and z are COLUMN vectors and
     B is a square matrix.  No warnings are checked here since they will show up in 
     the individual multiplications anyway
     Notice that the only difference between this method and the one above is that
     this method returns a full matrix and not just a double.
  */

  // Create a temporal instance for first multiplication result storage
  Matrix<H> temp(A.n, B.n);
  std::cout << "here bef. first mult" << std::endl;
  matrixMultiply(A, B, temp, 'T', 'N');
  std::cout << "here bef. second mult" << std::endl;
  matrixMultiply(temp, C, res, 'N', 'N');
  std::cout << "here final mult" << std::endl;
}


#endif // _NUMMATRIX_HPP_
