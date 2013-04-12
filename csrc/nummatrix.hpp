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
  Matrix(T *, int & , int & );
  ~Matrix();
  void initializeToZero();
  bool isSquare();
  friend void matrixMultiply(Matrix<T> *, Matrix<T> *, Matrix<T> *);
};

template<typename T>
Matrix<T>::Matrix():m(1), n(1){
  matrix = new T[1];
  matrix[0] = 0.0;
}

template<typename T>
Matrix<T>::Matrix(T* A, int& m0, int&n0):m(m0), n(n0){
  if (0 == m * n)
    std::cerr << "Impossible to have a dimension zero" << std::endl;
  matrix = new T[m * n];
}

template<typename T>
void Matrix<T>::initializeToZero(){
  for (int i = 0; i < m * n ; i++)
    matrix[i] = 0.0;
}

template<typename T>
bool Matrix<T>::isSquare(){
  (m == n) ? return true : return false;
}

template<typename T>
void matrixMultiply(Matrix<T>* A, Matrix<T> * B, Matrix<T> * C){
  /*  Multiplies A * B = C */
  
  // Perform some basic checks
  if(A->n != B->m)
    std::cerr << "Matrix Dimensions for the multiplicating matrix  do not agree" <<
      std::endl;
  if(A->m != C->m)
    std::cerr << "Size m of the result matrix is wrong" << std::endl;
  if(B->n != C->n)
    std::cerr << "Size n of the result matrix is wrong" << std::endl;
  
  

}

#endif // _NUMMATRIX_HPP_
