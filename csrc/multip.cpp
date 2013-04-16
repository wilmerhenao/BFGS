#include "nummatrix.hpp"
#include <iostream>

/* 
   This is just to test that my nummatrix.hpp is working correctly.  It contains a set
   of basic cases where everything should run smooth
*/

int main (){
  double * a;
  a = new double[6];
  for (int i = 0; i < 6; i++){
    a[i] = i;
  } 

  double * b;
  b = new double[3];
  for (int i = 0; i < 3; i++){
    b[i] = i;
  } 

  double * b2;
  b2 = new double[2];
  for (int i = 0; i < 2; i++){
    b2[i] = i + 1;
  } 

  double * b4;
  b4 = new double[3];
  for (int i = 0; i < 3; i++){
    b4[i] = i + 1;
  } 
 
  double * c = new double[2];
  double * c2 = new double[3];

  Matrix<double> A(a, 2, 3);
  for (int i = 0; i < 6; i++)
    std::cout << A(i) << std::endl;
  Matrix<double> B(b, 3, 1);
  Matrix<double> C(c, 2, 1);
  matrixMultiply(A, B, C, 'N', 'N');
  for (int i = 0; i < 2; i++)
    std::cout << C(i) << std::endl;
  std::cout << "-----------" << std::endl;

  std::cout << "B2: " << std::endl;
  Matrix<double> B2(b2, 2, 1);
  for (int i = 0; i < 2; i++)
    std::cout << B2(i) << std::endl;

  Matrix<double> C2(c2, 3, 1);
  std::cout << "A2: " << std::endl;
  for(int i = 0; i < 6; i++)
    std::cout << A(i) << std::endl;
  
  matrixMultiply(A, B2, C2, 'T', 'N');
  std::cout << "C2: " << std::endl;
  for (int i = 0; i < 3; i++)
    std::cout << C2(i) << std::endl;
  
  Matrix<double> B3(b2, 1, 2);
  matrixMultiply(A, B3, C2, 'T', 'T');
  for (int i = 0; i < 3; i++)
    std::cout << C2(i) << std::endl;
  
  std::cout << "ultimirris: " << std::endl;
  
  Matrix<double> B4(b4, 1, 3);
  matrixMultiply(A, B4, B2, 'N', 'T');
  for (int i = 0; i < 2; i++)
    std::cout << B2(i) << std::endl;

  std::cout << "Test the solution of a system with sgesv" << std::endl;
  double * amatrix = new double[100];
  for(int i = 0; i < 100; i++){
    amatrix[i] = i;
  }
  Matrix<double> AEq(amatrix, 10, 10);
  double * bmatrix = new double[10];
  for (int i = 0; i < 10; i++)
    bmatrix[i] = i;
  Matrix<double> BEq(bmatrix, 10, 1);

  double * x = new double[10];
  for (int i = 0; i < 10; i++)
    x[i] = i;
  Matrix<double> X(x, 10, 1);

  bfgssolver(A, BEq, X);
  for (int i = 0; i < 10; i++)
    std::cout << BEq(i) << std::endl;
  
  return 0;
}

