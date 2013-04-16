#include "nummatrix.hpp"
#include <iostream>

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

  return 0;
}

