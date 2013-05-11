#include "nummatrix.hpp"
#include <iostream>
#include <cmath>

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

  std::cout << "mat results" << std::endl;
  for (int i = 0; i < 3; i++)
    std::cout << C2(i) << std::endl;
  
  matrixMultiplywithPadding(A, B3, C2, 'T', 'T', 2, 1, 3);
  std::cout << "padding results" << std::endl;
  for (int i = 0; i < 3; i++)
    std::cout << C2(i) << std::endl;

  std::cout << "ultimirris: " << std::endl;
  
  Matrix<double> B4(b4, 1, 3);
  matrixMultiply(A, B4, B2, 'N', 'T');
  for (int i = 0; i < 2; i++)
    std::cout << B2(i) << std::endl;

  std::cout << "Test the solution of a system with dgesv" << std::endl;
  double * amatrix = new double[4];
  amatrix[0] = 1;
  amatrix[1] = 5;
  amatrix[2] = 8;
  amatrix[3] = 100;
  Matrix<double> AEq(amatrix, 2, 2);
  double * bmatrix = new double[2];
  bmatrix[0] = 35;
  bmatrix[1] = 22;
  Matrix<double> BEq(bmatrix, 2, 1);

  double * x = new double[2];
  for (int i = 0; i < 2; i++)
    x[i] = i;
  Matrix<double> X(x, 2, 1);

  bfgssolver(AEq, BEq, X);
  std::cout << "A: " << std::endl;
  for (int i = 0; i < 4; i++){
    std::cout << AEq(i) << std::endl;
  }

  for (int i = 0; i < 2; i++)
    std::cout << X(i) << " ----  " << B(i) << std::endl;

  std::cout << "testing the padded matrix code" << std::endl;

  int sizebig = 3;
  double * aa, *bb;
  aa = new double[sizebig * sizebig];
  bb = new double[sizebig * sizebig];

  aa[0] = bb[0] = 0;
  aa[1] = bb[1] = 1;
  aa[sizebig] = bb[sizebig] = 2;
  aa[(sizebig + 1)] = bb[(sizebig + 1)] = 3;

  Matrix<double> Aa(aa, 2, 2);
  Matrix<double> Aab(bb, 2, 2);
  
  double * cdd = new double[sizebig * sizebig];
  Matrix<double> Ccc(cdd, 2, 2);
  
  for(int i = 0; i < sizebig * sizebig; i++)
    std::cout << Aa(i) << std::endl;

  matrixMultiplywithPadding(Aa, Aab, Ccc, 'N', 'N', sizebig, sizebig, sizebig);
  for(int i = 0; i < sizebig * sizebig ; i++)
    std::cout << "nunmber: " << i << "  " << Ccc(i) << std::endl;
  
  return 0;
}

