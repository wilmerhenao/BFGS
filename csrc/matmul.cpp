#include <cstdio>
#include <iostream>
#include "nummatrix.hpp"
#include "libmatrix_template.hpp"
#include "../lib/qpspecial/lapackc.hpp"

int main(){
  
  double * Adouble = new double[9];
  double * Bdouble = new double[9];
  
  for(int i = 0; i < 9; i++){
    Adouble[i] = i;
    Bdouble[i] = 2 * i;
  }

  Matrix<double> A(Adouble, 3, 3), B(Bdouble, 3, 3), C(3, 3);
  matrixMultiplywithPadding(A, B, C, 'N', 'N', 3, 3, 3);
  
  for(int i = 0; i < 9; i++){
    std::cout << C(i) << std::endl;
  }
  
  return 0;
}
