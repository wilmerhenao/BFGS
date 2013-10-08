//#define NDEBUG //Uncomment if you are not debugging (FASTER)
/*
This is where I test my original matrix operations library
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "nummatrix.hpp"
#include <qd/dd_real.h>
/* ----- local macros  -------------------------------------------------- */

//Index must be defined in column major orderA
#define INDEXCM(i, j)     field_start + i * m + j
#define G(i, j)           G[INDEXCM(i, j)]
#define Q(i, j)           Q[INDEXCM(i, j)]
#define QD(i, j)          QD[INDEXCM(i, j)]
#define QL(i, j)          QL[INDEXCM(i, j)]
#define MAX(A,B)          ((A) > (B)) ? (A) : (B)
#define MIN(A,B)          ((A) > (B)) ? (B) : (A)
#define ABS(A)            ((A) >= 0) ? (A) : -(A)

int main(){
  std::cout << "Matrix tester " << std::endl;
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  int n = 4;
  double * mm = new double[n];

  // Creation of the matrix
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      mm[i + n * j] = i;
    }
  }

  Matrix<double> mdouble;
  fpu_fix_end(&old_cw);
  return 0;
}
