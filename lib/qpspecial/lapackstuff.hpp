#ifndef _LAPACKSTUFF_HPP_
#define _LAPACKSTUFF_HPP_

//Index must be defined in column major orderA
#define INDEXCM(i, j)     field_start + i * m + j
#define G(i, j)           G[INDEXCM(i, j)]
#define Q(i, j)           Q[INDEXCM(i, j)]
#define QD(i, j)          QD[INDEXCM(i, j)]
#define QL(i, j)          QL[INDEXCM(i, j)]
#define MAX(A,B)          ((A) > (B)) ? (A) : (B)
#define MIN(A,B)          ((A) > (B)) ? (B) : (A)
#define ABS(A)            ((A) >= 0) ? (A) : -(A)
// LAPACK function declarations
extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);
extern "C" void sgemm_(char *, char *, int*, int*,int*, float*, float*, int*, 
		       float*, int*, float*, float*, int*);
extern "C" double dlange_(char*, int*, int*, double*, int*, double* );
extern "C" float slange_(char*, int*, int*, float*, int*, float*);
extern "C" int dpotrf_(char *UPLO, int* N, double* A, int* LDA, int* INFO);
extern "C" int spotrf_(char *UPLO, int* N, float* A, int* LDA, int* INFO);
extern "C" int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" int sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);
// LAPACK function overloading
void mmul_(char transA, char transB, int M, int N,int K, double alpha, double*& A,
	   int LDA, std::vector<double>& B, int LDB, double beta, 
	   std::vector<double>& C, int LDC){
  double* pB = &*B.begin();
  double* pC = &*C.begin();
  dgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, pB, &LDB, &beta,
	 pC, &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, float alpha, float* A,
	   int LDA, std::vector<float>& B, int LDB, float beta, 
	   std::vector<float>& C, int LDC){
  sgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, & *B.begin(), &LDB, &beta, 
	 & *C.begin(), &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, double alpha, double*& A,
	   int LDA, double* B, int LDB, double beta, 
	   double* C, int LDC){
  dgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

void mmul_(char transA, char transB, int M, int N,int K, float alpha, float*& A,
	   int LDA, float* B, int LDB, float beta, 
	   float* C, int LDC){
  sgemm_(&transA, &transB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

double norm_(char* A, int* B, int* C, double*& D, int* E, double* F){
  return(dlange_(A, B, C, D, E, F));
}

float norm_(char* A, int* B, int* C, float*& D, int* E, float* F){
  return(slange_(A, B, C, D, E, F));
}

void cholesky_(char &UPLO, int* N, double*& A, int* LDA, int* INFO){
  dpotrf_(&UPLO, N, A, LDA, INFO);
  // fill the lower part with zeroes
  if ('U' == UPLO){
    for (int i = 1; i < *N; i++){
      for (int j = 0; j < i; j++)
	A[j * (*N) + i] = 0.0;
    }
  } else { // or the upper part if it's 'L'
    for (int i = 0; i < (*N - 1); i++){
      for (int j = i + 1; j < *N; j++)
	A[j * (*N) + i] = 0.0;
    }
  }
}

void cholesky_(char &UPLO, int* N, float*& A, int* LDA, int* INFO){
  spotrf_(&UPLO, N, A, LDA, INFO);
  // fill the lower part with zeroes
  if ('U' == UPLO){
    for (int i = 1; i < *N; i++){
      for (int j = 0; j < i; j++)
	A[j * (*N) + i] = 0.0;
    }
  } else { // or the upper part if it's 'L'
    for (int i = 0; i < (*N - 1); i++){
      for (int j = i + 1; j < *N; j++)
	A[j * (*N) + i] = 0.0;
    }
  }
}

void solve_(int B, int C, double*& D, int E, int* F, double* G, int H, int I){
  double* tD = new double[B * B];
  for(int i = 0; i < B; i++){
    for(int j = 0; j < B; j++){
      tD[j * B + i] = D[j * B + i];
    }
  }
  dgesv_(&B, &C, tD, &E, F, G, &H, &I);
  delete [] tD;
}

void solve_(int B, int C, float*& D, int E, int* F, float* G, int H, int I){
  float* tD = new float[B * B];
  for(int i = 0; i < B; i++){
    for(int j = 0; j < B; j++){
      tD[j * B + i] = D[j * B + i];
    }
  }
  sgesv_(&B, &C, tD, &E, F, G, &H, &I);
  delete [] tD;
}

#endif // _LAPACKSTUFF_HPP_
