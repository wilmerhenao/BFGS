#ifndef LAPACKC_HPP
#define LAPACKC_HPP

// LAPACK function declarations
extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);
extern "C" void sgemm_(char *, char *, int*, int*, int*, float*, float*, int*, 
		       float*, int*, float*, float*, int*);
extern "C" double dlange_(char*, int*, int*, double*, int*, double* );
extern "C" float slange_(char*, int*, int*, float*, int*, float*);
extern "C" int dpotrf_(char *UPLO, int* N, double* A, int* LDA, int* INFO);
extern "C" int spotrf_(char *UPLO, int* N, float* A, int* LDA, int* INFO);
extern "C" void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" void sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);
extern "C" void dgetri_(int*, double*, int*, int*, double*, int*, int*);

#endif // LAPACKC_HPP
