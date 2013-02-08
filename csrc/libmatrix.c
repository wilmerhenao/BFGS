#include <stdlib.h>
#include <math.h>

double vecip(double x[], double y[], int n) {
    double ip = 0;
    for(int i=0; i<n; i++) {
        ip += x[i] * y[i];
    }
    return ip;
}

double vecnorm(double v[], int n) {
    double nv = 0;
    double nv_squared = 0;
    for(int i=0; i<n; i++) {
        nv_squared += v[i] * v[i];
    }
    nv = sqrt(nv_squared);
    return nv;
}

void mxv(double y[], double A[], double x[], double alpha, double beta, int m, int n) {
    
     /* computes x = b*(A*v)
      *
      * where
      *
      * A = m x n
      * v = n x 1
      * x = m x 1
      *
     */

    double Ax[m];
    for(int i=0; i<m; i++) {
	double row = 0;
        for(int j=0; j<n; j++) {
	    row += A[i*n+j]*x[j];
        }
	Ax[i] = row;
    }

    for(int i=0; i<m; i++) {
	y[i] = alpha*Ax[i] + beta*y[i];
    }

}


//void mxm(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k) {
    
     /* computes C := alpha*A*B + beta*C
     cblas call is
      *(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                 const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
      */

//      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                  m,n,k,alpha,A,k,B,n,beta,C,n);
    
//}


/* BFGS VIRKER I HVERT FALD MED DENNE HER:
void mxm(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k) {
    
     /* computes C := alpha*A*B + beta*C
     cblas call is
      *(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                 const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  m,n,k,alpha,A,m,B,k,beta,C,m);
    
}
*/

//void mxmt(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k) {
    
     /* computes C := alpha*A*B' + beta*C
      *
      * where
      *
      * A = m x k
      * B = n x k    <-- NB: SO IT ONLY MAKES SENSE IF WE COMPUTE A*B'
      * C = m x n
      *

      cblas call is
      *(CblasRowMajor, CblasNoTrans, CblasTrans, 
                 const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
      */

//      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//                 m,n,k,alpha,A,k,B,k,beta,C,n);
    
//}


//void mtxm(double C[], double A[], double B[], double alpha, double beta, int m, int n, int k)
//{
    
    
     /* computes C := alpha*A'*B + beta*C
      *
      * where
      *
      * A = k x m   <-- NB: SO IT ONLY MAKES SENSE IF WE COMPUTE A'*B
      * B = k x n    
      * C = m x n
      *

      cblas call is
      *(CblasRowMajor, CblasNoTrans, CblasTrans, 
                 const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
      
    */
    
//    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
//                  m,n,k,alpha,A,m,B,n,beta,C,n);
    
    
    /* NAIVE CODE FOR THIS ROUTINE:
    int i,j,r;
    double tmp; 
    for (i = 0; i < M; i++) {
        int iM = i*M;
        for (j = 0; j < M; j++) {
            
            C[iM + j] = beta*C[iM + j];
            tmp       = 0.0;
            for (r = 0; r < K; r++) {
                tmp += A[r*M + i]*A[r*M + j];
            }
            C[iM + j] += alpha*tmp;
            
        }
    }
     */

//}
 
   
//int linsolve(double A[], double B[], int n, int m)
//{
    
    
     /* computes B := A \ B
      *
      * where
      *
      * A = n x n   
      * B = n x m    
      * X = n x m
      *
      * NB: A is overwritten with the LU factors of A

      clapack call is
      *int clapack_dgesv
             (const enum CBLAS_ORDER Order, const int N, const int NRHS,
              double *A, const int lda, int *ipiv,
              double *B, const int ldb);
      
    */
    
//    int *ipiv = malloc(n * sizeof(int));
//    int status = clapack_dgesv(CblasRowMajor, n, m, A, n, ipiv, B, n);  
//    free(ipiv);
//    return status; 
//}



//int minv(double A[], const int n) {
    
     /* computes A := A^{-1} 
      *
      * ALWAYS check the output integer.
      * if it is > 0, matrix was singular
      * if it is < 0, something else went wrong
      *
      */
    
//    int status = 0;
//    int *ipiv   = malloc(n * sizeof(int));  
//    int status1 = clapack_dgetrf(CblasRowMajor, n, n, A, n, ipiv);
//    int status2 = clapack_dgetri(CblasRowMajor, n, A, n, ipiv);
//    if ( (status1 > 0) || (status2 > 0) ) status = 1; 
//    if ( (status1 < 0) || (status2 < 0) ) status = -1; 
//    free(ipiv);    
//    return status; 
//}


void vpv(double y[], double x[], double a, int n) {
 /* add vectors (with a parameter):
  *  y = a*x + y
  */

    for(int i=0; i<n; i++) {
	y[i] = a*x[i] + y[i];
    }

}

void vcopy(double y[], double x[], int n) {
    /* copies x into y */

    for(int i=0; i<n; i++) {
	y[i] = x[i];
    }

}

void vcopyp(double y[], double x[], double a, int n) {
    /* copy with a parameter
     * y = a*x
     */

    for(int i=0; i<n; i++) {
	y[i] = a*x[i];
    }
    
}

void vscal(double x[], double a, int n) {
    /* scale x by parameter a */

    for(int i=0; i<n; i++) {
	x[i] = a*x[i];
    }

}


void mat_r1update(double A[], double x[], double y[], double alpha, int n) {
 /*  performs the rank one update A = A + alpha*x*y' using BLAS_dger
     General cblas call:
  cblas_dger(CblasRowMajor, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
  *
  * In the case of square matrices nxn, this becomes:
  */

    int dim = n*n;
    double xy[dim];
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
	    xy[i*n+j] = x[i]*y[j];
        }
    }

    for(int i=0; i<dim; i++) {
	A[i] = A[i] + alpha*xy[i];
    }

}


void mat_set_eye(double A[], int m, int n) {
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j] = 0;
            if (i==j) {
                A[i*n+j] = 1; /* if on the diagonal, set = 1*/
            }
        }
    }
}

/*
void changemajor(double A[], int m, int n) {
    /* Change the major of the mxn array A
     * i.e. from row major to column major
     * or   from column major to row major
     
    int i,j,in; 
    double dtmp; 
    for (i=0; i<m; i++) {
        in = i*n;
        for (j=0; j<n; j++) {
            dtmp     = A[in+j];
            A[in+j]  = A[j*n+i];
            A[j*n+i] = dtmp;
        }
    }
}
*/
