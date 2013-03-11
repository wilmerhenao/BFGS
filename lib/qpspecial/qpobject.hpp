
#ifndef _QPOBJECT_HPP_
#define _QPOBJECT_HPP_

#include "lapackstuff.hpp"
/*
  Notice that these functions solve the problem by using brute force.  I could have used
  one of the many lapack wrappers, but my idea here, was to go through thru the pain of
  solving a problem with basic use of LAPACK
*/
template<typename T>
void deepCopyvec(T* A, T* B, int n){
  //Makes a deep copy from vector A to vector B.  Assumes memory has been allocated
  for (int i = 0; i < n; i++){
    B[i] = A[i];
  }
}

// T can only be float or less than double though.  No higher prec. allowed.
template<typename T>
class qpclass{
private:
  int m, n, maxit, field_start;
  T eta, delta, mu0, tolmu, tolrs, kmu, nQ, krs, ap, ad, y, q;
  // G points to the matrix G that will be passed to us. e is a vector of ones.
  int info;
  T *G, *Q, *QD, *x, *e, *z, *d, *QL;
public:
  qpclass(int, int, T*&, int);
  ~qpclass();
  T dotprod(T* a, T* b);
  void copyQD();
  void fillGfromPointer(T*);
  void HessianfromG();
  void optimization();
  void postHessian();
  int iterativePhase();
  void printSolution();
  T* fetchSolution(T*);
};

template<typename T>
qpclass<T>::qpclass(int m0, int n0, T*& G0, int maxit0){
  m = m0;
  n = n0;
  maxit = maxit0;
  field_start = 0;
  eta = 0.9995;
  delta = 3.0;
  mu0 = 0.0;
  tolmu = 1e-5;
  tolrs = 1e-5;
  kmu = 0.0;
  nQ = 0.0;
  krs = 0.0;
  ap = 0.0;
  ad = 0.0;
  y = 0.0;
  q = 0.0;
  info = 11;
  assert(m * n);
  if(0 == m * n){
    std::cerr << "The matrix is empty"<< std::endl;
    exit(EXIT_FAILURE);
  }
  x = new T[n]; e = new T[n]; z = new T[n]; d = new T[n];

  for(int i = 0; i < n; i++)
    e[i] = z[i] = d[i] = x[i] = 1.0;

  fillGfromPointer(G0);
  Q = new T[n * n]; // initialize Q (needs to be zero at the beginning)
  QD = new T[n * n];
  QL = new T[n * n];
  for (int i = 0; i < n; i++){
    for(int j = i; j < n; j++){
      Q(i, j) = Q(j, i) = 0.0; //saving a pass
    }
  }
  mu0 = dotprod(x, z);
  mu0 /= T(n);
  kmu = tolmu * mu0;
  m = m0;
}

template<typename T>
qpclass<T>::~qpclass(){ 
  delete [] Q;
  delete [] QD;
  delete [] x;
  delete [] e;
  delete [] z;
  delete [] d;
}

template<typename T>
T qpclass<T>::dotprod(T* a, T* b){
  T res = 0.0;
  for(int i = 0; i < n; ++i){
    res += a[i] * b[i];
  }
  return(res);
}

template<typename T>
void qpclass<T>::fillGfromPointer(T *G0){
  G = G0;
}

template<typename T>
void qpclass<T>::HessianfromG(){
  char transA = 'T';
  char transB = 'N';
  T alpha = 1.0, beta = 0.0;
  mmul_(transA, transB, n, n, m, alpha, G, m, G, m, beta, Q, n);
}

template<typename T>
void qpclass<T>::postHessian(){
  // Norm of the Hessian
  char tnorm = 'I';
  T *WORK = new T[n];
  nQ = norm_(&tnorm, &n, &n, Q, &n, WORK) + 2;
  krs = tolrs * nQ;
}

template<typename T>
void qpclass<T>::copyQD(){
  for (int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      QL(i, j) = QD(i, j) = Q(i, j);
    }
  }
}

template<typename T>
int qpclass<T>::iterativePhase(){
  T rs, mu, r2, M, r5, r6, dy, sumx, alpha, beta, muaff, sig;
  //parameters for LAPACK functions
  char nTrans = 'N';
  char Uuplo, Luplo;
  int one, k, i, retval;
  
  //parameters for interior point method
  T* temp = new T[n];
  T *r1 = new T[n], *r3 = new T[n], *KT = new T[n];
  T *zdx = new T[n], *p = new T[n];

  T* r4 = new T[n];
  T* dx = new T[n]; //r1 just for the size of the vector, not for the values
  T* dz = new T[n];
  T* r7 = new T[n];
  
  for(k = 0; k < maxit; k++){
    rs = 0.0; mu = 0.0; r2 = -1.0; M = 0.0; sumx = 0.0;
    alpha = 1.0; beta = 0.0; muaff = 0.0; sig = 0.0;
    one = 1;
    mmul_(nTrans, nTrans, n, one, n, alpha, Q, n, x, n, beta, temp, n);
    for(i = 0; i < n; i++){
      r1[i] = -temp[i] + e[i] * y + z[i]; //residual
      rs = MAX(std::abs(r1[i]), rs);
      r2 += x[i];             //residual
      r3[i] = -(x[i] * z[i]); //slacks
      mu -= r3[i];            //current mu
    }
    rs = MAX(r2, rs); mu /= n; // residual norm
    
    /* Stopping if mu and the residual norm is small enough.  All went well */
    if(mu < kmu){
      if(rs < krs){
	retval = 0;
	break;
      }
    }
    
    copyQD();
    for(i = 0; i < n; i++){
      zdx[i] = z[i] / x[i]; // perturbation of the diagonal
      QL(i, i) = QD(i, i) += zdx[i];
    }
    
    // Perform Cholesky factorization on QD
    Uuplo = 'U';
    cholesky_(Uuplo, &n, QD, &n, &info);
    if (0 != info)
      return(2);
    Luplo = 'L';
    cholesky_(Luplo, &n, QL, &n, &info);
    if (0 != info)
      return(2);
    // Solve a system with QD and e (Notice that the solution will be stored in KT)
    int* ipiv = new int[n];
    deepCopyvec(e, KT, n); //If I solve directly over e in the next step 'e' will be 
                           //overwritten!
    solve_(n, one, QL, n, ipiv, KT, n, info);
    M = dotprod(KT, KT); // might need to make KT members of the class later?
    /* Compute approximate tangent direction using factorization from above */
    deepCopyvec(r1, r4, n); deepCopyvec(r1, dx, n); deepCopyvec(r1, dz, n);
    
    for(i = 0; i < n; i++)
      r4[i] += r3[i] / x[i];
    
    // r7 needs to be started here because r4 will be destroyed
    // next r4 keeps a temporary solution to a system.  So the original r4 ist kaput
    // This is a really bad practice but saves me a lot of memory
    deepCopyvec(r4, r7, n);
    solve_(n, one, QL, n, ipiv, r4, n, info);
    r5 = dotprod(KT, r4);
    r6 = r2 + r5;
    dy = -r6 / M;
    
    for(i = 0; i < n; i++)
      r7[i] += e[i] * dy;
    
    solve_(n, one, QL, n, ipiv, r7, n, info);
    solve_(n, one, QD, n, ipiv, r7, n, info);
    deepCopyvec(r7, dx, n);
    for(i = 0; i < n; i++)
      dz[i] = (r3[i] - z[i] * dx[i])/x[i];
    /*
      Determine maximal step possible in the approx. tangent direction here primal step
      size
    */
    ap = 1.0; ad = 1.0;
    for(i = 0; i < n; i++){
      p[i] = -x[i] / dx[i];
      if(p[i] > 0 )
	ap = MIN(p[i], ap);
      p[i] = -z[i] / dz[i];     /* Dual step size */
      if(p[i] > 0)  //Using different step sizes in primal and dual improves pfmnce a
	ad = MIN(p[i], ad); //bit
    }
    /* Heuristic for the centering paramater */

    for(i = 0; i < n; i++){
      muaff += (x[i] + (ap * dx[i])) * (z[i] + (ad * dz[i]));
    }
    muaff /= (T) n;
    /*
      Compute the new corrected search direction that now includes the appropriate
      amount of centering and mehrotras second order correction term (see r3).  We
      of course reuse the factorization from above
    */
    sig = std::pow((muaff / mu), delta);
    for(i = 0; i < n; i++){
      r3[i] += sig * mu - dx[i] * dz[i];
      r4[i] = r1[i] + r3[i] / x[i];
    }
    deepCopyvec(r4, r7, n);
    solve_(n, one, QL, n, ipiv, r4, n, info);
    r5 = dotprod(KT, r4);
    r6 = r2 + r5;
    dy = -r6 / M;

    for(i = 0; i < n; i++)
      r7[i] += e[i] * dy; //(r7 = r4 + e*dy) -- It is a little cryptic over here!
    solve_(n, one, QL, n, ipiv, r7, n, info);
    solve_(n, one, QD, n, ipiv, r7, n, info);
    deepCopyvec(r7, dx, n);
    for(i = 0; i < n; i++)
      dz[i] = (r3[i] - z[i] * dx[i])/x[i];
    /* Determine maximal step possible in the new direction here primal step size */
    ap = 1.0; ad = 1.0;
    for(i = 0; i < n; i++){
      p[i] = -(x[i]) / dx[i];
      if(p[i] > 0 )
	ap = MIN(p[i], ap);
      p[i] = -(z[i]) / (dz[i]);     /* Dual step size */
      if(p[i] > 0)  //Using different step sizes in primal and dual improves pfmnce a
	ad = MIN(p[i], ad); //bit
    }
    /* Update variables primal dual multipliers dual slacks */

    for(i = 0; i < n; i++){
      x[i] += eta * ap * dx[i];    
      z[i] += eta * ad * dz[i];
    }

    y += eta * ad * dy;
  }
  if (maxit == k){
    std::cout << "\nmax. num of iterations reached" << std::endl;
    retval = 1;
  }
  sumx = 0.0;
  for(i = 0; i < n; i++){
    x[i] = MAX(x[i], 0.0);
    sumx += x[i];
  }
  for(i = 0; i < n; i++)
    x[i] /= sumx;
  deepCopyvec(x, d, n);
  mmul_(nTrans, nTrans, m, one, n, alpha, G, n, x, n, beta, d, n);
  q = dotprod(d, d);
  return(retval);
}

template<typename T>
void qpclass<T>::printSolution(){
  std::cout << "Global value of the solution is " << q << std::endl;
}

template<typename T>
void qpclass<T>::optimization(){
  int ret = 0.0;
  HessianfromG();
  postHessian();
  ret = iterativePhase();
  switch(ret){
  case 0:
    std::cout << "Solution converged" << std::endl;
    break;
  case 1:
    std::cout << "Solution didn't converge to a global. More iterations are needed" <<
      std::endl;
    break;
  case 2:
    std::cerr << "A good solution was not found" << std::endl;
    break;
  default:
    std::cerr << "Output is undefined" << std::endl;
    break;
  }
  std::cout << ret << std::endl;
  printSolution();
}

template<typename T>
T* qpclass<T>::fetchSolution(T* ExternSol){
  // Assumes that the user supplies an empty vector with fully allocated memory
  for(int i = 0; i < n; i++)
    ExternSol[i] = x[i];
  return(ExternSol);
}

#endif // _QPOBJECT_HPP_
