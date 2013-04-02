#ifndef _BFGS_TEMPLATE_HPP_
#define _BFGS_TEMPLATE_HPP_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include "linesearch_template.hpp"
#include "quasinewt_updates_template.hpp"
#include "print_template.hpp"
#include "libmatrix_template.hpp"
#include "../lib/qpspecial/qpobject.hpp"
#include <qd/dd_real.h>
#include <qd/qd_real.h>

extern "C" void dgemm_(char *, char *, int*, int*,int*, double*, double*, int*, 
		       double*, int*, double*, double*, int*);
extern "C" int sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);

template<class T>
void bfgs(T*&, T*& fopt, size_t& n, short& lm, size_t& m, T& ftarget,  T& gnormtol,  
	  size_t& maxit,  long& J, T& taux,  T& taud, short& echo, 
	  int(*&testFunction)(T*, T*, T*, size_t),  std::string& datafilename, 
          double*&, size_t&, double*&, double*&, bool&);

/* Cast to double. */
inline double t_double(const dd_real &a) {
  return a.x[0];
}

inline double t_double(const qd_real &a) {
  return a[0];
}

inline double t_double(const double &a){
  return a;
}

// Quasinewton class declaration
template<typename T>
class quasinewton{
protected:
  const double C1 = 0.0001;
  const double C2 = 0.9;
  bool done;
  clock_t t1, t2;
  size_t n, n1, n2, nm, m1, tmp, maxit, m;
  size_t it = 0, ol = 1, cs = 0, nfevalval = 0;
  T *g, *p, *s, *y, *f, *qpoptvalptr, *x, *fopt;
  double *u, *l, *xcauchy;
  bool *freeVariable;
  T t, gnorm, gtp, fval, fprev, qpoptval, taud, ftarget, gnormtol;
  /* integer pointers: */  
  size_t* nfeval, gradientsamplingN;
  int *exitflag;
  int echo, jcur, lm, exitflagval = 0;
  int(*testFunction)(T*, T*, T*, size_t);
  std::ofstream* output;
  std::ofstream alloutput;
  const char * outputname;
  std::vector<T> breakpoints; //Contains the breakpoints to be ordered
  std::vector<T> breakpointsNOorder; 
  std::multimap<T, size_t> bpmemory; //Breakpoints automatically ordered
  
public:
  quasinewton(T [], T *, size_t ,  T,  int(*)(T*, T*, T*, size_t), std::ofstream&,  
	      T,  T,  size_t, short, short, const char *, size_t , size_t, double*, 
	      double*);
  ~quasinewton();
  void finish(){t2 = clock();};
  // The next functions prepare for the main loop
  void beforemainloop();
  virtual void befmainloopspecific() = 0; //To be implemented by each child
  void printbefmainloop(); // These parts are commong to both BFGS and LBFGS
  virtual void mainloop();
  virtual void mainloopspecific() = 0;
  void postmainloop();
  void runallsteps();
  void printallfinalinfo();
  void gettis();
  bool themin(double*, size_t);
  bool gradsamp(); // A few more iterations.  Defined only for double types
};

template<typename T>
quasinewton<T>::quasinewton(T x0[], T* fopt0, size_t n0,  T taud0,  
			    int(*tF)(T*, T*, T*, size_t), std::ofstream& output0,  
			    T ftarget0,  T gnormtol0,  size_t maxit0, short echo0, 
			    short lm0, const char * outputname0, size_t m0,
			    size_t gradientsamplingN0, double* u0, double* l0){
  x        = x0;
  fopt     = fopt0;
  n        = n0;
  taud     = taud0;
  ftarget  = ftarget0;
  gnormtol = gnormtol0;
  maxit    = maxit0;
  echo     = echo0;
  outputname = outputname0;
  output   = &output0;
  done     = false;
  quasinewton<T>::f = &(quasinewton<T>::fval);
  g        = new T[n];
  p        = new T[n];
  s        = new T[n];
  y        = new T[n];
  t1       = clock();  
  nfeval   = &(nfevalval);  
  exitflag = &(exitflagval);
  qpoptval = taud + 100;
  qpoptvalptr  = &(qpoptval);
  jcur     = 0;
  lm       = lm0;
  testFunction = tF;
  m        = m0;
  alloutput.open("../alloutput.txt", std::ios::app);
  gradientsamplingN = gradientsamplingN0;
  xcauchy = new double[n];
  freeVariable = new bool[n];
  u = u0;
  l = l0;
}

template<typename T>
quasinewton<T>::~quasinewton(){
  alloutput.close();
  delete [] g;
  delete [] p;
  delete [] s;
  delete [] y;
  delete [] xcauchy;
  delete [] freeVariable;
}

template<typename T>
void quasinewton<T>::beforemainloop(){
  // This function runs everything before the main loop.  
  // Puts all the other pieces together
  testFunction(f, g, x, n);
  *nfeval = *nfeval + 1;
  gnorm  = vecnorm<T>(g, n);
}

template<typename T>
void quasinewton<T>::printbefmainloop(){
  if (echo > 0) 
    print_init_info<T>(*output, n, ftarget, gnormtol, maxit, echo, lm, outputname);
  
  if (*f < ftarget) {
    done = 1;
    *exitflag = 3;
  }

  if (2 == echo) 
    print_iter_info<T>(*output, it, f, gnorm, jcur, qpoptvalptr, t);
}

template<typename T>
void quasinewton<T>::mainloop(){
  it++;
  gtp = vecip<T>(g, p, n);   // g'*p
  if(gtp > 0){
    done = true;
    *exitflag = -2;
  }
  vcopyp<T>(s, x, -1.0, n);
  vcopyp<T>(y, g, -1.0, n);
  /* Copy f before overwritten by linesearch in case of NaN */
  fprev = *f;
  /* line search:*/
  t = linesearch_ww<T>(x, f, g, p, C1, C2, n, testFunction, nfeval, 
		       ftarget, exitflag);
  /* If f is NaN, exit with best found f,x,g so far */
  if (::isnan(*f)) {
    *exitflag = -8;
    /* at this point, s=-x, so take x from s (same with y/g) */
    vcopyp<T>(x, s, -1.0, n);
    vcopyp<T>(g, y, -1.0, n);
    /* best f found is stored in fprev. Exit with that: */
    *f = fprev;
  }
  if (*f > fprev)
    *f = fprev;
  /* calculate s and y:
     before these calls, s and y contain
     previous -x and prev. -g, so e.g. s = x - xprev is s = s + x */
  vpv<T>(s, x, 1, n);
  vpv<T>(y, g, 1, n);
  gnorm = vecnorm<T>(g, n);
  
  //Quasinewton Update specific for each scheme
  mainloopspecific();
  /* print iteration info:*/
  if (2 == echo)
    print_iter_info<T>(*output, it, f, gnorm, jcur, qpoptvalptr, t);
  
  /* check convergence here*/

  if (it >= maxit)
    *exitflag = -1;
  if (*f < ftarget)
    *exitflag = 1;
  if (gnorm < gnormtol)
    *exitflag = 2;
  if (*qpoptvalptr < taud)
    *exitflag = 7;
  
  if(-1 == *exitflag){
    bool gradsampresult;
    gradsampresult = gradsamp();
    gradsampresult ? *exitflag = 1 : *exitflag = -1;
  }
  
  /* if exitflag was changed: exit main loop: */
  if (*exitflag != 0) done = true;
}

template<typename T>
void quasinewton<T>::postmainloop(){
  t2 = clock();
  double ttime = (double) ((double)(t2 - t1)/ (double)CLOCKS_PER_SEC );
    
  if (echo > 0) {
    print_final_info<T>(*output, it, f, gnorm, *nfeval, *exitflag, ttime);
    printallfinalinfo();
  }
  *fopt = *f;
  std::cout << "valor final: " << *f << std::endl;
  /* Gather rundata in the double array 'info' 
     info[0] = (double) (*nfeval);
     info[1] = (double) (*exitflag);
     info[2] = (double) (lm);
     info[3] = (double) (ttime);
  */
}

template<typename T>
void quasinewton<T>::printallfinalinfo(){
  alloutput << " f=" << *f << ", n=" << n << ", |g|=" << gnorm << std::endl;
}

template<typename T>
void quasinewton<T>::runallsteps(){
  beforemainloop();
  befmainloopspecific(); // Depending on whether BFGS or LBFGS is being implemented
  printbefmainloop();
  while(!done)
    mainloop();
  postmainloop();
}

template<typename T>
void quasinewton<T>::gettis(){
  // This function gets all the Ti points described in (4.1) of 8limited**
  // It also sorts them at the end
  for(size_t i = 0; i < n; i++){
    if(0 == g[i]){
      breakpoints.push_back(std::numeric_limits<T>::max()); // Assign \Infty if g == 0
      bpmemory.insert(std::pair<T, size_t>(std::numeric_limits<T>::max(), i));
    } else {
      if(g[i] < 0) {
	breakpoints.push_back( (x[i] - u[i]) / g[i]);
	bpmemory.insert(std::pair<T, size_t>((x[i] - u[i]) / g[i], i));
      } else {
	breakpoints.push_back((x[i] - l[i]) / g[i]);
	bpmemory.insert(std::pair<T, size_t>((x[i] - l[i]) / g[i], i));
      }
    }
  }
  breakpointsNOorder = breakpoints;
  std::make_heap(breakpoints.begin(), breakpoints.end());
  std::sort_heap(breakpoints.begin(), breakpoints.end());
}
/*
  template<typename T>
  size_t quasinewton<T>::findDimension(size_t i){
  //This function finds the correct dimension that we should work with and pops the
  // the corresponding element in the multimap container
  // receives the current cardinal dimension I'm working with.  Returns the real
  // dimension to work with
  T ti;
  size_t thisdimension;
  ti =  breakpoints[i];
  std::multimap<T, size_t>::iterator it = bpmemory.find(ti);
  thisdimension = it->second;
  bpmemory.erase(it);
  return(thisdimension);
  }
*/

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGS: public quasinewton<T>{
protected:
  T *q, *H;
  double* Hdouble;

public:
  BFGS(T*& x0, T*& fopt0, size_t&, T&, int(*&)(T*, T*, T*, size_t), 
       std::ofstream&, T&, T&, size_t&, short&, short&, const char *&, size_t&, size_t&,
       double*&, double*&);
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
  void createDoubleH();
};

template<typename T>
BFGS<T>::BFGS(T*& x0, T*& fopt0, size_t& n0,  T& taud0,  
	      int(*&tF)(T*, T*, T*, size_t), std::ofstream& output0,  T& ftarget0,  
	      T& gnormtol0,  size_t& maxit0, short& echo0, short& lm0, 
	      const char *& outputname0, size_t& m0, size_t& gradientsamplingN0, 
	      double*& u0, double*& l0):
  quasinewton<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
		 echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
  quasinewton<T>::n1 = quasinewton<T>::n;
  quasinewton<T>::n2 = quasinewton<T>::n * quasinewton<T>::n;
  quasinewton<T>::nm = 1;
  quasinewton<T>::m1 = 1; // on LBFGS this variable is m (the history)
  q = new T[quasinewton<T>::n1];
  H = new T[quasinewton<T>::n2];
  Hdouble = new double[quasinewton<T>::n * quasinewton<T>::n];
}

template<typename T>
void BFGS<T>::befmainloopspecific(){
  mat_set_eye(H, quasinewton<T>::n, quasinewton<T>::n);
  mxv<T>(quasinewton<T>::p, H, quasinewton<T>::g, -1.0, 0.0, quasinewton<T>::n, 
	 quasinewton<T>::n);
}

template<typename T>
void BFGS<T>::mainloopspecific(){
  update_bfgs<T>(H, quasinewton<T>::p, quasinewton<T>::g, quasinewton<T>::s, 
		 quasinewton<T>::y, q, quasinewton<T>::n);
}

template<typename T>
void BFGS<T>::createDoubleH(){
  for(size_t i = 0; i < quasinewton<T>::n; i++){
    for(size_t j = 0; j < quasinewton<T>::n; j++)
      Hdouble[i * quasinewton<T>::n + j] = t_double(H[i * quasinewton<T>::n + j]);
  }
}
/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGSB: public BFGS<T>{
protected:

public:
  BFGSB(T*& x0, T*& fopt0, size_t&, T&, int(*&)(T*, T*, T*, size_t), 
	std::ofstream&, T&, T&, size_t&, short&, short&, const char *&, size_t&, 
	size_t&, double*&, double*&);
  double findGeneralizedCauchyPoint();
  void findMinimum2ndApproximation();
  void mainloop();
};

template<typename T>
BFGSB<T>::BFGSB(T*& x0, T*& fopt0, size_t& n0,  T& taud0,  
		int(*&tF)(T*, T*, T*, size_t), std::ofstream& output0,  T& ftarget0,  
		T& gnormtol0,  size_t& maxit0, short& echo0, short& lm0, 
		const char *& outputname0, size_t& m0, size_t& gradientsamplingN0,
		double*& u0, double*&l0):
  BFGS<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
	   echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
}

template<typename T>
double BFGSB<T>::findGeneralizedCauchyPoint(){
  quasinewton<T>::gettis();
  BFGS<T>::createDoubleH();
  char yTrans = 'T';
  char nTrans = 'N';
  double alpha = 1.0, beta = 0.0;
  int ndouble = static_cast<int> (quasinewton<T>::n);
  int one = 1;
  double tj, fpj, fppj, deltatj, oldtj;
  double* di = new double[quasinewton<T>::n];
  double* z = new double[quasinewton<T>::n];
  double* C = new double[quasinewton<T>::n];
  double adouble;
  double dtstar, tstar;
  for(size_t i0 = 0; i0 < quasinewton<T>::n; i0++){
    di[i0] = z[i0] = 0.0;
    quasinewton<T>::freeVariable[i0] = true;
  }
  // bear in mind that the following multimap is already ordered
  typename std::multimap<T, size_t>::iterator iter = quasinewton<T>::bpmemory.begin();
  // First of all.  Run the zeroeth step from the multistep gradient projection
  iter = quasinewton<T>::bpmemory.begin();
  size_t b = (*iter).second;
  deltatj = t_double((*iter).first); // Change from zero
  oldtj = tj = deltatj;
  // Find the new x position
  for(size_t i0 = 0; i0 < quasinewton<T>::n; i0++){
    quasinewton<T>::xcauchy[i0] = (t_double(quasinewton<T>::x[i0]) - 
				  t_double((*iter).first) * 
				  t_double(quasinewton<T>::g[i0]));
    quasinewton<T>::xcauchy[i0]= MIN(quasinewton<T>::xcauchy[i0], quasinewton<T>::u[i0]);
    quasinewton<T>::xcauchy[i0]= MAX(quasinewton<T>::xcauchy[i0], quasinewton<T>::l[i0]);
    z[i0] = t_double(quasinewton<T>::xcauchy[i0] - quasinewton<T>::x[i0]);
  }
  // Update new d_i coordinate
  di[b] = -t_double(quasinewton<T>::g[b]);
  quasinewton<T>::freeVariable[b] = false;
  
  // z^T*B*z
  dgemm_(&nTrans, &nTrans, &ndouble, &one, &ndouble, &alpha, BFGS<T>::Hdouble, &ndouble,
	 z, &ndouble, &beta, C, &ndouble);
  dgemm_(&yTrans, &nTrans, &one, &ndouble, &ndouble, &alpha, di, &ndouble, C, &ndouble, 
	 &beta, &adouble, &one);
  fpj = t_double(veciptd<T>(quasinewton<T>::g, di, quasinewton<T>::n)) + 
    adouble; 
  //g^Td +  z^T*B*z
  
  // d^T*B*z
  dgemm_(&nTrans, &nTrans, &ndouble, &one, &ndouble, &alpha, BFGS<T>::Hdouble, &ndouble,
	 di, &ndouble, &beta, C, &ndouble);
  dgemm_(&yTrans, &nTrans, &one, &ndouble, &ndouble, &alpha, di, &ndouble, C, &ndouble, 
	 &beta, &adouble, &one);
  fppj = adouble;
  dtstar = -fpj / fppj;
  tstar = dtstar + oldtj;
  typename std::multimap<T, size_t>::iterator titer = quasinewton<T>::bpmemory.begin();
  titer++;
  tj = t_double((*titer).first);
  if (tstar < tj){
    if (tstar > 0)
      return(tstar);
  }
  
  for(iter++; iter != quasinewton<T>::bpmemory.end(); iter++){
    tj = t_double((*iter).first);
    b = (*iter).second;
    deltatj = t_double((*iter).first) - oldtj;
    di[b] = -t_double(quasinewton<T>::g[b]);  // This is equation 4.2 (minus?)
    quasinewton<T>::freeVariable[b] = false;
    for(size_t i = 0; i < quasinewton<T>::n; i++){
      quasinewton<T>::xcauchy[i] = t_double((quasinewton<T>::x[i]) - 
							(*iter).first * 
							quasinewton<T>::g[i]);
      quasinewton<T>::xcauchy[i] = MIN(quasinewton<T>::xcauchy[i], 
				       quasinewton<T>::u[i]);
      quasinewton<T>::xcauchy[i] = MAX(quasinewton<T>::xcauchy[i], 
				       quasinewton<T>::l[i]);
      z[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
    }
    dgemm_(&nTrans, &nTrans, &ndouble, &one, &ndouble, &alpha, BFGS<T>::Hdouble, 
	   &ndouble, z, &ndouble, &beta, C, &ndouble); //C = B^T * z
    
    fpj = fpj + deltatj * fppj + std::pow(t_double(quasinewton<T>::g[b]), 
					  2) + 
      t_double(quasinewton<T>::g[b]) * C[b];
    dgemm_(&nTrans, &nTrans, &ndouble, &one, &ndouble, &alpha, BFGS<T>::Hdouble, 
	   &ndouble, di, &ndouble, &beta, C, &ndouble); //C = B^T * d
    fppj = fppj + 2 * t_double(quasinewton<T>::g[b]) * C[b] + 
      std::pow(t_double(quasinewton<T>::g[b]), 2) * 
      t_double(BFGS<T>::H[b * quasinewton<T>::n + b]);
    dtstar = fpj / fppj;
    tstar = oldtj + dtstar;
    if (tstar >= oldtj){
      if (tstar <= tj)
	return(tstar);
    }
    if (fpj >= 0)
      return(oldtj);

    oldtj = tj; // Update the time to the new end of the time frame
  }

  // In case nothing was found.  Return the last point
  return(tj);
  // I still need to implement the last segment to locate xcauchy correctly.
}

template<typename T>
void BFGSB<T>::findMinimum2ndApproximation(){
  // Assuming xcauchy has been correctly found.  This function runs a minimization of
  // the quadratic approximation to the goal function
  int numfree = 0; // number of free variables
  size_t b = 0;
  double* Z, *r, *dx, *d, *dnsize;
  char yTrans = 'T', nTrans = 'N';
  int ndouble = static_cast<int>(quasinewton<T>::n);
  int one = 1;
  double alpha = 1.0, beta = 0.0;
  double* C = new double[quasinewton<T>::n];
  dx = new double[quasinewton<T>::n];
  r = new double[numfree];
  d = new double[numfree];
  dnsize = new double[quasinewton<T>::n];

  for(size_t i = 0; i < quasinewton<T>::n; i++){
    if(quasinewton<T>::freeVariable[i])
      numfree++;
  }
  Z = new double[static_cast<int>(quasinewton<T>::n) * numfree];
  
  typename std::multimap<T, size_t>::iterator titer = quasinewton<T>::bpmemory.begin();
  for(size_t i = 0; i < quasinewton<T>::n; i++, titer++){
    for(int j = 0; j < numfree; j++)
      Z[static_cast<int>(i) * numfree + j] = 0.0;
    b = (*titer).second; //position of the ith. crossed boundary
    Z[static_cast<int>(b) * numfree + static_cast<int>(i)] = 1.0;
  }
  
  // Now let's define the r vector.  r = Z(g + H(Xcauchy - X))
  // is it maybe worth representing r as in equation 5.4 instead?
  for(size_t i = 0; i < quasinewton<T>::n; i++)
    dx[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
  dgemm_(&nTrans, &nTrans, &ndouble, &one, &ndouble, &alpha, BFGS<T>::Hdouble,
	 &ndouble, dx, &ndouble, &beta, C, &ndouble);
  for(size_t i = 0; i < quasinewton<T>::n; i++)
    C[i] += t_double(quasinewton<T>::g[i]);
  dgemm_(&yTrans, &nTrans, &ndouble, &one, &numfree, &alpha, Z,
	 &ndouble, C, &ndouble, &beta, r, &ndouble);
  
  // Find Bhat = Z^TBZ
  double* BZ, *BHAT;
  BZ = new double[static_cast<int>(quasinewton<T>::n) * numfree];
  BHAT = new double[numfree * numfree];
  dgemm_(&nTrans, &nTrans, &ndouble, &numfree, &ndouble, &alpha, BFGS<T>::Hdouble,
	 &ndouble, Z, &numfree, &beta, BZ, &ndouble);
  dgemm_(&yTrans, &nTrans, &numfree, &ndouble, &ndouble, &alpha, Z,
	 &numfree, BZ, &ndouble, &beta, BHAT, &one); //Warning  Check that &one
  
  // Solve the system 5.5 and 5.6
  // Notice that this system could easily be solved by inverting the matrix *BHAT
  // Notice that BHAT will be completely overwritten with an L and U decomposition...
  int info = 11;
  int* ipiv = new int[numfree];
  dgesv_(&numfree, &one, BHAT, &numfree, ipiv, r, &numfree, &info);
  // The solution is now on variable r
  double alpha0 = 0.0;
  double alphacandidate = 0.0;
  // Define the new boundaries which appear on 5.6
  double* lbf = new double[numfree];
  double* ubf = new double[numfree];
  size_t ind = 0;
  titer = quasinewton<T>::bpmemory.begin();
  for(; titer != bpmemory.end(); titer++, ind++)
    {
      b = (*titer).second;
      lbf[ind] = quasinewton<T>::l[b] - xcauchy[b]; 
      ubf[ind] = quasinewton<T>::u[b] - xcauchy[b];
      alphacandidate = MAX(ubf[ind] / r[ind], lbf[ind] / r[ind]); //WARNING: r is d here
      alpha0 = MAX(alpha0, alphacandidate);
    }
  alpha0 = MIN(alpha0, 1.0);
  for(size_t i = 0; i < numfree; i++){
    d[i] = alpha0 * r[i];   // Move back to vector d since it is less confusing
  }
  // Find the new solution
  titer = quasinewton<T>::bpmemory.begin();
  // Z_k * d
  dgemm_(&nTrans, &nTrans, &ndouble, &one, &numfree, &alpha, Z,
	 &ndouble, d, &ndouble, &beta, dnsize, &ndouble);  
  for(size_t i = 0; i < quasinewton<T>::n; i++){
    quasinewton<T>::x[i] = quasinewton::xcauchy[i];
    if((*titer).second == i){
      quasinewton<T>::x[i] += dnsize[i];
    }
  }
 /* 
    calculate s and y:
    before these calls, s and y contain
    previous -x and prev. -g, so e.g. s = x - xprev is s = s + x 
 */
  vpv<T>(quasinewton<T>::s, quasinewton<T>::x, 1, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1, quasinewton<T>::n);
  // Here's a question.  Do I use the new g or the previous one before the update
  // My guess is that it's probably similar
  update_bfgs<T>(H, dnsize, quasinewton<T>::g, quasinewton<T>::s, 
		 quasinewton<T>::y, q, quasinewton<T>::n);
  if (2 == quasinewton<T>::echo)
    print_iter_info<T>(*output, quasinewton<T>::it, quasinewton<T>::f, 
		       quasinewton<T>::gnorm, quasinewton<T>::jcur, 
		       quasinewton<T>::qpoptvalptr, alpha);
  if (quasinewton<T>::it >= quasinewton<T>::maxit)
    quasinewton<T>::*exitflag = -1;
  if (quasinewton<T>::*f < quasinewton<T>::ftarget)
    quasinewton<T>::*exitflag = 1;
  if (quasinewton<T>::gnorm < quasinewton<T>::gnormtol)
    quasinewton<T>::*exitflag = 2;
  if (quasinewton<T>::*qoptvalptr < quasinewton<T>::taud)
    quasinewton<T>::*exitflag = 7;
  
  /*
  // Can I do gradient sampling here?
  if(-1 == quasinewton<T>::*exitflag){
    bool gradsampresult;
    gradsampresult = gradsamp();
    gradsampresult
  }
  */
  /* if exitflag was changed: exit main loop */
  if (0 != quasinewton<T>::*exitflag) done = true;
}

template<typename T>
void BFGSB<T>::mainloop(){
  quasinewton<T>::it++;
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  quasinewton<T>::fprev = quasinewton<T>::*f;
  findGeneralizedCauchyPoint();
  findMinimum2ndApproximation();
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class LBFGS: public quasinewton<T>{
protected:
  T *S, *Y, *rho, *a;
public:
  LBFGS(T*& x0, T*& fopt0, size_t&, T&, int(*&)(T*, T*, T*, size_t), 
	std::ofstream&, T&, T&, size_t&, short&, short&, const char *&, size_t&, 
	size_t&, double*&, double*&);
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
};

template<typename T>
LBFGS<T>::LBFGS(T*& x0, T*& fopt0, size_t& n0,  T& taud0,  
		int(*&tF)(T*, T*, T*, size_t), std::ofstream& output0,  T& ftarget0,
		T& gnormtol0, size_t& maxit0, short& echo0, short& lm0, 
		const char *& outputname0, size_t& m0, size_t& gradientsamplingN0, 
		double*& u0, double*& l0):
  quasinewton<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
		 echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
  quasinewton<T>::n1 = 1;
  quasinewton<T>::n2 = 1;
  quasinewton<T>::nm = quasinewton<T>::n * quasinewton<T>::m;
  quasinewton<T>::m1 = quasinewton<T>::m;
  S  = new T[quasinewton<T>::nm];
  Y  = new T[quasinewton<T>::nm];
  rho = new T[quasinewton<T>::m1];
  a  = new T[quasinewton<T>::m1];
}

template<typename T>
void LBFGS<T>::befmainloopspecific(){
  // p = -g (LBFGS)
  vcopyp<T>(quasinewton<T>::p, quasinewton<T>::g, -1.0, quasinewton<T>::n);
}

template<typename T>
void LBFGS<T>::mainloopspecific(){
  if(quasinewton<T>::it <= quasinewton<T>::m) quasinewton<T>::cs++;
  quasinewton<T>::tmp = (quasinewton<T>::ol - 1) * quasinewton<T>::n;
  vcopy<T>(S + quasinewton<T>::tmp, quasinewton<T>::s, quasinewton<T>::n);
  vcopy<T>(Y + quasinewton<T>::tmp, quasinewton<T>::y, quasinewton<T>::n);
  rho[quasinewton<T>::ol - 1] = 1.0 / vecip<T>(quasinewton<T>::y, quasinewton<T>::s, 
					       quasinewton<T>::n);
  update_lbfgs<T>(quasinewton<T>::p, S, Y, rho, a, quasinewton<T>::g, 
		  quasinewton<T>::cs, quasinewton<T>::ol, quasinewton<T>::m, 
		  quasinewton<T>::n);
  quasinewton<T>::ol = (quasinewton<T>::ol % quasinewton<T>::m) + 1;
}

/////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class LBFGSB: public LBFGS<T>{
protected:
public:
  LBFGSB(T*& x0, T*& fopt0, size_t&, T&, int(*&)(T*, T*, T*, size_t), 
	 std::ofstream&, T&, T&, size_t&, short&, short&, const char *&, size_t&, 
	 size_t&, double*&, double*&);
};

template<typename T>
LBFGSB<T>::LBFGSB(T*& x0, T*& fopt0, size_t& n0,  T& taud0,  
		  int(*&tF)(T*, T*, T*, size_t), std::ofstream& output0,  T& ftarget0,  
		  T& gnormtol0,  size_t& maxit0, short& echo0, short& lm0, 
		  const char *& outputname0, size_t& m0, size_t& gradientsamplingN0,
		  double*& u0, double*&l0):
  LBFGS<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
	  echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
}

/////////////////////////////////////////////////////////////////////////////////////
/* BFGS MAIN ALGORITHM: */
template<class T>
void bfgs(T*& x, T*& fopt,  size_t& n,  short& lm,  size_t& m, T& ftarget, 
	  T& gnormtol, size_t& maxit,  long& J, T& taux,  T& taud, short& echo, 
	  int(*&testFunction)(T*, T*, T*, size_t), std::string& datafilename, 
          double*& info, size_t& gradientsamplingN0, double*& u0, double*& l0, 
	  bool& boundedProblem){
  
  std::ofstream output;
  std::cout <<"echo is: "<< echo << " datafilename is: " << datafilename << std::endl;
  output.open(datafilename.c_str(), std::ios::app);
  const char * outputname = datafilename.c_str();
  if(!lm){
    if(boundedProblem){
      BFGS<T>* mybfgs;
      mybfgs = new BFGSB<T>(x, fopt, n, taud, testFunction, output, ftarget, gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      mybfgs->runallsteps();
    } else {
      BFGS<T>* mybfgs;
      mybfgs = new BFGS<T>(x, fopt, n, taud, testFunction, output, ftarget, gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      mybfgs->runallsteps();
    }
  } else {
    if(boundedProblem){
      LBFGS<T>* mylbfgs;
      mylbfgs = new LBFGSB<T>(x, fopt, n, taud, testFunction, output, ftarget,gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      mylbfgs->runallsteps();
    } else{
      LBFGS<T>* mylbfgs;
      mylbfgs = new LBFGS<T>(x, fopt, n, taud, testFunction, output, ftarget, gnormtol,
			     maxit, echo, lm, outputname, m, gradientsamplingN0, u0,l0);
      mylbfgs->runallsteps();
    }
  }
  taux = taux + 1; J++; info[0] = info[0] + 1;
  output.close();
}

template<typename T>
bool quasinewton<T>::themin(double *gradpoints, size_t i1){
  double mindist = 1.0;
  size_t mybase = i1 * n;
  size_t loc = 0;
  bool threshold = false;
  for(size_t i = 0; i < n; i++){
    loc = mybase + i;
    mindist = MIN((x[i] - gradpoints[loc]), mindist);
  }
  (mindist < 1e-14) ? threshold = true : threshold = false;
  return(threshold);
}

template<typename T>
bool quasinewton<T>::gradsamp(){
  std::cout << "gradsamp doesn't do anything for the current type.  Try doubles"<< 
    std::endl;
  //Do nothing unless T is a double type
  return(false);
}

template <> //template specialization for double (only doubles work with lapack)
bool quasinewton<double>::gradsamp(){
  // n is the size of x and gradientsamplingN is the number of gradient samples 
  // assigned in algoparameters.  Gradpoints contains all the points to be evaluated
  size_t j = 0;
  double * gradpoints;
  gradpoints = new double[static_cast<size_t>(gradientsamplingN) * 
			  static_cast<size_t>(n)];
  // Assign original x point to the gradpoints store
  for (size_t i = 0; i < n; i++)
    gradpoints[i] = x[i];
  
  // Assign random numbers to the rest of variables.
  for (size_t i = n; i < static_cast<size_t>(gradientsamplingN) * 
	 static_cast<size_t>(n); i++, j++){
    if(n == j)
      j = 0;
    gradpoints[i] = x[j] + (rand() / RAND_MAX) - 0.5;
  }
  
  // A distance of up to 0.5 is usually too big in any direction.  So I will divide all
  // coordinates by two until the smallest coordinate difference is < 1e-14
  for (size_t i = 1; i < gradientsamplingN; i++){
    do {
      for(size_t k = 0; k < n; k++)
	gradpoints[i + k] = (gradpoints[i + k] + x[i]) / 2;      
    } while (themin(gradpoints, i));
  }
  
  // found the x's now. Let's find the gradients
  double * gradgrads;
  // double f;
  
  gradgrads = new double[static_cast<size_t>(gradientsamplingN) * 
    static_cast<size_t>(n)];
for(size_t i = 0; i < static_cast<size_t>(gradientsamplingN); i++){
  testFunction(f, g +(i * n), gradpoints + (i * n), n);
 }
  
// Next step call qpspecial
 std::cout << "calling qpspecial" << std::endl;
 qpclass<double> * myopt = new qpclass<double>(static_cast<int>(n), 
					       static_cast<int>(gradientsamplingN), 
					       gradgrads, 100);
 myopt->optimization();
  
  double * solution = new double[n];
  myopt->fetchSolution(solution);
  return(true);
}


#endif // _BFGS_TEMPLATE_HPP_

// ** This is the paper: "A limited memory algorithm for bound constrained optimization"
