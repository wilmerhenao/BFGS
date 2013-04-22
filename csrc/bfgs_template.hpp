/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BFGS tools Copyright (C) 2013  Wilmer Henao
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#include "nummatrix.hpp"

template<class T>
void bfgs(T*&, T*& fopt, int& n, short& lm, int& m, T& ftarget,  T& gnormtol,  
	  int& maxit,  long& J, T& taux,  T& taud, short& echo, 
	  int(*&testFunction)(T*, T*, T*, int),  std::string& datafilename, 
          double*&, int&, double*&, double*&, bool&);

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
  int n, n1, n2, nm, m1, tmp, maxit, m;
  int it = 0, ol = 1, cs = 0, nfevalval = 0;
  T *g, *p, *s, *y, *f, *qpoptvalptr, *x, *fopt;
  double *u, *l, *xcauchy;
  bool *freeVariable;
  T t, gnorm, gtp, fval, fprev, qpoptval, taud, ftarget, gnormtol;
  /* integer pointers: */  
  int* nfeval, gradientsamplingN;
  int *exitflag;
  int echo, jcur, lm, exitflagval = 0;
  int(*testFunction)(T*, T*, T*, int);
  std::ofstream* output;
  std::ofstream alloutput;
  const char * outputname; 
  std::multimap<T, int> bpmemory; //Breakpoints automatically ordered
  
public:
  quasinewton(T [], T *, int ,  T,  int(*)(T*, T*, T*, int), std::ofstream&,  
	      T,  T,  int, short, short, const char *, int , int, double*, 
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
  void get_ti_s();
  bool themin(double*, int);
  bool gradsamp(); // A few more iterations.  Defined only for double types
  //int getn();
};
/*
template<typename T>
int quasinewton<T>::getn(){
  std::cout << "getn : " << this->n << std::endl; 
  return this->n;
}
*/
template<typename T>
quasinewton<T>::quasinewton(T x0[], T* fopt0, int n0,  T taud0,  
			    int(*tF)(T*, T*, T*, int), std::ofstream& output0,  
			    T ftarget0,  T gnormtol0,  int maxit0, short echo0, 
			    short lm0, const char * outputname0, int m0,
			    int gradientsamplingN0, double* u0, double* l0){
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
  for(int i0 = 0; i0 < n; i0++)
    freeVariable[i0] = true;
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
void quasinewton<T>::get_ti_s(){
  // This function gets all the Ti points described in (4.1) of 8limited**
  // It also sorts them at the end
  for(int i = 0; i < n; i++){
    if(0.0 == g[i]){
      // Assign \Infty if g == 0
      bpmemory.insert(std::pair<T, int>(std::numeric_limits<T>::max(), i));
    } else {
      if(g[i] < 0) {
	bpmemory.insert(std::pair<T, int>((x[i] - u[i]) / g[i], i));
      } else {
	bpmemory.insert(std::pair<T, int>((x[i] - l[i]) / g[i], i));
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGS: public quasinewton<T>{
protected:
  T *q, *H;
  double* Hdouble;
  Matrix<double> mHdouble;
public:
  BFGS(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
       std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, int&,
       double*&, double*&);
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
  void createDoubleH();
};

template<typename T>
BFGS<T>::BFGS(T*& x0, T*& fopt0, int& n0,  T& taud0,  
	      int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
	      T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
	      const char *& outputname0, int& m0, int& gradientsamplingN0, 
	      double*& u0, double*& l0):
  quasinewton<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
		 echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0),
  mHdouble(quasinewton<T>::n, quasinewton<T>::n){
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

// Hopefully this function will become obsolete if I include qd analysis
template<typename T>
void BFGS<T>::createDoubleH(){
  for(int i = 0; i < this->n; i++){ //left this->n here on purpose for whoever is
                                    //reading.  This way you don't know the scope!
    for(int j = 0; j < quasinewton<T>::n; j++)
      Hdouble[i * quasinewton<T>::n + j] = t_double(H[i * quasinewton<T>::n + j]);
  }
  Matrix<double> temp(Hdouble, quasinewton<T>::n, quasinewton<T>::n);
  mHdouble = temp;
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGSB: public BFGS<T>{
protected:
  double tstar;
  char yTrans, nTrans;
  double alpha, beta, tj, fpj, fppj, deltatj, oldtj, adouble, dtstar;
  int ndouble, one, b;
  double *di, *z, *C;
  Matrix<double> mZ, mdi;
  // bear in mind that the following multimap is already ordered
  typename std::multimap<T, int>::iterator iter;
public:
  BFGSB(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	int&, double*&, double*&);
  void zBz();
  void dBz();
  void zeroethstep();
  void lapackzerostep();
  void findXCauchymX(int);
  void lapackmanipulations();
  void tstarcalculation();
  void findGeneralizedCauchyPoint();
  void findMinimum2ndApproximation();
  void mainloop();
  void printFinalConditions();
};

// constructor
template<typename T>
BFGSB<T>::BFGSB(T*& x0, T*& fopt0, int& n0,  T& taud0,  
		int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
		T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
		const char *& outputname0, int& m0, int& gradientsamplingN0,
		double*& u0, double*&l0):
  BFGS<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
	  echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0), 
  mZ(quasinewton<T>::n, 1), mdi(quasinewton<T>::n, 1){
  tstar = 0.0;
  yTrans = 'T'; nTrans = 'N';
  alpha = 1.0; beta = 0.0;
  ndouble = quasinewton<T>::n;
  one = 1;
  di = new double[quasinewton<T>::n];
  z = new double[quasinewton<T>::n];
  for(int i0 = 0; i0 < quasinewton<T>::n; i0++)
    di[i0] = z[i0] = 0.0;
  C = new double[quasinewton<T>::n];
}

template<typename T>
void BFGSB<T>::zeroethstep(){
  // First of all.  Run the zeroeth step from the multistep gradient projection
  iter = quasinewton<T>::bpmemory.begin();
  b = iter->second;
  deltatj = t_double(iter->first); // Change from zero
  
  oldtj = tj = deltatj;
  // Find the new x position
  for(int i0 = 0; i0 < quasinewton<T>::n; i0++){
    quasinewton<T>::xcauchy[i0] = (t_double(quasinewton<T>::x[i0]) - tj *
				  t_double(quasinewton<T>::g[i0]));
    quasinewton<T>::xcauchy[i0]= MIN(quasinewton<T>::xcauchy[i0], quasinewton<T>::u[i0]);
    quasinewton<T>::xcauchy[i0]= MAX(quasinewton<T>::xcauchy[i0], quasinewton<T>::l[i0]);
    z[i0] = t_double(quasinewton<T>::xcauchy[i0] - quasinewton<T>::x[i0]);
  }
  Matrix<double> temp(z, quasinewton<T>::n, 1);
  mZ = temp;
  
  // Update new d_i coordinate
  di[b] = -t_double(quasinewton<T>::g[b]);
  Matrix<double> temp2(di, quasinewton<T>::n, 1);
  mdi = temp2;
  
  quasinewton<T>::freeVariable[b] = false;
}

template<typename T>
void BFGSB<T>::zBz(){
  adouble = squareForm(mZ, BFGS<T>::mHdouble, mZ);
}

template<typename T>
void BFGSB<T>::dBz(){
  adouble = squareForm(mdi, BFGS<T>::mHdouble, mZ);
}

template<typename T>
void BFGSB<T>::lapackzerostep(){
  /*
    This method includes all the lapack routines that have to be run at the beginning
    of the "find the cauchy point iteration"
  */
  
  double * grad;
  grad = new double[this->n];
  
  for(int i0 = 0; i0 < this->n; i0++){
    grad[i0] = t_double(this->g[i0]) * di[i0];
  }
  zBz();
  fpj = adouble;
  for(int i0 = 0; i0 < this->n; i0++){
    // veciptd<double>(quasinewton<T>::g, di, ndouble) + adouble;
    fpj = fpj + grad[i0];
  }
  
  std::cout << "checking existence after veciptd FPJ ->" << fpj << std::endl;
  //g^Td +  z^T*B*z
  
  // d^T*B*z
  dBz();
  fppj = adouble;
  dtstar = -fpj / fppj;
  tstar = dtstar + oldtj;
  typename std::multimap<T, int>::iterator titer = quasinewton<T>::bpmemory.begin();
  titer++;
  tj = t_double((*titer).first);
  if (tstar < tj){
    if (tstar > 0)
      exit(0);
  }
}

template<typename T>
void BFGSB<T>::findXCauchymX(int i){
  // Calculates the difference between the new "cauchy" point and X
  quasinewton<T>::xcauchy[i] = t_double((quasinewton<T>::x[i]) - 
					iter->first * 
					quasinewton<T>::g[i]);
  quasinewton<T>::xcauchy[i] = MIN(quasinewton<T>::xcauchy[i], 
				   quasinewton<T>::u[i]);
  quasinewton<T>::xcauchy[i] = MAX(quasinewton<T>::xcauchy[i], 
				   quasinewton<T>::l[i]);
  z[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
}

template<typename T>
void BFGSB<T>::lapackmanipulations(){
  // LAPACK manipulations for each of the loops in the xcauchy calculations

  Matrix<double> temp(1, quasinewton<T>::n);
  matrixMultiply(mZ, BFGS<T>::mHdouble, temp, 'T', 'N');
  fpj = fpj + deltatj * fppj + std::pow(t_double(quasinewton<T>::g[b]), 
					2) + 
    t_double(quasinewton<T>::g[b]) * temp(b);

  Matrix<double> temp2(quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mHdouble, mdi, temp2);
  // Reassigning values to C
  for(int i = 0; i < quasinewton<T>::n; i++)
    C[i] = temp2(i);
}

template<typename T>
void BFGSB<T>::tstarcalculation(){
  fppj = fppj + 2 * t_double(quasinewton<T>::g[b]) * C[b] + 
    std::pow(t_double(quasinewton<T>::g[b]), 2) * 
    t_double(BFGS<T>::H[b * quasinewton<T>::n + b]);
  dtstar = fpj / fppj;
  tstar = oldtj + dtstar;
}

template<typename T>
void BFGSB<T>::findGeneralizedCauchyPoint(){
  iter++;
  for(; iter != quasinewton<T>::bpmemory.end(); iter++){
    tj = t_double(iter->first);
    b = iter->second;
    deltatj = t_double(iter->first) - oldtj;
    di[b] = -t_double(quasinewton<T>::g[b]);  // This is equation 4.2 (minus?)
    quasinewton<T>::freeVariable[b] = false;

  for(int i = 0; i < quasinewton<T>::n; i++){
    std::cout << "freevariable: " << quasinewton<T>::freeVariable[i] << std::endl;
  }

    for(int i = 0; i < quasinewton<T>::n; i++){
      findXCauchymX(i);
    }

    // Organize matrices
    Matrix<double> temp(z, quasinewton<T>::n, 1);
    mZ = temp;
    Matrix<double> temp2(z, quasinewton<T>::n, 1);
    mdi = temp2;
    lapackmanipulations();
    tstarcalculation();
    if (tstar >= oldtj){
      if (tstar <= tj){
	std::cout << "found optimal cauchy point." << std::endl;
	return;
      }
    }
    if (fpj >= 0){
      tstar = oldtj;
      exit(0);
    }    
    oldtj = tj; // Update the time to the new end of the time frame
  }

  // In case nothing was found.  Return the last point
  tstar = tj;
  // I still need to implement the last segment to locate xcauchy correctly.
  std::cout << "exiting after finding generalized cauchy point (not optimal)"<<std::endl;
}

template<typename T>
void BFGSB<T>::findMinimum2ndApproximation(){
  // Assuming xcauchy has been correctly found.  This function runs a minimization of
  // the quadratic approximation to the goal function
  int numfree = 0; // number of free variables
  b = 0;
  double* ZfM2, *r, *dx, *dnsize;
  int myn;
  myn = quasinewton<T>::n;
  myn = myn + 1;
  dx = new double[quasinewton<T>::n];
  dnsize = new double[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    if(quasinewton<T>::freeVariable[i])
      numfree++;
  }
 
  ZfM2 = new double[(quasinewton<T>::n) * numfree];
  
  typename std::multimap<T, int>::iterator titer = quasinewton<T>::bpmemory.begin();
  for(int i = 0; i < quasinewton<T>::n; i++, titer++){
    for(int j = 0; j < numfree; j++)
      ZfM2[(i) * numfree + j] = 0.0;
    b = (*titer).second; //position of the ith. crossed boundary
    ZfM2[(b) * numfree + (i)] = 1.0;
  }
  std::cout << "Entered second approximation function" << std::endl;
  // Now let's define the r vector.  r = Z(g + H(Xcauchy - X))
  // is it maybe worth representing r as in equation 5.4 instead?
  for(int i = 0; i < quasinewton<T>::n; i++)
    dx[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);

  std::cout << "Before lapack computations" << std::endl;
  Matrix<double> mdx(dx, quasinewton<T>::n, 1), mC(C, quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mHdouble, mdx, mC); // Result kept in mC 
  std::cout << "After lapack computations" << std::endl;

  for(int i = 0; i < quasinewton<T>::n; i++){
    C[i] = mC(i);  //Warning!.  I need to correct for this double assignation
    C[i] += t_double(quasinewton<T>::g[i]);
  }

  r = new double[numfree];
  Matrix<double> mZfM2(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mmC(C, quasinewton<T>::n, 1);
  Matrix<double> mr(r, numfree, 1);
  std::cout << "numfree: " << numfree << "quasinewton<T>::n :  " << quasinewton<T>::n << std::endl;
  std::cout << "before the suspicious matrixMultiply" << std::endl;
  matrixMultiply(mZfM2, mmC, mr, 'T', 'N'); // the result is now on mr;

  // Find Bhat = Z^TBZ
  Matrix<double> mBHAT(numfree, numfree);
  std::cout << "Before more lapack computations" << std::endl;
  GensquareForm(mZfM2, BFGS<T>::mHdouble, mZfM2, mBHAT);
  std::cout << "After more lapack computations" << std::endl;


  // Solve the system 5.5 and 5.6
  // Notice that this system could easily be solved by inverting the matrix *BHAT
  // Notice that BHAT will be completely overwritten with an L and U decomposition...
  
  Matrix<double> md(numfree, 1);  // Where to put the solution 
  bfgssolver(mBHAT, mr, md);

  double alpha0 = 0.0;
  double alphacandidate = 0.0;
  // Define the new boundaries which appear on 5.6
  double* lbf = new double[numfree];
  double* ubf = new double[numfree];
  int ind = 0;
  titer = quasinewton<T>::bpmemory.begin();
  for(; titer != quasinewton<T>::bpmemory.end(); titer++, ind++)
    {
      b = (*titer).second;
      lbf[ind] = quasinewton<T>::l[b] - quasinewton<T>::xcauchy[b];
      ubf[ind] = quasinewton<T>::u[b] - quasinewton<T>::xcauchy[b];
      alphacandidate = MAX(ubf[ind] / md(ind), lbf[ind] / md(ind));
      alpha0 = MAX(alpha0, alphacandidate);
    }
  alpha0 = MIN(alpha0, 1.0);
  md *= alpha0;

  // Find the new solution
  titer = quasinewton<T>::bpmemory.begin();
  // Z_k * d
  Matrix<double> mdnsize(dnsize, quasinewton<T>::n, 1);
  matrixMultiply(mZfM2, md, mdnsize);

  for(int i = 0; i < quasinewton<T>::n; i++){
    quasinewton<T>::x[i] = quasinewton<T>::xcauchy[i];
    if((*titer).second == i){
      quasinewton<T>::x[i] += mdnsize(i);
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
  // ANSWER: USE THE NEW g!!!
  T* dnsizeT = new T[quasinewton<T>::n];
  for(int i = 0; i < quasinewton<T>::n; i++)
	dnsizeT[i] = mdnsize(i);
  update_bfgs<T>(BFGS<T>::H, dnsizeT, quasinewton<T>::g, 
		 quasinewton<T>::s, quasinewton<T>::y, BFGS<T>::q, quasinewton<T>::n);
  T alpha1 = static_cast<T>(alpha0);
  if (2 == quasinewton<T>::echo)
    print_iter_info<T>(*quasinewton<T>::output, quasinewton<T>::it, quasinewton<T>::f, 
		       quasinewton<T>::gnorm, quasinewton<T>::jcur, 
		       quasinewton<T>::qpoptvalptr, alpha1);
  if (quasinewton<T>::it >= quasinewton<T>::maxit)
    *quasinewton<T>::exitflag = -1;
  if (*quasinewton<T>::f < quasinewton<T>::ftarget)
    *quasinewton<T>::exitflag = 1;
  if (quasinewton<T>::gnorm < quasinewton<T>::gnormtol)
    *quasinewton<T>::exitflag = 2;
  if (*quasinewton<T>::qpoptvalptr < quasinewton<T>::taud)
    *quasinewton<T>::exitflag = 7;
  // Don't do Gradiend Sampling
  /* if exitflag was changed: exit main loop */
  if (0 != *quasinewton<T>::exitflag) quasinewton<T>::done = true;
}

template<typename T>
void BFGSB<T>::printFinalConditions(){
    std::cout << "values of vector g" << std::endl;

  for(int i = 0; i < quasinewton<T>::n; i++)
    std::cout << quasinewton<T>::g[i] << std::endl;

}

template<typename T>
void BFGSB<T>::mainloop(){
  quasinewton<T>::it++;
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  quasinewton<T>::fprev = *quasinewton<T>::f;
  quasinewton<T>::get_ti_s();
  BFGS<T>::createDoubleH();
  zeroethstep();
  lapackzerostep();
  findGeneralizedCauchyPoint(); // up to here all good
  findMinimum2ndApproximation();
  std::cout << "Print final conditions to see if we converged" << std::endl;
  printFinalConditions();
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class LBFGS: public quasinewton<T>{
protected:
  T *S, *Y, *rho, *a;
public:
  LBFGS(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	int&, double*&, double*&);
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
};

template<typename T>
LBFGS<T>::LBFGS(T*& x0, T*& fopt0, int& n0,  T& taud0,  
		int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,
		T& gnormtol0, int& maxit0, short& echo0, short& lm0, 
		const char *& outputname0, int& m0, int& gradientsamplingN0, 
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
  LBFGSB(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	 std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	 int&, double*&, double*&);
};

template<typename T>
LBFGSB<T>::LBFGSB(T*& x0, T*& fopt0, int& n0,  T& taud0,  
		  int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
		  T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
		  const char *& outputname0, int& m0, int& gradientsamplingN0,
		  double*& u0, double*&l0):
  LBFGS<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
	  echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
}

/////////////////////////////////////////////////////////////////////////////////////

template<typename T>
bool quasinewton<T>::themin(double *gradpoints, int i1){
  double mindist = 1.0;
  int mybase = i1 * n;
  int loc = 0;
  bool threshold = false;
  for(int i = 0; i < n; i++){
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
  int j = 0;
  double * gradpoints;
  gradpoints = new double[(gradientsamplingN) * 
			  (n)];
  // Assign original x point to the gradpoints store
  for (int i = 0; i < n; i++)
    gradpoints[i] = x[i];
  
  // Assign random numbers to the rest of variables.
  for (int i = n; i < (gradientsamplingN) * 
	 (n); i++, j++){
    if(n == j)
      j = 0;
    gradpoints[i] = x[j] + (rand() / RAND_MAX) - 0.5;
  }
  
  // A distance of up to 0.5 is usually too big in any direction.  So I will divide all
  // coordinates by two until the smallest coordinate difference is < 1e-14
  for (int i = 1; i < gradientsamplingN; i++){
    do {
      for(int k = 0; k < n; k++)
	gradpoints[i + k] = (gradpoints[i + k] + x[i]) / 2;      
    } while (themin(gradpoints, i));
  }
  
  // found the x's now. Let's find the gradients
  double * gradgrads;
  // double f;
  
  gradgrads = new double[(gradientsamplingN) * 
    (n)];
for(int i = 0; i < (gradientsamplingN); i++){
  testFunction(f, g +(i * n), gradpoints + (i * n), n);
 }
  
// Next step call qpspecial
 std::cout << "calling qpspecial" << std::endl;
 qpclass<double> * myopt = new qpclass<double>((n), 
					       (gradientsamplingN), 
					       gradgrads, 100);
 myopt->optimization();
  
  double * solution = new double[n];
  myopt->fetchSolution(solution);
  return(true);
}

/* BFGS MAIN ALGORITHM: */
template<class T>
void bfgs(T*& x, T*& fopt,  int& n,  short& lm,  int& m, T& ftarget, 
	  T& gnormtol, int& maxit,  long& J, T& taux,  T& taud, short& echo, 
	  int(*&testFunction)(T*, T*, T*, int), std::string& datafilename, 
          double*& info, int& gradientsamplingN0, double*& u0, double*& l0, 
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

#endif // _BFGS_TEMPLATE_HPP_

// ** This is the paper: "A limited memory algorithm for bound constrained optimization"
