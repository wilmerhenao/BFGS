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

#ifndef NDEBUG
  #define FLAG()								\
        std::cout  << "Checking position.  Reached line  " << __LINE__ << std::endl; \
        std::cout  << "in file " << __FILE__ << std::endl; 
  #define SHOW(x) \
        std::cout << "Value of variable " << #x << " is: " << x << std::endl;
  #define PRINTARRAY(Z,m,n) \
        std::cout << "Printing array "<<#Z << ": "<< std::endl;	\
        for(int i = 0; i < m; i++){ \
          for(int j = 0; j < n; j++){ \
            std::cout << Z[i + j * m] << " "; \
          } \
          std::cout << std::endl; \
        }
#else
  #define SHOW(x)
  #define FLAG()
  #define PRINTARRAY(Z,m,n)
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
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
  bool done, thisIterationConverged;
  clock_t t1, t2;
  int n, n1, n2, nm, m1, tmp, maxit, m;
  int it = 0, ol = 1, cs = 0, nfevalval = 0;
  T *g, *p, *s, *y, *f, *qpoptvalptr, *x, *fopt;
  double *u, *l, *xcauchy;
  bool *freeVariable;
  T t, gnorm, gtp, fprev, qpoptval, taud, ftarget, gnormtol;
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
  quasinewton(T [], T *, int,  T,  int(*)(T*, T*, T*, int), std::ofstream&,  
	      T,  T,  int, short, short, const char *, int , int, double*, 
	      double*);
  virtual ~quasinewton();
  void finish(){t2 = clock();};
  // The next functions prepare for the main loop
  void FirstEvaluation();
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
};

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
  thisIterationConverged = false;
  quasinewton<T>::f = &(quasinewton<T>::ftarget);
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
  for(int _i = 0; _i < n; _i++){
    freeVariable[_i] = true;
    y[_i] = s[_i] = p[_i] = g[_i] = 0.0;
  }
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
void quasinewton<T>::FirstEvaluation(){
  /*
    Initial evaluation of the function and its gradient
  */

  testFunction(f, g, x, n);
  *nfeval = *nfeval + 1;
  gnorm  = vecnorm<T>(g, n); // norm of the gradient (not used in BFGSB algos)
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
  gtp = vecip<T>(g, p, n);   // g'*p
  if(gtp > 0){
    done = true;
    *exitflag= -2;
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
  std::cout << "Final Value: " << *f << std::endl;
  
}

template<typename T>
void quasinewton<T>::printallfinalinfo(){
  alloutput << " f=" << *f << ", n=" << n << ", |g|=" << gnorm << std::endl;
}

template<typename T>
void quasinewton<T>::runallsteps(){
  FirstEvaluation(); // Evaluate function and gradient for the first time
  befmainloopspecific(); // Depending on whether BFGS or LBFGS is being implemented
                         // this is an abstract function and each class is forced to
                         // implement it its own way.
  printbefmainloop(); // No calculations performed.  Only output
  /*  
      In the next part run as many iterations as it is necessary.  This mainloop will
      be completely different for constrained problems.  It will also call the virtual
      mainloopspecific which is implemented different for BFGS and LBFGS problems
  */
  while(!done){
    it++;
    mainloop();
  }
  // Show the final results in a nice output
  postmainloop();
}

template<typename T>
void quasinewton<T>::get_ti_s(){
  // This function gets all the Ti points described in (4.1) of 8limited**
  // It also sorts them automatically
  T pp = 0.0;
  for(int i = 0; i < n; i++){
    if(0.0 == quasinewton<T>::g[i]){
      // Assign \Infty if g == 0
      
      bpmemory.insert(std::pair<T, int>(std::numeric_limits<T>::max(), i));
    } else {
      if(quasinewton<T>::g[i] < 0) {
	pp = (x[i] - u[i]) / quasinewton<T>::g[i];
	bpmemory.insert(std::pair<T, int>((x[i] - u[i]) / quasinewton<T>::g[i], i));
      } else {
	pp = (x[i] - l[i]) / quasinewton<T>::g[i];
	bpmemory.insert(std::pair<T, int>((x[i] - l[i]) / quasinewton<T>::g[i], i));
      }
    }
  }
  std::cout << pp << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGS: public virtual quasinewton<T>{
protected:
  T *q, *H, *B;
  double* Hdouble, *Bdouble;
  Matrix<double> mHdouble, mBdouble;
public:
  BFGS(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
       std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, int&,
       double*&, double*&);
  virtual ~BFGS();
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
  void createDoubleHandDoubleB();
};

template<typename T>
BFGS<T>::BFGS(T*& x0, T*& fopt0, int& n0,  T& taud0,  
	      int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
	      T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
	      const char *& outputname0, int& m0, int& gradientsamplingN0, 
	      double*& u0, double*& l0):
  quasinewton<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
		 echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0),
  mHdouble(quasinewton<T>::n, quasinewton<T>::n), mBdouble(quasinewton<T>::n, 
							   quasinewton<T>::n){
  quasinewton<T>::n1 = quasinewton<T>::n;
  quasinewton<T>::n2 = quasinewton<T>::n * quasinewton<T>::n;
  quasinewton<T>::nm = 1;
  quasinewton<T>::m1 = 1; // on LBFGS this variable is m (the history)
  q = new T[quasinewton<T>::n1];
  H = new T[quasinewton<T>::n2];
  B = new T[quasinewton<T>::n2];
  Hdouble = new double[quasinewton<T>::n * quasinewton<T>::n];
  Bdouble = new double[quasinewton<T>::n * quasinewton<T>::n];
}

template<typename T>
BFGS<T>::~BFGS(){
  delete [] q;
  delete [] H;
  delete [] Hdouble;
  delete [] B;
  delete [] Bdouble;
}

template<typename T>
void BFGS<T>::befmainloopspecific(){
  mat_set_eye(H, quasinewton<T>::n, quasinewton<T>::n);
  mat_set_eye(B, quasinewton<T>::n, quasinewton<T>::n);

  // This next variable doesn't need to be set during a run of constrained problems
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
void BFGS<T>::createDoubleHandDoubleB(){
  for(int i = 0; i < this->n; i++){ //left this->n here on purpose for whoever is
                                    //reading.  This way you don't know the scope!
    for(int j = 0; j < quasinewton<T>::n; j++){
      Hdouble[i * quasinewton<T>::n + j] = t_double(H[i * quasinewton<T>::n + j]);
      Bdouble[i * quasinewton<T>::n + j] = t_double(B[i * quasinewton<T>::n + j]);
    }
  }
  Matrix<double> temp(Hdouble, quasinewton<T>::n, quasinewton<T>::n);
  mHdouble = temp;
  Matrix<double> temp2(Bdouble, quasinewton<T>::n, quasinewton<T>::n);
  mBdouble = temp2;
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGSB: public virtual BFGS<T>{
protected:
  double tstar;
  char yTrans, nTrans;
  double alpha, beta, tj, fpj, fppj, deltatj, oldtj, adouble, dtstar, alpha0;
  int ndouble, one, b;
  double *di, *z, *C, *dnsize;
  Matrix<double> mZ, mdi;
  // bear in mind that the following multimap is already ordered
  typename std::multimap<T, int>::iterator iter;
public:
  BFGSB(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	int&, double*&, double*&);
  virtual ~BFGSB();
  void dBd();
  void dBz();
  void create_d();
  void update_d();
  void zeroethstep();
  virtual void nextIterationPrepare();
  virtual void lapackzerostep();
  void findXCauchymX(int);
  virtual void lapackmanipulations();
  void tstarcalculation();
  virtual void findGeneralizedCauchyPoint();
  void findMinimum2ndApproximation();
  void mainloop();
  void printFinalConditions();
  void prepareNextMainLoop();
  void exitSignalHandling();
  virtual void updatec(int);
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
  alpha0 = adouble = oldtj = deltatj = fppj = fpj = tj = dtstar = tstar = 0.0;
  yTrans = 'T'; nTrans = 'N';
  alpha = 1.0; beta = 0.0;
  ndouble = quasinewton<T>::n;
  one = 1;
  di = new double[quasinewton<T>::n];
  z = new double[quasinewton<T>::n];
  C = new double[quasinewton<T>::n];
  dnsize = new double[quasinewton<T>::n];
  for(int _i = 0; _i < quasinewton<T>::n; _i++){
    z[_i] = 0.0;
    di[_i] = 0.0;
    C[_i] = 0.0;
    dnsize[_i] = 0.0;
  }
}

template<typename T>
BFGSB<T>::~BFGSB(){
  delete [] di;
  delete [] z;
  delete [] C;
  delete [] dnsize;
}

template<typename T>
void BFGSB<T>::update_d(){
  //Do not move in the direction that reached the boundary from now on
    di[b] = 0.0;
    Matrix<double> temp2(di, quasinewton<T>::n, 1);
    mdi = temp2;
}

template<typename T>
void BFGSB<T>::create_d(){
  Matrix<double> temp(di, quasinewton<T>::n, 1);
  mdi = temp;
}

// friends of bfgs
template<typename T>
void BFGSB<T>::dBz(){
  adouble = squareForm(mdi, BFGS<T>::mBdouble, mZ);
}

template<typename T>
void BFGSB<T>::dBd(){
  adouble = squareForm(mdi, BFGS<T>::mBdouble, mdi);
}

template<typename T>
void BFGSB<T>::nextIterationPrepare(){
  update_d();
  quasinewton<T>::freeVariable[b] = false;
  oldtj = tj; // Because the next iteration tj will move one step to the front
}

template<typename T>
void BFGSB<T>::zeroethstep(){
  /*
    First of all.  Run the zeroeth step from the multistep gradient projection
    when you exit this function. Xcauchy will have the first value after the gradient
    hits a boundary. d will basically contain the same value as g and freeVariable
    will register the value of the dimension of the boundary that was hit the first 
    time.  Z will contain the variation between the two points (same as g for the 1st
    step.

    The only reason why mainloop has a zero step.  Is because after the first
    calculation you can save some time byusing the update formulae on (4.9) and (4.10)
  */
  
  // the first iter.  Will contain the information of the time when we first hit the 
  // boundary
  FLAG();
  iter = quasinewton<T>::bpmemory.begin();
  b = iter->second;
  deltatj = t_double(iter->first); // Change from zero.  "time" to first boundary

  // keep track of the ti_s that define the interval we are working on
  oldtj = 0.0;
  tj = deltatj;
  
  // Find the new x position.  Notice that in this case all the coordinates advance
  // (at least those with non-zero values in the gradient)
  // given that nothing will hit the boundary (until you hit the boundary corresponding
  // to dimension 'b' of course.
  for(int _i = 0; _i < quasinewton<T>::n; _i++){
    quasinewton<T>::xcauchy[_i] = (t_double(quasinewton<T>::x[_i]) - 
				   tj * t_double(quasinewton<T>::g[_i]));
    z[_i] = t_double(quasinewton<T>::xcauchy[_i] - quasinewton<T>::x[_i]);
  }
  
  Matrix<double> temp(z, quasinewton<T>::n, 1);
  mZ = temp;
  
  // Create di but Update new d_i coordinate.  this step it is basically just -gradient
  // This is from formula (4.2) in the paper
  for(int i = 0; i < quasinewton<T>::n; i++){
    di[i] = -1.0 * t_double(quasinewton<T>::g[i]);
  }
  
  create_d();
  
  // Run lapack intensive calculations
  lapackzerostep();
  
  // Prepare everything for the next iteration
  nextIterationPrepare();
  
  // check if I can finish now and graciously leave if that's the case.
  if (tstar < tj){
    if (tstar > 0){
      std::cout << "optimal value was found" << std::endl;
      this->thisIterationConverged = true; // you found the generalized cauchy point
      tj = tstar;
      for(int i = 0; i < quasinewton<T>::n; i++){
	  this->xcauchy[i] = this->xcauchy[i] - (tj - dtstar) * mdi(i); //stp bck a ltl
      }
    }
  }
}

template<typename T>
void BFGSB<T>::lapackzerostep(){
  /*
    This method includes all the lapack routines that have to be run at the beginning
    of the "find the cauchy point iteration"
  */

  // Calculation of variable fpj (page 6 of Nocedal's paper. Equation 4.4)
  // fpj =  g^Td + z^T*B*z
  dBz(); // this function modifies adouble (which is zero until now)
  fpj = adouble;
  for(int _i = 0; _i < this->n; _i++){
    fpj = fpj +  t_double(this->g[_i]) * di[_i]; //grad[_i];
  }
  
  // Calculation of variable fppj (page 6 of Nocedal's paper. Equation 4.5)
  // fppj = d^T*B*z
  dBd();
  fppj = adouble;
  
  tstarcalculation();  
  // garbage collection
  // delete [] grad;
}

template<typename T>
void BFGSB<T>::findXCauchymX(int i){
  // Calculates the difference between the new "cauchy" point and X one coordinate at a
  // time.  Of course it also calculates the new "cauchy" point
  quasinewton<T>::xcauchy[i] = quasinewton<T>::xcauchy[i] + deltatj * di[i];
  // this was the equation right before (4.2)
  
  // Just a revision of where I am (in reality. This should only work at the epsilon
  // level.  It might not be necessary at all. Consider erasing later).
  quasinewton<T>::xcauchy[i] = MIN(quasinewton<T>::xcauchy[i], 
				   quasinewton<T>::u[i]);
  quasinewton<T>::xcauchy[i] = MAX(quasinewton<T>::xcauchy[i], 
				   quasinewton<T>::l[i]);

  // zeta As in Equation (4.3) on the paper
  z[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
}
/*
template<typename T>
void BFGSB<T>::lapackmanipulations(){
  // LAPACK manipulations for each of the loops in the xcauchy calculations

  Matrix<double> temp(1, quasinewton<T>::n);
  matrixMultiply(mZ, BFGS<T>::mBdouble, temp, 'T', 'N');
  // Next step uses formula (4.9) page 7 of the paper.
  fpj = fpj + deltatj * fppj + std::pow(t_double(quasinewton<T>::g[b]), 
					2) + 
    t_double(quasinewton<T>::g[b]) * temp(b);

  Matrix<double> temp2(quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mBdouble, mdi, temp2);
  // Reassigning values to C
  for(int i = 0; i < quasinewton<T>::n; i++)
    C[i] = temp2(i);
}
*/

// I'm redefining lapack manipulations

template<typename T>
void BFGSB<T>::lapackmanipulations(){

  /*
   Literally just run the lapackzerostep method.  This is until I fix whether there
   is any change in methodology or not.  If anything... consider merging completely.

    This method includes all the lapack routines that have to be run at the beginning
    of the "find the cauchy point iteration"
  */

  /*
  double * grad;
  grad = new double[this->n];
  
  for(int _i = 0; _i < this->n; _i++){
    grad[_i] = t_double(this->g[_i]) * di[_i];
  }
  */

  // Calculation of variable fpj (page 6 of Nocedal's paper. Equation 4.4)
  // fpj =  g^Td + z^T*B*z
  dBz(); // this function modifies adouble (which is zero until now)
  fpj = adouble;
  for(int _i = 0; _i < this->n; _i++){
    fpj = fpj + t_double(this->g[_i]) * di[_i]; //grad[_i];
  }
  //delete [] grad;

  // Calculation of variable fppj (page 6 of Nocedal's paper. Equation 4.5)
  // fppj = d^T*B*z
  dBd();
  fppj = adouble;
  tstarcalculation();
}

template<typename T>
void BFGSB<T>::tstarcalculation(){
  // find optimal point dtstar and tstar.  From last paragraph on page 6 of the paper
  dtstar = -fpj / fppj;
  tstar = dtstar + oldtj;
}

template<typename T>
void BFGSB<T>::updatec(int){
  // do nothing
}

template<typename T>
void BFGSB<T>::findGeneralizedCauchyPoint(){
  iter++;  // this is a class member.  Starts at t_1 and the first time it moves to t_2
  for(; iter != quasinewton<T>::bpmemory.end(); iter++){
    tj = t_double(iter->first);
    b = iter->second;
    deltatj = tj - oldtj;

    SHOW(deltatj);

    // update xcauchy and update the new z (the array, not the Matrix<double>)
    T addDeviations = 0.0;
    FLAG();
    PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
    for(int i = 0; i < quasinewton<T>::n; i++){
      addDeviations = addDeviations + std::abs(deltatj * di[i]);
      findXCauchymX(i);
    }
    PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
    // Create Matrix<double> material
    Matrix<double> temp(z, quasinewton<T>::n, 1);
    mZ = temp;
    
    lapackmanipulations();
    
    if (tstar >= oldtj){
      if (tstar <= tj){
	std::cout << "found optimal cauchy point." << std::endl;
	this->thisIterationConverged = true;
	for(int i = 0; i < quasinewton<T>::n; i++){
	  this->xcauchy[i] = this->xcauchy[i] - (tj - dtstar) * mdi(i);//stpbck a lttl
	  updatec(i);
	}
	return;
      }
    }
    SHOW(addDeviations)
    if ((fpj >= 0) && (addDeviations > quasinewton<T>::taud)){
      std::cout << "the sum of deviations is: " << addDeviations << std::endl;
      std::cout << "found optimal cauchy point exactly at a breakpoint." << std::endl;
      tstar = oldtj;
      this->thisIterationConverged = true;
      for(int i = 0; i < quasinewton<T>::n; i++){
	this->xcauchy[i] = this->xcauchy[i] - (tj - oldtj) * mdi(i); // previous 
	                                                       // xcauchy was
	                                                       // the right point 
      }
      return;
    }
    nextIterationPrepare();
  }
  
  // In case nothing was found.  Return the last point
  tstar = tj;
  /*
  for(int i = 0; i < quasinewton<T>::n; i++){
    this->x[i] = this->xcauchy[i];// previous xcauchy was
    // the right point
    }*/
  // I still need to implement the last segment to locate xcauchy correctly.
  std::cout<< "exiting after finding generalized cauchy point (not optimal)"<<std::endl;
  std::cout << "You are in a corner at this point" << std::endl;
}



template<typename T>
void BFGSB<T>::prepareNextMainLoop(){
  // Set up everything for the next phase.  Calibration of B
  // These s and x will be substracting from xfinal later
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  //PRINTARRAY(quasinewton<T>::s, quasinewton<T>::n, 1);

  // See if we moved in this step.  If we didn't.  Finish with success
  double normminstep;
  double * difference = new double[quasinewton<T>::n];
  
  //PRINTARRAY(quasinewton<T>::x, quasinewton<T>::n, 1);
  PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
  PRINTARRAY(dnsize, quasinewton<T>::n, 1);
  for(int __i = 0; __i < quasinewton<T>::n; __i++){
    // Update the new position of x  this is the only point where this happens
    quasinewton<T>::xcauchy[__i] += dnsize[__i];
    difference[ __i ] = quasinewton<T>::xcauchy[ __i ] - 
      t_double(quasinewton<T>::x[ __i ]);
    quasinewton<T>::x[ __i ] = quasinewton<T>::xcauchy[ __i ];
  }
  
  PRINTARRAY(quasinewton<T>::x, quasinewton<T>::n, 1);
  PRINTARRAY(difference, quasinewton<T>::n, 1);
  normminstep = vecnorm(difference, quasinewton<T>::n);
  
  SHOW(normminstep);
  if(normminstep < quasinewton<T>::taud)
    *quasinewton<T>::exitflag = 1;
  
  delete [] difference;
  
  // Update the value of the function
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::x, quasinewton<T>::n);
  
  vpv<T>(quasinewton<T>::s, quasinewton<T>::x, 1, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1, quasinewton<T>::n);
  
  SHOW(quasinewton<T>::it);
  SHOW(*quasinewton<T>::f);
  
  PRINTARRAY(quasinewton<T>::x, quasinewton<T>::n, 1);
  
  // Call the new Hessian
  update_bfgs_B<T>(BFGS<T>::B, quasinewton<T>::s, quasinewton<T>::y, BFGS<T>::q,
		   quasinewton<T>::n);
}

template<typename T>
void BFGSB<T>::exitSignalHandling(){
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
  SHOW(*quasinewton<T>::exitflag);
  if (0 != *quasinewton<T>::exitflag) 
    quasinewton<T>::done = true;
}

template<typename T>
void BFGSB<T>::findMinimum2ndApproximation(){
  // Assuming xcauchy has been correctly found.  This function runs a minimization of
  // the quadratic approximation to the goal function
  int numfree = 0; // number of free variables
  b = 0;
  double * ZfM2, * r, * dx;
  
  dx = new double[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    // initialize the previous vectors just because...
    dnsize[i] = dx[i] = 0.0;
    if(quasinewton<T>::freeVariable[i])
      numfree++;
  }
  
  if(0 == numfree){
    // no free variables to work with.  This step can't be completed
    return;
  }
  // Construction of Matrix Z 
  ZfM2 = new double[quasinewton<T>::n * numfree];
  for(int i = 0; i < (quasinewton<T>::n * numfree); i++)
    ZfM2[i] = 0.0;
  
  int i_ = 0;
  PRINTARRAY(quasinewton<T>::freeVariable, quasinewton<T>::n, 1);
  for(iter = this->bpmemory.begin(); iter != quasinewton<T>::bpmemory.end(); iter++){
    b = (*iter).second; //position of the ith. crossed boundary
    // SHOW(b);
    if(quasinewton<T>::freeVariable[b]){
      // ZfM2 is a n x numfree matrix populated column-wise
      ZfM2[b + quasinewton<T>::n * i_] = 1.0; // fill with ones for free variables as 
                                              // explained on paragraph 2 of page 10 of 
                                              // the paper.
      ++i_;
    }
  }
  
  // Definition of r vector.  r = Z(g + B(Xcauchy - X))
  // is it maybe worth representing r as in equation 5.4 instead?
  for(int i = 0; i < quasinewton<T>::n; i++)
    dx[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
  
  Matrix<double> mdx(dx, quasinewton<T>::n, 1), mC(C, quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mBdouble, mdx, mC); // Result kept in mC 
  
  BFGS<T>::mBdouble.print('B');
  
  mC.print('c');
  PRINTARRAY(quasinewton<T>::g, quasinewton<T>::n, 1);
  for(int i = 0; i < quasinewton<T>::n; i++){
    C[i] = mC(i);  //Warning!.  I need to correct for this double assignation
    C[i] += t_double(quasinewton<T>::g[i]);
  }
  r = new double[numfree];
  // FLAG();
  // PRINTARRAY(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mZfM2(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mmC(C, quasinewton<T>::n, 1);
  Matrix<double> mr(r, numfree, 1);
  matrixMultiply(mZfM2, mmC, mr, 'T', 'N'); // the result is now on mr;
  mmC.print('C');
  mr.print('r');
  // Find Bhat = Z^TBZ
  Matrix<double> mBHAT(numfree, numfree);
  GensquareForm(mZfM2, BFGS<T>::mBdouble, mZfM2, mBHAT);
  
  /*
    Solve the system 5.5 and 5.6
    Notice that this system could easily be solved by inverting the matrix *BHAT
  */
  
  // First step. Solve the system with 
  Matrix<double> md(numfree, 1);  // Where to put the solution 
  bfgssolver(mBHAT, mr, md);
  
  double alphacandidate = 0.0;
  
  // Define the new boundaries which appear on 5.6
  double lbf, ubf;
  alpha0 = 1.0;
  
  Matrix<double> mdtemp(quasinewton<T>::n, 1);
  matrixMultiply(mZfM2, md, mdtemp); // the result is now on mr;
  mdtemp.print('D');
  for(iter = quasinewton<T>::bpmemory.begin(); 
      iter != quasinewton<T>::bpmemory.end(); iter++){
    b = (*iter).second;
    
    // only perform this analysis for free variables
    if(quasinewton<T>::freeVariable[b]){
      SHOW(b);
      lbf = quasinewton<T>::l[b] - quasinewton<T>::xcauchy[b];
      ubf = quasinewton<T>::u[b] - quasinewton<T>::xcauchy[b];
      
      if(mdtemp(b) > 0){
	alphacandidate = ubf / mdtemp(b);
      } else if(md(b) < 0){
	//both numbers are negative => div. is positive
	alphacandidate = lbf / mdtemp(b);
      } else{
	alphacandidate = 1.0;
      }
      SHOW(alphacandidate);
      alpha0 = MIN(alpha0, alphacandidate);
    }
  }
  
  md.print('d');
  mZfM2.print('Z');
  // if alpha == 1.0 that means the solution doesn't touch any constraint :)
  md *= alpha0;
  SHOW(alpha0);
  // Calculate the new solution
  iter = quasinewton<T>::bpmemory.begin();
  // Z_k * d
  
  Matrix<double> mdnsize(dnsize, quasinewton<T>::n, 1);
  //FLAG();
  matrixMultiply(mZfM2, md, mdnsize);
  
  for(int i = 0; i < quasinewton<T>::n; i++)
    dnsize[i] = mdnsize(i);
  
  PRINTARRAY(dnsize, quasinewton<T>::n, 1);
  
  // Delete memory
  // delete [] dnsize;
  delete [] dx;
  delete [] ZfM2;
  delete [] r;
}

template<typename T>
void BFGSB<T>::printFinalConditions(){
    std::cout << "values of vector g" << std::endl;

  for(int i = 0; i < quasinewton<T>::n; i++)
    std::cout << quasinewton<T>::g[i] << std::endl;

}

template<typename T>
void BFGSB<T>::mainloop(){
  
  /*
    This is the function that does all of the heavy lifting in BFGSB.  It calculates
    every single iteration.  Please notice that the projection step of BFGSB requires
    a few more iterations which here are run inside
  */
  
  quasinewton<T>::fprev = *quasinewton<T>::f;
  quasinewton<T>::get_ti_s();//calculate the ti_s from (4.1) in paper.Define breakpoints
  BFGS<T>::createDoubleHandDoubleB();//Convert stuff to double precision and create 
                                     // Matrix<double> container.
  this->thisIterationConverged = false;  //Once true leave cauchy iteration
  zeroethstep();
  
  // Only continue if this iteration has not converged
  if(!this->thisIterationConverged)
    findGeneralizedCauchyPoint(); // If this function converges inside it has to exit
  PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
  // Beginning of the optimization part
  findMinimum2ndApproximation();
  
  prepareNextMainLoop();
  exitSignalHandling();
  
  std::cout << "Print final conditions to see if we converged" << std::endl;
  printFinalConditions();
  
  // Clean class memory variables
  quasinewton<T>::bpmemory.clear();
  
  if (quasinewton<T>::it >= quasinewton<T>::maxit)
    *quasinewton<T>::exitflag = -1;
  if (*quasinewton<T>::f < quasinewton<T>::ftarget)
    *quasinewton<T>::exitflag = 1;
  if (quasinewton<T>::gnorm < quasinewton<T>::gnormtol)
    *quasinewton<T>::exitflag = 2;
  if (*quasinewton<T>::qpoptvalptr < quasinewton<T>::taud)
    *quasinewton<T>::exitflag = 7;
  
  /* if exitflag was changed: exit main loop: */
  if (*quasinewton<T>::exitflag != 0) quasinewton<T>::done = true;
  // deletion
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class LBFGS: public virtual quasinewton<T>{
protected:
  T *S, *Y, *rho, *a;
public:
  LBFGS(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	int&, double*&, double*&);
  virtual ~LBFGS();
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

  for(int i = 0; i < quasinewton<T>::nm; i++){
    S[i] = Y[i] = 0.0;
  }
  
  for(int i = 0; i < quasinewton<T>::m1; i++){
    rho[i] = 0.0;
    a[i] = 0.0;
  }
}

template<typename T>
LBFGS<T>::~LBFGS(){
  delete [] S;
  delete [] Y;
  delete [] rho;
  delete [] a;
}

template<typename T>
void LBFGS<T>::befmainloopspecific(){
  // p = -g (LBFGS)
  vcopyp<T>(quasinewton<T>::p, quasinewton<T>::g, -1.0, quasinewton<T>::n);
}

template<typename T>
void LBFGS<T>::mainloopspecific(){
  /*
    This function calculates the new ol which varies between 1 and m.  All the memory is
    stored on this arrays S and Y which contains the past s and y vectors (m of them).
    Everything is in Nocedal's book lbfgs chapter

    Notice that arrays S and Y are a huge array containing all the m 'y' and 's' vectors
  */
  
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
class LBFGSB: public BFGSB<T>: public LBFGS<T>{
protected:
  Matrix<double> mY(quasinewton<T>::n, quasinewton<T>::m);
  Matrix<double> mS(quasinewton<T>::n, quasinewton<T>::m);
  Matrix<double> Mmatrix(2 * quasinewton<T>::m, 2 * quasinewton<T>::m);
  std::list<std::vector<T>> Ycontainer;
  std::list<std::vector<T>> Scontainer; 
  double* c;
  Matrix<double> mc(2 * quasinewton<T>::m, 1);
  double * pvector, * pvectorbackup;
  double theta;
  int index, currentm;
public:
  LBFGSB(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	 std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	 int&, double*&, double*&);
  virtual ~LBFGSB();
  virtual void createDoubleHandDoubleB();
  virtual void lapackzerostep();
  virtual void nextIterationPrepare();
  virtual void lapackmanipulations();
  virtual void updatec(int);
};

// Constructor
template<typename T>
LBFGSB<T>::LBFGSB(T*& x0, T*& fopt0, int& n0,  T& taud0,  
		  int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
		  T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
		  const char *& outputname0, int& m0, int& gradientsamplingN0,
		  double*& u0, double*&l0):
  LBFGS<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
	  echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0){
  quasinewton<T>::nm = quasinewton<T>::n * quasinewton<T>::m;
  quasinewton<T>::m1 = quasinewton<T>::m;
  S  = new T[quasinewton<T>::nm];
  Y  = new T[quasinewton<T>::nm];
  rho = new T[quasinewton<T>::m1];
  a  = new T[quasinewton<T>::m1];
  c =  new T[2 * quasinewton<T>::n];
  pvector =  new T[2 * quasinewton<T>::m];
  pvectorbackup =  new T[2 * quasinewton<T>::m];
  theta = 1.0;
  index = 0;
  currentm = 0;
  for(int i = 0; i < quasinewton<T>::nm; i++){
    S[i] = Y[i] = 0.0;
  }
  
  for(int i = 0; i < quasinewton<T>::m1; i++){
    rho[i] = 0.0;
    a[i] = 0.0;
  }
  for(int i = 0; i < quasinewton<T>::n; i++){
    c[i] = 0.0;
  }
}

//destructor
template<typename T>
LBFGSB<T>::~LBFGSB(){
  // Most memory will be freed by the corresponding parent classes. 
  delete [] c;
  delete [] pvector;
  delete [] pvectorbackup;
}

// Overload function createdoubleHanddoubleB so that it does nothing (time-waster 
// otherwise
template<typename T>
LBFGSB<T>::createDoubleHandDoubleB(){
  // Do nothing since these matrices are not needed here
}

template<typename T>
void LBFGSB<T>::lapackzerostep(){
  /*
    This method will solve the same problems that the BFGSB version solves, but 
    without using the matrix B
  */
  
  //Prepare first part of the S and Y vectors.
  // While the matrices have dimension n x m.  I only pass here the first dimension
  // for the multiplication
  
  //p is initialized to -g (the paper is wrong saying p = W^Td.We still don't have W  
  vcopyp<T>(quasinewton<T>::p, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  
  fpj = vecip<T>(g, d, quasinewton<T>::n);
  fppj = -theta * fpj - vecip<T>(quasinewton<T>::p, quasinewton<T>::p, 
				 quasinewton<T>::n);
  BFGSB<T>::tstarcalculation();
}

template<typename T>
void LBFGSB<T>::updatec(int i){
  c[i] = c[i] - (tj - dtstar) * pvectorbackup[i];
}

template<typename T>
void LBFGSB<T>::lapackmanipulations(){
  for(int i = 0; i < 2 * quasinewton<T>::m; i++){
    c[i] = c[i] + deltatj * pvector[i];
    mc(i, 1) = c[i];
  }
  
  T * wbt = new T[2 * currentm];
  
  for (int i = 0; i < currentm; i++){
    wbt[i] = Ycontainer[i].at(b);
  }
  
  for (int i = currentm; i < 2 * currentm; i++){
    wbt[i] = theta * Scontainer[i].at(b);
  }
  
  Matrix<double> mpvector(pvector, 2 * quasinewton<T>::m, 1);
  
  Matrix<double> mwbt(wbt, 2 * currentm, 1);
  fpj = fpj + deltatj * fppj + quasinewton<T>::g[b] * quasinewton<T>::g[b] + theta *
    quasinewton<T>::g[b] * z[b] -
    squareFormwithPadding(mwbt, M, mc, 2 * quasinewton<T>::m);
  
  fppj = fppj - theta * quasinewton<T>::g[b] * quasinewton<T>::g[b] - 2 * 
    quasinewton<T>::g[b] * squareFormwithPadding(mwbt, M, mpvector, 2 * 
						 quasinewton<T>::m) - 
    quasinewton<T>::g[b] * quasinewton<T>::g[b] * 
    squareFormwithPadding(mwbt, M, mwbt, 2 * quasinewton<T>::m);
  
  for(int i = 0; i < 2 * quasinewton<T>::m; i++){
    pvectorbackup[i] = pvector[i];
    pvector[i] = pvector[i] + quasinewton<T>::g[b] * wbt[i];
  }
  BFGSB<T>::tstarcalculation();
}

template<typename T>
void LBFGSB<T>::nextIterationPrepare(){
  update_d();
  quasinewton<T>::freeVariable[b] = false;
  oldtj = tj; // Because the next iteration tj will move one step to the front
  // these were the regular steps for BFGSB  ----------------------------------


  // Update Sk, Yk and for Wk
  
  // updating s and y in this step
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  
  for(int __i = 0; __i < quasinewton<T>::n; __i++){
    // Update the new position of xcauchy.  Notice that we only need to do this
    // update in case that we haven't found an optimal tstar.  So we are allowed to
    // use the whole of deltatj as opposed to dtstar
    quasinewton<T>::xcauchy[__i] += deltatj * BFGSB<T>::d[__i];
  }
  
  // Update the value of the function
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::xcauchy, quasinewton<T>::n);
  
  vpv<T>(quasinewton<T>::s, quasinewton<T>::x, 1.0, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1.0, quasinewton<T>::n);
  
  // next step is to update s and y inside column index in the matrices
  std::vector<T> svector(quasinewton<T>::s, s + sizeof(quasinewton<T>::s) / 
			 sizeof(quasinewton<T>::s[0]));
  std::vector<T> yvector(quasinewton<T>::y, y + sizeof(quasinewton<T>::y) /
			 sizeof(quasinewton<T>::y[0]));
  Ycontainer.push_front(std::vector<T>());  // pass svector right away?
  Scontainer.push_front(std::vector<T>());
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    Ycontainer[0].push_back(quasinewton<T>::y[i]);
    Scontainer[0].push_back(quasinewton<T>::s[i]);
  }
  
  // if the number of elements is already larger than the limit.  Delete the oldest
  if(quasinewton<T>::m < Ycontainer.size()){
    Ycontainer.pop_back();
    Scontainer.pop_back();
  }
  
  currentm = Ycontainer.size();
  // index = (++index) % quasinewton<T>::m;
  
  // assign the D part of the matrix
  
  Mmatrix.setM(2 * currentm);
  Mmatrix.setN(2 * currentm);
  for(int i = currentm; i > 0; i--){
    T tempval = 0;
    int invi = currentm - i;
    //calculate the dot product Ycontainer[i] * Scontainer[i]
    for(int j = 0; j < quasinewton<T>::n; j++){
      tempval = tempval + Ycontainer[invi].at(j) * Scontainer[invi].at(j);
    }
    Mmatrix(i, i) = -tempval;
  }
  
  // Assign the L matrix
  Matrix<double> Lmatrix(currentm, currentm);
  for (int i = 0; i < currentm; i++){
    for(int j = 0; j < i; j++){ // only for j < i as stated in (3.5)
      T mytemp = 0.0;
      for(int z = 0; z < quasinewton<T>::n; z++){
	// WARNING! Review these.  what if there's not enough history?
        int index1 = currentm - i;  // try and review these 
	int index2 = currentm - j;
	mytemp = mytemp + Scontainer[index1].at(z) * Ycontainer[index2].at(z);
      }
      Lmatrix(i, j) = mytemp;
    }
  }
  
  // Assign Lmatrix to Mmatrix
  Mmatrix.insertMatrix(currentm, 0, 2 * currentm - 1, currentm - 1, Lmatrix);
  // Assign the S^TS matrix

  Matrix<double> Smatrix(currentm, currentm);
  for (int i = 0; i < currentm; i++){
    for(int j = 0; j <= i; j++){
      T mytemp = 0.0;
      for(int z = 0; z < quasinewton<T>::n; z++)
	mytemp = mytemp + Scontainer[i].at(z) * Scontainer[i].at(z);
      Smatrix(i, j) = theta * mytemp;
    }
  }
  
  Mmatrix.insertMatrix(currentm, currentm, 2 * currentm - 1, 2 * currentm - 1, Smatrix);
  
  // Reflect the matrix across the diagonal
  for(int i = 0; i < 2* currentm; i++){
    for(int j = i; j < 2 * currentm; j++){
      Mmatrix(i, j) = Mmatrix(j, i);
    }
  }
  
  // Calculate the inverse.
  Mmatrix.inverse();
}

template<typename T>
void LBFGSB<T>::findMinimum2ndApproximation(){
  // Assuming xcauchy has been correctly found.  This function runs a minimization of
  // the quadratic approximation to the goal function
  int numfree = 0; // number of free variables
  b = 0;
  double * ZfM2, * r, * dx;
  
  dx = new double[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    // initialize the previous vectors just because...
    dnsize[i] = dx[i] = 0.0;
    if(quasinewton<T>::freeVariable[i])
      numfree++;
  }
  
  if(0 == numfree){
    // no free variables to work with.  This step can't be completed
    return;
  }
  // Construction of Matrix Z 
  ZfM2 = new double[quasinewton<T>::n * numfree];
  for(int i = 0; i < (quasinewton<T>::n * numfree); i++)
    ZfM2[i] = 0.0;
  
  int i_ = 0;
  PRINTARRAY(quasinewton<T>::freeVariable, quasinewton<T>::n, 1);
  for(iter = this->bpmemory.begin(); iter != quasinewton<T>::bpmemory.end(); iter++){
    b = (*iter).second; //position of the ith. crossed boundary
    // SHOW(b);
    if(quasinewton<T>::freeVariable[b]){
      // ZfM2 is a n x numfree matrix populated column-wise
      ZfM2[b + quasinewton<T>::n * i_] = 1.0; // fill with ones for free variables as 
                                              // explained on paragraph 2 of page 10 of 
                                              // the paper.
      ++i_;
    }
  }
  // this function relies on matrix W.  So I will have to create it:
  int mnow = Ycontainer.size(); // so many calls to size must be expensive
  Matrix<double> Wmatrix(quasinewton<T>::n, 2 * mnow);
  // Assign values to Wmatrix
  for(int i = 0; i < quasinewton<T>::n; i++){
    for(int j = 0; j < mnow; j++){
      // Fill W two positions at a time
      Wmatrix(i, j) = Ycontainer[i].at(mnow - j - 1);
      Wmatrix(i, j + mnow) = theta * Scontainer[i].at(mnow - j - 1);
    }
  }
  // Form Matrix redgrad: eq. (5.4)
  // First of all calculate (reduced gradient helper) redgrad as WM:
  Matrix<double> redgrad(quasinewton<T>::n, 2 * mnow);
  matrixMultiplywithPadding(Wmatrix, Mmatrix, redgrad, 'N', 'N', quasinewton<T>::n,
			    quasinewton<T>::m, Wmatrix.genN());
  
  // Definition of r vector.  r = Z(g + B(Xcauchy - X))
  // notice that in the LBFGSB case I represent B as (\theta - WM) instead.
    ZfM2[i] = 0.0;
  for(int i = 0; i < quasinewton<T>::n; i++)
    dx[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
  
  Matrix<double> mdx(dx, quasinewton<T>::n, 1), mC(C, quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mBdouble, mdx, mC); // Result kept in mC 
  
  BFGS<T>::mBdouble.print('B');
  
  mC.print('c');
  PRINTARRAY(quasinewton<T>::g, quasinewton<T>::n, 1);
  for(int i = 0; i < quasinewton<T>::n; i++){
    C[i] = mC(i);  //Warning!.  I need to correct for this double assignation
    C[i] += t_double(quasinewton<T>::g[i]);
  }
  r = new double[numfree];
  // FLAG();
  // PRINTARRAY(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mZfM2(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mmC(C, quasinewton<T>::n, 1);
  Matrix<double> mr(r, numfree, 1);
  matrixMultiply(mZfM2, mmC, mr, 'T', 'N'); // the result is now on mr;
  mmC.print('C');
  mr.print('r');
  // Find Bhat = Z^TBZ
  Matrix<double> mBHAT(numfree, numfree);
  GensquareForm(mZfM2, BFGS<T>::mBdouble, mZfM2, mBHAT);
  
  /*
    Solve the system 5.5 and 5.6
    Notice that this system could easily be solved by inverting the matrix *BHAT
  */
  
  // First step. Solve the system with 
  Matrix<double> md(numfree, 1);  // Where to put the solution 
  bfgssolver(mBHAT, mr, md);
  
  double alphacandidate = 0.0;
  
  // Define the new boundaries which appear on 5.6
  double lbf, ubf;
  alpha0 = 1.0;
  
  Matrix<double> mdtemp(quasinewton<T>::n, 1);
  matrixMultiply(mZfM2, md, mdtemp); // the result is now on mr;
  mdtemp.print('D');
  for(iter = quasinewton<T>::bpmemory.begin(); 
      iter != quasinewton<T>::bpmemory.end(); iter++){
    b = (*iter).second;
    
    // only perform this analysis for free variables
    if(quasinewton<T>::freeVariable[b]){
      SHOW(b);
      lbf = quasinewton<T>::l[b] - quasinewton<T>::xcauchy[b];
      ubf = quasinewton<T>::u[b] - quasinewton<T>::xcauchy[b];
      
      if(mdtemp(b) > 0){
	alphacandidate = ubf / mdtemp(b);
      } else if(md(b) < 0){
	//both numbers are negative => div. is positive
	alphacandidate = lbf / mdtemp(b);
      } else{
	alphacandidate = 1.0;
      }
      SHOW(alphacandidate);
      alpha0 = MIN(alpha0, alphacandidate);
    }
  }
  
  md.print('d');
  mZfM2.print('Z');
  // if alpha == 1.0 that means the solution doesn't touch any constraint :)
  md *= alpha0;
  SHOW(alpha0);
  // Calculate the new solution
  iter = quasinewton<T>::bpmemory.begin();
  // Z_k * d
  // FLAG();
  Matrix<double> mdnsize(dnsize, quasinewton<T>::n, 1);
  //FLAG();
  matrixMultiply(mZfM2, md, mdnsize);
  
  for(int i = 0; i < quasinewton<T>::n; i++)
    dnsize[i] = mdnsize(i);
  
  PRINTARRAY(dnsize, quasinewton<T>::n, 1);  
  
  // Delete memory
  // delete [] dnsize;
  delete [] dx;
  delete [] ZfM2;
  delete [] r;
}

/////////////////////////////////////////////////////////////////////////////////////

template<typename T>
bool quasinewton<T>::themin(double *gradpoints, int i_){
  double mindist = 1.0;
  int mybase = i_ * n;
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

  //garbage collection
  delete [] gradpoints;
  delete [] gradgrads;
  delete [] solution;  //Warning.  Why delete so early?
  delete myopt;

  return(true);
}

/* BFGS MAIN ALGORITHM: */
template<class T>
void bfgs(T*& x, T*& fopt, int& n, short& lm, int& m, T& ftarget, 
	  T& gnormtol, int& maxit, long& J, T& taux, T& taud, short& echo, 
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

      delete mybfgs;
    } else {
      BFGS<T>* mybfgs;
      mybfgs = new BFGS<T>(x, fopt, n, taud, testFunction, output, ftarget, gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      mybfgs->runallsteps();
      delete mybfgs;
    }
  } else {
    if(boundedProblem){
      LBFGS<T>* mylbfgs;
      mylbfgs = new LBFGSB<T>(x, fopt, n, taud, testFunction, output, ftarget,gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      mylbfgs->runallsteps();
      delete mylbfgs;
    } else{
      LBFGS<T>* mylbfgs;
      mylbfgs = new LBFGS<T>(x, fopt, n, taud, testFunction, output, ftarget, gnormtol,
			     maxit, echo, lm, outputname, m, gradientsamplingN0, u0,l0);
      mylbfgs->runallsteps();
      delete mylbfgs;
    }
  }
  taux = taux + 1; J++; info[0] = info[0] + 1;
  output.close();
}

#endif // _BFGS_TEMPLATE_HPP_

// ** This is the paper: "A limited memory algorithm for bound constrained optimization"
