
/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#include "tdouble.hpp"

template<class T>
void bfgs(T*&, T*& fopt, int& n, short& lm, int& m, T& ftarget,  T& gnormtol,  
	  int& maxit,  long& J, T& taux,  T& taud, short& echo, 
	  int(*&testFunction)(T*, T*, T*, int),  std::string& datafilename, 
          double*&, int&, double*&, double*&, bool&);

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
  T * xtemp; // for when you might need a dd_real, qd_real copy of the double.
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
  t = 0;
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
  PRINTARRAY(x, n, 1);
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
class BFGS: public quasinewton<T>{
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
  mHdouble(n0, n0), mBdouble(n0, n0){
  quasinewton<T>::n1 = n0;
  quasinewton<T>::n2 = n0 * n0;
  quasinewton<T>::nm = 1;
  quasinewton<T>::m1 = 1; // on LBFGS this variable is m (the history)
  q = new T[quasinewton<T>::n1];
  H = new T[quasinewton<T>::n2];
  B = new T[quasinewton<T>::n2];
  Hdouble = new double[n0 * n0];
  Bdouble = new double[n0 * n0];
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
class BFGSB: public BFGS<T>{
protected:
  double tstar;
  char yTrans, nTrans;
  double alpha, beta, tj, fpj, fppj, deltatj, oldtj, adouble, dtstar, alpha0;
  int ndouble, one, b;
  double *di, *z, *C, *dnsize, *dibackup;
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
  virtual void zeroethstep();
  virtual void nextIterationPrepare();
  virtual void lapackzerostep();
  void findXCauchymX(int);
  virtual void lapackmanipulations();
  void tstarcalculation();
  virtual void findGeneralizedCauchyPoint();
  virtual void findMinimum2ndApproximation();
  virtual void mainloop();
  void printFinalConditions();
  virtual void prepareNextMainLoop();
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
  mZ(n0, 1), mdi(n0, 1){
  alpha0 = adouble = oldtj = deltatj = fppj = fpj = tj = dtstar = tstar = 0.0;
  yTrans = 'T'; nTrans = 'N';
  alpha = 1.0; beta = 0.0;
  ndouble = n0;
  one = 1;
  di = new double[n0];
  z = new double[n0];
  C = new double[n0];
  dibackup = new double[quasinewton<T>::n];
  dnsize = new double[n0];
  for(int _i = 0; _i < n0; _i++){
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
  T a = 1.0;
  a = a + 1;
  //Keep a backup of di in case we have to come back
  for(int i = 0; i < quasinewton<T>::n; i++)
    dibackup[i] = di[i];

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
  /*for(int _i = 0; _i < quasinewton<T>::n; _i++){
    quasinewton<T>::xcauchy[_i] = (t_double(quasinewton<T>::x[_i]) - 
				   0.0 * t_double(quasinewton<T>::g[_i]));
    z[_i] = t_double(quasinewton<T>::xcauchy[_i] - quasinewton<T>::x[_i]);
  }
  
  Matrix<double> temp(z, quasinewton<T>::n, 1);
  mZ = temp;
  */
  // Create di but Update new d_i coordinate.  this step it is basically just -gradient
  // This is from formula (4.2) in the paper
  for(int i = 0; i < quasinewton<T>::n; i++){
    di[i] = -1.0 * t_double(quasinewton<T>::g[i]);
  }
  
  create_d();
  
  // Run lapack intensive calculations
  lapackzerostep();
  FLAG();
  // Prepare everything for the next iteration
  //nextIterationPrepare();
  FLAG();  
  // check if I can finish now and graciously leave if that's the case.
  if (tstar < tj){
    if (tstar > 0){
      std::cout << "optimal value was found" << std::endl;
      this->thisIterationConverged = true; // you found the generalized cauchy point
      tj = tstar;
      for(int i = 0; i < quasinewton<T>::n; i++){
	this->xcauchy[i] = this->xcauchy[i] + (tstar) * mdi(i); //stp bck a ltl
      }
    }
  }
  FLAG();
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
void BFGSB<T>::updatec(int aa){
  // do nothing
  std::cout << aa << std::endl;
}

template<typename T>
void BFGSB<T>::findGeneralizedCauchyPoint(){
  // this is a class member.  Starts at t_1 and the first time it moves to t_2
  for(iter = quasinewton<T>::bpmemory.begin(); 
      iter != quasinewton<T>::bpmemory.end(); iter++){
    tj = t_double(iter->first);
    b = iter->second;
    deltatj = tj - oldtj;
    
    //SHOW(deltatj);
    
    // update xcauchy and update the new z (the array, not the Matrix<double>)
    T addDeviations = 0.0;
    for(int i = 0; i < quasinewton<T>::n; i++){
      addDeviations = addDeviations + std::abs(deltatj * di[i]);
      //findXCauchymX(i);
    }
    
    // Create Matrix<double> material
    //Matrix<double> temp(z, quasinewton<T>::n, 1);
    //mZ = temp;
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
    
    if ((fpj >= 0) && (addDeviations > quasinewton<T>::taud)){
      std::cout << "the sum of deviations is: " << addDeviations << std::endl;
      std::cout << "found optimal cauchy point exactly at a breakpoint." << std::endl;
      tstar = oldtj;
      FLAG();
      this->thisIterationConverged = true;
      for(int i = 0; i < quasinewton<T>::n; i++){
	this->xcauchy[i] = this->xcauchy[i] - (tj - oldtj) * mdi(i); // previous 
	                                                       // xcauchy was 
	                                                       // the right point 
      }
      FLAG();
      return;
    }
    FLAG();
    nextIterationPrepare();
  }
  FLAG();
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
  
  // See if we moved in this step.  If we didn't.  Finish with success
  double normminstep;
  double * difference = new double[quasinewton<T>::n];
  
  for(int __i = 0; __i < quasinewton<T>::n; __i++){
    // Update the new position of x  this is the only point where this happens
    quasinewton<T>::xcauchy[__i] += dnsize[__i];
    difference[ __i ] = quasinewton<T>::xcauchy[ __i ] - 
      t_double(quasinewton<T>::x[ __i ]);
    quasinewton<T>::x[ __i ] = quasinewton<T>::xcauchy[ __i ];
  }
  
  normminstep = vecnorm(difference, quasinewton<T>::n);
  
  if(normminstep < quasinewton<T>::taud)
    *quasinewton<T>::exitflag = 7;
  
  delete [] difference;
  
  // Update the value of the function
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::x, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::s, quasinewton<T>::x, 1, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1, quasinewton<T>::n);  
  SHOW(quasinewton<T>::it);
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

  for(int i = 0; i < quasinewton<T>::n; i++){
    C[i] = mC(i);  //Warning!.  I need to correct for this double assignation
    C[i] += t_double(quasinewton<T>::g[i]);
  }
  r = new double[numfree];
  // FLAG();

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
  
  // Delete memory
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
  std::cout << "xcauchy vector that is sent to the next step: " << std::endl; 
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
class LBFGS: public quasinewton<T>{
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
class LBFGSB: public BFGSB<T>{
protected:
  Matrix<double> mY;
  Matrix<double> mS;
  Matrix<double> Mmatrix;
  Matrix<double> mc;
  double* c;
  double * pvector, * pvectorbackup;
  double theta;
  int index, currentm;
  int cauchysteps;
  std::list<std::vector<T>> Ycontainer;
  std::list<std::vector<T>> Scontainer;
public:
  LBFGSB(T*& x0, T*& fopt0, int&, T&, int(*&)(T*, T*, T*, int), 
	 std::ofstream&, T&, T&, int&, short&, short&, const char *&, int&, 
	 int&, double*&, double*&);
  virtual ~LBFGSB();
  virtual void initializexcauchy();
  virtual void mainloop();
  virtual void dealwithFreeVariables(double);
  virtual void zeroethstep();
  virtual void findGeneralizedCauchyPoint();
  virtual void createDoubleHandDoubleB();
  virtual void lapackzerostep();
  virtual void nextIterationPrepare();
  virtual void prepareNextMainLoop();
  virtual void lapackmanipulations();
  virtual void updatec(int);
  virtual void findMinimum2ndApproximation();
  virtual void befmainloopspecific();
  virtual void mainloopspecific();
  virtual void updateYS();
  virtual void updatetheta();
  virtual void calculateMmatrix();
};

// Constructor
template<typename T>
LBFGSB<T>::LBFGSB(T*& x0, T*& fopt0, int& n0,  T& taud0,  
		  int(*&tF)(T*, T*, T*, int), std::ofstream& output0,  T& ftarget0,  
		  T& gnormtol0,  int& maxit0, short& echo0, short& lm0, 
		  const char *& outputname0, int& m0, int& gradientsamplingN0,
		  double*& u0, double*&l0):
  BFGSB<T>(x0, fopt0, n0, taud0, tF, output0, ftarget0, gnormtol0, maxit0, 
  	   echo0, lm0, outputname0, m0, gradientsamplingN0, u0, l0),
  mY(n0, m0), 
  mS(n0, m0), 
  Mmatrix(2 * m0, 2 * m0),
  mc(2 * m0, 1){
  //mc(2 * m0, BFGSB<T>::one){
  quasinewton<T>::nm = n0 * m0;
  quasinewton<T>::m1 = m0;
  c =  new double[2 * m0];
  //std::cout << "c has size: " << m0 << std::endl;
  pvector =  new double[2 * m0];
  pvectorbackup =  new double[2 * m0];
  theta = 1.0;
  index = 0;
  currentm = 0;
  int i;
  for(i = 0; i < (2 * n0); i++){
    pvector[i] = c[i] = 0.0;
  }
  cauchysteps = 0;
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
void LBFGSB<T>::createDoubleHandDoubleB(){
  // Do nothing since these matrices are not needed here
  std::cout << "" << std::endl;
}

template<typename T>
void LBFGSB<T>::mainloop(){
  /*
    This is the function that does all of the heavy lifting in BFGSB.  It calculates
    every single iteration.  Please notice that the projection step of BFGSB requires
    a few more iterations which here are run inside
  */
  
  quasinewton<T>::fprev = *quasinewton<T>::f;
  quasinewton<T>::get_ti_s();//calculate the ti_s from (4.1) in paper.Define breakpoints
  BFGS<T>::createDoubleHandDoubleB();//Convert stuff to double precision and create 
                                     // Matrix<double> container.
  initializexcauchy();
  this->thisIterationConverged = false;  //Once true leave cauchy iteration
  if(0 == Ycontainer.size())
    zeroethstep();
  
  // Only continue if this iteration has not converged
  if(!this->thisIterationConverged)
    findGeneralizedCauchyPoint(); // If this function converges inside it has to exit
  std::cout << "xcauchy vector that is sent to the next step: " << std::endl; 
  PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
  // Beginning of the optimization part
  findMinimum2ndApproximation();

  prepareNextMainLoop();
  BFGSB<T>::exitSignalHandling();
  
  std::cout << "Print final conditions to see if we converged" << std::endl;
  BFGSB<T>::printFinalConditions();
  
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


template<typename T>
void LBFGSB<T>::befmainloopspecific(){
  // p = -g (LBFGS)
  vcopyp<T>(quasinewton<T>::p, quasinewton<T>::g, -1.0, quasinewton<T>::n);
}

template<typename T>
void LBFGSB<T>::mainloopspecific(){
  std::cout << "" << std::endl;
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
  
  vcopyp<T>(quasinewton<T>::p, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  
  BFGSB<T>::fpj = t_double(veciptdd<T>(quasinewton<T>::g, BFGSB<T>::di, 
				    quasinewton<T>::n));
  BFGSB<T>::fppj = t_double(-theta * BFGSB<T>::fpj);// - vecip<T>(quasinewton<T>::p, 
						    //	      quasinewton<T>::p, 
						    //	      quasinewton<T>::n));
  BFGSB<T>::tstarcalculation();
}

template<typename T>
void LBFGSB<T>::updatec(int i){
  c[i] = c[i] - (BFGSB<T>::deltatj - BFGSB<T>::dtstar) * pvectorbackup[i];
}

template<typename T>
void LBFGSB<T>::dealwithFreeVariables(double thistime){
  double thistj;
  int thisb;
  for(BFGSB<T>::iter = quasinewton<T>::bpmemory.begin(); 
      BFGSB<T>::iter != quasinewton<T>::bpmemory.end(); BFGSB<T>::iter++){
    thistj = t_double(BFGSB<T>::iter->first);
    thisb = BFGSB<T>::iter->second;
    quasinewton<T>::freeVariable[thisb] = true;
    if(thistime >= thistj)
      quasinewton<T>::freeVariable[thisb] = false;
  }
}

template<typename T>
void LBFGSB<T>::findGeneralizedCauchyPoint(){
  // this is a class member.  Starts at t_1 and the first time it moves to t_2
  ++cauchysteps;
  for(BFGSB<T>::iter = quasinewton<T>::bpmemory.begin(); 
      BFGSB<T>::iter != quasinewton<T>::bpmemory.end(); BFGSB<T>::iter++){
    
    BFGSB<T>::tj = t_double(BFGSB<T>::iter->first);
    BFGSB<T>::b = BFGSB<T>::iter->second;
    BFGSB<T>::deltatj = BFGSB<T>::tj - BFGSB<T>::oldtj;
    if(this->deltatj < this->taud){
      quasinewton<T>::freeVariable[BFGSB<T>::b] = false;
      BFGSB<T>::di[BFGSB<T>::b] = 0.0;
      BFGSB<T>::update_d();
      continue;
    }
    //BFGSB<T>::oldtj = BFGSB<T>::tj; //Because the next iteration tj will move one step 
                                    // to the front
    
    // Don't do anything if deltatj is equal to zero
    if(this->oldtj < this->taud){
      // only do this in order to update p vector
      if(0 < currentm){
	Matrix<double> Wmatrix(quasinewton<T>::n, 2 * currentm);
	Matrix<double> mpvector(2 * currentm, 1);
	// Assign values to Wmatrix
	typename std::list<std::vector<T>>::iterator listityc, listitsc;
	double tempdob; int j__, k;
	
	for(j__ = 0, listityc = Ycontainer.begin(), listitsc = Scontainer.begin();
	    listityc != Ycontainer.end(); j__++, ++listityc, ++listitsc){
	  for(int i = 0; i < quasinewton<T>::n; i++){
	    // Fill W two positions at a time
	    unsigned tempor = static_cast<unsigned>(i);
	    tempdob = t_double((*listityc).at(tempor));
	    Wmatrix(i, j__) = tempdob;
	    k = currentm + j__;
	    Wmatrix(i, k) = theta * 
	      t_double((*listitsc).at(static_cast<unsigned>(i)));
	  }
	}
	double * neggi = new double[quasinewton<T>::n];
	for(int i = 0; i < quasinewton<T>::n; i++){
	  neggi[i] = -quasinewton<T>::g[i];
	}
	Matrix<double> mneggi(neggi, this->n, 1);
	matrixMultiply(Wmatrix, mneggi, mpvector, 'T', 'N');
	for(int i = 0; i < 2 * currentm; i++){
	  pvector[i] = mpvector(i);
	}
      }
    }
    // addDeviations is only for a checkup later to see that this actually moved
    T addDeviations = 0.0;
    for(int i = 0; i < quasinewton<T>::n; i++){
      addDeviations = addDeviations + std::abs(BFGSB<T>::deltatj * BFGSB<T>::di[i]);
    }
    
    nextIterationPrepare();
    lapackmanipulations();
    
    if ((BFGSB<T>::tstar >= BFGSB<T>::oldtj) && (BFGSB<T>::tstar <= BFGSB<T>::tj)){
      if (BFGSB<T>::dtstar <= BFGSB<T>::deltatj){
	std::cout << "found optimal cauchy point." << std::endl;
	this->thisIterationConverged = true;
	for(int i = 0; i < quasinewton<T>::n; i++){ 
          // we go back a little bit
	  this->xcauchy[i] =this->xcauchy[i] - (BFGSB<T>::deltatj - BFGSB<T>::dtstar) * 
	    BFGSB<T>::dibackup[i];//stpbck a lttl
	  updatec(i);
	}
	dealwithFreeVariables(BFGSB<T>::tstar);
	return;
      }
    }
    
    SHOW(addDeviations);
    
    if ((BFGSB<T>::fpj >= 0) && (addDeviations > quasinewton<T>::taud)){
      std::cout << "the sum of deviations is: " << addDeviations << std::endl;
      std::cout << "found optimal cauchy point exactly at a breakpoint." << std::endl;
      BFGSB<T>::tstar = BFGSB<T>::oldtj;
      FLAG();
      this->thisIterationConverged = true;
      for(int i = 0; i < quasinewton<T>::n; i++){
	this->xcauchy[i] = this->xcauchy[i] - (BFGSB<T>::deltatj) * 
	  BFGSB<T>::dibackup[i]; 
        // previous xcauchy was the right point 
      }
      dealwithFreeVariables(BFGSB<T>::oldtj); //look. Here it is tj as opposed to tstar
      
      double newtarget = *quasinewton<T>::f;
      T steplength;
      // Actual step length calculation
      std::cout << "cheating with line search study: " << std::endl;
      T* gtemp = new T[quasinewton<T>::n];
      PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
      steplength = linesearch_ww_lbfgsb<T>(quasinewton<T>::xcauchy, quasinewton<T>::f, 
					   gtemp, BFGSB<T>::dibackup, 0.0001, 0.9, 
					   quasinewton<T>::n, 
					   quasinewton<double>::testFunction, 
					   quasinewton<T>::nfeval, newtarget);
					   (steplength > 0) ? ++steplength: --steplength;
      return;
    }
    
    BFGSB<T>::oldtj = BFGSB<T>::tj;
    BFGSB<T>::update_d();
    quasinewton<T>::freeVariable[BFGSB<T>::b] = false;
  
    if(cauchysteps > 2000){
      return;
    }
  }/*
     for(int i = 0; i < quasinewton<T>::n; i++)
     quasinewton<T>::x[i] = quasinewton<T>::xcauchy[i];
     findGeneralizedCauchyPoint();
   */
  // In case nothing was found.  Return the last point
  BFGSB<T>::tstar = BFGSB<T>::tj;
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
void LBFGSB<T>::lapackmanipulations(){
  //int uno = 1;
  double xbcp;
  if(BFGSB<T>::dibackup[static_cast<unsigned>(BFGSB<T>::b)] > 0){
    xbcp = quasinewton<T>::u[static_cast<unsigned>(BFGSB<T>::b)];
  } else {
    xbcp = quasinewton<T>::l[static_cast<unsigned>(BFGSB<T>::b)];
  }
  double zb;
  zb = xbcp - quasinewton<T>::x[static_cast<unsigned>(BFGSB<T>::b)];
  if(currentm > 0){
    double * wbt = new double[2 * currentm];
  
    typename std::list<std::vector<T>>::iterator listit;
    int ind = 0;
    for (listit = Ycontainer.begin(); listit != Ycontainer.end(); ++listit){
      wbt[ind] = t_double(listit->at(static_cast<unsigned>(this->b)));
      ind++;
    }

    ind = 0;
    for (listit = Scontainer.begin(); listit != Scontainer.end(); ++listit){
      wbt[ind + currentm] = t_double(theta * 
				     (*listit).at(static_cast<unsigned>(BFGSB<T>::b)));
      ind++;
    }

    for(int i = 0; i < 2 * currentm; i++){
      pvectorbackup[i] = t_double(pvector[i]);
      pvector[i] = pvector[i] + t_double(quasinewton<T>::g[BFGSB<T>::b] * wbt[i]);
    }  
    Matrix<double> mpvector(pvector, 2 * currentm, 1);
    
    for(int i = 0; i < 2 * currentm; i++){
      //std::cout << "I changed!!!" << std::endl;
      c[i] = c[i] + BFGSB<T>::deltatj * pvector[i];
      mc.setPositionbyForce(i, c[i]);
    }
    mc.setM(2 * currentm);
    Matrix<double> mwbt(wbt, 2 * currentm, 1);
    double sqFormResult = squareFormwithPadding(mwbt, Mmatrix, mc);
    BFGSB<T>::fpj = t_double(BFGSB<T>::fpj + BFGSB<T>::deltatj * BFGSB<T>::fppj + 
			     quasinewton<T>::g[BFGSB<T>::b] * 
			     quasinewton<T>::g[BFGSB<T>::b] + theta * 
			     quasinewton<T>::g[BFGSB<T>::b] * zb - 
			     quasinewton<T>::g[BFGSB<T>::b] * sqFormResult);
  
    //squareFormwithPadding(mwbt, Mmatrix, mpvector, 2 * quasinewton<T>::m);
  
    BFGSB<T>::fppj = t_double(BFGSB<T>::fppj - theta * quasinewton<T>::g[BFGSB<T>::b] * 
			      quasinewton<T>::g[BFGSB<T>::b] - 2 * 
			      quasinewton<T>::g[BFGSB<T>::b] * 
			      squareFormwithPadding(mwbt, Mmatrix, mpvector) - 
			      quasinewton<T>::g[BFGSB<T>::b] * 
			      quasinewton<T>::g[BFGSB<T>::b] * 
			      squareFormwithPadding(mwbt, Mmatrix, mwbt));
  }

  else if(0 == currentm){
    BFGSB<T>::fpj = t_double(BFGSB<T>::fpj + BFGSB<T>::deltatj * BFGSB<T>::fppj +
			     BFGSB<T>::g[BFGSB<T>::b] * BFGSB<T>::g[BFGSB<T>::b] +
                             theta * BFGSB<T>::g[BFGSB<T>::b] * zb);
    BFGSB<T>::fppj = t_double(BFGSB<T>::fppj - theta * BFGSB<T>::g[BFGSB<T>::b] *
			     BFGSB<T>::g[BFGSB<T>::b]);
    
  }

  BFGSB<T>::tstarcalculation();
}

template<typename T>
void LBFGSB<T>::initializexcauchy(){
  // Create di but Update new d_i coordinate.  this step it is basically just -gradient
  // This is from formula (4.2) in the paper
  for(int i = 0; i < quasinewton<T>::n; i++){
    BFGSB<T>::di[i] = -1.0 * t_double(quasinewton<T>::g[i]);
  }
  BFGSB<T>::create_d();
  
  // initialize xcauchy
  for(int i = 0; i < quasinewton<T>::n; i++){
    this->xcauchy[i] = this->x[i];
  }
}

template<typename T>
void LBFGSB<T>::zeroethstep(){
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
  BFGSB<T>::iter = quasinewton<T>::bpmemory.begin();
  BFGSB<T>::b = BFGSB<T>::iter->second;
  BFGSB<T>::deltatj = t_double(BFGSB<T>::iter->first); 
  // Change from zero.  "time" to first boundary
  // keep track of the ti_s that define the interval we are working on
  BFGSB<T>::oldtj = 0.0;
  BFGSB<T>::tj = BFGSB<T>::deltatj;
  
  // Run lapack intensive calculations
  lapackzerostep();

  // Prepare everything for the next iteration
  //nextIterationPrepare();

  // check if I can finish now and graciously leave if that's the case.
  if (BFGSB<T>::tstar < BFGSB<T>::tj){
    if (BFGSB<T>::tstar > 0){
      std::cout << "optimal value was found" << std::endl;
      this->thisIterationConverged = true; // you found the generalized cauchy point
      BFGSB<T>::tj = BFGSB<T>::tstar;
      for(int i = 0; i < quasinewton<T>::n; i++){
	this->xcauchy[i] = this->xcauchy[i] + (BFGSB<T>::tstar) * BFGSB<T>::mdi(i); //stp bck a ltl
      }
    }
  }

}

template<typename T>
void LBFGSB<T>::updateYS(){
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);  
  // update the values of x
  // Update the value of the function for the new xcauchy position
  quasinewton<T>::xtemp = new T[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n ; i++){
    quasinewton<T>::xtemp[i] = (T)(quasinewton<T>::xcauchy[i]);
  }
  
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::xtemp, quasinewton<T>::n);
  (*quasinewton<T>::nfeval) = (*quasinewton<T>::nfeval) + 1;
  delete [] quasinewton<T>::xtemp;
  
  vpv<T>(quasinewton<T>::s, quasinewton<T>::xcauchy, 1.0, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1.0, quasinewton<T>::n);  

  typename std::list<std::vector<T>>::iterator listityc, listitsc, listitsc2;  
  // check the condition S^TY > 0  
  if(vecip(quasinewton<T>::s, quasinewton<T>::y, quasinewton<T>::n) > 0 ||
     0 == Ycontainer.size()){
    
    Ycontainer.push_front(std::vector<T>());  
    Scontainer.push_front(std::vector<T>());
  
    listityc = Ycontainer.begin();
    listitsc = Scontainer.begin();
    
    for(int i = 0; i < quasinewton<T>::n; i++){
      (*listityc).push_back(quasinewton<T>::y[i]);
      (*listitsc).push_back(quasinewton<T>::s[i]);
      //std::cout << "y: " << quasinewton<T>::y[i] << std::endl;
    }
  }
  
  // if the number of elements is already larger than the limit.  Delete the oldest
  if((unsigned)quasinewton<T>::m < Ycontainer.size()){
    Ycontainer.pop_back();
    Scontainer.pop_back();
  }
  calculateMmatrix();
  updatetheta();
}

template<typename T>
void LBFGSB<T>::updatetheta(){
  T numerator, denominator;

  denominator = vecip(quasinewton<T>::s, quasinewton<T>::y, quasinewton<T>::n);
  numerator = vecip(quasinewton<T>::y, quasinewton<T>::y, quasinewton<T>::n);

  this->theta = numerator / denominator;
}

template<typename T>
void LBFGSB<T>::calculateMmatrix(){
  currentm = static_cast<int>(Ycontainer.size());
  Mmatrix.setM(2 * currentm);
  Mmatrix.setN(2 * currentm);
  
  typename std::list<std::vector<T>>::reverse_iterator revlistityc, revlistitsc;
  revlistityc = Ycontainer.rbegin();
  revlistitsc = Scontainer.rbegin();
  
  for(int i = currentm; i > 0; i--){
    T tempval = 0.0; 
    //calculate the dot product Ycontainer[i] * Scontainer[i]
    for(int j = 0; j < quasinewton<T>::n; j++){
      tempval = tempval + (*revlistityc).at(static_cast<unsigned>(j)) * 
	(*revlistitsc).at(static_cast<unsigned>(j));
    } 
    if(i > 1){  // this is because otherwise you are running on inexistent memory
      // if(revlistityc != Ycontainer.rend()) if you wish
      ++revlistityc;
      ++revlistitsc;
    }
    int tempi = i - 1;
    Mmatrix(tempi, tempi) = -t_double(tempval);
  }

  // Assign the L matrix
  Matrix<double> Lmatrix(currentm, currentm);
  revlistitsc = Scontainer.rbegin();
  for (int i = 0; i < currentm; i++){
    for(int j = 0; j < i; j++){ // only for j < i as stated in (3.5)
      T mytemp = 0.0;
      revlistityc = Ycontainer.rbegin();
      for(int k = 0; k < quasinewton<T>::n; k++){
	// WARNING! Review these.  what if there's not enough history?
	mytemp = mytemp + (*revlistitsc).at(static_cast<unsigned>(k)) * 
	  (*revlistityc).at(static_cast<unsigned>(k));
      }
      if(j < (currentm - 1)) // check end of list
	++revlistityc;
      Lmatrix(i, j) = t_double(mytemp);
    }
    if(i < (currentm - 1)) //check end of list
      ++revlistitsc;
  }
  // Assign Lmatrix to Mmatrix
  Mmatrix.insertMatrix(currentm, 0, 2 * currentm - 1, currentm - 1, Lmatrix);
  
  // Assign the S^TS matrix
  typename std::list<std::vector<T>>::iterator listityc, listitsc, listitsc2;
  listitsc = Scontainer.begin();  
    
  Matrix<double> Smatrix(currentm, currentm);
  for (int i = 0; i < currentm; i++){
    listitsc2 = Scontainer.begin();
    for(int j = 0; j <= i; j++){ // avoid half the calculations
      T mytemp = 0.0;
      for(int k = 0; k < quasinewton<T>::n; k++)
	mytemp = mytemp + t_double((*listitsc).at(static_cast<unsigned>(k)) * 
				   (*listitsc2).at(static_cast<unsigned>(k)));
      Smatrix(i, j) = t_double(theta * mytemp);
      if(i != j)
	++listitsc2;
    }
    if((i + 1) < currentm)
      ++listitsc;
  }
  
  Mmatrix.insertMatrix(currentm, currentm, 2 * currentm - 1, 2 * currentm - 1, Smatrix);
  
  // Reflect the matrix across the diagonal
  for(int i = 0; i < 2* currentm; i++){
    for(int j = i + 1; j < 2 * currentm; j++){
      Mmatrix(i, j) = Mmatrix(j, i);
    }
  }
  
  // Calculate the inverse.
  Mmatrix.print();
  Mmatrix.matrixInverse();
  Mmatrix.print();
}

template<typename T>
void LBFGSB<T>::nextIterationPrepare(){
  
  /*
    This function does a couple of things.
    1) regular steps from BFGSB update vector d, find the new constrained variable
    2) update xcauchy
    3) find new s and y vectors and add them to the Y and S container lists
    4) Create the matrix 'M' from the paper
  */

  // Update Sk, Yk and for Wk
  
  // updating s and y in this step
  for(int __i = 0; __i < quasinewton<T>::n; __i++){
    // Update the new position of xcauchy.  Notice that we only need to do this
    // update in case that we haven't found an optimal tstar.  So we are allowed to
    // use the whole of deltatj as opposed to dtstar
    quasinewton<T>::xcauchy[__i] = quasinewton<T>::xcauchy[__i] + 
      BFGSB<T>::deltatj * BFGSB<T>::di[__i];
  }
  currentm = static_cast<int>(Ycontainer.size());
  // index = (++index) % quasinewton<T>::m;
  
  if(currentm > 0){
    calculateMmatrix();
  }
  BFGSB<T>::update_d();
  quasinewton<T>::freeVariable[BFGSB<T>::b] = false;
}

template<typename T>
void LBFGSB<T>::prepareNextMainLoop(){
  // Set up everything for the next phase.  Calibration of B
  // These s and x will be substracting from xfinal later
  
  // See if we moved in this step.  If we didn't.  Finish with success
  BFGSB<T>::deltatj = 0.0;
  double normminstep;
  double * difference = new double[quasinewton<T>::n];
  
  for(int __i = 0; __i < quasinewton<T>::n; __i++){
    // Update the new position of x  this is the only point where this happens
    difference[ __i ] = quasinewton<T>::xcauchy[ __i ] - 
      t_double(quasinewton<T>::x[ __i ]);
    quasinewton<T>::x[ __i ] = quasinewton<T>::xcauchy[ __i ];
    quasinewton<T>::freeVariable[__i] = true;
  }
  this->thisIterationConverged = false;
  cauchysteps = 0;
  normminstep = vecnorm(difference, quasinewton<T>::n);
  
  if(normminstep < quasinewton<T>::taud)
    *quasinewton<T>::exitflag = 7;
  
  delete [] difference;
  
  // Update the value of the function
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::x, quasinewton<T>::n);
  for(int __i = 0; __i < currentm; __i++){
    this->c[__i] = 0;
    if(0 == currentm)
      this->pvector[__i] = 0;
  }

  SHOW(quasinewton<T>::it);
}

template<typename T>
void LBFGSB<T>::findMinimum2ndApproximation(){
  // Assuming xcauchy has been correctly found.  This function runs a minimization of
  // the quadratic approximation to the goal function
  
  int numfree = 0; // number of free variables
  double * ZfM2, * r, * dx;
  dx = new double[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    // initialize the previous vectors just because...
    BFGSB<T>::dnsize[i] = dx[i] = 0.0;
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
  for(BFGSB<T>::iter = this->bpmemory.begin(); 
      BFGSB<T>::iter != quasinewton<T>::bpmemory.end(); BFGSB<T>::iter++){
    BFGSB<T>::b = (*BFGSB<T>::iter).second; //position of the ith. crossed boundary
    if(quasinewton<T>::freeVariable[BFGSB<T>::b]){
      // ZfM2 is a n x numfree matrix populated column-wise
      ZfM2[BFGSB<T>::b + quasinewton<T>::n * i_] = 1.0; // fill with ones for free variables as 
                                              // explained on paragraph 2 of page 10 of 
                                              // the paper.
      ++i_;
    }
  }
  
  currentm = static_cast<int>(Ycontainer.size()); // so many calls to size must be 
                                                  //expensive
  r = new double[numfree];  
  Matrix<double> mr(r, numfree, 1); //will hold the full reduced gradient
  for(int i = 0; i < quasinewton<T>::n; i++)
    dx[i] = theta * t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
  double * Cvec = new double[quasinewton<T>::n];  
  int wmatrixsize = MAX(currentm, 1);
  std::cout << "wmatrixsize: " <<wmatrixsize << std::endl;
  Matrix<double> Wmatrix(quasinewton<T>::n, 2 * wmatrixsize);
  
  // Enter only if we have some history to construct matrix W
  if(currentm > 0){
    int k = 0;
    // Assign values to Wmatrix
    typename std::list<std::vector<T>>::iterator listityc, listitsc;
    double tempdob; int j__;
    
    for(j__ = 0, listityc = Ycontainer.begin(), listitsc = Scontainer.begin();
	listityc != Ycontainer.end(); j__++, ++listityc, ++listitsc){
      for(int i = 0; i < quasinewton<T>::n; i++){
	// Fill W two positions at a time
	unsigned tempor = static_cast<unsigned>(i);
	tempdob = t_double((*listityc).at(tempor));
	Wmatrix(i, j__) = tempdob;
	k = currentm + j__;
	Wmatrix(i, k) = theta * 
	  t_double((*listitsc).at(static_cast<unsigned>(i)));
      }
    }
    
    // Form Matrix redgrad: eq. (5.4)
    // First of all calculate (reduced gradient helper) redgrad as WM:
    Matrix<double> redgrad(quasinewton<T>::n, 2 * currentm);
    matrixMultiplywithPadding(Wmatrix, Mmatrix, redgrad, 'N', 'N');
    Matrix<double> redgrad2(quasinewton<T>::n, 1);
    matrixMultiply(redgrad, mc, redgrad2); // Result kept in mC 
    for(int i = 0; i < quasinewton<T>::n; i++){
      Cvec[i] = t_double(quasinewton<T>::g[i]) + dx[i] - redgrad2(i);
    }
  }  
  
  // If we don't have any history to compute matrix W then work without it
  else if(0 == currentm){
    for(int i = 0; i < quasinewton<T>::n; i++){
      Cvec[i] = t_double(quasinewton<T>::g[i]) + dx[i];
    }
  } 

  Matrix<double> mZfM2(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mmC(Cvec, quasinewton<T>::n, 1);
  matrixMultiply(mZfM2, mmC, mr, 'T', 'N'); // the result is now on mr;
  
  // at this point I have r.  Now calculate the steps as they appear on
  // "Direct Primal Method". Page 12 of the paper.
  // Step 1: Zr^c
  Matrix<double> du(numfree, 1);
  if (currentm > 0){
    Matrix<double> zrc(quasinewton<T>::n, 1);
    matrixMultiply(mZfM2, mr, zrc, 'N', 'N');
    // Step 2: v = W^T Zr^c
    Matrix<double> vmatrix(2 * currentm, 1);
    matrixMultiply(Wmatrix, zrc, vmatrix, 'T', 'N');
    // Step 3: v = Mv;
    Matrix<double> vfinal(2 * currentm, 1);
    matrixMultiplywithPadding(Mmatrix, vmatrix, vfinal, 'N', 'N');
    // Step 4: Form N = (I - 1/theta * MW^TZZ^TW)
    // step 4.a Z^TW
    Matrix<double> ztw(numfree, 2 * currentm);
    matrixMultiply(mZfM2, Wmatrix, ztw, 'T', 'N');
    // step 4.b Z^TW
    Matrix<double> wzzw(2 * currentm, 2 * currentm);
    matrixMultiply(ztw, ztw, wzzw, 'T', 'N');
    // step 4.c divide by theta
    for(int i = 0; i < 2 * currentm; i++){
      for(int j = 0; j < 2 * currentm; j++){
	wzzw(i,j) = (-1.0 / theta) * wzzw(i, j);
      }
      wzzw(i, i) = wzzw(i, i) + 1.0;
    }
    // step 5: (Now wzzw is matrix N.  A poor choice of name but it just happened :(  )
    wzzw.matrixInverse();
    Matrix<double> vprev(2 * currentm, 1);
    matrixMultiply(wzzw, vmatrix, vprev);
    // step 6: du
    matrixMultiply(ztw, vprev, du);
    for(int i = 0; i < numfree; i++){
      du(i) = (1 / (theta * theta)) * du(i) + (1 / theta) * mr(i);
    }
  } 

  else if(0 == currentm){
    // this is just in case there is no Ycontainer and no W matrix
    for(int i = 0; i < numfree; i++){
      du(i) = (1 / theta) * mr(i);
    }
  }
  
  // step 7: Find \alpha^* satisfying 
  double alphacandidate = 0.0;
    // Define the new boundaries which appear on 5.6
  double lbf, ubf;
  BFGSB<T>::alpha0 = 1.0;
  
  Matrix<double> dutemp(quasinewton<T>::n, 1);
  matrixMultiply(mZfM2, du, dutemp); // the result is now on original n-dimensional sp
  for(BFGSB<T>::iter = quasinewton<T>::bpmemory.begin(); 
      BFGSB<T>::iter != quasinewton<T>::bpmemory.end(); BFGSB<T>::iter++){
    BFGSB<T>::b = (*BFGSB<T>::iter).second;
    
    // only perform this analysis for free variables
    if(quasinewton<T>::freeVariable[BFGSB<T>::b]){
      lbf = quasinewton<T>::l[BFGSB<T>::b] - quasinewton<T>::xcauchy[BFGSB<T>::b];
      ubf = quasinewton<T>::u[BFGSB<T>::b] - quasinewton<T>::xcauchy[BFGSB<T>::b];
      
      if(dutemp(BFGSB<T>::b) > 0){
	alphacandidate = ubf / dutemp(BFGSB<T>::b);
      } else if(dutemp(BFGSB<T>::b) < 0){
	//both numbers are negative => div. is positive
	alphacandidate = lbf / dutemp(BFGSB<T>::b);
      } else{
	alphacandidate = 1.0;
      }
      BFGSB<T>::alpha0 = MIN(BFGSB<T>::alpha0, alphacandidate);
    }
  }
  
  // if alpha == 1.0 that means the solution doesn't touch any constraint :)
  du *= -1 * BFGSB<T>::alpha0;
  
  // Calculate the new solution
  Matrix<double> dunsize(BFGSB<T>::dnsize, quasinewton<T>::n, 1);
  
  matrixMultiply(mZfM2, du, dunsize);
  T * vecdu = new T[quasinewton<T>::n];
  for(int i = 0; i < quasinewton<T>::n; i++){
    vecdu[i] = dunsize(i);
  }
  // Create a temporary holder for xtemp
  
  quasinewton<T>::xtemp = new T[quasinewton<T>::n];
  T* gtemp = new T[quasinewton<T>::n];
  for(int i = 0; i < quasinewton<T>::n ; i++){
    quasinewton<T>::xtemp[i] = (T)(quasinewton<T>::xcauchy[i]);
  }
  quasinewton<T>::testFunction(quasinewton<T>::f, gtemp, 
			       quasinewton<T>::xtemp, quasinewton<T>::n);
  (*quasinewton<T>::nfeval) = (*quasinewton<T>::nfeval) + 1;
  double newtarget = *quasinewton<T>::f;
  T steplength;
  // Actual step length calculation
  std::cout << "line search study: " << std::endl;
  PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
  steplength = linesearch_ww_lbfgsb<T>(quasinewton<T>::xcauchy, quasinewton<T>::f, 
				       gtemp, vecdu, 0.0001, 0.9, 
				       quasinewton<T>::n, 
				       quasinewton<double>::testFunction, 
				       quasinewton<T>::nfeval, newtarget);
  // Move back to xcauchy
  PRINTARRAY(quasinewton<T>::xcauchy, quasinewton<T>::n, 1);
  //delete [] quasinewton<T>::xtemp;
  SHOW(steplength);
  
  /*
  dunsize *= steplength;
    for(int i = 0; i < quasinewton<T>::n; i++){
    BFGSB<T>::dnsize[i] = dunsize(i);
    quasinewton<T>::xcauchy[i] = quasinewton<T>::xcauchy[i] + dunsize(i);
    }
  */
  
  updateYS();
  // Delete memory
  delete [] dx;
  delete [] ZfM2;
  delete [] r;
  delete [] Cvec;
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
      BFGSB<T>* mylbfgs;
      mylbfgs = new LBFGSB<T>(x, fopt, n, taud, testFunction, output, ftarget,gnormtol, 
			   maxit, echo, lm, outputname, m, gradientsamplingN0, u0, l0);
      //std::cout << "Hola LBFGSB" << std::endl;
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
