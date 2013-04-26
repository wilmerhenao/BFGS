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
  quasinewton(T [], T *, int ,  T,  int(*)(T*, T*, T*, int), std::ofstream&,  
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
    it++
    mainloop();
  }
  // Show the final results in a nice output
  postmainloop();
}

template<typename T>
void quasinewton<T>::get_ti_s(){
  // This function gets all the Ti points described in (4.1) of 8limited**
  // It also sorts them automatically
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
  Matrix<double> temp(Bdouble, quasinewton<T>::n, quasinewton<T>::n);
  mBdouble = temp;
}

/////////////////////////////////////////////////////////////////////////////////////
template<typename T>
class BFGSB: public BFGS<T>{
protected:
  double tstar;
  char yTrans, nTrans;
  double alpha, beta, tj, fpj, fppj, deltatj, oldtj, adouble, dtstar;
  int ndouble, one, b, numberOfIterations;
  double *di, *z, *C;
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
  void zeroethstep();
  void lapackzerostep();
  void findXCauchymX(int);
  void lapackmanipulations();
  void tstarcalculation();
  void findGeneralizedCauchyPoint();
  void findMinimum2ndApproximation();
  void mainloop();
  void printFinalConditions();
  void update_d();
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
  for(int _i = 0; _i < quasinewton<T>::n; _i++){
    z[_i] = 0.0;
    di[_i] = 0.0;
  }
  C = new double[quasinewton<T>::n];
  numberOfIterations = 0;
}

template<typename T>
BFGSB<T>::~BFGSB(){
  delete [] di;
  delete [] z;
  delete [] C;
}

template<typename T>
void BFGSB<T>::update_d(){
    di[b] = 0.0;  // Do not move in the direction that reached the boundary from now on
    Matrix<double> temp2(di, quasinewton<T>::n, 1);
    mdi = temp2;
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
    calculation you can save some time by using the update formulae on (4.9) and (4.10)
  */
  
  // the first iter.  Will contain the information of the time when we first hit the 
  // boundary
  iter = quasinewton<T>::bpmemory.begin();
  b = iter->second;
  deltatj = t_double(iter->first); // Change from zero.  "time" to first boundary

  // keep track of the ti_s that the define the intarval we are working on
  oldtj = 0.0;
  tj = deltatj;

  // Find the new x position.  Notice that in this case all the coordinates advance
  // (at least those with non-zero values in the gradient)
  // given that nothing will hit the boundary (until you hit the boundary corresponding
  // to dimension 'b' of course.
  for(int _i = 0; _i < quasinewton<T>::n; _i++){
    quasinewton<T>::xcauchy[_i] = (t_double(quasinewton<T>::x[_i]) - tj *
				  t_double(quasinewton<T>::g[_i]));
    z[_i] = t_double(quasinewton<T>::xcauchy[_i] - quasinewton<T>::x[_i]);
  }


  Matrix<double> temp(z, quasinewton<T>::n, 1);
  mZ = temp;
  
  // Create di but Update new d_i coordinate.  this step it is basically just -gradient
  // This is from formula (4.2) in the paper
  for(int i = 0; i < quasinewton<T>::n; i++){
    di[i] = -1.0 * t_double(quasinewton<T>::g[i]);
  }
  
  // Run lapack intensive calculations
  lapackzerostep();
}

template<typename T>
void BFGSB<T>::dBz(){
  adouble = squareForm(mdi, BFGS<T>::mBdouble, mZ);
}

template<typename T>
void BFGSB<T>::dBd(){
  adouble = squareForm(mdi, BFGS<T>::mBdouble, mdi);
}

template<typename T>
void BFGSB<T>::lapackzerostep(){
  /*
    This method includes all the lapack routines that have to be run at the beginning
    of the "find the cauchy point iteration"
  */
  
  double * grad;
  grad = new double[this->n];
  
  for(int _i = 0; _i < this->n; _i++){
    grad[_i] = t_double(this->g[_i]) * di[_i];
  }
  
  // Calculation of variable fpj (page 6 of Nocedal's paper. Equation 4.4)
  // fpj =  g^Td + z^T*B*z
  dBz(); // this function modifies adouble (which is zero until now)
  fpj = adouble;
  for(int _i = 0; _i < this->n; _i++){
    fpj = fpj + grad[_i];
  }
  
  // Calculation of variable fppj (page 6 of Nocedal's paper. Equation 4.5)
  // fppj = d^T*B*z
  dBd();
  fppj = adouble;

  tstarcalculation();  
  update_d();
  quasinewton<T>::freeVariable[b] = false;
  oldtj = tj; // Because the next iteration tj will move one step to the front
  // garbage collection
  delete [] grad;
}

template<typename T>
void BFGSB<T>::findXCauchymX(int i){
  // Calculates the difference between the new "cauchy" point and X one coordinate at a
  // time.  Of course it also calculates the new "cauchy" point
  
  quasinewton<T>::xcauchy[i] = quasinewton<T>::xcauchy[i] + deltatj * di[i];

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
  // Literally just run the lapackzerostep method.  This is until I fix whether there
  // is any change in methodology or not.  If anything... consider merging completely
  /*
    This method includes all the lapack routines that have to be run at the beginning
    of the "find the cauchy point iteration"
  */
  
  double * grad;
  grad = new double[this->n];
  
  for(int _i = 0; _i < this->n; _i++){
    grad[_i] = t_double(this->g[_i]) * di[_i];
  }
  
  // Calculation of variable fpj (page 6 of Nocedal's paper. Equation 4.4)
  // fpj =  g^Td + z^T*B*z
  zBz(); // this function modifies adouble (which is zero until now)
  fpj = adouble;
  for(int _i = 0; _i < this->n; _i++){
    fpj = fpj + grad[_i];
  }
  
  delete [] grad;

  // Calculation of variable fppj (page 6 of Nocedal's paper. Equation 4.5)
  // fppj = d^T*B*z
  dBz();
  fppj = adouble;

  tstarcalculation();
}

template<typename T>
void BFGSB<T>::tstarcalculation(){
  // find optimal point dtstar and tstar.  From last paragraph on page 6 of the paper
  dtstar = -fpj / fppj;
  tstar = dtstar + oldtj;
}

/*
template<typename T>
void BFGSB<T>::tstarcalculation(){
  fppj = fppj + 2 * t_double(quasinewton<T>::g[b]) * C[b] + 
    std::pow(t_double(quasinewton<T>::g[b]), 2) * 
    t_double(BFGS<T>::B[b * quasinewton<T>::n + b]);
  dtstar = fpj / fppj;
  tstar = oldtj + dtstar;
}
*/
template<typename T>
void BFGSB<T>::findGeneralizedCauchyPoint(){
  iter++;  // this is a class member.  Starts at t_1 and the first time it moves to t_2
  for(; iter != quasinewton<T>::bpmemory.end(); iter++){
    tj = t_double(iter->first);
    b = iter->second;
    deltatj = tj - oldtj;
    
    // update xcauchy and update the new z (the array, not the Matrix<double>)
    for(int i = 0; i < quasinewton<T>::n; i++){
      findXCauchymX(i);
    }
    
    // Create Matrix<double> material
    Matrix<double> temp(z, quasinewton<T>::n, 1);
    mZ = temp;
    
    lapackmanipulations();
    tstarcalculation();
    if (tstar >= oldtj){
      if (tstar <= tj){
	std::cout << "found optimal cauchy point." << std::endl;
	this->thisIterationConverged = true;
	for(int i = 0; i < quasinewton<T>::n; i++){
	  this->xcauchy[i] = this->xcauchy[i] - (tj - dtstar) * mdi(i);//stpbck a lttl
	}
	return;
      }
    }
    if (fpj >= 0){
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
    update_d();
    oldtj = tj;
    quasinewton<T>::freeVariable[b] = false; 
  }
  
  // In case nothing was found.  Return the last point
  tstar = tj;
  for(int i = 0; i < quasinewton<T>::n; i++){
    this->x[i] = this->xcauchy[i];// previous xcauchy was
    // the right point
  }
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
  dx = new double[quasinewton<T>::n];
  dnsize = new double[quasinewton<T>::n];
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    if(quasinewton<T>::freeVariable[i])
      numfree++;
  }
  
  ZfM2 = new double[quasinewton<T>::n * numfree];
  typename std::multimap<T, int>::iterator titer = quasinewton<T>::bpmemory.begin();
  for(int i = 0; i < quasinewton<T>::n; i++, titer++){
    for(int j = 0; j < numfree; j++)
      ZfM2[i * numfree + j] = 0.0;
  }
  
  int i1 = 0;
  PRINTARRAY(quasinewton<T>::freeVariable, quasinewton<T>::n, 1);
  for(titer = this->bpmemory.begin(); titer != quasinewton<T>::bpmemory.end(); titer++){
    b = (*titer).second; //position of the ith. crossed boundary
    // SHOW(b);
    if(quasinewton<T>::freeVariable[b]){
      // ZfM2 is a n x numfree matrix populated column-wise
      ZfM2[b + quasinewton<T>::n * i1] = 1.0; // fill with ones for free variables as 
                                              // explained on paragraph 2 of page 10 of 
                                              // the paper.
      i1++;
    }
  }

  // Now let's define the r vector.  r = Z(g + B(Xcauchy - X))
  // is it maybe worth representing r as in equation 5.4 instead?
  for(int i = 0; i < quasinewton<T>::n; i++)
    dx[i] = t_double(quasinewton<T>::xcauchy[i] - quasinewton<T>::x[i]);
  
  Matrix<double> mdx(dx, quasinewton<T>::n, 1), mC(C, quasinewton<T>::n, 1);
  matrixMultiply(BFGS<T>::mBdouble, mdx, mC); // Result kept in mC 
  
  for(int i = 0; i < quasinewton<T>::n; i++){
    C[i] = mC(i);  //Warning!.  I need to correct for this double assignation
    C[i] += t_double(quasinewton<T>::g[i]);
  }
  r = new double[numfree];
  // FLAG();
  PRINTARRAY(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mZfM2(ZfM2, quasinewton<T>::n, numfree);
  Matrix<double> mmC(C, quasinewton<T>::n, 1);
  Matrix<double> mr(r, numfree, 1);
  matrixMultiply(mZfM2, mmC, mr, 'T', 'N'); // the result is now on mr;
  // Find Bhat = Z^TBZ
  Matrix<double> mBHAT(numfree, numfree);
  GensquareForm(mZfM2, BFGS<T>::mBdouble, mZfM2, mBHAT);
  
  // Solve the system 5.5 and 5.6
  // Notice that this system could easily be solved by inverting the matrix *BHAT
  
  Matrix<double> md(numfree, 1);  // Where to put the solution 
  bfgssolver(mBHAT, mr, md);
  
  double alpha0 = 0.0;
  double alphacandidate = 0.0;
  
  // Define the new boundaries which appear on 5.6
  double lbf;
  double ubf;
  int ind = 0;
  titer = quasinewton<T>::bpmemory.begin();
  for(; titer != quasinewton<T>::bpmemory.end(); titer++, ind++){
      b = (*titer).second;
      lbf = quasinewton<T>::l[b] - quasinewton<T>::xcauchy[b];
      ubf = quasinewton<T>::u[b] - quasinewton<T>::xcauchy[b];
      alphacandidate = MAX(ubf / md(ind), lbf / md(ind));
      alpha0 = MAX(alpha0, alphacandidate);
    }
  alpha0 = MIN(alpha0, 1.0);
  md *= alpha0;
  // Find the new solution
  
  titer = quasinewton<T>::bpmemory.begin();
  // Z_k * d
  // FLAG();
  Matrix<double> mdnsize(dnsize, quasinewton<T>::n, 1);
  //FLAG();
  matrixMultiply(mZfM2, md, mdnsize);
  
  // Set up everything for the next phase.  Calibration of B
  // These s and x will be substracting from xfinal later
  vcopyp<T>(quasinewton<T>::s, quasinewton<T>::x, -1.0, quasinewton<T>::n);
  vcopyp<T>(quasinewton<T>::y, quasinewton<T>::g, -1.0, quasinewton<T>::n);
  //PRINTARRAY(quasinewton<T>::s, quasinewton<T>::n, 1);

  for(int i = 0; i < quasinewton<T>::n; i++){
    // Update the new position of x
    quasinewton<T>::x[i] = quasinewton<T>::xcauchy[i];
    if((*titer).second == i){
      quasinewton<T>::x[i] += mdnsize(i);
    }
  }
  quasinewton<T>::testFunction(quasinewton<T>::f, quasinewton<T>::g, 
			       quasinewton<T>::x, quasinewton<T>::n);

  vpv<T>(quasinewton<T>::s, quasinewton<T>::x, 1, quasinewton<T>::n);
  vpv<T>(quasinewton<T>::y, quasinewton<T>::g, 1, quasinewton<T>::n);
  numberOfIterations++;
  SHOW(numberOfIterations);

  //PRINTARRAY(quasinewton<T>::s, quasinewton<T>::n, 1);
  //PRINTARRAY(quasinewton<T>::x, quasinewton<T>::n, 1);
  FLAG();
  update_bfgs_B<T>(BFGS<T>::B, quasinewton<T>::s, quasinewton<T>::y, BFGS<T>::q,
		   quasinewton<T>::n);
  FLAG();
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
  if (0 != *quasinewton<T>::exitflag) quasinewton<T>::done = true;
  
  // Delete memory
  delete [] dnsize;
  delete [] dx;
  delete [] ZfM2;
  delete [] r;
  
  //SHOW(numfree);
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
  // check if I can finish now and graciously leave if that's the case.
  if (tstar < tj){
    if (tstar > 0){
      std::cout << "optimal value was found" << std::endl;
      this->thisIterationConverged = true; // you found the generalized cauchy point
      tj = tj + tstar;
      for(int i = 0; i < quasinewton<T>::n; i++){
	  this->xcauchy[i] = this->xcauchy[i] - (tj - dtstar) * mdi(i); //stp bck a ltl
      }
    }
  }
    
  // Only continue if this iteration has not not converged
  if(!this->thisIterationConverged)
    findGeneralizedCauchyPoint(); // If this function converges inside it has to exit
  findMinimum2ndApproximation();
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
  virtual ~LBFGSB();
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

template<typename T>
LBFGSB<T>::~LBFGSB(){
  
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
