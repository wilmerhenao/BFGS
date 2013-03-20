#ifndef _TEST_HPP_
#define _TEST_HPP_

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include "container.hpp"
#include "../testfunctions/functions.hpp"

#define MIN(x, y) ((x) < (y) ? (x) : (y));

void printhowtodoit(); // implemented on the cpp

// This class contains all of the parameters.  It was created as a template class
// in order to initialize all the parameters with the correct type at once

template<typename T>
class algoparameters{
private:
  friend void bfgs<T>(T *, T * , size_t, int, size_t, T, T, 
		       size_t,  long int,  T,  T,  int, 
		      int(*)(T*, T*, T*, size_t),  std::string , double*, size_t);
  T ftarget, gnormtol, taux, taud;
  short echo, lm;
  size_t n, m, maxit, gradientsamplingN;
  std::string datafilename;
  long J;
  T * fopt, * u, *l;
  double * info;
  T * xf;
  double * x0;    
  int(*fun_ptr)(T*, T*, T*, size_t);
  allfunctions<T> * pFunctions;
public:
  algoparameters(size_t, std::string, std::string, short);
  algoparameters(size_t, std::string, std::string, short, double*, double*);
  initializationcode(size_t, std::string, std::string, short);
  ~algoparameters();
  void generateXF();
  void BFGSfunction();
  bool boundedProblem;
};

template<typename T>
algoparameters<T>::initializationcode(size_t k, std::string locfunc, 
				      std::string outstr, short lmprm){
  // All the parameters initialized with the same values
  // It is very possible that we may want to do this depending on the type
  try{
    J = (long)MIN(15, ceil(2.0 * static_cast<double> (n) / 3.0));
  } catch(std::exception ex){
    std::cerr << "issue with casting? " << ex.what() << std::endl;
  }
  try{
    fopt = new T;
    u = new T[n];
    l = new T[n];
    info = new double[4];
    xf = new T[n];
    x0 = new double[n];
  } catch(std::bad_alloc& ex){
    std::cerr << "Issues allocating inside algoparameters constructor: " << 
      ex.what() << std::endl;
    assert(false);
  }
  try{
    pFunctions = new allfunctions<T>;
  } catch(std::bad_alloc& ex){
    std::cerr << "Problem allocating function map in algoparameters constructor" << 
      ex.what() << std::endl;
    assert(false);  
  }
  int asserter = 1;
  asserter = pFunctions->fillMap();
  assert(0 == asserter); ++asserter;
  // The next line declares and initializes an iterator to the map container
  // the keys are the names of the non-smooth functions (std::strings) and the
  // elements pointed to are function pointers corresponding to the function

  typename std::map<std::string, int(*)(T*, T*, T*, size_t), 
		    StringComparerForMap>::iterator it = 
    pFunctions->tMap.find(locfunc);

  if(it != pFunctions->tMap.end()) 
    fun_ptr = it->second;
  else
    std::cerr << "This function was not found!" << std::endl;
  boundedProblem = false;
  generateXF();
}

template<typename T>
algoparameters<T>::algoparameters(size_t k, std::string locfunc, std::string outstr, 
				  short lmprm, double * u0, double * l0):
  ftarget(1e-100), gnormtol(0.0), taux(1e-16), taud(1e-14), echo(2), lm(lmprm), n(k), 
  m(7), maxit(10e8), gradientsamplingN(1000), datafilename(outstr){
  initializationcode(k, locfunc, outstr, lmprm);
  for(size_t i; i < n; i++){
    u[i] = u0[i]; l[i] = l0[i];
  }
  // turn on the boundedProblem flag variable
  boundedProblem = true;
}

template<typename T>
algoparameters<T>::algoparameters(size_t k, std::string locfunc, std::string outstr, 
				  short lmprm):
  ftarget(1e-100), gnormtol(0.0), taux(1e-16), taud(1e-14), echo(2), lm(lmprm), n(k), 
  m(7), maxit(10e8), gradientsamplingN(1000), datafilename(outstr){
  initializationcode(k, locfunc, outstr, lmprm);
}

template<typename T>
algoparameters<T>::~algoparameters(){
  try{
    delete fopt;
    delete [] u;
    delete [] l;
    delete [] info;
    delete [] xf;
    delete [] x0;
  } catch(std::exception ex){
    std::cerr << ex.what() << std::endl;
  }
}

template <typename T> 
void algoparameters<T>::generateXF(){
  srand(static_cast<unsigned> (time(NULL)));
  rand_real_vec<double>(x0, n, -1, 1);
  for (size_t j = 0 ; j < n; j++)
    xf[j] = (T) x0[j];
}

// The next function will be explicitly called by the user.  But if we want
// we can call it from the algoparameters constructor instead, so that it
// runs automatically

template<typename T>
void algoparameters<T>::BFGSfunction(){
  bfgs<T>(xf, fopt, n, lm, m, ftarget, gnormtol, maxit, J, taux, taud, 
	  echo, fun_ptr, datafilename, info, gradientsamplingN);
}

#endif // _TEST_HPP_
