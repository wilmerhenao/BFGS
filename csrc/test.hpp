#ifndef _TEST_HPP_
#define _TEST_HPP_

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "container.hpp"
#include "../testfunctions/functions.hpp"

#define MIN(x, y) ((x) < (y) ? (x) : (y));

// This class contains all of the parameters.  It was created as a template class
// in order to initialize all the parameters with the correct type at once

template<typename T>
class algoparameters{
private:
  friend void bfgs<T>(T *, T * , size_t, int, size_t, T, T, 
		       size_t,  long int,  T,  T,  int, 
		      void(*)(T*, T*, T*, size_t),  std::string , double*);
  T ftarget, gnormtol, taux, taud;
  int echo, lm;
  size_t n, m, maxit;
  long J;
  std::string datafilename;
  T * fopt;
  double * info;
  T * xf;
  double * x0;    
  void(*fun_ptr)(T*, T*, T*, size_t);
  allfunctions<T> * pFunctions;
public:
  algoparameters(size_t, std::string);
  ~algoparameters();
  void generateXF();
  void BFGSfunction();
};

template<typename T>
algoparameters<T>::algoparameters(size_t k, std::string locfunc){

  // All the parameters initialized with the same values
  // It is very possible that we may want to do this depending on the type

  n = k;
  lm = 0;
  m = 7;
  taud = 1e-14;
  taux = 1e-16;
  ftarget = -1e100;
  gnormtol = 0.0;
  maxit = 10e8;
  echo = 2;
  datafilename = "stddump.txt";
  J = (long)MIN(15, ceil(2.0 * static_cast<double> (n) / 3.0));
  fopt = new T;
  info = new double[4];
  xf = new T[n];
  x0 = new double[n];
  pFunctions = new allfunctions<T>;
  pFunctions->fillMap();
  
  // The next line declares and initializes an iterator to the map container
  // the keys are the names of the non-smooth functions (std::strings) and the
  // elements pointed to are function pointers corresponding to the function

  typename std::map<std::string, void(*)(T*, T*, T*, size_t), 
		    StringComparerForMap>::iterator it = pFunctions->tMap.find(locfunc);

  if(it != pFunctions->tMap.end()) 
    fun_ptr = it->second;
  else
    std::cerr << "This function was not found!" << std::endl;

  generateXF();
}

template<typename T>
algoparameters<T>::~algoparameters(){
  // Destructor is probably not necessary, but is implemented here for the future
  std::cout << "Destructor called" << std::endl;
  delete fopt;
  delete [] info;
  delete [] xf;
  delete [] x0;
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
	  echo, fun_ptr, datafilename, info);
}

#endif // _TEST_HPP_