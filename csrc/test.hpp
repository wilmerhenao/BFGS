
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

void printhowtodoit(); // implemented on the cpp

// This class contains all of the parameters.  It was created as a template class
// in order to initialize all the parameters with the correct type at once

template<typename T>
class algoparameters{
private:
  friend void bfgs<T>(T*&, T *& , int&, short&, int&, T&, T&, 
		       int&,  long&, T&, T&,  short&, 
		      int(*&)(T*, T*, T*, int),  std::string& , double*&, int&, 
		      double*&, double*&, bool&);
  T ftarget, gnormtol, taux, taud;
  short echo, lm;
  int n, m, maxit, gradientsamplingN;
  std::string datafilename;
  long J;
  T * fopt; 
  double * u, *l;
  double * info;
  T * xf;
  double * x0;
  int(*fun_ptr)(T*, T*, T*, int);
  allfunctions<T> * pFunctions;
  bool boundedProblem;
public:
  algoparameters(int, std::string, std::string, short);
  algoparameters(int, std::string, std::string, short, double*, double*);
  void initializationcode(std::string);
  ~algoparameters();
  void generateXF();
  void BFGSfunction();
};

template<typename T>
void algoparameters<T>::initializationcode(std::string locfunc){
  // All the parameters initialized with the same values
  // It is very possible that we may want to do this depending on the type
  try{
    J = static_cast<long>(MIN(15, ceil(2.0 * static_cast<double>(n) / 3.0)));
  } catch(std::exception ex){
    std::cerr << "issue with casting? " << ex.what() << std::endl;
  }
  try{
    fopt = new T;
    u = new double[n];
    l = new double[n];
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
    //*pFunctions = 0;
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

  typename std::map<std::string, int(*)(T*, T*, T*, int), 
		    StringComparerForMap>::iterator it = 
    pFunctions->tMap.find(locfunc);

  if(it != pFunctions->tMap.end()) 
    fun_ptr = it->second;
  else
    std::cerr << "This function was not found!" << std::endl;

  fopt[0] = 0.0;
  for(int i = 0; i < 4; i++){
    info[i] = 0.0;
  }
  
  for(int i = 0; i < n; i++){
    u[i] = l[i] = 0.0;
    xf[i] = x0[i] = 0.0;
  }

  boundedProblem = false;
}


// There are two initializers.  The first one initializes when there are u, l vectors
template<typename T>
algoparameters<T>::algoparameters(int k, std::string locfunc, std::string outstr, 
				  short lmprm, double * u0, double * l0):
  ftarget(1e-100), gnormtol(0.0), taux(1e-16), taud(1e-14), echo(2), lm(lmprm), n(k), 
  m(7), maxit(10e8), gradientsamplingN(1000), datafilename(outstr){
  initializationcode(locfunc);
  for(int i =0; i < n; i++){
    u[i] = u0[i]; l[i] = l0[i];
  }
  // turn on the boundedProblem flag variable
  boundedProblem = true;
  generateXF();
}

// The second initializer initializes when there are no u, l vectors
template<typename T>
algoparameters<T>::algoparameters(int k, std::string locfunc, std::string outstr, 
				  short lmprm):
  ftarget(1e-100), gnormtol(0.0), taux(1e-16), taud(1e-14), echo(2), lm(lmprm), n(k), 
  m(7), maxit(10e8), gradientsamplingN(1000), datafilename(outstr){
  initializationcode(locfunc);
  // This for loop is just here for consistency and won't be implemented
  for(int i = 0; i < n; i++){
    u[i] = l[i] = 0.0;  // 
  }
  boundedProblem = false;
  generateXF();
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
    delete pFunctions;
  } catch(std::exception ex){
    std::cerr << ex.what() << std::endl;
  }
}

template <typename T> 
void algoparameters<T>::generateXF(){
  /*
    Generate a starting point.  If the problem is bounded, starting point is exactly
    in the middle of the box.  If the problem is unbounded a random starting point is
    applied
  */

  if(boundedProblem){
    for(int j = 0; j < n; j++)
      xf[j] = (u[j] + l[j]) / 2;
  } else {
    srand(static_cast<unsigned> (time(NULL)));
    rand_real_vec<double>(x0, n, -1, 1);
    for (int j = 0 ; j < n; j++)
      xf[j] = (T) x0[j];
  }
}

// The next function will be explicitly called by the user.  But if we want
// we can call it from the algoparameters constructor instead, so that it
// runs automatically

template<typename T>
void algoparameters<T>::BFGSfunction(){
  bfgs<T>(xf, fopt, n, lm, m, ftarget, gnormtol, maxit, J, taux, taud, 
	  echo, fun_ptr, datafilename, info, gradientsamplingN, u, l, boundedProblem);
}

#endif // _TEST_HPP_
