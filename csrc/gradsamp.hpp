#ifndef _GRADSAMP_HPP_
#define _GRADSAMP_HPP_

#include <cmath>
#include <iostream>
#include "../lib/qpspecial/qpobject.hpp"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

using std::endl;

bool themin(double *, double *, size_t, const size_t);

template<typename T>
bool gradsamp(){
  std::cout << "gradsamp doesn't do anything for the current type.  Try doubles"<< endl;
  //Do nothing unless T is a double type
  return(false);
};


template <> //template specialization for double (only doubles work with lapack)
bool gradsamp<double>(){
  // n is the size of x and gradientsamplingN is the number of gradient samples 
  // assigned in algoparameters.  Gradpoints contains all the points to be evaluated
  size_t j = 0;
  double * gradpoints;
  gradpoints = new[gradientsamplingN * n];
  // Assign original x point to the gradpoints store
  for (size_t i = 0; i < n; i++)
    gradpoints[i] = x[i];

  // Assign random numbers to the rest of variables.
  for (size_t i = n; i < gradientsamplingN * n; i++, j++){
    if(n == j)
      j = 0;
    gradpoints[i] = x[j] + (rand() / RAND_MAX) - 0.5;
  }
  
  // A distance of up to 0.5 is usually too big in any direction.  So I will divide all
  // coordinates by two until the smallest coordinate difference is < 1e-14
  for (size_t i = 1; i < gradientsamplingN; i++){
    do{
      for(size_t k = 0; k < n; k++)
	gradpoints[i + k] = (gradpoints[i + k] + x[i]) /2;
      
    } while (themin(x, gradpoints, i, n))
  }
  
  // found the x's now. Let's find the gradients
  double * gradgrads;
  double f;
  
  gradgrads = new[gradientsamplingN * n];
  for(size_t i = 0; i < gradientsamplingN; i++){
    testFunction(*f, &g[i * n], &gradpoints[i * n], n)
  }
  
  // Next step call qpspecial
  std::cout << "calling qpspecial" << std::endl;
  qpclass<double> * myopt = new qpclass<double>(n, gradientsamplingN, gradgrads, 100);
  myopt->optimization();
  
  double * solution = new double[n];
  myopt->fetchSolution(solution);
  return(true);
}

#endif // _GRADSAMP_HPP_
