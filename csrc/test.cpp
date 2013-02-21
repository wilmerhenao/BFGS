#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <string>
#include "randnums_template.hpp"
#include "bfgs_template.hpp"
#include "../testfunctions/functions.hpp"
#include "libmatrix_template.hpp"
#include "test.hpp"
#include <qd/dd_real.h>

int main(int argc, char *argv[]){    
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
    
  std::string str;
  for(int i = 0; argv[1][i] != 0; i++)
    str += argv[1][i];

  // If you want to add a new type, just:
  // 1) declare it.  
  // 2) Set precision.
  // 3) Run and delete

  algoparameters<double> * doubleparameters = 
    new algoparameters<double>(atoi(argv[3]), str);
  algoparameters<dd_real> * dd_realparameters =
    new algoparameters<dd_real>(atoi(argv[3]), str);

  cout.precision(16);	
  doubleparameters->BFGSfunction();
  cout.precision(32);
  dd_realparameters->BFGSfunction();
    
  delete doubleparameters;
  delete dd_realparameters;
    
  fpu_fix_end(&old_cw);
  return 0;
}
