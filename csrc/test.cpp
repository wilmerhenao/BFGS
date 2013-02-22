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

// Executable will return 1 if there are wrong set of arguments

int main(int argc, char *argv[]){
  std::cout << "BFGS called with " << (argc - 1) << "arguments" << std::endl;

  // If the user uses no arguments.  Teach him how to use this tool
  bool callflag = 0;
  (1 == argc) ? callflag = true: callflag = false;
  if (callflag){
    std::cerr << "You have called this function without any parameters." << std::endl;
    std::cerr << "Please call this function again with the name of the " << std::endl;
    std::cerr << "function (don't use quotes) and an integer.  Example:" << std::endl;
    std::cerr << "./mytest chained_mifflin2 10" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  
  std::string str;
  for(int i = 0; argv[1][i] != 0; i++)
    str += argv[1][i];

  // If you want to add a new type, just:
  // 1) declare it.  
  // 2) Set precision.
  // 3) Run and delete
  
  std::cout << "Running phase started" << std::endl;
  unsigned int ninitial = static_cast<unsigned>(strtoul(argv[2], NULL, 0));
  algoparameters<double> * doubleparameters = 
    new algoparameters<double>(ninitial, str);
  algoparameters<dd_real> * dd_realparameters =
    new algoparameters<dd_real>(ninitial, str);
  
  std::cout << "Objects were created safely.  Running BFGS next" << std::endl;
  std::cout.precision(16);	
  doubleparameters->BFGSfunction();
  std::cout.precision(32);
  dd_realparameters->BFGSfunction();

  delete doubleparameters;
  delete dd_realparameters;
    
  fpu_fix_end(&old_cw);
  return 0;
}
