//#define NDEBUG //Uncomment if you are not debugging (FASTER)

#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <string>
#include <cassert>
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
    std::cerr << "function (don't use quotes) an int and output. Example:" << std::endl;
    std::cerr << "./mytest chained_mifflin2 10 chained_mifflin2_10.txt" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  
  std::string str;
  for(size_t i = 0; 0 != argv[1][i]; i++)
    str += argv[1][i];
  assert(str.length());
  
  std::string strout;
  for(size_t i = 0; 0 != argv[3][i]; i++)
    strout += argv[3][i];
  assert(str.length());

  // If you want to add a new type, just:
  // 1) declare it.  
  // 2) Set precision.
  // 3) Run and delete
  
  std::cout << "output file" << strout << std::endl;
  unsigned int ninitial = static_cast<unsigned>(strtoul(argv[2], NULL, 0)); 
  algoparameters<double> * doubleparameters;
  algoparameters<dd_real> * dd_realparameters;
  try{
     doubleparameters =  new algoparameters<double>(ninitial, str, strout);
     dd_realparameters =  new algoparameters<dd_real>(ninitial, str, strout);
  } catch(std::bad_alloc& ex){
    std::cerr << "Problem allocating memory in the stack" << ex.what() << std::endl;
    assert(false);
  }
  assert(sizeof(doubleparameters));
  assert(sizeof(dd_realparameters));

  std::cout << "Objects were created safely.  Running BFGS next" << std::endl;
  try{
    std::cout.precision(16);	
    doubleparameters->BFGSfunction();
    std::cout.precision(32);
    dd_realparameters->BFGSfunction();
  } catch(std::exception ex){
    std::cerr << "General Problem during execution. " << ex.what() << std::endl;
    assert(false);
  }

  try{
    delete doubleparameters;
    delete dd_realparameters;
  } catch(std::exception ex){
    std::cerr << "Issues with delete: " << ex.what() << std::endl;
    assert(false);
  }

  fpu_fix_end(&old_cw);
  return 0;
}
