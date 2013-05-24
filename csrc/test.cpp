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

//#define NDEBUG //Comment if you are debugging

#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <string>
#include <cassert>
#include "../lib/qpspecial/qpobject.hpp"
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
  bool callflag = false;
  (1 == argc) ? callflag = true: callflag = false;
  if (callflag){
    std::cerr << "You have called this function without any parameters." << std::endl;
    printhowtodoit();  
  }
  (1 == argc) ? callflag = true: callflag = false;
  if (callflag) {
    std::cerr<<"You have called this function without enough parameters "<<std::endl;
    printhowtodoit();
  }

  (argc > 4) ? callflag = true: callflag = false;
  if (callflag) {
    std::cerr<<"You have called this function with too many parameters "<<std::endl;
    printhowtodoit();
  }
  
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  
  std::string str;
  for(int i = 0; 0 != argv[1][i]; i++)
    str += argv[1][i];
  assert(str.length());
  
  std::string strout;
  for(int i = 0; 0 != argv[3][i]; i++)
    strout += argv[3][i];
  assert(str.length());

  // If you want to add a new type, just:
  // 1) declare it.  
  // 2) Set precision.
  // 3) Run and delete
  
  std::cout << "output file" << strout << std::endl;
  int ninitial = atoi(argv[2]); 
  algoparameters<double> * doubleparameters;
  algoparameters<dd_real> * dd_realparameters;
  algoparameters<dd_real> * dd_realparLBFGS;
  try{
    doubleparameters = new algoparameters<double>(ninitial, str, strout, 0);
    dd_realparameters = new algoparameters<dd_real>(ninitial, str, strout, 0);
    dd_realparLBFGS = new algoparameters<dd_real>(ninitial, str, strout, 1);
  } catch(std::bad_alloc& ex){
    std::cerr << "Problem allocating memory in the stack" << ex.what() << std::endl;
    assert(false);
  }
  assert(sizeof(doubleparameters));
  assert(sizeof(dd_realparameters));
  assert(sizeof(dd_realparLBFGS));
  
  std::cout << "Objects were created safely.  Running BFGS next" << std::endl;
  try{
    std::cout.precision(16);	
    //doubleparameters->BFGSfunction();
    std::cout.precision(32);
    //dd_realparameters->BFGSfunction();
    //dd_realparLBFGS->BFGSfunction();  // This one redirects to the right LBFGS
  } catch(std::exception ex){
    std::cerr << "General Problem during execution. " << ex.what() << std::endl;
    assert(false);
  }
  
  try{
    delete doubleparameters;
    delete dd_realparameters;
    delete dd_realparLBFGS;
  } catch(std::exception ex){
    std::cerr << "Issues with delete: " << ex.what() << std::endl;
    assert(false);
  }
  
  // This is for simple bounded problems.  Definition of the border
  
  double * u = new double[static_cast<unsigned>(strtoul(argv[2], NULL, 0)) ];
  double * l = new double[static_cast<unsigned>(strtoul(argv[2], NULL, 0)) ];
  // dummy initialization just to have something to work with...
  for(int counter = 0; counter < static_cast<int>(strtoul(argv[2], NULL, 0)) ;
      counter++){
    u[counter] = 3.0; l[counter] = 1.0;
  }
  l[0] = -3.0;
  l[3] = -3.0;
  
  // Declaration of the problem
  
  algoparameters<double> * constrainedproblem;
  try{
    constrainedproblem = new algoparameters<double>(ninitial, str, strout, 1, u, l);
  } catch(std::exception ex) {
    std::cerr << "Problem assigning memory to constrained" << ex.what() << std::endl;
    assert(false);
  }
  assert(sizeof(constrainedproblem));
  
  std::cout  << "Constrained problem was created safely" << std::endl;
  try{
    std::cout.precision(16);
    constrainedproblem->BFGSfunction();
  } catch(std::exception ex){
    std::cerr << "General Problem during execution. " << ex.what() << std::endl;
    assert(false);
  }
  
  
  try{
    delete [] u;
    delete [] l;
  } catch(std::exception ex){
    std::cerr << "Issues deleting mem. constrd. Problem" << ex.what() << std::endl;
    assert(false);
  }
  fpu_fix_end(&old_cw);
  return 0;
}

void printhowtodoit(){
  std::cerr << "Please call this function again with the name of the " << std::endl;
  std::cerr <<"function (don't use quotes) an int and output. Example:"<< std::endl;
  std::cerr << "./mytest chained_mifflin2 10 chained_mifflin2_10.txt" << std::endl;
  exit(EXIT_FAILURE);
}
