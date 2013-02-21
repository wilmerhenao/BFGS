#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <dlfcn.h>
#include <string>
#include "randnums_template.hpp"
#include "bfgs_template.hpp"
#include "F10.hpp"
#include "libmatrix_template.hpp"
#include "T10.hpp"
#include "yurirosen.hpp"
#include "test.hpp"
#include <qd/dd_real.h>

int main(int argc, char *argv[]){    
    unsigned int old_cw;
    fpu_fix_start(&old_cw);//necessary for double-double (see Sherry's paper 
    
    std::string str;
    for(int i = 0; argv[1][i] != 0; i++)
	str += array[i];

    algoparameters<double> * doubleparameters = 
	new algoparameters<double>(atoi(argv[3]), str);
    algoparameters<dd_real> * dd_realparameters =
	new algoparameters<dd_real>(atoi(argv[3]), str);

    doubleparameters->BFGSfunction();
    dd_realparameters->BFGSfunction();
    
    cout.precision(16);	
    delete doubleparameters;
    cout.precision(32);
    delete dd_realparameters;
    
    fpu_fix_end(&old_cw);
    return 0;
}
