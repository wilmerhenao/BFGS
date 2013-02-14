using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <dlfcn.h>
#include "randnums_template.hpp"
#include "libmatrix_dd.h"
#include "bfgs_dd.h"
#include "F10_dd.h"
#include "T10_dd.h"
#include "yurirosen_dd.h" 
#include "libmatrix.h"
#include "bfgs.h"
#include "F10.h"
#include "T10.h"
#include "yurirosen.h"
#include <qd/dd_real.h>

int main(int argc, char *argv[]) {
    unsigned int old_cw;
    fpu_fix_start(&old_cw);		//necessary for double-double
    cout.precision(32);			//set precision
    
    /* Declare needed variables (input) */
    // Don't forget that this is just a pointer to a function
    // we might want to replace this with templates
    typedef void(*test_fun_ptr_type)(double*, double*, double*, int);
    test_fun_ptr_type testFunction;
    typedef void(*test_fun_ptr_type_dd)(dd_real*, dd_real*, dd_real*, int);
    test_fun_ptr_type_dd testFunction_dd;
    int n;
    int lm;
    int m;
    
    // All these are candidates for template definitions
    double ftarget;
    double gnormtol;
    double taux, taud;
    dd_real ftarget_dd;
    dd_real gnormtol_dd;
    dd_real taux_dd, taud_dd;
    
    int J;
    int maxit;
    int echo;
    const char * datafilename;
    const char * datafilename_dd;
    
    // Perhaps this RTLD_NOW here is wasting some time, this is something that could be
    // investigated further
    // Catch error if dlopen returns NULL (see wikipedia)
    // It would be nice to close this library
    void *prg_handle = dlopen("../lib/libtestfun.so", RTLD_NOW);    

    /* Declare needed variables (output) */

    // more candidates for templating
    double * fopt = new double;
    dd_real * fopt_dd = new dd_real;
    double * info = new double[4];
    double * info_dd = new double[4];

    /* Assign pointers: */
    /* input: */
    testFunction = (test_fun_ptr_type) dlsym(prg_handle, argv[1]);
    testFunction_dd = (test_fun_ptr_type_dd) dlsym(prg_handle, argv[2]);

    n            = atoi(argv[3]);	//# variables
    lm           = 0;			//lbfgs or not
    m            = 7;			//lbfgs paramete
    J            = fmin(15,ceil(2*n/3));
    taux         = 1e-16;		//stpsize tolerance
    taud         = 1e-4;
    ftarget      = -1e100;		//fvalquit
    gnormtol     = 0;
    maxit        = 10e8;		//maximum iterations
    echo         = 2;			//print level
    datafilename = "stddump.txt";
    datafilename_dd = "stddump_dd.txt";

    taux_dd         = 1e-16;		//stepsize tolerance
    taud_dd         = 1e-4;
    ftarget_dd      = -1e100;		//fvalquit
    gnormtol_dd     = 0.0;

    double *x0 = new double[n];        
    double *xf = new double[n];
    dd_real *xf_dd = new dd_real[n];

    struct timeval time;
    long useconds;
    gettimeofday(&time, NULL);
    useconds = time.tv_usec;	//pseudo-random number generator seed

    srand(useconds);
    rand_real_vec<double
>(x0, n, -1, 1);

/*
    x0[0] = -1;
    for (int i=1; i<n; i++) {
	x0[i] = 1;
    }
*/

    /* output: */    
    /* Create matrices for the return arguments. */

    /* Copy xf <-- x0 so x0 is still available after BFGS run */ 	
    for (int j=0; j<n; j++) {
	xf[j] = x0[j];
    }

    for (int j=0; j<n; j++) {
	xf_dd[j] = x0[j];
    }
    
// Call the BFGS D and DD subroutine.

    bfgs(xf, fopt, n, lm, m, ftarget, gnormtol, maxit, J, taux, taud, echo, testFunction, datafilename, info);

    bfgs_dd(xf_dd, fopt_dd, n, lm, m, ftarget_dd, gnormtol_dd, maxit, J, taux_dd, taud_dd, echo, testFunction_dd, datafilename_dd, info_dd);

    fpu_fix_end(&old_cw);
}
