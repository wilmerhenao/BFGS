#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <dlfcn.h>
#include "randnums_template.hpp"
#include "bfgs_template.hpp"
#include "F10_dd.h"
#include "T10_dd.h"
#include "yurirosen_dd.h" 
#include "libmatrix_template.hpp"
#include "F10.h"
#include "T10.h"
#include "yurirosen.h"
#include <qd/dd_real.h>

template<typename T>
struct algoparameters{
    T ftarget, gnormtol, taux, taud;
    int J, maxit, echo, n, lm, m;
    const char * datafilename;
    T * fopt;
    double * info;
    T * xf;
    double * x0;    
    algoparameters(int);
    void generateXF();
};

template<typename T>
algoparameters<T>::algoparameters(int k){
    n = k;
    lm=0;
    m = 7;
    taud = 1e-14;
    taux = 1e-16;
    ftarget = -1e100;
    gnormtol = 0.0;
    maxit = 10e8;
    echo = 2;
    datafilename = "stddump.txt";
    J = fmin(15, ceil(2 * n / 3));
    fopt = new T;
    info = new double[4];
    xf = new T[n];
    x0 = new double[n];
}

template <typename T> 
void algoparameters<T>::generateXF(){
    struct timeval time;
    long useconds;
    gettimeofday(&time, NULL);
    useconds = time.tv_usec; //pseudo-randon number generator seed
    srand(useconds);
    rand_real_vec<double>(x0, n, -1, 1);
    for (int j = 0 ; j < n; j++)
	xf[j] = (T)x0[j];
}

int main(int argc, char *argv[]){    
    unsigned int old_cw;
    fpu_fix_start(&old_cw);//necessary for double-double (see Sherry's paper 
    cout.precision(32);	
    
    typedef void(*test_fun_ptr_type)(double*, double*, double*, int);
    test_fun_ptr_type testFunction;
    typedef void(*test_fun_ptr_type_dd)(dd_real*, dd_real*, dd_real*, int);
    test_fun_ptr_type_dd testFunction_dd;
	
    algoparameters<double> doubleparameters(atoi(argv[3]));
    algoparameters<dd_real> dd_realparameters(atoi(argv[3]));
    /*
    int n;
    int lm;
    int m;
    
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
    */
    // Perhaps this RTLD_NOW here is wasting some time, 
    //this is something that could be
    // investigated further
    // Catch error if dlopen returns NULL (see wikipedia)
    // It would be nice to close this library
    void *prg_handle = dlopen("../lib/libtestfun.so", RTLD_NOW);    
    
    /* Declare needed variables (output) */
    /*
    double * fopt = new double;
    dd_real * fopt_dd = new dd_real;
    double * info = new double[4];
    double * info_dd = new double[4];
    */
    /* Assign pointers: */
    /* input: */
    testFunction = (test_fun_ptr_type) dlsym(prg_handle, argv[1]);
    testFunction_dd = (test_fun_ptr_type_dd) dlsym(prg_handle, argv[2]);
    /*
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
    */
    /*
    double *x0 = new double[n];        
    double *xf = new double[n];
    dd_real *xf_dd = new dd_real[n];

    struct timeval time;
    long useconds;
    gettimeofday(&time, NULL);
    useconds = time.tv_usec;	//pseudo-random number generator seed

    srand(useconds);
    rand_real_vec<double>(x0, n, -1, 1);*/
    /*
      x0[0] = -1;
      for (int i=1; i<n; i++) {
      x0[i] = 1;
      }
    */

    /* output: */    
    /* Create matrices for the return arguments. */

    /* Copy xf <-- x0 so x0 is still available after BFGS run */ 	
    /*for (int j=0; j<n; j++) {
	xf[j] = xf_dd[j] = x0[j];
    }*/
    // Call the BFGS D and DD subroutine.

    doubleparameters.generateXF();
    dd_realparameters.generateXF();

    bfgs<double>(doubleparameters.xf, doubleparameters.fopt, 
		 doubleparameters.n, doubleparameters.lm, doubleparameters.m,
		 doubleparameters.ftarget, doubleparameters.gnormtol, 
		 doubleparameters.maxit, doubleparameters.J, 
		 doubleparameters.taux, doubleparameters.taud, 
		 doubleparameters.echo, testFunction, 
		 doubleparameters.datafilename, doubleparameters.info);

    bfgs<dd_real>(dd_realparameters.xf, dd_realparameters.fopt, 
		  dd_realparameters.n, dd_realparameters.lm, 
		  dd_realparameters.m, dd_realparameters.ftarget, 
		  dd_realparameters.gnormtol, dd_realparameters.maxit, 
		  dd_realparameters.J, dd_realparameters.taux, 
		  dd_realparameters.taud, dd_realparameters.echo, 
		  testFunction_dd, dd_realparameters.datafilename, 
		  dd_realparameters.info);
    
    fpu_fix_end(&old_cw);
}
