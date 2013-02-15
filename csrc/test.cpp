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
class algoparameters{
public:
    T ftarget, gnormtol, taux, taud;
    int J, maxit, echo, n, lm, m;
    const char * datafilename;
    T * fopt;
    double * info;
    T * xf;
    double * x0;    
    void(*test_fun_ptr)(T*, T*, T*, int);
    algoparameters(int, char *);
    void generateXF();
  void BFGSfunction();
};

template<typename T>
algoparameters<T>::algoparameters(int k, char * locfunc){
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
    void *prg_handle = dlopen("../lib/libtestfun.so", RTLD_NOW); 
    test_fun_ptr = (void(*FPTR)(T*, T*, T*, int)) dlsym(prg_handle, locfunc);
    generateXF();
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

template<typename T>
void algoparameters<T>::BFGSfunction(){
    bfgs<T>(xf, fopt, n, lm, m, ftarget, gnormtol, maxit, J, taux, taud, 
		 echo, test_fun_ptr, datafilename, info);
}

int main(int argc, char *argv[]){    
    unsigned int old_cw;
    fpu_fix_start(&old_cw);//necessary for double-double (see Sherry's paper 
    cout.precision(32);	
    
    /*typedef void(*test_fun_ptr_type)(double*, double*, double*, int);
    test_fun_ptr_type testFunction;
    typedef void(*test_fun_ptr_type_dd)(dd_real*, dd_real*, dd_real*, int);
    test_fun_ptr_type_dd testFunction_dd;*/
	
    algoparameters<double> doubleparameters(atoi(argv[3]), argv[1]);
    algoparameters<dd_real> dd_realparameters(atoi(argv[3]), argv[1]);

    //void *prg_handle = dlopen("../lib/libtestfun.so", RTLD_NOW);    
    
    //testFunction = (test_fun_ptr_type) dlsym(prg_handle, argv[1]);
    //testFunction_dd = (test_fun_ptr_type_dd) dlsym(prg_handle, argv[2]);

    doubleparameters.BFGSfunction();
    dd_realparameters.BFGSfunction();
    
    fpu_fix_end(&old_cw);
}
