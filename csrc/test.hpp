#ifndef _TEST_HPP_
#define _TEST_HPP_

#include <string>
#include "../testfunctions/F10.hpp"
#include "../testfunctions/T10.hpp"
#include "../testfunctions/yurirosen.hpp"

//This class contains an instantiation of all functions of the given type
template<typename T>
class allfunctions: public allfunctionsF10: public allfunctionsT10: 
    public allfunctionsyurirosen {
public:
    allfunctions();
}

template<typename T>
allfunctions<T>::allfunctions(){
    allfunctionsF10.fillMap();
    allfunctionsT10.fillMap();
    allfunctionsyurirosen.fillMap();
}

template<typename T>
class algoparameters{
private:
    friend void bfgs<T>(T [], T * , const int, const int, const int, const T, const T, 
			const int, const int, const T, const T, const int, 
			void(*)(T*, T*, T*, int), const char * , double[]);
    T ftarget, gnormtol, taux, taud;
    int J, maxit, echo, n, lm, m;
    const char * datafilename;
    T * fopt;
    double * info;
    T * xf;
    double * x0;    
    void(*fun_ptr)(T*, T*, T*, int);
    allfunctions<T> * pFunctions;
public:
    algoparameters(int, string);
    ~algoparameters();
    //function that gets a string and returns a pointer to that function
    void find_function(string);
    void generateXF();
    void BFGSfunction();
};

template<typename T>
algoparameters<T>::algoparameters(int k, string locfunc){
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
    pFunctions = new allfunctions<T>;
    //void *prg_handle = dlopen("../lib/libtestfun.so", RTLD_NOW);
    fun_ptr = find_function(locfunc);
    generateXF();
}

template<typename T>
algoparameters<T>::~algoparameters(){
    delete fopt;
    delete [] info;
    delete [] xf;
    delete [] x0;
    delete fun_ptr;
}

template<typename T>
void find_function(char *){
    fun_ptr = 
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
	    echo, fun_ptr, datafilename, info);
}

#endif // _TEST_HPP_
