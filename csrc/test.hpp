#ifndef _TEST_HPP_
#define _TEST_HPP_

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
public:
    algoparameters(int, char *);
    ~algoparameters();
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
    fun_ptr = (void(*)(T*, T*, T*, int)) dlsym(prg_handle, locfunc);
    generateXF();
}

template<typename T>
algoparameters<T>::~algoparameters(){
    delete fopt;
    delete [] info;
    delete [] xf;
    delete [] x0;
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
