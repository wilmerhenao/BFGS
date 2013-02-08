#include <qd/dd_real.h>

void print_iter_info(FILE *output, int it, dd_real * f, dd_real gnorm, int j, dd_real * q, dd_real * x, dd_real t,  int n);
void print_init_info(FILE *output, int n, dd_real ftarget, dd_real gnormtol, int maxit, int echo, int lm, const char*outputname, void(*testFunction)(dd_real*, dd_real*, dd_real*, int));
void print_final_info(FILE *output, int it, dd_real f, dd_real gnorm, int nfeval, int exitflag, double ttime);

void print_mat(dd_real A[], int m, int n, char * str);
void print_vec(dd_real * v, int n, char * str);
void print_int_vec(int * v, int n, char * str);
void print_str(char * str);
void print_double(char * str, dd_real num);
void print_int(char * str, int num);

void print_gs0(dd_real r, int k, dd_real f, dd_real gnorm, dd_real t, dd_real qpinfo[]);
void print_gs1(dd_real r);
void print_gs2(int k, dd_real f);
void print_gs4(double info[]);
void print_gs5(dd_real f, double info[]);


