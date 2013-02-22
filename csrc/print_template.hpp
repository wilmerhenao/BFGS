#ifndef _PRINT_TEMPLATE_HPP_
#define _PRINT_TEMPLATE_HPP_

#include <ios>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <qd/dd_real.h>

template<class T> void print_iter_info(std::ofstream&, size_t& it, T*& f, 
				       T& gnorm, int& j, 
				       T *& q, T& t);
template<class T> void print_init_info(std::ofstream&,const size_t& n, const T& ftarget,
				       const T& gnormtol, const size_t& maxit, 
				       const int& echo, const int& lm, 
				       const char *& outputname);
template<class T> void print_final_info(std::ofstream&, size_t& it, T*& f,
					 T& gnorm,  int& nfeval,
					 int& exitflag,  double& ttime);

template<class T> void print_mat(T A[], int m, int n, char * str);
template<class T> void print_vec(T * v, int n, char * str);
void print_str(char * str);
template<class T> void print_double(char * str, T num);
void print_int(char * str, int num);

template<class T> void print_gs0(T r, int k, T f, T gnorm, T t, T qpinfo[]);
template<class T> void print_gs1(T r);
template<class T> void print_gs2(int k, T f);
void print_gs4(double info[]);
template<class T> void print_gs5(T f, double info[]);

template<class T>
void print_iter_info(std::ofstream& output, size_t& it,  T *& f,  T& gnorm,
		     int& j,  T *& q, T& t) {
  output << "it=" << it << ", f=" << *f << ", |g|=" << gnorm << ", t=" << t << 
    ", j=" << j << ", q=" << *q << std::endl << std::endl;
}

template<class T>
void print_init_info(std::ofstream& output, const size_t& n, const T& ftarget, 
		     const T& gnormtol, const size_t& maxit, const int& echo,
		     const int& lm, const char*& outputname) {
    
  output << "================================\n" << std::endl;
  if (!lm) {
    output << "BFGS DD on Test Function with" << std::endl;
  }
  else {
    output << "LM-BFGS DD on Test function with " << std::endl;
  }
  output << "n        = " << n << std::endl;
  output << "ftarget  = " << ftarget << std::endl;
  output << "gnormtol = " << gnormtol << std::endl;
  output << "maxit    = " << maxit << std::endl;
  output << "echo     = " << echo << std::endl;
  output << "output   = " << outputname << std::endl;
  output << "--------------------------------" << std::endl;
}

template <class T>
void print_final_info(std::ofstream& output,  size_t& it, T*& f, T& gnorm,
		      int& nfeval,  int& exitflag,  double& ttime) {
  const char * exitstr;
  switch (exitflag) {
  case 1:
    exitstr = "f < ftarget";
    break;
            
  case 2:
    exitstr = "gnorm < gnormtol";
    break;
            
  case 3:
    exitstr = "f(x0) < ftarget";
    break;
            
  case 7:
    exitstr = "optcertval < taud";
    // exitstr = "x in conv{g_1,g_2,...,g_J} with ||x|| < taud";
    break;
            
  case -1:
    exitstr = "maximal # iterations reached";
    break;
            
  case -2:
    exitstr = "gtp > 0, i.e. non-descent dir. encountered (rounding)";
    break;
            
  case -4:
    exitstr = "linesearch: nbisect >= nbisectmax. Wolfe conds not satisfied.";
    break;

  case -5:  
    exitstr = "linesearch: nexpand >= nexpandmax. Wolfe conds not satisfied.";
    break;
            
  case -8:
    exitstr = "NaN encountered. Exiting with best found f,x and g.";
    break;
            
  case 0:
    exitstr = "exitflag was never changed";
    break;
            
  default:
    exitstr = "exited for some other reason";
    printf("exitflag = %i \n",exitflag);
    break;
  }
    
  output << "--------------------------------" << std::endl;
  output << "Stopped after" << std::endl;
  output << it << "iters and " << std::endl;
  output << nfeval << "function evaluations" << std::endl; 
  std::streamsize ss = output.precision();
  output.precision(3);
  output << std::scientific << ttime << std::fixed << " seconds " << std::endl;
  output.precision(ss);
  output <<"Because " << std::endl << exitstr << std::endl;
  output <<"Found: \n" << std::endl << "f       =" << *f << std::endl;
  output <<"norm(g) = " << gnorm << std::endl;
  output <<"================================" << std::endl;
}

template<class T>
void print_mat(T A[], int m, int n, char * str) {
  // prints matrix A
  std::cout << str << std::endl;
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      std::cout << A[i * n + j] << "  ";
    std::cout << std::endl;
  }
}

template <class T>
void print_vec(T * v, int n, char * str) {
  // prints vector v
  std::cout << str << std::endl;
  for (int j = 0; j < n; j++)
    std::cout << "     " << v[j] << std::endl;
}

void print_str(char * str) {
  std::cout << str << std::endl;
}

template <class T>
void print_double(char * str, T num) {
  std::cout << str << " = " << num << std::endl;
}

void print_int(char * str, int num) {
  std::cout << str << " = " << num << std::endl;
}

template <class T>
void print_gs0(T r, int k, T f, T gnorm, T t, T qpinfo[]) {
  std::cout << "r=" << r << ", k=" << k << ", f=" << f << ", |g|=" << gnorm << 
    ", t=" << t << ", qps=" << qpinfo[0] << ", qpk=" << qpinfo[1] << ", qpP=" << 
    qpinfo[2] << std::endl;
}

template <class T>
void print_gs1(T r) {
  std::cout << "rad=" << r << ", ";
}

template <class T>
void print_gs2(int k, T f) {
  std::cout << "k=" << k << ", f=" << f << std::endl;
}

template <class T>
void print_gs4(double info[]) {
  std::cout << "Stopped because:" << std::endl;
  if (info[0] >= 7) std::cout << "f < ftarget. \n" << std::endl;
  if (info[0] == -9.0) std::cout << "f apparantly unbounded below. \n" << std::endl;
  if (info[0] == -2.0) std::cout << "max. iters reached. " << std::endl;
  else exit(EXIT_FAILURE);
}

template <class T>
void print_gs5(T f, double info[]) {
  std::cout << "GradSamp found: f=" << f << ", info(1)=" << info[0] << std::endl;
}

#endif // _PRINT_TEMPLATE_HPP_
