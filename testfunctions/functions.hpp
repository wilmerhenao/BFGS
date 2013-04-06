#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <qd/dd_real.h>
#include <string>
#include <iostream>
#include <map>
#include "../csrc/container.hpp"

template<class T> int chained_CB3v1(T *, T *, T *, int);
template<class T> int chained_CB3v2(T *, T *, T *, int);
template<class T> int chained_crescent1(T *, T *, T *, int);
template<class T> int chained_crescent2(T *, T *, T *, int);
template<class T> int chained_LQ(T *, T *, T *, int);
template<class T> int chained_mifflin2(T *, T *, T *, int);
template<class T> int gen_brownfunc2(T *, T *, T *, int);
template<class T> int gen_maxhilbert(T *, T *, T *, int);
template<class T> int gen_maxq(T *, T *, T *, int);
template<class T> int nactfaces(T *, T *, T *, int);
template<class T> int test29f02(T *, T *, T *, int);
template<class T> int test29f05(T *, T *, T *, int);
template<class T> int test29f06(T *, T *, T *, int);
template<class T> int test29f11(T *, T *, T *, int);
template<class T> int test29f22(T *, T *, T *, int);
template<class T> int yurirosen(T *, T *, T *, int);
template<class T> int yurirosen_ns1(T *, T *, T *, int);
template<class T> int yurirosen_ns2(T *, T *, T *, int);

template<typename T>
class allfunctions{
public:
  int (*pchained_CB3v1)(T *f, T *g, T *x, int n)  = &chained_CB3v1<T>;
  int (*pchained_CB3v2)(T *f, T *g, T *x, int n) = &chained_CB3v2<T>;
  int (*pchained_crescent1)(T *f, T *g, T *x, int n) = &chained_crescent1<T>;
  int (*pchained_crescent2)(T *f, T *g, T *x, int n) = &chained_crescent2<T>;
  int (*pchained_LQ)(T *f, T *g, T *x, int n) = &chained_LQ<T>;
  int (*pchained_mifflin2)(T *f, T *g, T *x, int n) = &chained_mifflin2<T>;
  int (*pgen_brownfunc2)(T *f, T *g, T *x, int n) = &gen_brownfunc2<T>;
  int (*pgen_maxhilbert)(T *f, T *g, T *x, int n) = &gen_maxhilbert<T>;
  int (*pgen_maxq)(T *f, T *g, T *x, int n) = &gen_maxq<T>;
  int (*pnactfaces)(T *f, T *g, T  *x, int n) = &nactfaces<T>;
  int (*ptest29f02)(T *f, T *g, T *x, int n) = &test29f02<T>;
  int (*ptest29f05)(T *f, T *g, T *x, int n) = &test29f05<T>;
  int (*ptest29f06)(T *f, T *g, T *x, int n) = &test29f06<T>;
  int (*ptest29f11)(T *f, T *g, T *x, int n) = &test29f11<T>;
  int (*ptest29f22)(T *f, T *g, T *x, int n) = &test29f22<T>;
  int (*pyurirosen)(T *f, T *g, T *x, int n) = &yurirosen<T>;
  int (*pyurirosen_ns1)(T *f, T *g, T *x, int n) = &yurirosen_ns1<T>;
  int (*pyurirosen_ns2)(T *f, T *g, T *x, int n) = &yurirosen_ns2<T>;
  std::map<std::string, int(*)(T*, T*, T*, int), StringComparerForMap> tMap;
  int fillMap();
};

template<typename T>
int allfunctions<T>::fillMap(){
  tMap["chained_CB3v1"] = pchained_CB3v1;
  tMap["chained_CB3v2"] = pchained_CB3v2;
  tMap["chained_crescent1"] = pchained_crescent1;
  tMap["chained_crescent2"] = pchained_crescent2;
  tMap["chained_LQ"] = pchained_LQ;
  tMap["chained_mifflin2"] = pchained_mifflin2;
  tMap["gen_brownfunc2"] = pgen_brownfunc2;
  tMap["gen_maxhilbert"] = pgen_maxhilbert;
  tMap["gen_maxq"] = pgen_maxq;
  tMap["nactfaces"] = pnactfaces;
  tMap["test29f02"] = ptest29f02;
  tMap["test29f05"] = ptest29f05;
  tMap["test29f06"] = ptest29f06;
  tMap["test29f11"] = ptest29f11;
  tMap["test29f22"] = ptest29f22;
  tMap["yurirosen"] = pyurirosen;
  tMap["yurirosen_ns1"] = pyurirosen_ns1;
  tMap["yurirosen_ns2"] = pyurirosen_ns2;
  return 0;
}

// Implementations:
template<class T> 
int chained_CB3v1(T *f, T *g, T *x, int n){
  /* printthis:
     % Chained CB3 I (convex, non-smooth func):
     % f = 
     % 0 parameters: w, w
     % show: n, n
  */
  T y;
  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;

  for (int i = 0; i < (n-1); i++) {
    T a = pow(x[i], 4.0) + pow(x[i + 1], 2.0);
    T b = pow((2 - x[i]), 2.0) + pow((2 - x[i + 1]), 2.0);
    T c = 2 * exp(-x[i] + x[i + 1]);
    (a >= b) ? y = a : y = b;
    if (y < c)
      y = c;
    if (y == a) {
      g[i]   = g[i] + 4 * pow(x[i], 3.0);
      g[i + 1] = 2 * x[i + 1];
    }
    else if (y == b) {
      g[i] = g[i] + 2 * x[i] - 4;
      g[i + 1] = 2 * x[i + 1] - 4;
    }
    else {
      g[i] = g[i] - c;
      g[i + 1] = c;
    }
    *f = *f + y;
  }    
  return 0;
}

template<class T> 
int chained_CB3v2(T *f, T *g, T *x, int n){
  /* printthis:
     % Chained CB3 II (convex, non-smooth func):
     % f = 
     % 0 parameters: w, w
     % show: n, n0
  */
  T a = 0.0; 
  T b = 0.0; 
  T c = 0.0;

  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;
  for (int i = 0; i < (n-1); i++) {
    a = a + pow(x[i], 4.0) + pow(x[i + 1], 2.0);
    b = b + pow((2 - x[i]), 2.0) + pow((2 - x[i + 1]), 2.0);
    c = c + 2 * exp(-x[i] + x[i + 1]);
  }    
  (a >= b) ? *f = a : *f = b;
  if (*f < c)
    *f = c;

  if (*f == a) {
    for (int i = 0; i < (n-1); i++) {
      g[i]   = g[i] + 4 * pow(x[i], 3.0);
      g[i + 1] = 2 * x[i + 1];
    }
  }
  else if (*f == b) {
    for (int i = 0; i < (n-1); i++) {
      g[i] = g[i] + 2 * x[i] - 4;
      g[i + 1] = 2 * x[i + 1] - 4;
    }
  }
  else {
    for (int i = 0; i < (n-1); i++) {
      g[i] = g[i] - 2 * exp(-x[i] + x[i + 1]);
      g[i + 1] = 2 * exp(-x[i] + x[i + 1]);
    }
  }
  return 0;
}

template<class T> 
int chained_crescent1(T *f, T *g, T *x, int n)
{

  /* printthis:
     % Chained Crescent 1:
     % f = 
     % 0 parameters: w, w
     % show: n, n
  */

  T t2 = 0.0;
  T t3 = 0.0;

  *f = 0.0;
  for (int i = 0; i < n ; i++)
    g[i] = 0.0;

  for (int i = 0; i < (n - 1); i++) { 
    t2 = t2 + pow(x[i], 2.0) + pow((x[i + 1] - 1), 2.0) + x[i + 1] - 1;
    t3 = t3 - pow(x[i], 2.0) - pow((x[i + 1] - 1), 2.0) + x[i + 1] + 1;
  }
  (t2>=t3) ? *f = t2 : *f = t3;
  if (t2 >= t3) {
    for (int i = 0; i < (n-1); i++) {
      g[i] = g[i] + 2 * x[i];
      g[i + 1] = 2 * (x[i + 1] - 1) + 1;
    }
  }
  else {
    for (int i = 0; i < (n-1); i++) {
      g[i] = g[i] - 2 * x[i];
      g[i + 1] = -2 * (x[i + 1] - 1) + 1;
    }
  }
  return 0;
}

template<class T> 
int chained_crescent2(T *f, T *g, T *x, int n){
  /* printthis:
     % Chained Crescent 2:
     % f = 
     % 0 parameters: w, w
     % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i<n; i++) 
    g[i] = 0.0;

  for (int i = 0; i < (n-1); i++) { 
    T t2 = pow(x[i], 2.0) + pow((x[i + 1] - 1), 2.0) + x[i + 1] - 1;
    T t3 = -pow(x[i], 2.0) - pow((x[i + 1] - 1), 2.0) + x[i + 1] + 1;
    if (t2 >= t3) {
      *f = *f + t2;
      g[i]   = g[i] + 2 * x[i];
      g[i + 1] = 2 * (x[i + 1] - 1) + 1;
    }
    else {
      *f = *f + t3;
      g[i]   = g[i] - 2 * x[i];
      g[i + 1] = -2 * (x[i + 1] - 1) + 1;
    }
  }
  return 0;
}

template<class T> 
int chained_LQ(T *f, T *g, T *x, int n)
{
  /* printthis:
     % Chained LQ (convex func):
     % f = sum_i( max( -x_i-x_{i+1},-x_i-x_{i+1}+x_i^2+x_{i+1}^2-1 ) )
     % 0 parameters: w, w
     % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;

  for (int i = 0; i < (n-1); i++) {
    T a = -x[i] - x[i + 1];
    T b = a + ( pow(x[i], 2.0) + pow(x[i + 1], 2.0) - 1 );
    if (a >= b) {
      *f = *f + a;
      g[i]   = g[i] - 1;
      g[i + 1] = -1;
    }
    else {
      *f = *f + b;
      g[i]   = g[i] - 1 + 2 * x[i];
      g[i + 1] = -1 + 2 * x[i + 1];
    }
  }
  return 0;
}

template<class T> 
int chained_mifflin2(T *f, T *g, T *x, int n)
{

  /* printthis:
     % Chained Mifflin II (non-convex, non-smooth func):
     % f = 
     % 0 parameters: w, w
     % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i < n; i++) 
    g[i] = 0.0;

  for (int i = 0; i < (n-1); i++) {
    T y = x[i] * x[i] + x[i + 1] * x[i + 1] - 1;
    *f = *f - x[i] + 2 * y + 1.75 * fabs(y);
    T s;
    if (y >= 0) {
      s = 3.5;
    }
    else {
      s = -3.5;
    }
    y = s + 4;
    g[i] = g[i] + y * x[i] - 1;
    g[i + 1] = y * x[i + 1];
    
  }
  return 0;
}

template<class T> 
int gen_brownfunc2(T *f, T *g, T *x, int n)
{

  /* printthis:
     % Non-smooth generalization of brown function 2 (non-convex, non-smooth):
     % f = sum_i ( abs(x_i)^{x_{i+1}^2+1} + abs(x_{i+1})^{x_i^2+1} )
     % 0 parameters: w, w
     % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;

  for (int i = 0; i < (n-1); i++) {
    T a = fabs(x[i]);
    T b = fabs(x[i+1]);
    T c = pow(x[i], 2.0) + 1;
    T d = pow(x[i+1], 2.0) + 1;
    *f = *f + pow(b, c) + pow(a, d);
    
    T p = 0.0;
    T q = 0.0;
    if (x[i] < 0) {
      if (b > p)
	p = log(b);
      g[i] = g[i] - d * pow(a, (d - 1)) + 2 * x[i] * p * pow(b, c);
    }
    else {
      if (b > p)
	p = log(b);
      g[i] = g[i] + d * pow(a, (d - 1)) + 2 * x[i] * p * pow(b, c);
    }
    if (x[i + 1] == 0)
      g[i + 1] = 0.0;
    else if (x[i + 1] < 0) {
      if (a > q) {
	q = log(a);
      }
      g[i + 1] = -c * pow(b, (c - 1)) + 2 * x[i + 1] * q * pow(a, d);
    }
    else {
      if (a > q) {
	q = log(a);
      }
      g[i + 1] = c * pow(b, (c - 1)) + 2 * x[i + 1]* q * pow(a, d);
    }
  }
  return 0;
}

template<class T> 
int gen_maxhilbert(T *f, T *g, T *x, int n){
  /*
    % printthis:
    % Max of Hilberts (convex func):
    % f = max_i sum_j abs( x_j / (i+j-1) )
    % 0 parameters: w, w
    % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;

  int k = 0;
  for (int j = 0; j < n; j++)	    
    *f = *f + x[j] / static_cast<double> (j + 1);

  if (*f > 0) {
    g[0] = 1.0;
  }
  else if (*f < 0) {
    g[0] = -1.0;
  }
  else {
    g[0] = 0.0;
  }

  *f = fabs(*f);

  for (int i = 1; i < n; i++) {
    T t2 = 0;
    for (int j = 0; j < n; j++)
      t2 = t2 + x[j] / static_cast<double> (i + j + 1);
	    
    if (t2 > 0) {
      g[i] = 1.0;
    }
    else if (t2 < 0) {
      g[i] = -1.0;
    }
    else {
      g[i] = 0.0;
    }
    
    t2 = fabs(t2);

    if (t2 > *f) {
      *f = t2;
      k = i;
    }
  }

  T t3 = g[k];

  for (int j = 0; j < n; j++) {
    g[j] = t3 / static_cast<double> (k + j + 1);
  }
  return 0;
}

template<class T> 
int gen_maxq(T *f, T *g, T *x, int n)
{

  /* printthis:
     % Max of quadratics (convex func):
     % f = max (x_i)^2, i = 1,..,n
     % 0 parameters: w, w
     % show: n, n
  */

  /* Need to multiply to get 2d matrix, and 
     take max first value along first row
  */

  for (int i = 0; i<n; i++) {
    g[i] = 0.0;
  }

  *f   = x[0] * x[0];
  int k = 0;

  for (int j = 1; j < n; j++) {
    T y = x[j] * x[j];
    if (y > *f) {
      *f = y;
      k = j;
    }
  }
  g[k] = 2 * x[k];
  return 0;
}

template<class T> 
int nactfaces(T *f, T *g, T *x, int n)
{

  /* printthis:
     % Number of active faces (non-convex, non-smooth func):
     % f = max_i {g(x_i),g(-sum_j x_j)}, g(y) = ln(abs(y)+1)
     % 0 parameters: w, w
     % show: n, n
  */

  *f = 0.0;
  for (int i = 0; i < n; i++)
    g[i] = 0.0;
  int k = 0;
	
  T t3 = 1.0;
  T y  = -x[0];
  *f  = log( fabs(x[0]) + 1);
  T  t2 = *f;
	
  for (int i = 1; i < n; i++) {
    
    y    = y - x[i];
    g[i] = 0.0;
    if (*f <= log(fabs(x[i]) + 1))
      *f = log(fabs(x[i]) + 1);
    if (*f > t2) {
      k  = i;
      t2 = *f;
    }
  }
	
  if (*f <= log(fabs(y) + 1))
    *f = log(fabs(y) + 1);

  if (*f > t2) {
    if (y > 0)
      t3 = -1;
    for (int i = 0; i < n; i++)
      g[i] = t3 * ( 1 / (fabs(y) + 1) );
  }
  else {
    if (x[k] < 0) {
      t3 = -1;
    }
    g[k] = t3 * ( 1/ (fabs(x[k]) + 1));
  }
  return 0;
}

template<class T> 
int test29f02(T *f, T *g, T *x, int n) {
    
  T fval = 0.0;
  T tmp  = 0.0;
  T tmp2 = 1.0;
  int idx = 0;
    
  for (int j = 0; j < n; j++) {
    g[j] = 0.0;
    tmp  = fabs(x[j]);        
    if (tmp > fval) {
      fval = tmp;
      idx  = j;            
    }
  }
    
  if (x[idx] < 0) tmp2 = -1.0;
  g[idx] = tmp2;
  *f     = fval;
  return 0;
}

template<class T> 
int test29f05(T *f, T *g, T *x, int n) {
  // This function can't really handle the whole range of ints
    
  T fval = 0.0;
  T fa   = 0.0;
  T tmp3 = 1.0;
    
  for (int j = 0; j < n; j++) 
    g[j] = 0.0;	
    
  for (int j = 0; j < n; j++) {
    fa = 0.0;
    for (int i = 0; i < n; i++) {
      fa = fa + x[i] / ( ((T) static_cast<int> (i + j + 1)) );
    }
    fval = fval + fabs(fa);
    tmp3 = 1.0;
    if (fa < 0) tmp3 = -1.0;
    for (int i = 0; i < n; i++) {
      g[i] = g[i] + tmp3 / ( ((T)static_cast<int> ( i + j + 1)) );
    }
  }
  *f = fval;
  return 0;  
}

template<class T> 
int test29f06(T *f, T *g, T *x, int n) {
    
  T fval = 0.0;
  T tmp  = 0.0;
  T tmp2 = 0.0;
  T tmp3 = 1.0;
  int k = 0;
    
  for (int j = 0; j < n; j++) {
    g[j] = 0.0;
    tmp  = (3.0 - 2.0 * x[j]) * x[j] + 1;
    if (j > 0) tmp = tmp - x[j - 1];
    if (j < (n - 1)) tmp = tmp - x[j + 1];
    tmp2 = fabs(tmp);
    if (fval <= tmp2) {
      fval = tmp2;
      k = j;
      tmp3 = 1.0;
      if (tmp < 0.0) tmp3 = -1.0;
    }
  }
  g[k] = tmp3 * (3.0 - 4.0 * x[k]);
  if (k > 0) g[k - 1] = -tmp3;
  if (k < (n - 1)) g[k + 1] = -tmp3;
  *f = fval;
  return 0;
}

template<class T> 
int test29f11(T *f, T *g, T *x, int n) {
    
  T fval = 0.0;
  T fa   = 0.0;
  T tmp3 = 1.0;
  int i  = 0;
    
  for (int j = 0; j<n; j++) 
    g[j] = 0.0;
    
  for (int ka = 0; ka < (2 * n - 2); ka++) {
    i = ka / 2;
    if ((ka % 2)==0) {
      fa   = x[i] + x[i + 1] * ((5 - x[i + 1]) * x[i + 1] - 2) - 13;
      (fa < 0) ?  tmp3 = -1.0 : tmp3 = 1.0;
      g[i] = g[i] + tmp3;
      g[i + 1] = g[i + 1] + (10 * x[i + 1]- 3 * x[i + 1] * x[i + 1] - 2) * 
	tmp3;
    }
    else {
      fa   = x[i] + x[i+1]*((1+x[i+1])*x[i+1]-14) - 29;
      (fa < 0) ?  tmp3 = -1.0 : tmp3 = 1.0;
      g[i] = g[i] + tmp3;
      g[i + 1] = g[i + 1] + (2 * x[i + 1] + 3 * x[i + 1] * x[i + 1] - 14) *
	tmp3;
    }
    fval = fval + fabs(fa);
  }    
  *f = fval;
  return 0;    
}

template<class T> 
int test29f22(T *f, T *g, T *x, int n) {
    
  T fval = 0.0;
  T fa   = 0.0;
  T U    = 1.0 / ((T) (static_cast<int>(n)) + 1.0);
  T V    = 0.0;
  T tmp  = 0.0;
  T tmp3 = 1.0;
  int k  = 1;
    
  for (int j = 0; j<n; j++)
    g[j] = 0.0;
    
  for (int ka = 0; ka < n; ka++) {
    V   = ( ((T) static_cast<int> (ka)) + 1.0) * U;
    tmp = (x[ka] + V + 1.0);
    fa  = 2 * x[ka] + 0.50 * U * U * tmp * tmp * tmp;
    if (ka > 0) fa = fa - x[ka - 1];
    if (ka < (n-1)) fa = fa - x[ka + 1];
    if (fval < fabs(fa) ){
      k  = ka;
      tmp3 = 1.0;
      if (fa < 0) tmp3 = -1.0;
      fval = fabs(fa);
    }
  }
  V    = ((T) static_cast<int>(k)) * U;
  tmp  = (x[k] + V + 1);
  g[k] = (2.0 + 1.5 * U * U * tmp * tmp) * tmp3;
  if (k > 0) g[k - 1] = -tmp3;
  if (k < (n - 1)) g[k + 1] = -tmp3;
  *f = fval;
  return 0;
}

template<class T>
int yurirosen(T *f, T *g, T *x, int n)
{
  // based on Chebyshev polynomials
  // x_{i+1} = 2x_i^2 - 1 = T_2(x_i) = T_{2^i}(x_1) = cos(2^i arccos(x_1))

  for (int i = 0; i < n; i++) 
    g[i] = 0.0;	

  *f = pow((1 - x[0]), 2.0)/4; // the 1/4 gives initial value 1 and prevents method
  // skipping over the nasty prescribed path
  g[0] = (x[0] - 1) / 2;

  for (int i = 0; i < (n - 1); i++) {
    *f = *f + pow((1 + x[i + 1] - 2 * pow(x[i], 2.0)), 2.0);
    T r = 1 + x[i + 1] - 2 * pow(x[i], 2.0);
    g[i+1] = g[i + 1] + 2 * r;
    g[i] = g[i] - 8 * x[i] * r;
  }
  return 0;
}

/*
  % Nesterov's Chebyshev-Rosenbrock SMOOTH version
  pars.nvar = n;
  pars.fgname = 'yurirosen';
  pars.title = sprintf('Nesterov-Chebyshev-Rosenbrock Smooth, n=%d', n);
  pars.optval = 0;
  options.x0 = ones(n,1);
  options.x0(1) = -1;
  options.maxit = 10^(n-1);  % note!
*/

template<class T>
int yurirosen_ns1(T *f, T *g, T *x, int n)
{
  /*n = pars.nvar;
    % this is version 1 of Yuri's nonsmooth Chebyshev Rosenbrock, but there are
    % two versions of this, according to the choice of the initial term
    % f = abs(1-x(1))/4;       % in this case the dimension of the V space is 2
    % g(1) = sign(x(1)-1)/4;   % (not sure if f is regular at the solution)
    % BFGS does well for n=4 50% of time, and occasionally for n=5
    % we go with the 2nd variant of version 1 in the paper, because it
    % illustrates the difficulty BFGS can have when there is a nontrivial
    % U space and it takes a long time to traverse along it
  */
	
  for (int i = 0; i < n; i++) 
    g[i] = 0.0;	
	
  *f = pow((1 - x[0]), 2.0) / 4;  // in this case the U and V spaces both have dim 1
  g[0] = -(1 - x[0]) / 2;      // problem is much more difficult: cannot solve n=4
  // accurately, probabaly because of rounding
	
  for (int i = 0; i < (n - 1); i++) {
    T y = 1 + x[i + 1] - 2 * pow(x[i], 2.0);
    *f = *f + fabs(y);
    T r;
    if (y > 0) {
      r = 1.0;
    }
    else if (y < 0) {
      r = -1.0;
    }
    else {
      r = 0.0;
    }
    g[i+1] = g[i + 1] + r;
    g[i] = g[i] - 4.0 * x[i] * r;
  }
  return 0;
}
    
/*
  % Nesterov's Chebyshev-Rosenbrock nonsmooth version 1
  pars.fgname = 'yurirosen_ns1';
  pars.title = 'Nesterov-Chebyshev-Rosenbrock Nonsmooth#1'; 
  pars.nvar = n;
  pars.optval = 0;
  pars.title = 'Nesterov-Chebyshev-Rosenbrock 1';
  if nargin >= 2 % if not, don't set options.x0
  options.x0 = scalex0*ones(n,1);  % BFGS fails immediately if scalex0 is 1, 
  but GS does not
  options.x0(1) = -scalex0;
  end
  options.H0 = eye(n); % so does not do initial scaling for BFGS //turns off scaling
  options.fvalquit = 1e-15;
  options.maxit = 10^(n+1); % !!!  for BFGS
  % options.phasenum = [1 10 6]; % for HANSO
  % options.phasemaxit = [10^n 10^n 1000]; % for HANSO
*/

template<class T>
int yurirosen_ns2(T *f, T *g, T *x, int n){
  /*%%%% NOTE THE 2, MORE INTERESING OF THE TWO,
    % because for this version, there are 2^(n-1) Clarke stationary
    % points, all of which are attractors for BFGS
  */
	
  for (int i = 0; i < n; i++) 
    g[i] = 0.0;	
	
  *f = fabs(1 - x[0]) / 4;	
  T k;
  if ((x[0] - 1) > 0) 
    k = 1;
  else if ((x[0] - 1) < 0)
    k = -1;
  else
    k = 0.0;
  g[0] = k / 4;

  for (int i = 0; i < (n - 1); i++) {
    T y = 1 + x[i + 1] - 2 * fabs(x[i]);
    *f = *f + fabs(y);
    T r;
    if (y > 0)
      r = 1.0;	    
    else if (y < 0)
      r = -1.0;
    else
      r = 0.0;

    g[i+1] = g[i+1] + r;

    T l;
    if (x[i] > 0) {
      l = 1.0;
    }
    else if (x[i] < 0) {
      l = -1.0;
    }
    else {
      l = 0.0;
    }

    g[i] = g[i] - 2 * l * r;
  }
  return 0;
}

/*
  % Nesterov's Chebyshev-Rosenbrock nonsmooth version 2
  pars.fgname = 'yurirosen_ns2';
  pars.title= 'Nesterov-Chebyshev-Rosenbrock Nonsmooth#2';
  %%%% NOTE THE 2, MORE INTERESING OF THE TWO,
  % because for this version, there are 2^(n-1) Clarke stationary
  % points, all of which are attractors for BFGS
  pars.nvar = n;
  pars.optval = 0;
  if nargin >= 2 % otherwise don't set options.x0
  options.x0 = scalex0*ones(n,1);  % BFGS fails immediately if scalex0 is 1, 
  but GS does not
  options.x0(1) = -scalex0;
  end
  options.H0 = eye(n); % so does not do initial scaling for BFGS
  options.maxit = 10^(n+1); % !!!  for BFGS
  options.fvalquit = 1e-15;
  % options.phasenum = [1 10 6]; % for HANSO
  % options.phasemaxit = [10^n 10^n 1000]; % for HANSO
  pars.title = 'Nesterov-Chebyshev-Rosenbrock'; % omit the 2
  pars.varytitle = pars.title;
*/

#endif // _FUNCTIONS_HPP_
