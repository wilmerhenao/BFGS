#ifndef _F10_HPP_
#define _F10_HPP_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <qd/dd_real.h>

template<class T> void chained_CB3v1(T *f, T *g, T *x, int n);
template<class T> void chained_CB3v2(T *f, T *g, T *x, int n);
template<class T> void chained_crescent1(T *f, T *g, T *x, int n);
template<class T> void chained_crescent2(T *f, T *g, T *x, int n);
template<class T> void chained_LQ(T *f, T *g, T *x, int n);
template<class T> void chained_mifflin2(T *f, T *g, T *x, int n);
template<class T> void gen_brownfunc2(T *f, T *g, T *x, int n);
template<class T> void gen_maxhilbert(T *f, T *g, T *x, int n);
template<class T> void gen_maxq(T *f, T *g, T *x, int n);
template<class T> void nactfaces(T *f, T *g, T  *x, int n);

template<class T> 
void chained_CB3v1(T *f, T *g, T *x, int n){
    /* printthis:
       % Chained CB3 I (convex, non-smooth func):
       % f = 
       % 0 parameters: w, w
       % show: n, n
    */
    T y;
    *f = 0.0;
    for (int i = 0; i<n; i++)
	g[i] = 0.0;

    for (int i = 0; i < (n-1); i++) {
	T a = pow(x[i], 4) + pow(x[i + 1], 2);
	T b = pow((2 - x[i]), 2) + pow((2 - x[i + 1]), 2);
	T c = 2 * exp(-x[i] + x[i + 1]);
	(a >= b) ? y = a : y = b;
	if (y < c)
	    y = c;
	if (y == a) {
	    g[i]   = g[i] + 4 * pow(x[i], 3);
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

}

template<class T> 
void chained_CB3v2(T *f, T *g, T *x, int n){
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
	a = a + pow(x[i], 4) + pow(x[i + 1], 2);
	b = b + pow((2 - x[i]),2) + pow((2 - x[i + 1]), 2);
	c = c + 2 * exp(-x[i] + x[i + 1]);
    }    
    (a >= b) ? *f = a : *f = b;
    if (*f < c)
	*f = c;

    if (*f == a) {
	for (int i = 0; i < (n-1); i++) {
	    g[i]   = g[i] + 4 * pow(x[i], 3);
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

}

template<class T> 
void chained_crescent1(T *f, T *g, T *x, int n)
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
	t2 = t2 + pow(x[i], 2) + pow((x[i + 1] - 1), 2) + x[i + 1] - 1;
	t3 = t3 - pow(x[i], 2) - pow((x[i + 1] - 1), 2) + x[i + 1] + 1;
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

}

template<class T> 
void chained_crescent2(T *f, T *g, T *x, int n){
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
	T t2 = pow(x[i], 2) + pow((x[i + 1] - 1),2) + x[i + 1] - 1;
	T t3 = -pow(x[i], 2) - pow((x[i + 1] - 1),2) + x[i + 1] + 1;
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

}

template<class T> 
void chained_LQ(T *f, T *g, T *x, int n)
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
	T b = a + ( pow(x[i], 2) + pow(x[i + 1], 2) - 1 );
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

}

template<class T> 
void chained_mifflin2(T *f, T *g, T *x, int n)
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

}

template<class T> 
void gen_brownfunc2(T *f, T *g, T *x, int n)
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
	T c = pow(x[i], 2) + 1;
	T d = pow(x[i+1], 2) + 1;
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

}

template<class T> 
void gen_maxhilbert(T *f, T *g, T *x, int n){
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
	*f = *f + x[j] / (j + 1);

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
	    t2 = t2 + x[j] / (i + j + 1);
	    
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
	g[j] = t3 / (k + j + 1);
    }

}

template<class T> 
void gen_maxq(T *f, T *g, T *x, int n)
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
}

template<class T> 
void nactfaces(T *f, T *g, T *x, int n)
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

}

#endif // _F10_HPP_
