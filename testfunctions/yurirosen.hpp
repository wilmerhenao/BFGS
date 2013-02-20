#ifndef _YURIROSEN_HPP_
#define _YURIROSEN_HPP_

#include<cmath>
#include <qd/dd_real.h>

template<class T> void yurirosen(T *f, T *g, T *x, int n);
template<class T> void yurirosen_ns1(T *f, T *g, T *x, int n);
template<class T> void yurirosen_ns2(T *f, T *g, T *x, int n);

template<class T>
    void yurirosen(T *f, T *g, T *x, int n)
    {
    // based on Chebyshev polynomials
    // x_{i+1} = 2x_i^2 - 1 = T_2(x_i) = T_{2^i}(x_1) = cos(2^i arccos(x_1))

	for (int i = 0; i < n; i++) 
	    g[i] = 0.0;	

	*f = pow((1 - x[0]),2)/4; // the 1/4 gives initial value 1 and prevents method
	// skipping over the nasty prescribed path
	g[0] = (x[0] - 1) / 2;

	for (int i = 0; i < (n - 1); i++) {
	    *f = *f + pow((1 + x[i + 1] - 2 * pow(x[i], 2)), 2);
	    T r = 1 + x[i + 1] - 2 * pow(x[i], 2);
	    g[i+1] = g[i + 1] + 2 * r;
	    g[i] = g[i] - 8 * x[i] * r;
	}
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
    void yurirosen_ns1(T *f, T *g, T *x, int n)
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
	
	*f = pow((1 - x[0]), 2) / 4;  // in this case the U and V spaces both have dim 1
	g[0] = -(1 - x[0]) / 2;      // problem is much more difficult: cannot solve n=4
	// accurately, probabaly because of rounding
	
	for (int i = 0; i < (n - 1); i++) {
	    T y = 1 + x[i + 1] - 2 * pow(x[i], 2);
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
    void yurirosen_ns2(T *f, T *g, T *x, int n){
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

#endif // _YURIROSEN_HPP_
