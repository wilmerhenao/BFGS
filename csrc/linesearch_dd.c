using namespace std;
#include <math.h>
#include "libmatrix_dd.h"
#include <stdio.h>
#include <iostream>

dd_real linesearch_ww 
(dd_real x[], dd_real *f, dd_real g[], dd_real d[], dd_real C1, dd_real C2,
int n, void(*testFunction)(dd_real*, dd_real*, dd_real*, int), int*nfeval, dd_real ftarget, int*exitflag) 
{
    dd_real alpha = 0;		/* lower bound on step length*/
    dd_real beta_max = 1e100;	/* upper bound on step length*/
    dd_real beta  = beta_max;	
    
    dd_real t = 1;               /* try step length 1 first*/
    dd_real tprev = 0;           /* variable to store prev. t*/
    
    double done = 0;
    int j;
    
    dd_real nbisect = 0;
    dd_real nexpand = 0;
    dd_real dnorm = vecnorm(d,n);
    
    /* MATLAB:
    nbisectmax = max(30, round(log2(1e5*dnorm))); % allows more if ||d|| big
    nexpandmax = max(10, round(log2(1e5/dnorm))); % allows more if ||d|| small
    */
    dd_real nbisectmax = ceil( log( 100000*dnorm )/log(2.0) );
    if (nbisectmax < 100) nbisectmax = 100;

    dd_real nexpandmax = ceil( log( 100000 / dnorm )/log(2.0) );    
    if (nexpandmax < 100) nexpandmax = 100;

    dd_real gtd  = 0;
    dd_real g0td;
    dd_real armijo_rhs_p, wwolfe_rhs, f0;
    f0 = *f;
    
    g0td = vecip(g,d,n);	/* first gradient * search direction d (g0 in overton matlab)*/
    
    armijo_rhs_p = g0td*C1;
    wwolfe_rhs   = g0td*C2;
    while (!done) {
        /* x = x0 + t*d   (next x to check at):
        /* this is the same as x = x + (t-tprev)*d
        /* because x is overwritten each time*/
        vpv(x,d,(t-tprev),n);
        
        testFunction(f,g,x,n);    /* evaluate f and g at new x*/
        *nfeval = *nfeval + 1;    /* counting number of function evaluations*/
        
        /* break out if f below target:*/
        if (*f < ftarget) {
            done = 1;
            *exitflag = 1;
            return t;
        }
        
        /* inner product of current g and d (result is gtd):*/
        gtd = vecip(g,d,n);
        if (*f > f0 + t*armijo_rhs_p) {	/* armijo fails, gone too far*/
            beta = t;
        }
        else if (gtd < wwolfe_rhs) {	/* weak wolfe fails, not gone far enough*/
            alpha = t;
        }
        else {							/* both conditions are ok, so quit:*/
            alpha = t;
            beta = t;
            done = 1;
        }
        
        tprev = t;					/* store t before changing it*/
        /* next t:*/
        if (beta < beta_max) {		/* bisection:*/
            if (nbisect < nbisectmax) {
                t = (alpha + beta) / 2; 
                nbisect+=1;
            }
            else {
                done = 1;
                *exitflag = -4;
            }
        }							
        else {						/* expansion:*/
            if (nexpand < nexpandmax) {
                t = 2*alpha;
                nexpand+=1;
            }
            else {
                done = 1;
                *exitflag = -5;
            }
        }
    } /* end while*/

    return t;
}
