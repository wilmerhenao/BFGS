#include <cmath>
#include "libmatrix_template.hpp"
#include <cstdio>

double linesearch_ww 
(double x[], double *f, double g[], double d[], double C1, double C2,
int n, void(*testFunction)(double*, double*, double*, int), int*nfeval, double ftarget, int*exitflag) 
{    
    double alpha = 0;		/* lower bound on step length*/
    double beta_max = 1e100; 	/* upper bound on step length*/
    double beta  = beta_max;	
    
    double t = 1;               /* try step length 1 first*/
    double tprev = 0;           /* variable to store prev. t*/
    
    double done = 0;
    int j;
    
    int nbisect = 0;
    int nexpand = 0;
    double dnorm = vecnorm<double>(d,n);
    
    /* MATLAB:
    nbisectmax = max(30, round(log2(1e5*dnorm))); % allows more if ||d|| big
    nexpandmax = max(10, round(log2(1e5/dnorm))); % allows more if ||d|| small
    */
    int nbisectmax = (int) ceil( log2( 100000*dnorm ) );
    if (nbisectmax < 100) nbisectmax = 100;

    int nexpandmax = ceil( log2( 100000 / dnorm ) );    
    if (nexpandmax < 100) nexpandmax = 100;

    double gtd  = 0;
    double g0td;
    double armijo_rhs_p, wwolfe_rhs, f0;
    f0 = *f;
    
    g0td = vecip<double>(g,d,n);/* first gradient * search direction d (g0 in overton matlab)*/
    
    armijo_rhs_p = g0td*C1;
    wwolfe_rhs   = g0td*C2;

    while (!done) {
        /* x = x0 + t*d   (next x to check at):
        /* this is the same as x = x + (t-tprev)*d
        /* because x is overwritten each time*/
        vpv<double>(x,d,(t-tprev),n);
        
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
                nbisect++;
            }
            else {
                done = 1;
                *exitflag = -4;
            }
        }							
        else {						/* expansion:*/
            if (nexpand < nexpandmax) {
                t = 2*alpha;
                nexpand++;
            }
            else {
                done = 1;
                *exitflag = -5;
            }
        }
    } /* end while*/

    return t;
}
