#ifndef _BFGS_TEMPLATE_HPP_
#define _BFGS_TEMPLATE_HPP_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include "linesearch_template.hpp"
#include "quasinewt_updates_template.hpp"
#include "print_template.hpp"
#include "libmatrix_template.hpp"
#include <qd/dd_real.h>

#define debug 0

template<class T>
void bfgs(T x[], T * fopt, const int n, const int lm, const int m, 
	  const T ftarget, const T gnormtol, const int maxit, const int J,
          const T taux, const T taud, const int echo, 
          void(*testFunction_dd)(T*, T*, T*, int), const char * datafilename, 
          double[]);


/* BFGS MAIN ALGORITHM: */
template<class T>
void bfgs(T x[], T * fopt, const int n, const int lm, const int m, 
          const T ftarget, const T gnormtol, const int maxit, const int J,
          const T taux, const T taud, const int echo, 
          void(*testFunction_dd)(T*, T*, T*, int), const char * datafilename, 
          double info[]){
    
  /* ============ INITIALIZATION =================== */
    
  const double C1 = 0.0001;
  const double C2 = 0.9;
    
  clock_t t1, t2;
  t1 = clock();
    
  int absecho = abs(echo);
    
  /*
    echo   = 0: no output, no file dump
    = 1: print init and final data to +
    = 2: print iter info to +
    = -k: as k, but to - (only for above 0,1,2)
    = -9: print final data (minimal) to -
    where + is the screen and - is the datafile
  */
    
  FILE *output = stdout;
  const char *outputname = "stdout";

  if (echo < 0) {
    output = fopen(datafilename, "a" );
    outputname = datafilename;
  }
  else {
    output = stdout;
    outputname = "stdout";
  }

  /* ============= ALLOCATE memory: ================= */
    
  int n1, n2, nm, m1;
    
  /* Common to both BFGS and LBFGS */
  T *g = new T[n];
  T *p = new T[n];
  T *s = new T[n];
  T *y = new T[n];
  T *d = new T[n];
    
  /* Only allocate what is actually needed: */
  if (!lm) {
    n1 = n;
    n2 = n * n;
    nm = 1;
    m1 = 1;
  }
  else {
    n1 = 1;
    n2 = 1;
    nm = n*m;
    m1 = m;
  }
    
  /* BFGS specific arrays: */
  T *q = new T[n1];
  T *H = new T[n2];
    
  /* LBFGS specific arrays: */
  T *S   = new T[nm];
  T *Y   = new T[nm];
  T *rho = new T[m1];
  T *a   = new T[m1];
    
  /* doubles:*/
  T t, gnorm, bgnorm, gtp, fval, fprev, dtmp;
  /* double pointers: */
  T *f;
  f = &(fval);
    
  /* integers:*/
  int it = 0, done = 0, j, ol = 1, cs = 0, tmp, bgnormidx = 1;
  /* integer pointers: */
  int nfevalval = 0, exitflagval = 0;
  int *nfeval    = &(nfevalval);
  int *exitflag  = &(exitflagval);
    
  /* VARS FOR QP STOPPING CRITERION */
  const T R = 10;
  int qpmaxit    = 100;
  int jcur       = 0;
  int oldestg    = 2;
  T snorm   = 0;

  T qpoptval = taud + 100;
  T * qpoptvalptr;
  qpoptvalptr = &(qpoptval);
    
  T * G      = new T[J * n];
  T * qpx0   = new T[J];
  T * qpinfo = new T[3];
    
  /* ============ INITIALIZATION END ===============*/
    
  /* ============ BEFORE MAIN LOOP: ================ */
  testFunction_dd(f, g, x, n);

  *nfeval = *nfeval + 1;
    
  /* calculate gnorm:*/
  gnorm  = vecnorm<T>(g, n);
  bgnorm = gnorm;

  /* first search direction  p = -H*g (BFGS) and p = -g (LBFGS) */
  if(!lm) {
    /* initialize H to the identity matrix:*/
    mat_set_eye(H, n, n);
    mxv<T>(p, H, g, -1.0, 0.0, n, n);
  }
  else {
    vcopyp<T>(p, g, -1.0, n);
  }
    
  if (echo > 0) 
    print_init_info<T>(output, n, ftarget, gnormtol, maxit, echo, lm, 
		       outputname, testFunction_dd);
    
  if (*f < ftarget) {
    done = 1;
    *exitflag = 3;
  }
    
  if (echo==2) 
    print_iter_info<T>(output, it, f, gnorm, jcur, qpoptvalptr, x, t, n);
    
  /* ================= MAIN LOOP: =================== */
  while (!done) {
    /* increase iteration counter:*/
    it++;        
    /* gtp = g'*p*/
    gtp = vecip<T>(g,p,n);        
    if (gtp > 0) {
      done = 1;
      *exitflag = -2;
    }
        
    /* need previous x and g for update later (for s and y), */
    /* so copy them before the linsearch overwrites them:*/
    vcopyp<T>(s,x,-1.0,n);
    vcopyp<T>(y,g,-1.0,n);
        
    /* Copy f before overwritten by linesearch in case of NaN */
    fprev = *f;
        
    /* line search:*/
    t = linesearch_ww<T>(x, f, g, p, C1, C2, n, testFunction_dd, nfeval, 
			 ftarget, exitflag);

    /* If f is NaN, exit with best found f,x,g so far */
    if (isnan(*f)) {
      *exitflag = -8;
      /* at this point, s=-x, so take x from s (same with y/g) */
      vcopyp<T>(x,s,-1.0,n);
      vcopyp<T>(g,y,-1.0,n);
      /* best f found is stored in fprev. Exit with that: */
      *f = fprev;
    }

    if (*f > fprev)
      *f = fprev;
       
    /* calculate s and y:
       /* before these calls, s and y contain
       /* previous -x and prev. -g, so e.g. s = x - xprev is s = s + x */
    vpv<T>(s,x,1,n);
    vpv<T>(y,g,1,n);
    gnorm = vecnorm<T>(g,n);

    /* ============= QP STOPPING CRITERING =========== */
    //        snorm = vecnorm(s,n);

    //        if (snorm > taux) {
    //            jcur = 1;
    //            /* set first row in G equal to g: */
    //            vcopy(G, g, n);
            
    /* for initial point: */
    //            bgnorm    = gnorm;
    //            bgnormidx = 1;
    //        }
    //        else {
    //            jcur = jcur + 1;
    //            if (jcur > J) jcur = J;
            
    /* write new g in row of oldest g here */
    //            vcopy(G+(oldestg-1)*n, g, n);
            
    /* for initial point: */
    //            if (gnorm < bgnorm) {
    //                bgnorm    = gnorm; 
    //                bgnormidx = oldestg;
    //            }
            
    /* change oldestg to new location of oldest g: */
    //            oldestg = (oldestg % J) + 1;            
    //        }
    //        if (jcur > 1) {
    /* set initial point: */
    //            dtmp = 1 / (jcur-1+R);
    //            for(j=0;j<jcur;j++) qpx0[j] = dtmp;
    //            qpx0[bgnormidx-1] = R*dtmp;
            
    /* call qpsolv here. General call is
     * qpsubprob(G[jcur*n], x[jcur], double * q, info[2], maxit, jcur, n)
     * so in this case: */
    //            qpsubprob(G, qpx0, d, qpoptvalptr, qpinfo, qpmaxit, jcur, n);
            
    //            if (qpinfo[0] > 0) *qpoptvalptr = taud+100;
             
    //        }
    /* ========== END QP STOPPING CRITERION ========== */
        
        
    /* QUASINEWTON update:*/
    if (!lm) {
      /* BFGS update: */
      update_bfgs<T>(H, p, g, s, y, q, n);
      /* BFGS update end */
    }
    else {
      /* LBFGS update:*/
      if(it <= m) cs++;
      tmp = (ol-1) * n;
      vcopy<T>(S + tmp, s, n);
      vcopy<T>(Y + tmp, y, n);
      rho[ol-1] = 1.0 / vecip<T>(y, s, n);
      update_lbfgs<T>(p, S, Y, rho, a, g, cs, ol, m, n);
      ol = (ol % m) + 1;
      /* LBFGS update end */
    }
        
    /* print iteration info:*/
    if (echo==2) print_iter_info<T>(output, it, f, gnorm, jcur, qpoptvalptr,
				    x,t,n);
        
    /* check convergence here*/
    if (it >= maxit) {
      *exitflag = -1;
    }
    if (*f < ftarget) {
      *exitflag = 1;
    }
    if (gnorm < gnormtol) {
      *exitflag = 2;
    }
    if (*qpoptvalptr < taud) {
      *exitflag = 7;
    }
        
    /* if exitflag was changed: exit main loop: */
    if (*exitflag != 0) done = 1;

  } /*end while*/
    
  t2 = clock();
  double ttime = (double) ((double)(t2 - t1)/ (double)CLOCKS_PER_SEC );
    
  if (echo>0) 
    print_final_info<T>(output, it, *f, gnorm, *nfeval, *exitflag, ttime);
  if (echo<0) 
    fclose(output);

  *fopt = *f;
    
  /* Gather rundata in the double array 'info' */
  info[0] = (double) (*nfeval);
  info[1] = (double) (*exitflag);
  info[2] = (double) (lm);
  info[3] = (double) (ttime);
}

#endif // _BFGS_TEMPLATE_HPP_
