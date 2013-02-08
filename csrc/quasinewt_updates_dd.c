#include "libmatrix_dd.h"

/* BFGS update */
void update_bfgs (dd_real H[], dd_real p[], dd_real g[], dd_real s[], dd_real y[], dd_real q[], int n)  {    
	
    dd_real rho;
    int j, in;
    
    rho = 1 / vecip(s,y,n);
    mxv(q,H,y,1.0,0.0,n,n);
    
    dd_real a = rho*vecip(y,q,n) + 1.0;
    
    /* reuse y to save memory and precompute a*s-q to save flops: */
    vcopyp(y, q, -1.0, n); /* y := -q */
    vpv(y, s, a, n);       /* y := a*s + y */
    /* so now: y = a*s - q */
    
    /* save flops by multiplying s by rho before the two rank one updates */
    vscal(s, rho, n);
    
    /* the two rank one updates (in total a rank two update): */
    /* NB: THERE IS ALSO A DIRECT RANK 2 UPDATE IN BLAS!! */
    mat_r1update(H, s, y, 1.0, n);
    mat_r1update(H, q, s, -1.0, n);
    
    /* compute next search dir p = -H*g: */
    mxv(p,H,g,-1.0,0.0,n,n);
}


/* LBFGS update */
void update_lbfgs (dd_real p[], dd_real S[],dd_real Y[], dd_real rho[], dd_real a[], dd_real g[], 
int cs, int ne, int m, int n) {
    /* newest s/y/rho is in col ne of S/Y/rho */    
    dd_real b, gam;
    int ci,i,j;
    
    gam = (1.0 / rho[ne-1]) / vecip(Y+n*(ne-1),Y+n*(ne-1),n);
    vcopy(p,g,n);
    
    /* FIRST RECURSION: */
    for (j=1; j<=cs; j++) 
    {
        i  = (ne + m - j) % m + 1;
        ci = (i-1)*n;
        
        /* a(i) = rho(i)*S(:,i)*p */
        a[i-1] = rho[i-1]*vecip(S+ci,p,n);
        
        /* p = p - a(i)*Y(:,i) */
        vpv(p, Y+ci, -a[i-1], n);
        
    } 
    /* FIRST RECURSION END */
    
    /* SCALING OF P: */
    vscal(p,gam,n);
    
    /* SECOND RECURSION: */
    for (j=cs; j>0; j--) {
        i = (ne - j + m) % m + 1;
        ci = (i-1)*n;
        
        /* b = rho(i)*Y(:,i)*p */
        b = rho[i-1]*vecip(Y+ci,p,n);
        
        /* p = p + (a(i)-b)*S(:,i) */
        vpv(p,S+ci,(a[i-1]-b),n);
        
    } /* SECOND RECURSION END */
    
    vscal(p,-1.0,n);
}
