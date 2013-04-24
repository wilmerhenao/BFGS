/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BFGS tools Copyright (C) 2013  Wilmer Henao
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _QUASINEWT_UPDATES_TEMPLATE_HPP_
#define _QUASINEWT_UPDATES_TEMPLATE_HPP_

#include "libmatrix_template.hpp"
templte <typename T> void  update_bfgs_B(T [], T [], T [], T [], T [], T [], int );

template <class T> void update_bfgs (T [], T p[], T g[], T s[], T y[], T q[], 
				     int n);
template <class T> void update_lbfgs (T p[], T S[],T Y[], T rho[], T a[], T g[], 
				      int cs, int ne, int m, int n);
// BFGS update of the actual matrix

template<typename T>
void update_bfgs (T H[], T p[], T g[], T s[], T y[], T q[], int n){
  
}


/* BFGS update of the inverse matrix */
template<class T>
void update_bfgs (T H[], T p[], T g[], T s[], T y[], T q[], int n)  {    
  
  T rho;
  rho = 1 / vecip<T>(s, y, n);
  mxv<T>(q, H, y, 1.0, 0.0, n, n);
  
  T a = rho * vecip<T>(y, q, n) + 1.0;
  
  /* reuse y to save memory and precompute a*s-q to save flops: */
  vcopyp<T>(y, q, -1.0, n); /* y := -q */
  vpv<T>(y, s, a, n);       /* y := a*s + y */
  /* so now: y = a*s - q */
  
  /* save flops by multiplying s by rho before the two rank one updates */
  vscal<T>(s, rho, n);
  
  /* the two rank one updates (in total a rank two update): */
  /* NB: THERE IS ALSO A DIRECT RANK 2 UPDATE IN BLAS!! */
  mat_r1update<T>(H, s, y, 1.0, n);
  mat_r1update<T>(H, q, s, -1.0, n);
  
  /* compute next search dir p = -H*g: */
  mxv<T>(p, H, g, -1.0, 0.0, n, n);
}

/* LBFGS update */
template <class T>
void update_lbfgs (T p[], T S[],T Y[], T rho[], T a[], T g[], 
		   int cs, int ne, int m, int n) {
  /* newest s/y/rho is in col ne of S/Y/rho */    
  T b, gam;
  int ci, i, j;
    
  gam = (1.0 / rho[ne - 1]) / vecip<T>(Y + n * (ne - 1), Y + n * (ne - 1), n);
  vcopy<T>(p, g, n);
    
  /* FIRST RECURSION: */
  for (j = 1; j <= cs; j++) {
    i  = (ne + m - j) % m + 1;
    ci = (i - 1) * n;
    
    /* a(i) = rho(i)*S(:,i)*p */
    a[i - 1] = rho[i - 1] * vecip<T>(S + ci, p, n);
    
    /* p = p - a(i)*Y(:,i) */
    vpv<T>(p, Y + ci, -a[i - 1], n);
        
  } 
  /* FIRST RECURSION END */
    
  /* SCALING OF P: */
  vscal<T>(p, gam, n);
    
  /* SECOND RECURSION: */
  for (j = cs; j > 0; j--) {
    i = (ne - j + m) % m + 1;
    ci = (i - 1) * n;
        
    /* b = rho(i)*Y(:,i)*p */
    b = rho[i - 1] * vecip<T>(Y + ci, p, n);
        
    /* p = p + (a(i)-b)*S(:,i) */
    vpv<T>(p, S + ci,(a[i - 1] - b), n);
        
  } /* SECOND RECURSION END */
    
  vscal<T>(p, -1.0, n);
}

#endif // _QUASINEWT_UPDATES_TEMPLATE_HPP_
