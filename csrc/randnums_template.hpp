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

#ifndef _RANDNUMS_TEMPLATE_HPP_
#define _RANDNUMS_TEMPLATE_HPP_

#include <cstdlib>

// The first template function produces a uniform random number between 'a' and 'b' 
// The second template function produces a uniform distribution.  It is a whole
// vector of dimension n

// Notice that the original implementation called these function as "double"
template <class T> T rand_real(T a, T b);
template <class T> void rand_real_vec(T v[], int n, T a, T b);

template <class T>
T rand_real(T a, T b){

  return(a + (b - a) * ((T)rand() / (T)RAND_MAX));

}

template <class T>
void rand_real_vec(T v[], int n, T a, T b){
  int j;
  T tmp = b - a;
  for (j = 0; j < n; j++){
    v[j] = (a + tmp * ((T)rand() / (T)RAND_MAX));
  }
}

#endif // _RANDNUMS_TEMPLATE_HPP__
