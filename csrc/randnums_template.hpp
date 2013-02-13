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
