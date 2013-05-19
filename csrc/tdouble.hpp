
#ifndef _TDOUBLE_HPP_
#define _TDOUBLE_HPP_

/* Cast to double. */
inline double t_double(const dd_real &a) {
  return a.x[0];
}

inline double t_double(const qd_real &a) {
  return a[0];
}

inline double t_double(const double &a){
  return a;
}

#endif // _TDOUBLE_HPP_
