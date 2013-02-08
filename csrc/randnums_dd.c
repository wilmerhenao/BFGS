#include <stdlib.h>
#include <qd/dd_real.h>

/* return (uniform) random double between a and b: */
dd_real rand_double(dd_real a, dd_real b) {
	
	return ( a+(b-a)*( (dd_real)rand() / (dd_real)RAND_MAX) ); 
}

void rand_double_vec_dd(dd_real v[], int n, dd_real a, dd_real b) {
    int j;
    dd_real tmp = b - a;
    for(j=0; j<n; j++) {
        v[j] = ( a+tmp*( (dd_real)rand() / (dd_real)RAND_MAX) );
    }
    
}
