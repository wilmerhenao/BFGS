#include <stdlib.h>
/* return (uniform) random double between a and b: */
double rand_double(double a, double b) {
	
	return ( a+(b-a)*( (double)rand() / (double)RAND_MAX) ); 
}

void rand_double_vec(double v[], int n, double a, double b) {
    int j;
    double tmp = b - a;
    for(j=0; j<n; j++) {
        v[j] = ( a+tmp*( (double)rand() / (double)RAND_MAX) );
    }
    
}
