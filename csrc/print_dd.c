using namespace std;
#include <iostream>
#include <cstring>
#include <cstdio>
#include <qd/dd_real.h>

void print_iter_info(FILE *output, int it, dd_real * f, dd_real gnorm, int j, dd_real * q, dd_real * x, dd_real t,  int n) {

    cout << "it=" << it << ", f=" << *f << ", |g|=" << gnorm << ", t=" << t << ", j=" << j << ", q=" << *q << endl << endl;

//cout << it << " " << *f << endl;

}

void print_init_info(FILE *output, int n, dd_real ftarget, dd_real gnormtol, int maxit, int echo, int lm, const char*outputname, void(*testFunction)(dd_real*, dd_real*, dd_real*, int)) {
    
    printf("================================\n");
    if (!lm) {
		printf("BFGS DD on %s with \n", "Test Function");
	}
	else {
		printf("LM-BFGS DD on %s with \n", "Test Function");
	}
    printf("n        = %i \n",n);
    cout << "ftarget  = " << ftarget << endl;//printf("ftarget  = %f \n",ftarget);
    cout << "gnormtol = " << gnormtol << endl;//printf("gnormtol = %f \n",gnormtol);
    printf("maxit    = %i \n",maxit);
    printf("echo     = %i \n",echo);
    printf("output   = %s \n",outputname);
    printf("--------------------------------\n");

}

void print_final_info(FILE *output, int it, dd_real f, dd_real gnorm, int nfeval, int exitflag, double ttime) {
    
    const char * exitstr;
    
    switch (exitflag) {
        case 1:
            exitstr = "f < ftarget";
            break;
            
        case 2:
            exitstr = "gnorm < gnormtol";
            break;
            
        case 3:
            exitstr = "f(x0) < ftarget";
            break;
            
        case 7:
            exitstr = "optcertval < taud";
            // exitstr = "x in conv{g_1,g_2,...,g_J} with ||x|| < taud";
            break;
            
        case -1:
            exitstr = "maximal # iterations reached";
            break;
            
        case -2:
            exitstr = "gtp > 0, i.e. non-descent dir. encountered (rounding)";
            break;
            
        case -4:
            exitstr = "linesearch: nbisect >= nbisectmax. Wolfe conds not satisfied.";
            break;

        case -5:  
            exitstr = "linesearch: nexpand >= nexpandmax. Wolfe conds not satisfied.";
            break;
            
        case -8:
			exitstr = "NaN encountered. Exiting with best found f,x and g.";
			break;
            
        case 0:
            exitstr = "exitflag was never changed";
            break;
            
        default:
            exitstr = "exited for some other reason";
            printf("exitflag = %i \n",exitflag);
            break;
    }
    
    printf("--------------------------------\n");
    printf("Stopped after\n",it);
    printf("%i iters and \n",it);
    printf("%i function evaluations\n",nfeval);  
    printf("%1.3e seconds \n",ttime);
    printf("Because \n",nfeval);
    printf("   %s",exitstr);
    printf("\n");
    printf("Found: \n");
    cout << "f       = " << f << endl;
    cout <<"norm(g) = " << gnorm << endl;
    printf("================================\n");

//cout << it << " " << f << " " << ttime << endl;
//cout << "gnorm=" << gnorm << endl;
}

void print_mat(dd_real A[], int m, int n, char * str) {
    printf("%s = \n",str);
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
	    cout << A[i*n+j] << "  ";
        }
        printf("\n");
    }
}

void print_vec(dd_real * v, int n, char * str) {
    printf("%s = \n",str);
    for (int j=0; j<n; j++) {
        cout << "     " << v[j] << endl;
    }
}

void print_int_vec(int * v, int n, char * str) {
    printf("%s = \n",str);
    for (int j=0; j<n; j++) {
        printf("     %i \n",v[j]);
    }
}

void print_str(char * str) {
    printf("%s = \n",str);
}

void print_double(char * str, dd_real num) {
    cout << str << " = " << num << endl;
}

void print_int(char * str, int num) {
    printf("%s = %i \n",str,num);
}

void print_gs0(dd_real r, int k, dd_real f, dd_real gnorm, dd_real t, dd_real qpinfo[]) {
    cout << "r=" << r << ", k=" << k << ", f=" << f << ", |g|=" << gnorm << ", t=" << t << ", qps=" << qpinfo[0] << ", qpk=" << qpinfo[1] << ", qpP=" << qpinfo[2] << endl;
}

void print_gs1(dd_real r) {
    cout << "rad=" << r << ", ";
}

void print_gs2(int k, dd_real f) {
    cout << "k=" << k << ", f=" << f << endl;
}

void print_gs4(double info[]) {
    printf("Stopped because: \n");
    if (info[0] >= 7) printf("f < ftarget. \n");
    if (info[0] == -9.0) printf("f apparantly unbounded below. \n");
    if (info[0] == -2.0) printf("max iters reached. \n");
}

void print_gs5(dd_real f, double info[]) {
    cout << "GradSamp found: f=" << f << ", info(1)=" << info[0] << endl;
}
