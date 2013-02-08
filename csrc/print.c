#include <string.h>
#include <stdio.h>

void print_iter_info(FILE *output, int it, double * f, double gnorm, int j, double * q, double * x, double t,  int n) {

    printf("it=%i, f=%2.2e, |g|=%2.2e, t=%2.2e, j=%i, q=%2.2e \n",it,*f,gnorm,t,j,*q);

//printf("%i %2.16e\n", it, *f);

}

void print_init_info(FILE *output, int n, double ftarget, double gnormtol, int maxit, int echo, int lm, const char*outputname, void(*testFunction)(double*, double*, double*, int)) {
    
    printf("================================\n");
    if (!lm) {
		printf("BFGS D on %s with \n", "Test Function");
	}
	else {
		printf("LM-BFGS DD on %s with \n", "Test Function");
	}
    printf("n        = %i \n",n);
    printf("ftarget  = %f \n",ftarget);
    printf("gnormtol = %f \n",gnormtol);
    printf("maxit    = %i \n",maxit);
    printf("echo     = %i \n",echo);
    printf("output   = %s \n",outputname);
    printf("--------------------------------\n");

}

void print_final_info(FILE *output, int it, double f, double gnorm, int nfeval, int exitflag, double ttime) {
    
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
    printf("f       = %1.6e \n",f);
    printf("norm(g) = %1.6e \n",gnorm);
    printf("================================\n");

//printf("%i %1.16e %1.3e\n", it, f, ttime);
//printf("gnorm=%e\n", gnorm);
}

void print_mat(double A[], int m, int n, char * str) {
    printf("%s = \n",str);
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            printf("%1.7e  ",A[i*n+j]);
        }
        printf("\n");
    }
}

void print_vec(double * v, int n, char * str) {
    printf("%s = \n",str);
    for (int j=0; j<n; j++) {
        printf("     %1.3e \n",v[j]);
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

void print_double(char * str, double num) {
    printf("%s = %1.8e \n",str,num);
}

void print_int(char * str, int num) {
    printf("%s = %i \n",str,num);
}

void print_gs0(double r, int k, double f, double gnorm, double t, double qpinfo[]) {
    printf("r=%1.2e, k=%i, f=%1.2e, |g|=%1.2e, t=%1.2e, qps=%f, qpk=%f, qpP=%f \n",r,k,f,gnorm,t,qpinfo[0],qpinfo[1],qpinfo[2]);
}

void print_gs1(double r) {
    printf("rad=%1.3e, ",r);
}

void print_gs2(int k, double f) {
    printf("k=%i, f=%1.3e \n",k,f);
}

void print_gs4(double info[]) {
    printf("Stopped because: \n");
    if (info[0] >= 7) printf("f < ftarget. \n");
    if (info[0] == -9.0) printf("f apparantly unbounded below. \n");
    if (info[0] == -2.0) printf("max iters reached. \n");
}

void print_gs5(double f, double info[]) {
    printf("GradSamp found: f=%1.3e, info(1)=%f \n",f,info[0]);
}
