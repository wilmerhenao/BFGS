void print_iter_info(FILE *output, int it, double * f, double gnorm, int j, double * q, double * x, double t,  int n);
void print_init_info(FILE *output, int n, double ftarget, double gnormtol, int maxit, int echo, int lm, const char*outputname, void(*testFunction)(double*, double*, double*, int));
void print_final_info(FILE *output, int it, double f, double gnorm, int nfeval, int exitflag, double ttime);

void print_mat(double A[], int m, int n, char * str);
void print_vec(double * v, int n, char * str);
void print_int_vec(int * v, int n, char * str);
void print_str(char * str);
void print_double(char * str, double num);
void print_int(char * str, int num);

void print_gs0(double r, int k, double f, double gnorm, double t, double qpinfo[]);
void print_gs1(double r);
void print_gs2(int k, double f);
void print_gs4(double info[]);
void print_gs5(double f, double info[]);


