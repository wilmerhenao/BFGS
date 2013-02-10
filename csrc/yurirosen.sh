#set library paths
BFGS_HOME=~/Documents/thesis/BFGS
export LD_LIBRARY_PATH=$(BFGS_HOME)/BFGS/csrc:$(BFGS_HOME)/BFGS/lib

make clean
make
make mytest

#x0 = {-1, 1, 1..1}

for n in {2..10}
do

for i in {1..1000}
do
    ./mytest yurirosen yurirosen_dd $n >> ../testresults/yurirosen/yurirosen_$n.txt
    echo $i
done

done
	
exit 0

#    n            = atoi(argv[3]);		//# variables
#    lm           = 0;				//lbfgs or not
#    m            = 7;//atoi(argv[5]);
#    J            = fmin(15,ceil(2*n/3));    
#    taux         = 1e-16;			//stpsize tolerance
#    taud         = 1e-4;
#    ftarget      = -1e100;//atoi(argv[4]);//	//fvalquit
#    gnormtol     = 0;
#    maxit        = 10e8;//atoi(argv[4]);//
#    echo         = 1;				//prtlevel
#    datafilename = "stddump.txt";
#    datafilename_dd = "stddump_dd.txt";

#    taux_dd         = 1e-16;			//stpsize tolerance
#    taud_dd         = 1e-4;
#    ftarget_dd      = -1e100;//atoi(argv[5]);//	//fvalquit
#    gnormtol_dd     = 0.0;
