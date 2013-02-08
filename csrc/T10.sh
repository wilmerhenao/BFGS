#set library paths
export LD_LIBRARY_PATH=/home/ubuntu/Desktop/BFGS/csrc:/home/ubuntu/Desktop/BFGS/lib

make clean
make
make mytest

for n in 10 50 200
do

for i in {1..10}
do
    ./mytest test29f02 test29f02_dd $n >> ../testresults/T10/test29f02_$n.txt
    ./mytest test29f05 test29f05_dd $n >> ../testresults/T10/test29f05_$n.txt
    ./mytest test29f06 test29f06_dd $n >> ../testresults/T10/test29f06_$n.txt
    ./mytest test29f22 test29f22_dd $n >> ../testresults/T10/test29f22_$n.txt
    ./mytest test29f24 test29f24_dd $n >> ../testresults/T10/test29f24_$n.txt
    ./mytest test29f11 test29f11_dd $n >> ../testresults/T10/test29f11_$n.txt
    echo $i
done
echo $n
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
