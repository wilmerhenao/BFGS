
export LD_LIBRARY_PATH=${PWD}:${PWD/'/csrc'/'/lib'}

#set library paths

make clean
make
make mytest

# Remove this exit 0 line later when you are running the whole thing

./mytest yurirosen 5 ../testresults/yurirosen/yurirosen_5.txt

exit 0

./mytest chained_mifflin2 10 ../testresults/F10/chained_mifflin2_10.txt 
./mytest nactfaces 50 ../testresults/F10/nactfaces_50.txt

exit 0

for n in 10 #50 200
  do
    
  for i in {1..1}
    do
    ./mytest chained_crescent1 $n ../testresults/F10/chained_crescent1_$n.txt
    ./mytest gen_brownfunc2 $n ../testresults/F10/gen_brownfunc2_$n.txt
    ./mytest nactfaces $n ../testresults/F10/nactfaces_$n.txt
    ./mytest chained_CB3v2 $n ../testresults/F10/chained_CB3v2_$n.txt
    ./mytest gen_maxq $n ../testresults/F10/gen_maxq_$n.txt
    ./mytest chained_CB3v1 $n ../testresults/F10/chained_CB3v1_$n.txt
    ./mytest chained_LQ $n ../testresults/F10/chained_LQ_$n.txt
    ./mytest gen_maxhilbert $n ../testresults/F10/gen_maxhilbert_$n.txt
    ./mytest chained_mifflin2 $n ../testresults/F10/chained_mifflin2_$n.txt
    echo $i
  done
  echo $n
done

for n in {2..10}
  do

  for i in {1..1000}
    do
    ./mytest yurirosen $n ../testresults/yurirosen/yurirosen_$n.txt
    echo $i
  done
done

for n in {2..8}
  do

  #MAXIT=$[10**($n-1)] for NS1, n=6,7

  for i in {1..1000}
    do
    ./mytest yurirosen_ns1 yurirosen_ns1_dd $n ../testresults/yurirosen/yurirosen_ns1_$n.txt
    ./mytest yurirosen_ns2 yurirosen_ns2_dd $n ../testresults/yurirosen/yurirosen_ns2_$n.txt
    echo $i
  done
done

for n in 10 50 200
  do

  for i in {1..10}
    do
    ./mytest test29f02 $n ../testresults/T10/test29f02_$n.txt
    ./mytest test29f05 $n ../testresults/T10/test29f05_$n.txt
    ./mytest test29f06 $n ../testresults/T10/test29f06_$n.txt
    ./mytest test29f22 $n ../testresults/T10/test29f22_$n.txt
    ./mytest test29f24 $n ../testresults/T10/test29f24_$n.txt
    ./mytest test29f11 $n ../testresults/T10/test29f11_$n.txt
    echo $i
  done
  echo $n
done

exit 0

