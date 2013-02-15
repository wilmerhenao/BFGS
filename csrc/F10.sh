export LD_LIBRARY_PATH=${PWD}:${PWD/'/csrc'/'/lib'}

#set library paths

make clean
make
make mytest

# Remove this exit 0 line later when you are compiling the whole thing
# exit 0
for n in 10 50 200
do

for i in {1..10}
do
    ./mytest chained_crescent1 chained_crescent1_dd $n >> ../testresults/F10/chained_crescent1_$n.txt
    ./mytest gen_brownfunc2 gen_brownfunc2_dd $n >> ../testresults/F10/gen_brownfunc2_$n.txt
    ./mytest nactfaces nactfaces_dd $n >> ../testresults/F10/nactfaces_$n.txt
    ./mytest chained_CB3v2 chained_CB3v2_dd $n >> ../testresults/F10/chained_CB3v2_$n.txt

    ./mytest gen_maxq gen_maxq_dd $n >> ../testresults/F10/gen_maxq_$n.txt
    ./mytest chained_CB3v1 chained_CB3v1_dd $n >> ../testresults/F10/chained_CB3v1_$n.txt
    ./mytest chained_LQ chained_LQ_dd $n >> ../testresults/F10/chained_LQ_$n.txt
    ./mytest gen_maxhilbert gen_maxhilbert_dd $n >> ../testresults/F10/gen_maxhilbert_$n.txt
    ./mytest chained_mifflin2 chained_mifflin2_dd $n >> ../testresults/F10/chained_mifflin2_$n.txt
    echo $i
done
echo $n
done

exit 0

