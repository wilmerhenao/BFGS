
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%  BFGS tools Copyright (C) 2013  Wilmer Henao
#%%  This program is free software: you can redistribute it and/or modify
#%%  it under the terms of the GNU General Public License as published by
#%%  the Free Software Foundation, either version 3 of the License, or
#%%  (at your option) any later version.
#%%
#%%  This program is distributed in the hope that it will be useful,
#%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
#%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#%%  GNU General Public License for more details.
#%%
#%%  You should have received a copy of the GNU General Public License
#%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


export LD_LIBRARY_PATH=${PWD}:${PWD/'/csrc'/'/lib'}

#set library paths

make clean
make
make mytest

# Remove this exit 0 line later when you are running the whole thing

#./mytest randomsq 3 ../testresults/randomsq/randomsq_10.txt
./mytest parabola 5 ../testresults/parabola/parabola_5.txt

exit 0

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

