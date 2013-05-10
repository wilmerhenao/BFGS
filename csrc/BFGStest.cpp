#include<iostream>
#include"quasinewt_updates_template.hpp"

int main(){
  double B[9] = {0,1,2,3,4,5,6,7,8};
  double s[3] = {4,5,6};
  double y[3] = {3,4,5};
  double q[3] = {1,2,3};
  int n = 3; 
  update_bfgs_B(B,s,y,q,n);
  for(int i = 0; i < n*n; i++)
    std::cout<< B[i] << " " << std::endl;
  return 0;
}
