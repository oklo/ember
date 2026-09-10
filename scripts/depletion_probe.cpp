// SPDX-License-Identifier: GPL-3.0-or-later
// Independent offline driver for the optional Fortran depletion interface.
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
extern "C" void ember_condense_(const double*,const double*,const double*,double*,const int*,int*);
int main() {
  std::array<double,92> bulk{},gas{};std::ifstream input("bulk.dat");
  for(auto& value:bulk)if(!(input>>value))return 1;
  const int n=92;double T,P;std::cout<<std::setprecision(17);
  while(std::cin>>T>>P) {
    int error;ember_condense_(&T,&P,bulk.data(),gas.data(),&n,&error);
    if(error)return 1;
    std::cout<<T<<' '<<P;for(double value:gas)std::cout<<' '<<value;std::cout<<'\n';
  }
  return std::cin.eof()?0:1;
}
