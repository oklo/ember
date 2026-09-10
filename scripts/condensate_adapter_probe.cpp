// Offline regression probe for the depletion adapter's bulk-input contract.
// Input: 92 bulk element/H number ratios, then T[K], P[dyn/cm2], reset[0/1].
// Default matches MOLEQ: fixed cbase input, distinct mutable gas output.
// --mutable-bulk deliberately violates that contract as a robustness test.
#include <array>
#include <iomanip>
#include <iostream>
#include <string>

extern "C" void ember_condense_(const double*,const double*,const double*,
                                double*,const int*,int*);
int main(int argc,char** argv) {
  if(argc>2 || (argc==2 && std::string(argv[1])!="--mutable-bulk")) return 2;
  const bool mutable_bulk=argc==2;
  constexpr int n=92;
  std::array<double,n> original{},current{},gas{};
  for(auto& x:original) if(!(std::cin>>x)) return 2;
  current=original;
  double T,P; int reset;
  std::cout<<std::setprecision(17);
  while(std::cin>>T>>P>>reset) {
    if(reset!=0 && reset!=1) return 2;
    if(reset) current=original;
    int error;
    ember_condense_(&T,&P,current.data(),gas.data(),&n,&error);
    if(error) return 1;
    for(int j=0;j<n;++j) std::cout<<(j?" ":"")<<gas[j];
    std::cout<<'\n'; if(mutable_bulk) current=gas;
  }
}
