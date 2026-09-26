// Ideal-electron contribution for offline pure-He phase comparisons.
// The separate EOSFI22 probe supplies ideal ions and interaction terms.
#include "ember/eos_component.hpp"
#include <cmath>
#include <cstdio>
#include <exception>
int main() {
  try {
    ember::Composition c{};c.basis=ember::AbundanceBasis::baryon_mass;
    c[ember::Species::He4]=1;
    ember::ElectronGas electrons;
    double rho,T;int fields;
    while((fields=std::scanf("%lf %lf",&rho,&T))==2) {
      const auto e=electrons.eval(T,rho,c);
      std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
          e.E-T*e.S,e.E,e.S,e.P,e.dE_dlnT/T,e.dP_dlnT,e.dP_dlnRho);
      std::fflush(stdout);
    }
    if(fields!=EOF)return 2;
  } catch(const std::exception& e) {std::fprintf(stderr,"%s\n",e.what());return 1;}
}
