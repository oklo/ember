// Offline electron pressure response for ion-screening comparisons.
#include "ember/eos_component.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
int main() {
  try {
    ember::ElectronGas electrons;
    double rho,T,X,Y3,Y4,Z;
    int fields;
    while((fields=std::scanf("%lf %lf %lf %lf %lf %lf",&rho,&T,&X,&Y3,&Y4,&Z))==6) {
      if(!std::isfinite(rho+T+X+Y3+Y4+Z) || rho<=0 || T<=0
         || X<0 || Y3<0 || Y4<0 || Z<0 || std::abs(X+Y3+Y4+Z-1)>1e-12)
        throw std::domain_error("invalid baryonic state");
      ember::Composition c{};
      c.basis=ember::AbundanceBasis::baryon_mass;
      c.metal_inventory=ember::MetalInventory::gs98;
      c[ember::Species::H1]=X;c[ember::Species::He3]=Y3;
      c[ember::Species::He4]=Y4;c[ember::Species::Zrest]=Z;
      const auto e=electrons.eval(T,rho,c);
      const double ne=c.mu_elec_inv()*ember::constants::NA*rho;
      std::printf("%.17g %.17g %.17g %.17g\n",ne,e.P,e.dP_dlnRho,e.dP_dlnT);
    }
    if(fields!=EOF)return 2;
  } catch(const std::exception& e) {std::fprintf(stderr,"%s\n",e.what());return 1;}
}
