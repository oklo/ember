#include "ember/eos_mixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  try {
    if(argc!=2)throw std::invalid_argument("metal_eos_probe family.dat");
    ember::MetalHelmholtzEos eos(argv[1],ember::HelmholtzTableEos::Mixture::allow_documented_proxy);
    double x,y,T,rho;std::cout<<std::setprecision(17);
    while(std::cin>>x>>y>>T>>rho) {
      auto c=ember::solar_scaled(x,.02);c.basis=ember::AbundanceBasis::baryon_mass;c.X[1]=y;c.X[2]-=y;
      c.metal_inventory=ember::MetalInventory::gs98;
      const auto e=eos.eval_with_derivatives(T,rho,c);const auto d=eos.composition_response(T,rho,c);
      const auto& s=e.state;
      std::cout<<x<<' '<<y<<' '<<T<<' '<<rho<<' '<<s.P<<' '<<s.E<<' '<<s.S<<' '<<s.chiRho<<' '<<s.chiT<<' '
               <<s.cp<<' '<<s.cv<<' '<<s.grad_ad<<' '<<s.delta<<' '<<eos.rho_from_PT(T,s.P,c,rho*1.01)
               <<' '<<d.dP[0]<<' '<<d.dP[1]<<' '<<d.dE[0]<<' '<<d.dE[1]<<'\n';
    }
    if(!std::cin.eof())throw std::invalid_argument("malformed input");
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
