// Evaluate a selected physical metal-atmosphere response at independent
// source coordinates, without relaxing a stellar model.
#include "ember/atmosphere_metal_response.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/eos_variable_metal.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace ember;
int main(int argc,char** argv) {
  if(argc!=4)return 2;
  try {
    VariableMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    CompositionAtmosphereGrid reference(eos,argv[2],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    MetalResponseAtmosphere atmosphere(eos,reference,std::filesystem::path(argv[3]),
        MetalResponseAtmosphere::Approximation::separable_gs98_response,{1e-5});
    std::cout<<std::setprecision(17);
    double x,y,z,t,lg;
    while(std::cin>>x>>y>>z>>t>>lg)try {
      auto c=solar_scaled(x,z);c.X[1]=y;c.X[2]-=y;
      c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      const auto s=atmosphere.eval(t,std::pow(10.,lg),c);
      std::cout<<"{\"T\":"<<s.T<<",\"Pgas\":"<<s.Pgas<<",\"P\":"<<s.P
        <<",\"rho\":"<<s.rho<<",\"dlnT_dlnTeff\":"<<s.dlnT_dlnTeff
        <<",\"dlnT_dlng\":"<<s.dlnT_dlng<<",\"dlnP_dlnTeff\":"<<s.dlnP_dlnTeff
        <<",\"dlnP_dlng\":"<<s.dlnP_dlng<<"}\n";
    }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
