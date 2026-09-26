#include "ember/material_flux.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
  using namespace ember;
  std::cout<<std::setprecision(17);
  MaterialMatrix reduced{};
  while(std::cin>>reduced[0][0]) {
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)if(i || j)std::cin>>reduced[i][j];
    std::array<double,2> h{},lo{},hi{};double el,e0,s0,Ta,Tb,rho,area,dm;
    std::cin>>h[0]>>h[1]>>el>>e0>>s0>>Ta>>Tb>>rho>>area>>dm>>lo[0]>>lo[1]>>hi[0]>>hi[1];
    if(!std::cin)return 2;
    try {
      const auto full=full_material_mobility(reduced,h,el,e0,s0);
      const auto split=split_material_heat(full,std::sqrt(Ta)*std::sqrt(Tb),e0,s0);
      const auto f=material_face_flux(full,lo,hi,Ta,Tb,rho,area,dm,e0,s0);
      std::cout<<"{\"ok\":true,\"enthalpy\":["<<split.transported_enthalpy[0]<<','<<split.transported_enthalpy[1]
        <<"],\"conductivity\":"<<split.conductivity<<",\"species_rate\":["<<f.species_rate[0]<<','<<f.species_rate[1]
        <<"],\"material_luminosity\":"<<f.material_luminosity<<",\"carried_luminosity\":"<<f.carried_luminosity
        <<",\"conductive_luminosity\":"<<f.conductive_luminosity<<",\"entropy_production\":"<<f.entropy_production<<"}\n";
    } catch(const std::exception& e){std::cout<<"{\"ok\":false,\"error\":"<<std::quoted(e.what())<<"}\n";}
  }
  return std::cin.eof()?0:2;
}
