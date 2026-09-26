// Fixed-state CN and pp comparison, in the specified abundance basis.
#include "ember/nuclear.hpp"
#include <iostream>
#include <iomanip>
int main() {
  using namespace ember;
  int basis;double carbon,T,rho,X,Y3,Y4;
  std::cout<<std::setprecision(17);
  while(std::cin>>basis>>carbon>>T>>rho>>X>>Y3>>Y4) {
    try {
      auto c=solar_scaled(X,.02);c.metal_inventory=MetalInventory::gs98;
      c.basis=basis?AbundanceBasis::baryon_mass:AbundanceBasis::atomic_mass;
      c.X[1]=Y3;c.X[2]=Y4;
      const auto cn=CNCycle(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,carbon).eval(T,rho,c);
      const auto pp=PPChains(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn).eval(T,rho,c);
      const auto rate=cn_bare_rate(T,PPRates::solar_fusion_iii);
      const auto screen=cn_screening(T,rho,c,PPScreening::salpeter_van_horn);
      std::cout<<"{\"eps_cn\":"<<cn.eps<<",\"eps_pp\":"<<pp.eps<<",\"neutrino_cn\":"<<cn.eps_neutrino
        <<",\"dXH_cn\":"<<cn.dXdt[0]<<",\"dXHe4_cn\":"<<cn.dXdt[2]
        <<",\"T_slope\":"<<cn.dlneps_dlnT<<",\"rho_slope\":"<<cn.dlneps_dlnRho
        <<",\"bare_rate\":"<<rate.molar_rate<<",\"bare_slope\":"<<rate.dlnrate_dlnT
        <<",\"screen\":"<<screen.log_factor<<",\"zeta\":"<<screen.zeta<<"}\n";
    }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
  }
  return std::cin.eof()?0:1;
}
