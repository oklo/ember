#include "ember/nuclear_cn.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace ember;
template<class T> void array(const T& a) {
  std::cout<<'[';bool first=true;for(auto x:a){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';
}
int main() {
  std::cout<<std::setprecision(17);std::string line;
  while(std::getline(std::cin,line))try {
    std::istringstream in(line);std::string mode;in>>mode;
    if(mode=="implicit") {
      double dt;std::array<double,3> rate;CNAbundances old;
      if(!(in>>dt>>rate[0]>>rate[1]>>rate[2]>>old[0]>>old[1]>>old[2]))throw std::invalid_argument("bad implicit query");
      std::cout<<"{\"molality\":";array(cn_backward_euler(old,rate,dt));std::cout<<"}\n";
    }else if(mode=="source") {
      double T,rho,X,X3;CNAbundances y;int rates,screen;
      if(!(in>>T>>rho>>X>>X3>>y[0]>>y[1]>>y[2]>>rates>>screen))throw std::invalid_argument("bad source query");
      auto c=solar_scaled(X,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=X3;c.X[2]=1-X-X3-.02;
      const auto choice=rates==2?PPRates::solar_fusion_ii:PPRates::solar_fusion_iii;
      const auto screening=screen==0?PPScreening::debye_fermi:PPScreening::salpeter_van_horn;
      const auto result=CNNetwork(choice,screening).response(T,rho,c,y);const auto& n=result.physical.state;
      std::array<double,3> bare,slope,scr,zeta;
      for(std::size_t k=0;k<3;++k) {
        const auto r=cn_bare_rate(T,static_cast<CNReaction>(k),choice);
        const auto s=cn_screening(T,rho,c,static_cast<CNReaction>(k),screening);
        bare[k]=r.molar_rate;slope[k]=r.dlnrate_dlnT;scr[k]=s.log_factor;zeta[k]=s.zeta;
      }
      std::cout<<"{\"bare\":";array(bare);std::cout<<",\"slope\":";array(slope);
      std::cout<<",\"screen\":";array(scr);std::cout<<",\"zeta\":";array(zeta);
      std::cout<<",\"frequency\":";array(result.frequency);std::cout<<",\"reaction_rate\":";array(result.reaction_rate);
      std::cout<<",\"dXdt\":";array(n.dXdt);std::cout<<",\"eps\":"<<n.eps<<",\"eps_neutrino\":"<<n.eps_neutrino<<"}\n";
    }else throw std::invalid_argument("unknown query");
  }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
}
