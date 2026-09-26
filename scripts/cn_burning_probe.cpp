#include "ember/cn_burning.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace ember;
int main() {
  std::cout<<std::setprecision(17);std::string line;
  PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);CNNetwork cn;
  while(std::getline(std::cin,line))try {
    std::istringstream in(line);std::size_t count,nr;double dt,tolerance;
    if(!(in>>count>>nr>>dt>>tolerance) || count<2 || count>2048)throw std::invalid_argument("invalid query header");
    Model m;std::vector<CNAbundances> catalysts;MixingRegions regions;
    for(std::size_t i=0;i<count;++i) {
      double mass,r,rho,T,L,X,X3;CNAbundances y;
      if(!(in>>mass>>r>>rho>>T>>L>>X>>X3>>y[0]>>y[1]>>y[2]))throw std::invalid_argument("invalid model row");
      auto c=solar_scaled(X,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=X3;c.X[2]=1-X-X3-.02;m.comp.push_back(c);m.m.push_back(mass);
      m.y.push_back({std::log(r),std::log(rho),std::log(T),L});catalysts.push_back(y);
    }
    m.M=m.m.back();
    for(std::size_t i=0;i<nr;++i) {std::size_t begin,end;if(!(in>>begin>>end))throw std::invalid_argument("invalid region");regions.emplace_back(begin,end);}
    const auto start=std::chrono::steady_clock::now();const auto result=burn_cn_and_mix(m,m,catalysts,pp,cn,regions,dt,tolerance);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"seconds\":"<<seconds<<",\"iterations\":"<<result.maximum_iterations
      <<",\"residual\":"<<result.maximum_equation_residual<<",\"nuclear_mass_balance\":"<<result.nuclear_mass_balance
      <<",\"nuclear_luminosity\":"<<result.nuclear_luminosity<<",\"neutrino_luminosity\":"<<result.neutrino_luminosity
      <<",\"composition\":[";
    for(std::size_t i=0;i<count;++i) {
      if(i)std::cout<<',';const auto& c=result.lookup[i];const auto& y=result.catalysts[i];
      std::cout<<'['<<c.X[0]<<','<<c.X[1]<<','<<c.X[2]<<','<<y[0]<<','<<y[1]<<','<<y[2]<<']';
    }
    std::cout<<"],\"sources\":[";
    for(std::size_t i=0;i<count;++i) {
      if(i)std::cout<<',';const auto& c=result.lookup[i];const auto& y=result.catalysts[i];
      const auto p=pp.eval(m.T(i),m.rho(i),c);const auto n=cn.response(m.T(i),m.rho(i),c,y).physical.state;
      std::cout<<'['<<p.dXdt[0]+n.dXdt[0]<<','<<p.dXdt[1]<<','<<p.dXdt[2]+n.dXdt[2]<<','<<n.dXdt[3]<<','<<n.dXdt[4]<<','<<n.dXdt[5]<<']';
    }
    std::cout<<"]}\n";
  }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
}
