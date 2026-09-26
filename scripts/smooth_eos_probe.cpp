#include "ember/eos_smooth_mixture.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {
    const bool smooth=std::string(argv[2])=="smooth";
    if(!smooth && std::string(argv[2])!="linear")return 2;
    const auto start=std::chrono::steady_clock::now();
    std::unique_ptr<ember::Eos> eos;
    const auto mixture=ember::HelmholtzTableEos::Mixture::allow_documented_proxy;
    if(smooth)eos=std::make_unique<ember::SmoothMetalHelmholtzEos>(argv[1],mixture);
    else eos=std::make_unique<ember::MetalHelmholtzEos>(argv[1],mixture);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);
    double x,y,T,rho;int chemical;
    while(std::cin>>x>>y>>T>>rho>>chemical) {
      try {
        auto c=ember::solar_scaled(x,.02);c.basis=ember::AbundanceBasis::baryon_mass;
        c.metal_inventory=ember::MetalInventory::gs98;c.X[1]=y;c.X[2]-=y;
        const auto a=eos->eval_with_derivatives(T,rho,c);const auto& e=a.state;
        const auto d=eos->composition_response(T,rho,c);
        std::vector<double> v{e.P,e.E,e.S,e.cv,e.cp,e.chiT,e.chiRho,e.delta,e.grad_ad,e.Gamma1,
          a.dE_dlnRho,a.dcp_dlnT,a.dcp_dlnRho,a.ddelta_dlnT,a.ddelta_dlnRho,
          a.dgrad_ad_dlnT,a.dgrad_ad_dlnRho,d.dP[0],d.dP[1],d.dE[0],d.dE[1]};
        if(chemical) {
          if(!smooth)throw std::invalid_argument("linear EOS has no chemical derivative API");
          const auto g=static_cast<const ember::SmoothMetalHelmholtzEos&>(*eos).composition_potential(T,rho,c);
          v.insert(v.end(),{g.phi,g.gradient[0],g.gradient[1],g.hessian[0][0],g.hessian[0][1],
              g.hessian[1][0],g.hessian[1][1],g.dgradient_dlnT[0],g.dgradient_dlnT[1],
              g.dgradient_dlnRho[0],g.dgradient_dlnRho[1]});
        }
        const auto range=eos->density_range(T,c).value();v.push_back(range.min);v.push_back(range.max);
        std::cout<<"{\"ok\":true,\"values\":[";
        for(std::size_t i=0;i<v.size();++i)std::cout<<(i?",":"")<<v[i];
        std::cout<<"]}\n";
      } catch(const std::exception& error) {
        std::cout<<"{\"ok\":false,\"error\":"<<std::quoted(error.what())<<"}\n";
      }
    }
    return std::cin.eof()?0:2;
  } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
