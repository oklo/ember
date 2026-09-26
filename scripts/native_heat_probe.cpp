#include "ember/eos_smooth_mixture.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {
void number(double v) {if(std::isfinite(v))std::cout<<v;else std::cout<<"null";}
template<std::size_t N> void array(const std::array<double,N>& a) {
  std::cout<<'[';for(std::size_t i=0;i<N;++i) {if(i)std::cout<<',';number(a[i]);}std::cout<<']';
}
}
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    ember::SmoothMetalHelmholtzEos eos(argv[1],ember::HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::string line;std::cout<<std::setprecision(17);
    while(std::getline(std::cin,line)) {
      std::istringstream in(line);double x,y,T,rho;int a,b;
      if(!(in>>x>>y>>T>>rho>>a>>b))return 2;
      try {
        auto c=ember::solar_scaled(x,.02);c.basis=ember::AbundanceBasis::baryon_mass;
        c.metal_inventory=ember::MetalInventory::gs98;c.X[1]=y;c.X[2]-=y;
        const std::array<bool,2> active{a!=0,b!=0};
        const auto begin=std::chrono::steady_clock::now();
        const auto h=eos.composition_heat(T,rho,c,active);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        const auto p=eos.active_composition_potential(T,rho,c,active);
        const auto e=eos.eval_with_derivatives(T,rho,c);const auto& s=e.state;
        std::cout<<"{\"delta\":"<<h.material_delta<<",\"enthalpy\":";array(h.exchange_enthalpy);
        std::cout<<",\"delta_partials\":";array(h.delta_partials);
        std::cout<<",\"enthalpy_partials\":[";array(h.enthalpy_partials[0]);std::cout<<',';array(h.enthalpy_partials[1]);std::cout<<']';
        std::cout<<",\"phi\":"<<p.phi<<",\"gradient\":";array(p.gradient);
        std::cout<<",\"gradient_T\":";array(p.dgradient_dlnT);
        std::cout<<",\"gradient_rho\":";array(p.dgradient_dlnRho);
        std::cout<<",\"hessian\":[";array(p.hessian[0]);std::cout<<',';array(p.hessian[1]);std::cout<<']';
        std::cout<<",\"state\":["<<s.P<<','<<s.E<<','<<s.S<<','<<s.cv<<','<<s.cp<<','<<s.chiT<<','<<s.chiRho<<','
          <<s.delta<<','<<s.grad_ad<<','<<s.Gamma1<<"],\"seconds\":"<<seconds<<"}\n";
      }catch(const std::exception& e) {std::cout<<"{\"error\":\""<<e.what()<<"\"}\n";}
      std::cout.flush();
    }
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
