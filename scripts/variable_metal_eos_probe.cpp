#include "ember/eos_variable_metal.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {
void number(double v) {if(std::isfinite(v))std::cout<<v;else std::cout<<"null";}
template<class T,std::size_t N> void array(const std::array<T,N>& a) {
  std::cout<<'[';
  for(std::size_t i=0;i<N;++i) {
    if(i)std::cout<<',';
    if constexpr(std::is_same_v<T,double>)number(a[i]);else array(a[i]);
  }
  std::cout<<']';
}
}
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    ember::VariableMetalHelmholtzEos eos(argv[1],ember::HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::string line;std::cout<<std::setprecision(17);
    while(std::getline(std::cin,line)) {
      std::istringstream in(line);double x,y,z,T,rho;int a,b,d;
      if(!(in>>x>>y>>z>>T>>rho>>a>>b>>d))return 2;
      try {
        auto c=ember::solar_scaled(x,z);c.basis=ember::AbundanceBasis::baryon_mass;
        c.metal_inventory=ember::MetalInventory::gs98;c.X[1]=y;c.X[2]=1-(x+y+z);
        const std::array<bool,3> active{a!=0,b!=0,d!=0};
        const auto begin=std::chrono::steady_clock::now();
        const auto h=eos.composition_heat(T,rho,c,active);
        const auto p=eos.composition_potential(T,rho,c,active);
        const auto e=eos.eval_with_derivatives(T,rho,c);const auto& s=e.state;
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        std::cout<<"{\"delta\":"<<h.material_delta<<",\"enthalpy\":";array(h.exchange_enthalpy);
        std::cout<<",\"radiation_enthalpy\":";array(h.radiation_enthalpy);
        std::cout<<",\"delta_partials\":";array(h.delta_partials);
        std::cout<<",\"enthalpy_partials\":";array(h.enthalpy_partials);
        std::cout<<",\"radiation_enthalpy_partials\":";array(h.radiation_enthalpy_partials);
        std::cout<<",\"phi\":"<<p.phi<<",\"gradient\":";array(p.gradient);
        std::cout<<",\"gradient_T\":";array(p.dgradient_dlnT);
        std::cout<<",\"gradient_rho\":";array(p.dgradient_dlnRho);
        std::cout<<",\"hessian\":";array(p.hessian);
        std::cout<<",\"state\":["<<s.P<<','<<s.E<<','<<s.S<<','<<s.cv<<','<<s.cp<<','<<s.chiT<<','<<s.chiRho<<','
          <<s.delta<<','<<s.grad_ad<<','<<s.Gamma1<<"],\"seconds\":"<<seconds<<"}\n";
      }catch(const std::exception& e) {std::cout<<"{\"error\":\""<<e.what()<<"\"}\n";}
      std::cout.flush();
    }
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
