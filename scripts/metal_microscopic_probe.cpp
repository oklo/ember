#include "ember/eos_variable_metal.hpp"
#include "ember/metal_microscopic_transport.hpp"
#include "conditional_metal_envelope_heat.hpp"
#include "ember/conduction.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace ember;
namespace {
template<class T>void array(const T& a) {
  std::cout<<'[';bool first=true;for(const auto& v:a){if(!first)std::cout<<',';first=false;std::cout<<v;}std::cout<<']';
}
template<class T>void matrix(const T& a) {
  std::cout<<'[';bool first=true;for(const auto& v:a){if(!first)std::cout<<',';first=false;array(v);}std::cout<<']';
}
Composition composition(double x,double y,double z) {
  auto c=solar_scaled(x,z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=y;c.X[2]-=y;return c;
}
struct ConstantConduction final:Conduction {
  OpacityState eval(double,double,const Composition&)const override{return {1e2,0,0};}
  const char* name()const override{return "analytic conductivity control";}
};
}
int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {
    VariableMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);ConstantConduction conduction;
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line))try {
      std::istringstream in(line);std::string mode;int ions,derivatives,mask,radiation;
      double ma,mb,xa,ya,za,xb,yb,zb;Point a,b;
      if(!(in>>mode>>ions>>derivatives>>mask>>radiation>>ma>>mb>>a.lnr>>a.lnrho>>a.lnT>>a.L>>xa>>ya>>za
          >>b.lnr>>b.lnrho>>b.lnT>>b.L>>xb>>yb>>zb))return 2;
      const auto ca=composition(xa,ya,za),cb=composition(xb,yb,zb);
      ScreenedMetalMicroscopicTransport hot(eos,collisions,ions,2e6,{bool(mask&1),bool(mask&2),bool(mask&4)},radiation);
      ConditionalMetalEnvelopeHeat envelope(eos,collisions,conduction,ions,2e6,3e6,radiation);
      const MetalMicroscopicTransport& provider=mode=="envelope"?static_cast<const MetalMicroscopicTransport&>(envelope):hot;
      const auto start=std::chrono::steady_clock::now();
      MicroscopicHeatResponse h;MetalMicroscopicFaceResponse f;
      if(mode=="total" || mode=="envelope") {
        MetalSpeciesVector total;if(!(in>>total[0]>>total[1]>>total[2]))return 2;
        h=microscopic_heat_with_total_metal_rate(provider,0,ma,mb,a,ca,b,cb,total,derivatives);
      }else if(mode=="heat")h=microscopic_heat(provider,0,ma,mb,a,ca,b,cb,derivatives);
      else if(mode=="face"){f=metal_microscopic_face(provider,0,ma,mb,a,ca,b,cb,derivatives);h=f;}
      else if(mode=="wrong"){h=microscopic_face(provider,0,ma,mb,a,ca,b,cb,derivatives);}
      else return 2;
      std::cout<<"{\"carried\":"<<h.carried_luminosity<<",\"conductivity\":"<<h.conductivity<<",\"rate\":";array(f.species.rate);
      std::cout<<",\"dleft\":";matrix(f.species.dleft);std::cout<<",\"dright\":";matrix(f.species.dright);
      std::cout<<",\"dcarried_lo\":";array(h.dcarried_lo);std::cout<<",\"dcarried_hi\":";array(h.dcarried_hi);
      std::cout<<",\"dconductivity_lo\":";array(h.dconductivity_lo);std::cout<<",\"dconductivity_hi\":";array(h.dconductivity_hi);
      std::cout<<",\"enthalpy\":";array(MetalSpeciesVector{h.total_rate_enthalpy[0],h.total_rate_enthalpy[1],h.total_metal_rate_enthalpy});
      std::cout<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
