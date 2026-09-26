// Read-only diagnosis of composition buoyancy on saved variable-metal models.
// Uses the actual temperature gradient and common-pressure EOS comparisons.
// No diffusivity or physical mixing is selected by this probe.
#include "ember/eos_variable_metal.hpp"
#include "ember/convection.hpp"
#include "ember/nuclear_cn.hpp"
#include "ember/model.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ember;
int main(int argc, char** argv) {
  if(argc!=2)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    VariableMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::string label;std::size_t count;
    while(std::cin>>label>>count) {
      if(count<2 || count>20000)throw std::invalid_argument("invalid profile size");
      Model model;std::vector<EosState> states;
      for(std::size_t i=0;i<count;++i) {
        double mass,radius,rho,T,L,H,Y,Z;CNAbundances cn;
        if(!(std::cin>>mass>>radius>>rho>>T>>L>>H>>Y>>cn[0]>>cn[1]>>cn[2]>>Z))
          throw std::invalid_argument("invalid physical profile row");
        if(!(mass>0 && radius>0 && rho>0 && T>0) || (i && mass<=model.m.back()))
          throw std::invalid_argument("invalid physical geometry");
        auto c=solar_scaled(H,Z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
        c.X[1]=Y;c.X[2]-=Y;c.cn_molality=cn;c.cn_mass_convention=CNMassConvention::explicit_metal_mass;
        (void)cn_physical_ledger(c,cn);
        model.m.push_back(mass);model.y.push_back({std::log(radius),std::log(rho),std::log(T),L});
        model.comp.push_back(c);states.push_back(eos.eval(T,rho,c));
      }
      model.M=model.m.back();
      std::ostringstream out;out<<std::setprecision(17)<<"{\"label\":"<<std::quoted(label)<<",\"pressure\":[";
      for(std::size_t i=0;i<count;++i)out<<(i?",":"")<<states[i].P;
      out<<"],\"faces\":[";
      for(std::size_t i=0;i+1<count;++i) {
        const double T=.5*(model.T(i)+model.T(i+1)),rho=.5*(model.rho(i)+model.rho(i+1));
        const double P=.5*(states[i].P+states[i+1].P),delta=.5*(states[i].delta+states[i+1].delta);
        const double dlnP=std::log(states[i+1].P/states[i].P);
        if(!(dlnP<0))throw std::domain_error("pressure must decrease outward");
        const auto b=composition_buoyancy(eos,T,P,delta,dlnP,model.comp[i],model.comp[i+1],rho);
        const double grad=(model.y[i+1].lnT-model.y[i].lnT)/dlnP;
        const double ad=.5*(states[i].grad_ad+states[i+1].grad_ad);
        if(!std::isfinite(b.B+grad+ad))throw std::domain_error("nonfinite buoyancy");
        if(i)out<<',';
        out<<"{\"face\":"<<i<<",\"q\":"<<.5*(model.m[i]+model.m[i+1])/model.M
           <<",\"T\":"<<T<<",\"rho\":"<<rho<<",\"P\":"<<P<<",\"delta\":"<<delta
           <<",\"cp\":"<<.5*(states[i].cp+states[i+1].cp)<<",\"grad\":"<<grad
           <<",\"grad_ad\":"<<ad<<",\"B\":"<<b.B<<",\"R0\":";
        if(b.B<0 && grad<ad)out<<(grad-ad)/b.B;else out<<"null";
        out<<",\"XH\":"<<.5*(model.comp[i].X[0]+model.comp[i+1].X[0])
           <<",\"X3\":"<<.5*(model.comp[i].X[1]+model.comp[i+1].X[1])
           <<",\"Z\":"<<.5*(model.comp[i].Z()+model.comp[i+1].Z())<<'}';
      }
      out<<"]}\n";std::cout<<out.str()<<std::flush;
    }
    if(!std::cin.eof())throw std::invalid_argument("invalid profile header");
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
