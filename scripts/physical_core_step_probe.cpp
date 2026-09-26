// Coupled hot-interior control with fixed outer pressure and temperature.
// This boundary represents a held envelope. It is not a photosphere or a
// complete stellar trajectory; species flux is sealed at the truncated edge.
#include "ember/screened_microscopic_transport.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/evolution.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace ember;
namespace {
std::string scalar(double value) {
  if(!std::isfinite(value))return "null";
  std::ostringstream out;out<<std::setprecision(17)<<value;return out.str();
}
struct HeldEnvelope final:Atmosphere {
  double temperature,pressure;
  HeldEnvelope(double T,double P):temperature(T),pressure(P){}
  AtmosphereState eval(double,double,const Composition&)const override {
    AtmosphereState out;out.T=temperature;out.P=pressure;return out;
  }
  const char* name()const override{return "fixed pressure and temperature at truncated hot-interior edge";}
};
}
int main(int argc,char** argv) {
  if(argc!=4)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);
    MixtureOpacity atomic(std::filesystem::path(argv[3])/"tops_gs98_mixture_high.dat");
    ElementalOpacity opacity(atomic);
    const PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line)) {
      try {
        std::istringstream in(line);int ions,absent;double dt;std::size_t n;
        if(!(in>>ions>>absent>>dt>>n) || n<2 || n>512)return 2;
        Model old;
        for(std::size_t i=0;i<n;++i) {
          double mass,r,rho,T,L,x,y;if(!(in>>mass>>r>>rho>>T>>L>>x>>y))return 2;
          auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
          c.X[1]=absent?0:y;c.X[2]-=c.X[1];
          old.m.push_back(mass);old.y.push_back({std::log(r),std::log(rho),std::log(T),L});old.comp.push_back(c);
        }
        old.M=old.m.back();const auto unchanged=old;
        ScreenedMicroscopicTransport micro(eos,collisions,ions,2e6);
        Physics physics{&eos,&opacity,&nuclear,1.9};physics.microscopic=&micro;
        physics.criterion=ConvectiveCriterion::ledoux;
        const auto boundary=eos.eval(old.T(n-1),old.rho(n-1),old.comp.back());
        HeldEnvelope envelope(old.T(n-1),boundary.P);
        EvolutionOptions options;options.max_abundance_change=.01;options.max_coupling_iterations=12;
        const auto begin=std::chrono::steady_clock::now();
        const auto regions=convective_mixing_regions(old,physics);
        const auto result=evolve_step(old,physics,envelope,dt,options);
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        bool preserved=old.m==unchanged.m && old.age==unchanged.age;
        for(std::size_t i=0;i<n;++i)preserved &= old.comp[i].X==unchanged.comp[i].X
          && old.y[i].lnr==unchanged.y[i].lnr && old.y[i].lnrho==unchanged.y[i].lnrho
          && old.y[i].lnT==unchanged.y[i].lnT && old.y[i].L==unchanged.y[i].L;
        double max_temperature_change=0,max_density_change=0,max_composition_change=0;
        for(std::size_t i=0;i<n;++i) {
          max_temperature_change=std::max(max_temperature_change,std::abs(result.model.y[i].lnT-old.y[i].lnT));
          max_density_change=std::max(max_density_change,std::abs(result.model.y[i].lnrho-old.y[i].lnrho));
          for(std::size_t k=0;k<NSPEC;++k)max_composition_change=std::max(max_composition_change,std::abs(result.model.comp[i].X[k]-old.comp[i].X[k]));
        }
        std::cout<<"{\"converged\":"<<result.converged<<",\"message\":"<<std::quoted(result.message)
          <<",\"initial_regions\":"<<regions.size()<<",\"iterations\":"<<result.coupling_iterations
          <<",\"residual\":"<<result.residual<<",\"correction\":"<<scalar(result.correction)
          <<",\"abundance_residual\":"<<result.abundance_residual<<",\"luminosity_balance\":"<<result.luminosity_balance
          <<",\"nuclear_mass_balance\":"<<result.nuclear_mass_balance<<",\"age\":"<<result.model.age
          <<",\"input_preserved\":"<<preserved<<",\"seconds\":"<<seconds
          <<",\"max_dlnT\":"<<max_temperature_change<<",\"max_dlnrho\":"<<max_density_change
          <<",\"max_dX\":"<<max_composition_change<<",\"model\":[";
        for(std::size_t i=0;i<n;++i) {
          if(i)std::cout<<',';const auto& m=result.model;
          std::cout<<'['<<m.m[i]<<','<<m.r(i)<<','<<m.rho(i)<<','<<m.T(i)<<','<<m.y[i].L<<','<<m.comp[i].X[0]<<','<<m.comp[i].X[1]<<']';
        }
        std::cout<<"]}\n"<<std::flush;
      }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
