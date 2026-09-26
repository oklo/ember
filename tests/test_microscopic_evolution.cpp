#include "microscopic_test_transport.hpp"
#include "ember/stellar_seed.hpp"
#include "ember/eos_composite.hpp"
#include "ember/evolution.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace ember;
// Analytic opacity isolates evolution integration from external table coverage.
// Together with the fully ionized EOS and grey atmosphere this is a numerical
// control, not a predicted 0.5-solar-mass stellar model.
struct AnalyticOpacity final:Opacity {
  OpacityState eval(double T,double rho,const Composition& c)const override {
    const double absorption=1e23*rho*std::pow(T,-3.5),scattering=.2*(1+c.X[0]);
    const double base=absorption+scattering;
    // A smooth opacity enhancement creates a convective region, so the full
    // evolution test exercises constrained mixing as well as radiative zones.
    const double z=std::log(T/3e5)/.5,bump=20*std::exp(-z*z);
    return {base*(1+bump),-3.5*absorption/base-4*z*bump/(1+bump),absorption/base};
  }
  const char* name()const override{return "analytic opacity integration control";}
};
int main(int argc,char** argv) {
  try {
    const auto start=std::chrono::steady_clock::now();
    CompositeEos eos;AnalyticOpacity opacity;
    PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
    GreyAtmosphere atmosphere(eos,opacity);
    Physics physics{&eos,&opacity,&nuclear,1.9};
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
    auto seed=example::stellar_seed(256,.5*constants::Msun,.6*constants::Rsun,c,nuclear,atmosphere);
    auto equilibrium=relax(seed,physics,atmosphere);
    if(!equilibrium.converged) {
      std::cerr<<std::setprecision(17)<<"initial residual "<<equilibrium.residual<<", correction "<<equilibrium.correction
        <<", radius "<<equilibrium.model.r(255)<<", Tc "<<equilibrium.model.T(0)<<", L "<<equilibrium.model.y.back().L<<'\n';
      for(std::size_t i=0;i<equilibrium.history.size();++i) {
        const auto& h=equilibrium.history[i];std::cerr<<i<<' '<<h.residual<<' '<<h.correction<<' '<<h.damping<<'\n';
      }
      throw std::runtime_error("starting equilibrium: "+equilibrium.message);
    }
    seed=equilibrium.model;
    for(std::size_t i=0;i<seed.size();++i) {
      const double change=.003*(1-seed.m[i]/seed.M);
      seed.comp[i].X[0]-=change;seed.comp[i].X[2]+=change;
    }
    ThermalSubspaceControl micro;micro.value.diffusion=100;micro.value.thermal_force=0;micro.value.heat_conductivity=1e7;
    micro.total_species_heat=argc==2 && std::string(argv[1])=="total";
    physics.microscopic=&micro;physics.criterion=ConvectiveCriterion::ledoux;
    equilibrium=relax(seed,physics,atmosphere);
    if(!equilibrium.converged)throw std::runtime_error("microscopic equilibrium: "+equilibrium.message);
    const auto initial=equilibrium.model;
    const auto regions=convective_mixing_regions(initial,physics);
    if(regions.size()<2 || regions.size()==initial.size())
      throw std::runtime_error("test needs mixed and unmixed regions; found "+std::to_string(regions.size()));
    EvolutionOptions options;options.max_abundance_change=.01;
    const double dt=1e6*365.25*86400;
    const auto advanced=evolve_step(initial,physics,atmosphere,dt,options);
    if(!advanced.converged)throw std::runtime_error("coupled microscopic evolution: "+advanced.message);
    if(micro.total_species_heat && (micro.total_heat_calls==0 || advanced.total_species_rates.size()!=initial.size()-1
        || advanced.material_heat_residual>options.material_heat_tolerance))
      throw std::runtime_error("total-species heat was not reconstructed and converged");
    if(std::abs(advanced.luminosity_balance)>2e-8 || std::abs(advanced.nuclear_mass_balance)>2e-6
        || advanced.abundance_residual>options.abundance_tolerance)
      throw std::runtime_error("coupled energy, nuclear mass or abundance balance failed");
    auto no_species=micro;no_species.value.diffusion=0;auto control_physics=physics;control_physics.microscopic=&no_species;
    const auto control=evolve_step(initial,control_physics,atmosphere,dt,options);
    if(!control.converged)throw std::runtime_error("zero-species control: "+control.message);
    double difference=0,spread=0;
    for(std::size_t i=0;i<initial.size();++i) {
      const auto& comp=advanced.model.comp[i];
      if(*std::min_element(comp.X.begin(),comp.X.end())<0 || std::abs(comp.sum()-1)>1e-12)
        throw std::runtime_error("negative or non-normalized evolved composition");
      for(std::size_t j=0;j<NSPEC;++j)difference=std::max(difference,std::abs(comp.X[j]-control.model.comp[i].X[j]));
    }
    for(auto [begin,end]:convective_mixing_regions(advanced.model,physics))
      for(std::size_t i=begin;i<end;++i)for(std::size_t j=0;j<NSPEC;++j)
        spread=std::max(spread,std::abs(advanced.model.comp[i].X[j]-advanced.model.comp[begin].X[j]));
    if(!(difference>1e-10) || spread>1e-13 || advanced.model.age!=initial.age+dt
        || initial.age!=0 || initial.comp[0].X[0]!=equilibrium.model.comp[0].X[0])
      throw std::runtime_error("microscopic flux, convective homogeneity or age integration failed");
    if(micro.absent_heat_calls==0 || micro.species_calls==0 || initial.comp[0].X[1]!=0
        || advanced.model.comp[0].X[1]<=0)
      throw std::runtime_error("initially absent He3 did not exercise the thermal/species split");
    options.max_abundance_change=1e-12;
    const auto rejected=evolve_step(initial,physics,atmosphere,dt,options);
    if(rejected.converged || rejected.model.age!=initial.age || rejected.model.comp[0].X!=initial.comp[0].X
        || rejected.model.y[0].lnT!=initial.y[0].lnT)
      throw std::runtime_error("rejected microscopic step changed the original model");
    std::cout<<std::setprecision(17)<<"{\"outcome\":\"passed\",\"scope\":\"synthetic transport closure; grey 0.5 solar mass integration test\","
      <<"\"initial_regions\":"<<regions.size()<<",\"coupling_iterations\":"<<advanced.coupling_iterations
      <<",\"luminosity_balance\":"<<advanced.luminosity_balance<<",\"nuclear_mass_balance\":"<<advanced.nuclear_mass_balance
      <<",\"abundance_residual\":"<<advanced.abundance_residual<<",\"maximum_species_control_difference\":"<<difference
      <<",\"total_species_heat\":"<<micro.total_species_heat<<",\"material_heat_residual\":"<<advanced.material_heat_residual
      <<",\"convective_abundance_spread\":"<<spread<<",\"elapsed_seconds\":"
      <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
