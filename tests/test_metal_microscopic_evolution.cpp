#include "ember/metal_microscopic_transport.hpp"
#include "ember/eos_composite.hpp"
#include "ember/stellar_seed.hpp"
#include "ember/detail/differential.hpp"
#include "ember/evolution_checkpoint.hpp"
#include "ember/convective_evolution_checks.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>
using namespace ember;
namespace {
struct AnalyticOpacity final:Opacity {
  OpacityState eval(double T,double rho,const Composition& c)const override {
    const double absorption=1e23*rho*std::pow(T,-3.5),base=absorption+.2*(1+c.X[0]);
    const double z=std::log(T/3e5)/.5,bump=20*std::exp(-z*z);
    return {base*(1+bump),-3.5*absorption/base-4*z*bump/(1+bump),absorption/base};
  }
  const char* name()const override{return "analytic opacity integration control";}
};
// Synthetic smooth transport isolates conservative abundance and heat
// integration from source-table coverage and physical collision assumptions.
struct MetalControl final:MetalMicroscopicTransport {
  double diffusion{100},metal_heat{3e14};
  mutable std::size_t total_calls{};
  MetalMicroscopicFaceResponse response(double ma,double mb,const Point& a,const Composition& ca,
      const Point& b,const Composition& cb,const MetalSpeciesVector* total,bool derivatives)const {
    using D=detail::Differential<8>;using detail::exp;
    const D r=.5*(exp(D::variable(a.lnr,0))+exp(D::variable(b.lnr,4)));
    const D rho=.5*(exp(D::variable(a.lnrho,1))+exp(D::variable(b.lnrho,5)));
    const D T=.5*(exp(D::variable(a.lnT,2))+exp(D::variable(b.lnT,6)));
    const D area=4*M_PI*r*r,g=diffusion*area*area*rho*rho/(mb-ma);
    const MetalSpeciesVector left{ca.X[0],ca.X[1],ca.Z()},right{cb.X[0],cb.X[1],cb.Z()};
    const MetalSpeciesVector hbase{3e14,-1.2e14,metal_heat};
    D carried;MetalMicroscopicFaceResponse out;out.conductivity=1e7;
    for(std::size_t k=0;k<3;++k) {
      const D rate=g*(left[k]-right[k]),h=hbase[k]*T/1e7;
      out.species.rate[k]=rate.value;
      if(derivatives){out.species.dleft[k][k]=g.value;out.species.dright[k][k]=-g.value;}
      carried+=h*(total?D((*total)[k]):rate);
      if(total){if(k<2)out.total_rate_enthalpy[k]=h.value;else out.total_metal_rate_enthalpy=h.value;}
    }
    out.carried_luminosity=carried.value;
    if(derivatives)for(std::size_t j=0;j<NVAR;++j){out.dcarried_lo[j]=carried.d[j];out.dcarried_hi[j]=carried.d[4+j];}
    return out;
  }
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double ma,double mb,const Point& a,
      const Composition& ca,const Point& b,const Composition& cb,bool derivatives)const override {
    return response(ma,mb,a,ca,b,cb,nullptr,derivatives);
  }
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double ma,double mb,const Point& a,
      const Composition& ca,const Point& b,const Composition& cb,const MetalSpeciesVector& total,bool derivatives)const override {
    std::atomic_ref(total_calls).fetch_add(1,std::memory_order_relaxed);return response(ma,mb,a,ca,b,cb,&total,derivatives);
  }
  const char* name()const override{return "analytic three-mass integration control";}
};
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
}
int main(int argc,char** argv) {
  try {
    CompositeEos eos;AnalyticOpacity opacity;PPCNNetwork nuclear;GreyAtmosphere atmosphere(eos,opacity);
    Physics physics{&eos,&opacity,&nuclear,1.9};
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.001;c.X[2]-=.001;c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);
    const auto grid=argc>1?LuminosityGrid::volume_faces:LuminosityGrid::mass_nodes;
    auto seed=example::stellar_seed(256,.5*constants::Msun,.6*constants::Rsun,c,nuclear,atmosphere,3.,0.,grid);
    auto equilibrium=relax(seed,physics,atmosphere);
    if(!equilibrium.converged)throw std::runtime_error("initial equilibrium: "+equilibrium.message);
    seed=equilibrium.model;
    for(std::size_t i=0;i<seed.size();++i) {
      auto v=metal_cn_abundances(seed.comp[i]);const double fraction=1-seed.m[i]/seed.M;
      v[0]-=.003*fraction;v[5]+=.003*fraction;
      seed.comp[i]=metal_cn_composition(seed.comp[i],v);
    }
    MetalControl transport;physics.microscopic=&transport;physics.criterion=ConvectiveCriterion::ledoux;
    equilibrium=relax(seed,physics,atmosphere);
    if(!equilibrium.converged)throw std::runtime_error("transport equilibrium: "+equilibrium.message);
    const auto initial=equilibrium.model;const auto regions=convective_mixing_regions(initial,physics);
    require(regions.size()>1 && regions.size()<initial.size(),"mixed and radiative regions required");
    EvolutionOptions options;options.max_abundance_change=.01;options.abundance_tolerance=1e-13;
    const double dt=1e6*365.25*86400;
    const bool finite=argc>2;
    std::vector<MetalSpeciesVector> initial_heat_rates;
    double full_half_error=0;
    if(finite) {
      options.convective_mixing=std::string(argv[2])=="lagged"?
          ConvectiveMixing::finite_lagged:ConvectiveMixing::finite_implicit;
      // This initial static structure was solved with microscopic heat only.
      // Record exactly that face heat, rather than inventing a previous step.
      for(std::size_t i=0;i+1<initial.size();++i)
        initial_heat_rates.push_back(transport.metal_eval(i,initial.m[i],initial.m[i+1],initial.y[i],
            initial.comp[i],initial.y[i+1],initial.comp[i+1],false).species.rate);
      options.previous_metal_heat_rates=initial_heat_rates;
      auto missing=options;missing.convective_mixing=ConvectiveMixing::finite_lagged;
      missing.previous_metal_heat_rates={};
      const auto rejected=evolve_step(initial,physics,atmosphere,dt,missing);
      require(!rejected.converged && rejected.message.find("previous total species heat rates")!=std::string::npos,
          "stratified lagged convection silently omitted previous material heat");
    }
    const auto step=evolve_step(initial,physics,atmosphere,dt,options);
    if(!step.converged)throw std::runtime_error("physical metal evolution: "+step.message);
    require(transport.total_calls>0 && step.total_metal_species_rates.size()==initial.size()-1 && step.total_species_rates.empty(),"three-mass total rates not used");
    require(std::abs(step.luminosity_balance)<2e-8 && std::abs(step.nuclear_mass_balance)<2e-6,"discrete stellar energy or nuclear mass balance");
    require(step.abundance_residual<=options.abundance_tolerance && step.material_heat_residual<=options.material_heat_tolerance,"abundance or heat convergence");
    require(driver::check_interval(initial,step,dt,nuclear,1e-14).pass,
        "returned composition must conserve integrated nuclear sources and energy");
    auto fast=options;
    fast.relaxation.zone_threads=4;
    fast.coupling_stop_tolerance=1e-10;
    fast.material_heat_tolerance=1e-7;
    fast.verification_residual_tolerance=1e-8;
    fast.verification_correction_tolerance=1e-7;
    const auto parallel=evolve_step(initial,physics,atmosphere,dt,fast);
    if(!parallel.converged)throw std::runtime_error("parallel coupling: "+parallel.message);
    require(driver::check_interval(initial,parallel,dt,nuclear,1e-14).pass,
        "outer stopping rule changed integrated conservation");
    require(convective_mixing_regions(parallel.model,physics)==convective_mixing_regions(step.model,physics),
        "parallel coupling changed convection boundaries");
    auto response_options=fast;response_options.linearized_burning=true;
    const auto responsive=evolve_step(initial,physics,atmosphere,dt,response_options);
    if(!responsive.converged)throw std::runtime_error("linearized burning: "+responsive.message);
    require(driver::check_interval(initial,responsive,dt,nuclear,1e-14).pass,
        "linearized burning changed isotope or energy conservation");
    require(convective_mixing_regions(responsive.model,physics)==convective_mixing_regions(step.model,physics),
        "linearized burning changed convection boundaries");
    for(std::size_t i=0;i<initial.size();++i) {
      for(auto v:{Var::lnr,Var::lnrho,Var::lnT})
        require(std::abs(responsive.model.y[i][v]-parallel.model.y[i][v])<1e-6,
            "linearized burning changed the thermal solution");
      for(std::size_t k=0;k<NSPEC;++k)
        require(std::abs(responsive.model.comp[i].X[k]-parallel.model.comp[i].X[k])<1e-8,
            "linearized burning changed the species solution");
    }
    for(std::size_t i=0;i<initial.size();++i) {
      for(auto v:{Var::lnr,Var::lnrho,Var::lnT})
        require(std::abs(parallel.model.y[i][v]-step.model.y[i][v])<1e-6,
            "outer stopping rule changed the thermal solution");
      for(std::size_t k=0;k<NSPEC;++k)
        require(std::abs(parallel.model.comp[i].X[k]-step.model.comp[i].X[k])<1e-8,
            "outer stopping rule changed the species solution");
    }
    const auto serial_faces=convective_mixing_faces(initial,physics);
    const auto parallel_faces=convective_mixing_faces(initial,physics,4);
    for(std::size_t i=0;i<serial_faces.size();++i) {
      const auto& a=serial_faces[i];const auto& b=parallel_faces[i];
      require(a.diffusivity==b.diffusivity && a.velocity==b.velocity && a.length==b.length
          && a.buoyancy_contrast==b.buoyancy_contrast && a.composition_term==b.composition_term,
          "parallel face test changed a transport coefficient");
    }
    for(int k=0;k<3;++k) {
      auto bad=options;
      if(k==0)bad.coupling_stop_tolerance=-1;
      if(k==1)bad.verification_residual_tolerance=std::numeric_limits<double>::quiet_NaN();
      if(k==2)bad.verification_correction_tolerance=-1;
      bool refused=false;
      try{(void)evolve_step(initial,physics,atmosphere,dt,bad);}
      catch(const std::invalid_argument&){refused=true;}
      require(refused,"invalid coupling setting accepted");
    }
    const auto weights=nodal_mass_weights(initial);double inert=0,maximum_rate=0;
    for(std::size_t i=0;i<initial.size();++i) {
      const auto before=metal_cn_abundances(initial.comp[i]),after=metal_cn_abundances(step.model.comp[i]);
      inert+=weights[i]/initial.M*(after[5]-before[5]);
      require(step.model.comp[i].X[2]>=0 && std::abs(step.model.comp[i].sum()-1)<1e-12,"physical composition positivity or normalization");
      for(double v:after)require(v>=0,"negative independent physical abundance");
    }
    for(const auto& rate:step.total_metal_species_rates)maximum_rate=std::max(maximum_rate,std::abs(rate[2]));
    require(std::abs(inert)<1e-13 && maximum_rate>0,"inert metal conservation and nonzero redistribution");
    if(!finite)for(auto [begin,end]:convective_mixing_regions(step.model,physics))
      for(std::size_t i=begin;i<end;++i)require(step.model.comp[i]==step.model.comp[begin],"mixed region must be homogeneous");
    if(finite) {
      require(step.mixed_regions>0 && step.convective_mass_fraction>0,"finite convection lost physical region diagnostics");
      const auto half=evolve_step(initial,physics,atmosphere,dt/2,options);
      if(!half.converged)throw std::runtime_error("finite half step: "+half.message);
      auto follow=options;follow.previous_metal_heat_rates=half.total_metal_species_rates;
      const auto twice=evolve_step(half.model,physics,atmosphere,dt/2,follow);
      if(!twice.converged)throw std::runtime_error("finite second half step: "+twice.message);
      require(driver::check_interval(initial,half,dt/2,nuclear,1e-14).pass
          && driver::check_interval(half.model,twice,dt/2,nuclear,1e-14).pass,
          "each half interval must conserve its own nuclear sources and energy");
      const auto stamp=std::chrono::high_resolution_clock::now().time_since_epoch().count();
      const auto file=std::filesystem::temp_directory_path()/("ember-finite-"+std::to_string(stamp)+".checkpoint");
      const driver::Selections selection{"synthetic","CN","physical metals",argv[2],"total material heat"};
      const driver::Identities identities{{"executable","finite-test-only"}};
      driver::write_checkpoint(file,{half.model,dt/2,1,0,half.total_metal_species_rates},selection,1e-12,identities);
      bool refused=false;
      try{(void)driver::read_checkpoint(file,half.model.size(),half.model.M,half.model.comp[0],selection,1e-12,identities,grid);}
      catch(const std::runtime_error&){refused=true;}
      require(refused,"reader silently discarded saved material heat");
      const auto saved=driver::read_checkpoint(file,half.model.size(),half.model.M,half.model.comp[0],selection,1e-12,identities,grid,true);
      require(saved.metal_heat_rates==half.total_metal_species_rates,"checkpoint changed material heat rates");
      follow.previous_metal_heat_rates=saved.metal_heat_rates;
      const auto replay=evolve_step(saved.model,physics,atmosphere,dt/2,follow);
      require(replay.converged && replay.model.comp==twice.model.comp
          && replay.total_metal_species_rates==twice.total_metal_species_rates,"checkpoint changed finite-mixing continuation");
      for(std::size_t i=0;i<initial.size();++i)for(std::size_t k=0;k<NVAR;++k)
        require(replay.model.y[i][static_cast<Var>(k)]==twice.model.y[i][static_cast<Var>(k)],"restart changed thermal state");
      std::filesystem::remove(file);
      for(std::size_t i=0;i<initial.size();++i)
        full_half_error=std::max(full_half_error,std::abs(step.model.y[i].lnT-twice.model.y[i].lnT));
      require(full_half_error<2e-4,"finite full/two-half thermal difference");
    }
    auto noheat=transport;noheat.metal_heat=0;auto p=physics;p.microscopic=&noheat;
    const auto control=evolve_step(initial,p,atmosphere,dt,options);
    if(!control.converged)throw std::runtime_error("metal enthalpy control: "+control.message);
    double thermal_difference=0;
    for(std::size_t i=0;i<initial.size();++i)thermal_difference=std::max(thermal_difference,std::abs(control.model.y[i].lnT-step.model.y[i].lnT));
    require(thermal_difference>1e-12,"metal enthalpy did not affect thermal solution");
    options.max_abundance_change=1e-12;const auto rejected=evolve_step(initial,physics,atmosphere,dt,options);
    require(!rejected.converged && rejected.model.comp==initial.comp && rejected.model.age==initial.age,"failed step changed input");
    auto wrong=physics;wrong.microscopic=nullptr;const auto incompatible=evolve_step(initial,wrong,atmosphere,dt);
    require(!incompatible.converged && incompatible.model.comp==initial.comp,"missing physical metal transport not rejected");
    require(step.model.age==initial.age+dt,"physical metal step age");
    require(step.model.luminosity_grid==grid,"luminosity grid changed");
    std::cout<<std::setprecision(17)<<"{\"outcome\":\"passed\",\"scope\":\"synthetic transport and grey 0.5 solar mass integration\",\"regions\":"<<regions.size()
      <<",\"face_luminosities\":"<<face_luminosities(initial)
      <<",\"finite_convection\":"<<finite<<",\"lagged\":"<<(options.convective_mixing==ConvectiveMixing::finite_lagged)
      <<",\"convection_composition_ratio\":"<<step.convection_composition_ratio
      <<",\"full_two_half_logT_difference\":"<<full_half_error
      <<",\"coupling_iterations\":"<<step.coupling_iterations<<",\"luminosity_balance\":"<<step.luminosity_balance
      <<",\"nuclear_mass_balance\":"<<step.nuclear_mass_balance<<",\"inert_balance\":"<<inert
      <<",\"abundance_residual\":"<<step.abundance_residual<<",\"material_heat_residual\":"<<step.material_heat_residual
      <<",\"maximum_metal_rate\":"<<maximum_rate<<",\"metal_heat_thermal_difference\":"<<thermal_difference<<"}\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
