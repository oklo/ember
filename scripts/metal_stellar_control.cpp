// Finite stellar control of the physical metal-mass and heat integration.
// Opacity composition is fixed by default. Explicit selectors allow changing
// interior opacity and a measured atmospheric metal response. The caller
// records the selected approximations and supplies complete physical profiles.
#include "conditional_metal_envelope_heat.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_hydrogen.hpp"
#include "ember/opacity_hydrogen_continuation.hpp"
#include "ember/opacity_hydrogen_share.hpp"
#include "ember/opacity_metal_extension.hpp"
#include "ember/opacity_metal_lower_continuation.hpp"
#include "ember/conduction_table.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/atmosphere_trace_helium.hpp"
#include "ember/atmosphere_metal_interval.hpp"
#include "ember/atmosphere_metal_response.hpp"
#include "ember/atmosphere_metal_piecewise.hpp"
#include "ember/atmosphere_metal_chain.hpp"
#include "ember/atmosphere_hydrogen_interval.hpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <sstream>
using namespace ember;
namespace {
Composition fixed_metal(const Composition& c) {
  auto p=solar_scaled(c.X[0],.02);p.basis=AbundanceBasis::baryon_mass;p.metal_inventory=MetalInventory::gs98;
  p.X[1]=c.X[1];p.X[2]-=c.X[1];return p;
}
struct FixedMetalOpacity final:Opacity {
  const Opacity& source;explicit FixedMetalOpacity(const Opacity& s):source(s){}
  OpacityState eval(double T,double rho,const Composition& c)const override{return source.eval(T,rho,fixed_metal(c));}
  std::optional<DensityRange> density_range(double T,const Composition& c)const override{return source.density_range(T,fixed_metal(c));}
  const char* name()const override{return "control: opacity composition held at Z=.02";}
};
struct FixedMetalAtmosphere final:Atmosphere {
  const Atmosphere& source;const Eos& eos;
  FixedMetalAtmosphere(const Atmosphere& s,const Eos& e):source(s),eos(e){}
  AtmosphereState eval(double teff,double gravity,const Composition& c)const override {
    auto state=source.eval(teff,gravity,fixed_metal(c));
    state.rho=eos.rho_from_PT(state.T,state.P,c,state.rho);return state;
  }
  const char* name()const override{return "control: atmospheric opacity composition held at Z=.02";}
};
template<class T>void array(const T& a) {
  std::cout<<'[';bool first=true;for(const auto& v:a){if(!first)std::cout<<',';first=false;std::cout<<v;}std::cout<<']';
}
}
int main(int argc,char** argv) {
  if(argc<11 || argc>24 || argc==14 || argc==16)return 2;
  try {
    const bool hydrogen_share_opacity=argc>=12 && std::string(argv[11])=="variable-opacity-hydrogen-share";
    const bool variable_opacity=hydrogen_share_opacity || (argc>=12 && std::string(argv[11])=="variable-opacity");
    if(argc>=12 && !variable_opacity)throw std::invalid_argument("unknown interior opacity selection");
    VariableMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::unique_ptr<Opacity> molecular_source,hydrogen_source;
    if(hydrogen_share_opacity) {
      molecular_source=std::make_unique<HydrogenShareMixtureOpacity>(argv[2]);
      hydrogen_source=std::make_unique<HydrogenShareMixtureOpacity>(argv[8]);
    } else {
      molecular_source=std::make_unique<MixtureOpacity>(argv[2]);
      hydrogen_source=std::make_unique<MixtureOpacity>(argv[8]);
    }
    const Opacity& molecular=*molecular_source;
    const Opacity& dependence=*hydrogen_source;
    MixtureOpacity warm(argv[3]),bridge(argv[4]),hot(argv[5]);
    const double maximum_warm_hydrogen=argc>=21?std::stod(argv[20]):.90;
    const double allowed_hydrogen=hydrogen_share_opacity?.99999:.97;
    if(!std::isfinite(maximum_warm_hydrogen) || maximum_warm_hydrogen<.90 || maximum_warm_hydrogen>allowed_hydrogen)
      throw std::invalid_argument(hydrogen_share_opacity
          ?"warm hydrogen continuation must be limited to [.90,.99999]"
          :"warm hydrogen continuation must be limited to [.90,.97]");
    // The hydrogen-share selection keeps source-plane helium nonnegative.
    // Existing selections retain their limits; actual table bounds still apply.
    const double maximum_atomic_hydrogen=std::max(.95,maximum_warm_hydrogen);
    HydrogenOpacityContinuation warm_x(warm,.75,.05,maximum_warm_hydrogen);
    const std::string bridge_hydrogen_method=argc>=22?argv[21]:"oplib-ratio";
    if(bridge_hydrogen_method!="oplib-ratio" && bridge_hydrogen_method!="source-log")
      throw std::invalid_argument("unknown bridge hydrogen opacity method");
    HydrogenOpacityExtension bridge_ratio(bridge,dependence,.75,maximum_atomic_hydrogen),hot_x(hot,dependence,.75,maximum_atomic_hydrogen);
    HydrogenOpacityContinuation bridge_log(bridge,.75,.05,maximum_atomic_hydrogen);
    const Opacity& bridge_x=bridge_hydrogen_method=="source-log"
        ?static_cast<const Opacity&>(bridge_log):static_cast<const Opacity&>(bridge_ratio);
    BlendedOpacity warm_bridge(warm_x,bridge_x,5.05,5.10),upper(warm_bridge,hot_x,5.6,5.7);
    std::unique_ptr<Opacity> lower_upper;
    const double maximum_metal=argc>=18?std::stod(argv[17]):.05;
    const double minimum_metal=argc>=19?std::stod(argv[18]):.004;
    if(argc>=17) {
      const std::string method=argv[16];
      if(method!="linear-kappa" && method!="log-kappa")
        throw std::invalid_argument("unknown lower-metal opacity method");
      const auto selected_method=method=="linear-kappa"
          ?LowerMetalOpacityContinuation::Method::linear_kappa
          :LowerMetalOpacityContinuation::Method::logarithmic_kappa;
      if(hydrogen_share_opacity)
        lower_upper=std::make_unique<LowerMetalHydrogenShareOpacity>(upper,.01,minimum_metal,selected_method);
      else
        lower_upper=std::make_unique<LowerMetalOpacityContinuation>(upper,.01,minimum_metal,selected_method);
    }
    BlendedOpacity atomic(molecular,lower_upper?static_cast<const Opacity&>(*lower_upper):upper,4.4,4.47);
    std::unique_ptr<MixtureOpacity> metal_dependence;
    std::unique_ptr<MetalOpacityExtension> extended_atomic;
    if(argc>=15) {
      const std::string method=argv[14];
      if(method!="oplib-ratio" && method!="linear-kappa")
        throw std::invalid_argument("unknown upper-metal opacity method");
      metal_dependence=std::make_unique<MixtureOpacity>(argv[13]);
      extended_atomic=std::make_unique<MetalOpacityExtension>(atomic,*metal_dependence,.03,maximum_metal,
          method=="oplib-ratio"?MetalOpacityExtension::Method::source_ratio:
                               MetalOpacityExtension::Method::linear_kappa);
    }
    ElementalOpacity radiation(extended_atomic?static_cast<const Opacity&>(*extended_atomic):atomic);
    FixedMetalOpacity opacity(radiation);
    TabulatedConduction table(argv[6]);HotConduction conduction(table);
    const std::string atmosphere_reference_mode=argc>=23?argv[22]:"nested-reference";
    if(atmosphere_reference_mode!="nested-reference" && atmosphere_reference_mode!="low-z-reference"
        && atmosphere_reference_mode!="continuous-reference" && atmosphere_reference_mode!="trace-metal-interval")
      throw std::invalid_argument("unknown atmosphere reference mode");
    const bool trace_metal_interval=atmosphere_reference_mode=="trace-metal-interval";
    const bool continuous_reference=atmosphere_reference_mode=="continuous-reference";
    if((continuous_reference||trace_metal_interval)!=(argc==24))
      throw std::invalid_argument("continuous atmosphere mode requires exactly one interval manifest");
    const bool low_z_reference=atmosphere_reference_mode=="low-z-reference"||trace_metal_interval;
    CompositionAtmosphereGrid atmosphere_source(eos,argv[7],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    FixedMetalAtmosphere atmosphere(atmosphere_source,eos);
    std::unique_ptr<MetalResponseAtmosphere> metal_atmosphere;
    if(argc>=13 && !low_z_reference)metal_atmosphere=std::make_unique<MetalResponseAtmosphere>(
        eos,atmosphere_source,std::filesystem::path(argv[12]),
        MetalResponseAtmosphere::Approximation::separable_gs98_response,
        MetalResponseAtmosphere::Options{1e-5});
    std::unique_ptr<MetalResponseAtmosphere> lower_atmosphere;
    std::unique_ptr<PiecewiseMetalAtmosphere> joined_atmosphere;
    if(argc>=17 && !low_z_reference) {
      lower_atmosphere=std::make_unique<MetalResponseAtmosphere>(eos,*metal_atmosphere,
          std::filesystem::path(argv[15]),MetalResponseAtmosphere::Approximation::separable_gs98_response,
          MetalResponseAtmosphere::Options{0});
      joined_atmosphere=std::make_unique<PiecewiseMetalAtmosphere>(*metal_atmosphere,*lower_atmosphere,.01);
    }
    std::unique_ptr<MetalAtmosphereChain> depleted_join;
    if(argc>=20) {
      const Atmosphere& chain_reference=low_z_reference
          ?static_cast<const Atmosphere&>(atmosphere_source):static_cast<const Atmosphere&>(*joined_atmosphere);
      depleted_join=std::make_unique<MetalAtmosphereChain>(eos,chain_reference,
          std::filesystem::path(argv[19]),low_z_reference?atmosphere_source.reference_metallicity():.005);
    }
    const Atmosphere& original_boundary=depleted_join?static_cast<const Atmosphere&>(*depleted_join):
        joined_atmosphere?static_cast<const Atmosphere&>(*joined_atmosphere):
        metal_atmosphere?static_cast<const Atmosphere&>(*metal_atmosphere):atmosphere;
    std::unique_ptr<HydrogenIntervalAtmosphere> interval_atmosphere;
    if(continuous_reference)interval_atmosphere=std::make_unique<HydrogenIntervalAtmosphere>(
        eos,original_boundary,std::filesystem::path(argv[23]));
    std::unique_ptr<TraceHeliumAtmosphereGrid> trace_atmosphere;
    std::unique_ptr<MetalIntervalAtmosphere> trace_join;
    if(trace_metal_interval) {
      trace_atmosphere=std::make_unique<TraceHeliumAtmosphereGrid>(eos,argv[23],
          TraceHeliumAtmosphereGrid::Approximation::neglect_atmospheric_helium3,.0003);
      trace_join=std::make_unique<MetalIntervalAtmosphere>(eos,*trace_atmosphere,original_boundary,1e-6,1e-5);
    }
    const Atmosphere& boundary=trace_join?static_cast<const Atmosphere&>(*trace_join):
        interval_atmosphere?static_cast<const Atmosphere&>(*interval_atmosphere):original_boundary;
    ScreenedCollisionTransport collision(argv[9]);ConditionalMetalEnvelopeHeat transport(eos,collision,conduction,true,2e6,3e6,true);
    PPCNNetwork nuclear;PlasmaNeutrinoLosses losses;
    Physics physics{&eos,variable_opacity?static_cast<const Opacity*>(&radiation):&opacity,
        &nuclear,1.9,ConvectiveCriterion::ledoux,0.,0.,&losses};physics.microscopic=&transport;
    std::string line;std::cout<<std::setprecision(17);
    while(std::getline(std::cin,line))try {
      std::istringstream in(line);double age,dt,tolerance,heat_tolerance;std::size_t count;
      if(!(in>>age>>dt>>tolerance>>heat_tolerance>>count) || count!=512)throw std::invalid_argument("invalid control header");
      Model old;old.age=age;transport.prescribed={};
      for(std::size_t i=0;i<count;++i) {
        double mass,radius,rho,T,L,X,Y,Z;CNAbundances cn;
        if(!(in>>mass>>radius>>rho>>T>>L>>X>>Y>>cn[0]>>cn[1]>>cn[2]>>Z))throw std::invalid_argument("invalid physical composition row");
        auto c=solar_scaled(X,Z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
        c.X[1]=Y;c.X[2]-=Y;c.cn_molality=cn;c.cn_mass_convention=CNMassConvention::explicit_metal_mass;
        (void)cn_physical_ledger(c,cn);
        old.m.push_back(mass);old.y.push_back({std::log(radius),std::log(rho),std::log(T),L});old.comp.push_back(c);
      }
      in>>std::ws;if(!in.eof())throw std::invalid_argument("extra query values");old.M=old.m.back();
      const auto preserved=old;
      EvolutionOptions options;options.abundance_tolerance=tolerance;options.material_heat_tolerance=heat_tolerance;
      options.max_abundance_change=.005;options.relaxation.zone_threads=std::stoul(argv[10]);
      const auto start=std::chrono::steady_clock::now();const auto cpu_start=std::clock();
      const auto step=evolve_step(old,physics,boundary,dt,options);const auto& model=step.model;
      // Build the entire response before emitting it, including diagnostics.
      std::ostringstream output;auto* saved=std::cout.rdbuf(output.rdbuf());
      try {
        std::cout<<"{\"converged\":"<<step.converged<<",\"message\":"<<std::quoted(step.message)
          <<",\"variable_interior_opacity\":"<<variable_opacity
          <<",\"hydrogen_share_opacity\":"<<hydrogen_share_opacity
          <<",\"upper_metal_opacity_extension\":"<<bool(extended_atomic)
          <<",\"lower_metal_opacity_continuation\":"<<bool(lower_upper)
          <<",\"lower_metal_atmospheric_response\":"<<(bool(joined_atmosphere)||low_z_reference)
          <<",\"depleted_metal_atmospheric_response\":"<<bool(depleted_join)
          <<",\"maximum_opacity_metallicity\":"<<maximum_metal
          <<",\"minimum_opacity_metallicity\":"<<minimum_metal
          <<",\"maximum_warm_hydrogen\":"<<maximum_warm_hydrogen
          <<",\"bridge_hydrogen_method\":"<<std::quoted(bridge_hydrogen_method)
          <<",\"atmosphere_reference_mode\":"<<std::quoted(atmosphere_reference_mode)
          <<",\"variable_atmospheric_response\":"<<(bool(metal_atmosphere)||low_z_reference)
          <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
          <<",\"cpu_seconds\":"<<static_cast<double>(std::clock()-cpu_start)/CLOCKS_PER_SEC
          <<",\"iterations\":"<<step.coupling_iterations<<",\"residual\":"<<step.residual
          <<",\"abundance_residual\":"<<step.abundance_residual<<",\"material_heat_residual\":"<<step.material_heat_residual
          <<",\"nuclear_mass_balance\":"<<step.nuclear_mass_balance<<",\"luminosity_balance\":"<<step.luminosity_balance
          <<",\"nuclear_luminosity\":"<<step.nuclear_luminosity<<",\"neutrino_luminosity\":"<<step.neutrino_luminosity
          <<",\"elapsed_age_seconds\":"<<model.age-age
          <<",\"input_preserved\":"<<(old.comp==preserved.comp && old.age==preserved.age && old.m==preserved.m &&
              std::equal(old.y.begin(),old.y.end(),preserved.y.begin(),[](const Point& a,const Point& b){
                return a.lnr==b.lnr && a.lnrho==b.lnrho && a.lnT==b.lnT && a.L==b.L;}))
          <<",\"age_seconds\":"<<model.age<<",\"Teff\":"<<std::pow(model.y.back().L/(4*M_PI*constants::sigma_SB*std::pow(model.r(count-1),2)),.25)
          <<",\"convective_mass_fraction\":"<<step.convective_mass_fraction<<",\"model\":[";
        for(std::size_t i=0;i<count;++i){if(i)std::cout<<',';const auto& c=model.comp[i];const auto y=*c.cn_molality;
          array(std::array{model.m[i],model.r(i),model.rho(i),model.T(i),model.y[i].L,c.X[0],c.X[1],y[0],y[1],y[2],c.Z()});}
        std::cout<<"],\"total_metal_species_rates\":[";
        for(std::size_t i=0;i<step.total_metal_species_rates.size();++i){if(i)std::cout<<',';array(step.total_metal_species_rates[i]);}
        std::cout<<"],\"inert_metal_mass_fraction_change\":";
        const auto weights=nodal_mass_weights(model);double inert=0;
        for(std::size_t i=0;i<count;++i)inert+=weights[i]/model.M*(metal_cn_abundances(model.comp[i])[5]-metal_cn_abundances(old.comp[i])[5]);
        std::cout<<inert;
        if(step.converged) {
          transport.prescribed=step.total_metal_species_rates;
          const auto diffusivity=secular_mixing_diffusivities(model,physics);
          const auto regions=convective_mixing_regions(model,physics);
          std::cout<<",\"secular_diffusivity\":";array(diffusivity);
          std::cout<<",\"mixing_regions\":[";
          for(std::size_t i=0;i<regions.size();++i){if(i)std::cout<<',';array(std::array{regions[i].first,regions[i].second});}
          std::cout<<"],\"metal_boundary_fluxes\":[";
          for(std::size_t j=0;j+1<regions.size();++j) {
            const auto i=regions[j].second-1;
            auto f=common_metal_cn_flux(metal_microscopic_face(transport,i,model.m[i],model.m[i+1],
                model.y[i],model.comp[i],model.y[i+1],model.comp[i+1],false).species,
                model.comp[i],model.comp[i+1],false);
            const auto a=metal_cn_abundances(model.comp[i]),b=metal_cn_abundances(model.comp[i+1]);
            const double r=.5*(model.r(i)+model.r(i+1)),rho=.5*(model.rho(i)+model.rho(i+1));
            const double rate=std::pow(4*std::numbers::pi*r*r*rho,2)*diffusivity[i]/(model.m[i+1]-model.m[i]);
            for(std::size_t k=0;k<6;++k)f.rate[k]+=rate*(a[k]-b[k]);
            if(j)std::cout<<',';std::cout<<'['<<i<<',';array(f.rate);std::cout<<']';
          }
          std::cout<<"],\"physical_sources\":[";
          for(std::size_t i=0;i<count;++i) {
            if(i)std::cout<<',';const auto& c=model.comp[i];
            const auto p=nuclear.pp().eval(model.T(i),model.rho(i),c);
            const auto n=nuclear.cn().response(model.T(i),model.rho(i),c,*c.cn_molality).physical.state;
            array(std::array{p.dXdt[0]+n.dXdt[0],p.dXdt[1],p.dXdt[2]+n.dXdt[2],n.dXdt[3],n.dXdt[4],n.dXdt[5],0.});
          }
          std::cout<<"],\"energy_cells\":[";
          for(std::size_t i=0;i<count;++i) {
            if(i)std::cout<<',';
            const auto e=eos.eval(model.T(i),model.rho(i),model.comp[i]);
            const auto prev=eos.eval(old.T(i),old.rho(i),old.comp[i]);
            const auto n=nuclear.eval(model.T(i),model.rho(i),model.comp[i]);
            const auto loss=losses.eval(model.T(i),model.rho(i),model.comp[i]);
            array(std::array{e.E,e.P,prev.E,n.eps,n.eps_neutrino,loss.eps});
          }
          std::cout<<']';
        }
        std::cout<<"}\n";
      }catch(...){std::cout.rdbuf(saved);throw;}
      std::cout.rdbuf(saved);std::cout<<output.str()<<std::flush;
    }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
