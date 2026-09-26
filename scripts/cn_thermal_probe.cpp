// Coupled stellar control with explicit CN abundances. Fixed-GS98 material
// tables; optional secular mixing and explicitly chosen trace CN transport.
#include "ember/eos_smooth_mixture.hpp"
#include "conditional_envelope_heat.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/conduction_table.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/cn_burning.hpp"
#include "ember/cn_transport.hpp"
#include "ember/constants.hpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
using namespace ember;
namespace {
// Build each reply atomically, including failures in optional diagnostics.
struct BufferedReply {
  std::ostringstream text;std::streambuf* original;
  BufferedReply():original(std::cout.rdbuf(text.rdbuf())) {}
  ~BufferedReply(){if(original)std::cout.rdbuf(original);}
  void finish(){std::cout.rdbuf(original);original=nullptr;std::cout<<text.str()<<std::flush;}
};
double teff(const Model& m) {
  return std::pow(m.y.back().L/(4*std::numbers::pi*constants::sigma_SB*std::pow(m.r(m.size()-1),2)),.25);
}
template<class T> void array(const T& a) {
  std::cout<<'[';bool first=true;for(const auto& x:a){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';
}
class PairedAtmosphere final:public Atmosphere {
public:
  const Atmosphere& source;const Eos& eos;double temperature,pressure;
  PairedAtmosphere(const Atmosphere& a,const Eos& e,double t,double p):source(a),eos(e),temperature(t),pressure(p) {
    if(!(t>0 && p>0) || !std::isfinite(t+p))throw std::invalid_argument("invalid boundary factors");
  }
  AtmosphereState eval(double Teff,double g,const Composition& c)const override {
    auto a=source.eval(Teff,g,c);if(temperature==1 && pressure==1)return a;
    const double oldP=a.P,rad=constants::a_rad*std::pow(a.T,4)/3;
    const double newrad=rad*std::pow(temperature,4);
    const double newP=pressure*(oldP-rad)+newrad;
    auto derivative=[&](double dP,double dT){return (pressure*(oldP*dP-4*rad*dT)+4*newrad*dT)/newP;};
    a.dlnP_dlnTeff=derivative(a.dlnP_dlnTeff,a.dlnT_dlnTeff);
    a.dlnP_dlng=derivative(a.dlnP_dlng,a.dlnT_dlng);
    a.T*=temperature;a.P=newP;a.Pgas=newP-newrad;
    a.rho=eos.rho_from_PT(a.T,a.P,c,a.rho);
    return a;
  }
  const char* name()const override{return "paired matching-layer T and gas-pressure sensitivity";}
};
}
int main(int argc,char** argv) {
  if(argc!=9 && argc!=11 && argc!=12 && argc!=13)return 2;
  try {
    const auto begin=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    MixtureOpacity molecular(argv[2]),warm(argv[3]),bridge(argv[4]),hot(argv[5]),retained_warm(argv[8]);
    TabulatedConduction conduction(argv[6]);
    CompositionAtmosphereGrid source_atmosphere(eos,argv[7],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    BlendedOpacity warm_bridge(warm,bridge,5.05,5.10),joined_hot(warm_bridge,hot,5.6,5.7);
    BlendedOpacity atomic(molecular,joined_hot,4.4,4.47);
    CombinedOpacity opacity(std::make_shared<ElementalOpacity>(atomic),std::make_shared<HotConduction>(conduction));
    BlendedOpacity retained_bridge(retained_warm,bridge,5.05,5.10),retained_hot(retained_bridge,hot,5.6,5.7);
    BlendedOpacity retained_atomic(molecular,retained_hot,4.4,4.47);
    CombinedOpacity retained_opacity(std::make_shared<ElementalOpacity>(retained_atomic),std::make_shared<HotConduction>(conduction));
    ElementalOpacity radiative(atomic),retained_radiative(retained_atomic);
    HotConduction material(conduction);
    std::unique_ptr<ScreenedCollisionTransport> collisions;
    std::unique_ptr<ConditionalEnvelopeHeat> transport;
    if(argc>=12) {
      collisions=std::make_unique<ScreenedCollisionTransport>(argv[11]);
      transport=std::make_unique<ConditionalEnvelopeHeat>(eos,*collisions,material,true,2e6,3e6,true);
    }
    PPCNNetwork nuclear;PlasmaNeutrinoLosses losses;
    Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux,0.,0.,&losses};
    if(argc>=11){physics.alpha_semiconvection=std::stod(argv[9]);physics.alpha_thermohaline=std::stod(argv[10]);}
    physics.microscopic=transport.get();
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<'\n';
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line))try {
      std::istringstream input(line);int warm_choice;double tf,pf,dt,tolerance;std::size_t count;
      if(!(input>>warm_choice>>tf>>pf>>dt>>tolerance>>count) || count!=512 || (warm_choice!=0 && warm_choice!=1))
        throw std::invalid_argument("invalid query header");
      physics.opacity=warm_choice?&retained_opacity:&opacity;
      double heat_tolerance=1e-6;
      if(transport) {
        int mode;
        if(!(input>>mode>>heat_tolerance) || mode<0 || mode>2)throw std::invalid_argument("invalid CN transport selection");
        physics.cn_microscopic=mode?CNMicroscopicApproximation::helium_velocity:CNMicroscopicApproximation::zero_catalyst_drift;
        transport->fully_convective_only=mode==2;
        physics.opacity=warm_choice?&retained_radiative:&radiative;
        transport->prescribed={};
      }
      Model old;
      for(std::size_t i=0;i<count;++i) {
        double m,r,rho,T,L,X,X3;CNAbundances y;
        if(!(input>>m>>r>>rho>>T>>L>>X>>X3>>y[0]>>y[1]>>y[2]))throw std::invalid_argument("invalid query model");
        auto c=solar_scaled(X,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
        c.X[1]=X3;c.X[2]=.98-X-X3;c.cn_molality=y;(void)cn_physical_ledger(c,y);
        old.m.push_back(m);old.y.push_back({std::log(r),std::log(rho),std::log(T),L});old.comp.push_back(c);
      }
      input>>std::ws;if(!input.eof())throw std::invalid_argument("extra query values");
      old.M=old.m.back();const auto preserved=old;
      PairedAtmosphere atmosphere(source_atmosphere,eos,tf,pf);
      EvolutionOptions options;options.abundance_tolerance=tolerance;options.max_abundance_change=.001;
      if(argc==13)options.relaxation.zone_threads=std::stoul(argv[12]);
      options.max_coupling_iterations=30;options.material_heat_tolerance=heat_tolerance;
      const auto start=std::chrono::steady_clock::now();
      const auto cpu_start=std::clock();
      const auto step=evolve_step(old,physics,atmosphere,dt,options);
      const double cpu_seconds=static_cast<double>(std::clock()-cpu_start)/CLOCKS_PER_SEC;
      const auto& model=step.model;
      BufferedReply reply;
      std::cout<<"{\"converged\":"<<step.converged<<",\"message\":"<<std::quoted(step.message)
        <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
        <<",\"cpu_seconds\":"<<cpu_seconds
        <<",\"iterations\":"<<step.coupling_iterations<<",\"residual\":"<<step.residual
        <<",\"abundance_residual\":"<<step.abundance_residual<<",\"nuclear_mass_balance\":"<<step.nuclear_mass_balance
        <<",\"material_heat_residual\":"<<step.material_heat_residual
        <<",\"luminosity_balance\":"<<step.luminosity_balance<<",\"nuclear_luminosity\":"<<step.nuclear_luminosity
        <<",\"neutrino_luminosity\":"<<step.neutrino_luminosity<<",\"Teff\":"<<teff(model)
        <<",\"convective_mass_fraction\":"<<step.convective_mass_fraction
        <<",\"elapsed_age_seconds\":"<<model.age<<",\"input_preserved\":"<<(old.comp==preserved.comp && old.age==preserved.age)
        <<",\"model\":[";
      for(std::size_t i=0;i<count;++i) {
        if(i)std::cout<<',';const auto& c=model.comp[i];const auto& y=*c.cn_molality;
        array(std::array{model.m[i],model.r(i),model.rho(i),model.T(i),model.y[i].L,c.X[0],c.X[1],y[0],y[1],y[2]});
      }
      if(!step.converged){std::cout<<"]}\n";reply.finish();continue;}
      if(transport)transport->prescribed=step.total_species_rates;
      const auto diffusivity=secular_mixing_diffusivities(model,physics);
      std::cout<<"],\"secular_diffusivity\":";array(diffusivity);
      std::cout<<",\"mixing_regions\":[";
      const auto regions=convective_mixing_regions(model,physics);
      for(std::size_t i=0;i<regions.size();++i){if(i)std::cout<<',';array(std::array{regions[i].first,regions[i].second});}
      std::cout<<"],\"cn_boundary_fluxes\":[";
      if(transport)for(std::size_t j=0;j+1<regions.size();++j) {
        const auto i=regions[j].second-1;
        auto f=trace_cn_flux(microscopic_face(*transport,i,model.m[i],model.m[i+1],
          model.y[i],model.comp[i],model.y[i+1],model.comp[i+1],false).species,
          model.comp[i],model.comp[i+1],physics.cn_microscopic,false);
        const auto a=cn_transport_abundances(model.comp[i]),b=cn_transport_abundances(model.comp[i+1]);
        const double r=.5*(model.r(i)+model.r(i+1)),rho=.5*(model.rho(i)+model.rho(i+1));
        const double g=std::pow(4*std::numbers::pi*r*r*rho,2)*diffusivity[i]/(model.m[i+1]-model.m[i]);
        for(std::size_t k=0;k<5;++k)f.rate[k]+=g*(a[k]-b[k]);
        if(j)std::cout<<',';std::cout<<'['<<i<<',';array(f.rate);std::cout<<']';
      }
      std::cout<<"],\"total_species_rates\":[";
      for(std::size_t i=0;i<step.total_species_rates.size();++i){if(i)std::cout<<',';array(step.total_species_rates[i]);}
      std::cout<<"],\"physical_sources\":[";
      for(std::size_t i=0;i<count;++i) {
        if(i)std::cout<<',';const auto& c=model.comp[i];
        const auto p=nuclear.pp().eval(model.T(i),model.rho(i),c);
        const auto n=nuclear.cn().response(model.T(i),model.rho(i),c,*c.cn_molality).physical.state;
        array(std::array{p.dXdt[0]+n.dXdt[0],p.dXdt[1],p.dXdt[2]+n.dXdt[2],n.dXdt[3],n.dXdt[4],n.dXdt[5]});
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
      std::cout<<"]}\n";reply.finish();
    }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
