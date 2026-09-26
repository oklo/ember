// Coupled 0.5-Msun regression with analytic EOS/opacity and a grey
// boundary. This tests the integration, not a physical atmosphere or track.
#include "../examples/stellar_seed.hpp"
#include "ember/eos_composite.hpp"
#include "ember/metal_cn_transport.hpp"
#include "ember/deuterium_burning.hpp"
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace ember;
namespace {
class ControlOpacity final:public Opacity {
 public:
  OpacityState eval(double T,double rho,const Composition& c)const override {
    const double absorption=1e23*rho*std::pow(T,-3.5),base=absorption+.2*(1+c.X[0]);
    const double z=std::log(T/3e5)/.5,bump=20*std::exp(-z*z);
    return {base*(1+bump),-3.5*absorption/base-4*z*bump/(1+bump),absorption/base};
  }
  const char* name()const override{return "analytic opacity for coupled regression";}
};
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
}
int main(int argc,char**) {
  const auto start=std::clock();
  try {
    CompositeEos eos;ControlOpacity opacity;GreyAtmosphere atmosphere(eos,opacity);
    PPCNNetwork nuclear;
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.001;c.X[2]-=.001;
    c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);
    Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux};
    physics.explicit_metal_mixing_only=true;
    auto seed_physics=physics;
    const auto grid=argc>1?LuminosityGrid::volume_faces:LuminosityGrid::mass_nodes;
    auto initial=example::stellar_seed(256,.5*constants::Msun,.6*constants::Rsun,c,nuclear,atmosphere,3.,0.,grid);
    RelaxationOptions ro;ro.max_iterations=200;
    const auto relaxed=relax(initial,seed_physics,atmosphere,ro);
    require(relaxed.converged,"seed: "+relaxed.message);initial=relaxed.model;
    // Declared numerical coupling test: inject trace D into the relaxed
    // thermal control. This is not a formation model or a production restart.
    for(auto& q:initial.comp){q.X[0]-=2e-8;q[Species::H2]=2e-8;}
    EvolutionOptions options;options.relaxation=ro;options.abundance_tolerance=1e-14;
    options.max_coupling_iterations=80;
    constexpr double dt=1000*31557600.;
    const auto full=evolve_step(initial,physics,atmosphere,dt,options);
    require(full.converged,"combined step: "+full.message);
    const auto half=evolve_step(initial,physics,atmosphere,dt/2,options);
    require(half.converged,"first half: "+half.message);
    const auto end=evolve_step(half.model,physics,atmosphere,dt/2,options);
    require(end.converged,"second half: "+end.message);
    double time_error=0,first_law=0,mass_error=0,initial_D=0,final_D=0;
    const auto weights=nodal_mass_weights(initial);
    for(std::size_t i=0;i<initial.size();++i) {
      const auto& a=full.model.y[i];const auto& q=end.model.y[i];
      time_error=std::max({time_error,std::abs(a.lnr-q.lnr),std::abs(a.lnrho-q.lnrho),std::abs(a.lnT-q.lnT)});
      require(full.model.comp[i].cn_molality.has_value(),"CN inventory lost");
      require(full.model.comp[i][Species::H2]<=initial.comp[i][Species::H2]*(1+1e-12),"D increased");
      initial_D+=weights[i]/initial.M*initial.comp[i][Species::H2];
      final_D+=weights[i]/initial.M*end.model.comp[i][Species::H2];
    }
    for(const auto* s:{&full,&half,&end}) {
      first_law=std::max(first_law,std::abs(s->luminosity_balance));
      mass_error=std::max(mass_error,std::abs(s->nuclear_mass_balance));
    }
    require(time_error<1e-4,"unresolved control interval");
    require(final_D<.99*initial_D,"total D did not burn");
    require(first_law<2e-7 && mass_error<2e-5,"coupled energy accounting");
    require(full.model.age==initial.age+dt && end.model.age==initial.age+dt,"clock mismatch");
    require(full.model.luminosity_grid==grid && end.model.luminosity_grid==grid,"luminosity grid changed");
    std::cout<<std::setprecision(17)<<"{\"outcome\":\"passed\",\"scope\":\"analytic EOS/grey coupled 0.5 solar mass test with trace-D injection; not a physical track\""
      <<",\"face_luminosities\":"<<face_luminosities(initial)
      <<",\"full_half_structure_error\":"<<time_error<<",\"maximum_first_law_error\":"<<first_law
      <<",\"maximum_mass_defect_error\":"<<mass_error<<",\"mean_initial_D\":"<<initial_D
      <<",\"mean_final_D\":"<<final_D<<",\"cpu_seconds\":"<<double(std::clock()-start)/CLOCKS_PER_SEC<<"}\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
