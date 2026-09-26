// Read-only transport diagnosis of an evolved CN checkpoint. Prescribed
// total species rates reproduce the heat flux used by its converged step.
#include "ember/eos_smooth_mixture.hpp"
#include "conditional_envelope_heat.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/conduction_table.hpp"
#include "ember/cn_transport.hpp"
#include "ember/nuclear_cn.hpp"
#include "ember/convection.hpp"
#include "thermal_transport.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace ember;
int main(int argc,char** argv) {
  if(argc!=9)return 2;
  try {
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    MixtureOpacity molecular(argv[2]),warm(argv[3]),bridge(argv[4]),hot(argv[5]);
    BlendedOpacity wb(warm,bridge,5.05,5.10),wh(wb,hot,5.6,5.7),atomic(molecular,wh,4.4,4.47);
    ElementalOpacity radiation(atomic);
    TabulatedConduction conduction(argv[6]);HotConduction material(conduction);
    ScreenedCollisionTransport collisions(argv[7]);
    ConditionalEnvelopeHeat transport(eos,collisions,material,true,2e6,3e6,true);
    PPCNNetwork nuclear;
    Physics physics{&eos,&radiation,&nuclear,1.9,ConvectiveCriterion::ledoux,.1,1.};
    physics.microscopic=&transport;physics.cn_microscopic=CNMicroscopicApproximation::helium_velocity;
    std::ifstream input(argv[8]);std::size_t count;input>>count;
    if(count!=512)throw std::invalid_argument("512-point checkpoint required");
    Model model;
    for(std::size_t i=0;i<count;++i) {
      double m,r,rho,T,L,H,He3;CNAbundances y;
      if(!(input>>m>>r>>rho>>T>>L>>H>>He3>>y[0]>>y[1]>>y[2]))throw std::invalid_argument("invalid model");
      auto c=solar_scaled(H,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=He3;c.X[2]=.98-H-He3;c.cn_molality=y;(void)cn_physical_ledger(c,y);
      model.m.push_back(m);model.y.push_back({std::log(r),std::log(rho),std::log(T),L});model.comp.push_back(c);
    }
    model.M=model.m.back();std::vector<SpeciesVector> rates(count-1);
    for(auto& rate:rates)for(double& value:rate)if(!(input>>value))throw std::invalid_argument("missing total species rate");
    input>>std::ws;if(!input.eof())throw std::invalid_argument("extra input");transport.prescribed=rates;
    const auto regions=convective_mixing_regions(model,physics);
    const auto secular=secular_mixing_diffusivities(model,physics);
    std::vector<EosState> eos_states;std::vector<double> opacity;
    for(std::size_t i=0;i<count;++i) {
      eos_states.push_back(eos.eval(model.T(i),model.rho(i),model.comp[i]));
      opacity.push_back(radiation.eval(model.T(i),model.rho(i),model.comp[i]).kappa);
    }
    std::ostringstream out;out<<std::setprecision(17)<<"{\"mixing_regions\":[";
    for(std::size_t i=0;i<regions.size();++i){if(i)out<<',';out<<'['<<regions[i].first<<','<<regions[i].second<<']';}
    out<<"],\"faces\":[";
    for(std::size_t i=0;i+1<count;++i) {
      const double T=.5*(model.T(i)+model.T(i+1)),rho=.5*(model.rho(i)+model.rho(i+1));
      const double P=.5*(eos_states[i].P+eos_states[i+1].P),k=.5*(opacity[i]+opacity[i+1]);
      const double m=.5*(model.m[i]+model.m[i+1]),L=.5*(model.y[i].L+model.y[i+1].L);
      const auto heat=microscopic_heat(transport,i,model.m[i],model.m[i+1],model.y[i],model.comp[i],model.y[i+1],model.comp[i+1],false);
      const auto thermal=detail::thermal_transport<0>(T,rho,P,m,k,L,true,heat.carried_luminosity,heat.conductivity,model.y[i].lnT,model.y[i+1].lnT);
      const auto buoyancy=composition_buoyancy(eos,T,P,.5*(eos_states[i].delta+eos_states[i+1].delta),std::log(eos_states[i+1].P/eos_states[i].P),model.comp[i],model.comp[i+1],rho);
      const double ad=.5*(eos_states[i].grad_ad+eos_states[i+1].grad_ad);
      const double Krad=4*constants::a_rad*constants::c*T*T*T/(3*rho*k);
      if(i)out<<',';
      out<<"{\"face\":"<<i<<",\"q\":"<<m/model.M<<",\"T\":"<<T<<",\"rho\":"<<rho
         <<",\"radiative_opacity\":"<<k<<",\"required_gradient\":"<<thermal.gradient.value
         <<",\"adiabatic_gradient\":"<<ad<<",\"composition_buoyancy\":"<<buoyancy.B
         <<",\"conductive_fraction_of_diffusive_heat\":"<<heat.conductivity/(heat.conductivity+Krad)
         <<",\"carried_heat_over_local_L\":"<<heat.carried_luminosity/L<<",\"secular_D\":"<<secular[i]
         <<",\"convective\":"<<(thermal.gradient.value>ad+buoyancy.B?"true":"false")<<'}';
    }
    out<<"]}\n";std::cout<<out.str();
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
