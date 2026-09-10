#include "ember/eos_mixture.hpp"
#include "ember/conduction_table.hpp"
#include "ember/convection.hpp"
#include "ember/evolution.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* label,double error=0) {
  failures+=!ok;std::printf("[%s] %s (%.8g)\n",ok?"PASS":"FAIL",label,error);
}
template<class F>bool throws(F f){try{f();}catch(const std::exception&){return true;}return false;}
Composition composition(double x,double y) {
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=y;c.X[2]-=y;return c;
}
}
int main() {
  const std::string data=EMBER_DATA_DIR;
  MetalHelmholtzEos eos(data+"/eos/freeeos300_gs98_z020.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
  auto c=composition(.55,.06);constexpr double h=1e-6;
  double firstlaw=0,derivatives=0,transport=0,buoyancy_response=0;bool stabilizing=true;
  StellarMixtureOpacity opacity(data+"/opacity");TabulatedConduction conduction(data+"/conduction/condtab21wd_metals.dat");
  PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
  for(auto [T,rho]:{std::pair{4700.,1e-4},std::pair{32000.,.02},std::pair{2.31e6,17.},std::pair{5.11e6,267.}}) {
    const auto e=eos.eval_with_derivatives(T,rho,c);const auto d=eos.composition_response(T,rho,c);const auto& s=e.state;
    const auto a=eos.eval(T,rho*std::exp(h),c),b=eos.eval(T,rho*std::exp(-h),c);
    firstlaw=std::max({firstlaw,std::abs((a.E-b.E)/(2*h)-s.P/rho*(1-s.chiT))/(T*s.cv),
        std::abs((a.S-b.S)/(2*h)+s.P*s.chiT/(rho*T))/s.cv,std::abs(s.P*s.delta/(rho*T*s.cp*s.grad_ad)-1)});
    for(int k=0;k<2;++k) {
      auto plus=c,minus=c;plus.X[k]+=h;plus.X[2]-=h;minus.X[k]-=h;minus.X[2]+=h;
      const auto p=eos.eval(T,rho,plus),m=eos.eval(T,rho,minus);
      derivatives=std::max({derivatives,std::abs((p.P-m.P)/(2*h)-d.dP[k])/s.P,
          std::abs((p.E-m.E)/(2*h)-d.dE[k])/(T*s.cv)});
      if(T>1e6) {
        auto verify=[&](const auto& op) {
          const auto v=op.eval(T,rho,c);
          const double measured=(std::log(op.eval(T,rho,plus).kappa)-std::log(op.eval(T,rho,minus).kappa))/(2*h);
          transport=std::max(transport,std::abs(measured-(k==0?v.dlnk_dX:v.dlnk_dY3)));
        };verify(opacity);verify(conduction);
      }
    }
    const auto p=eos.eval_with_derivatives(T*std::exp(h),rho,c),m=eos.eval_with_derivatives(T*std::exp(-h),rho,c);
    derivatives=std::max({derivatives,std::abs((p.state.cp-m.state.cp)/(2*h)-e.dcp_dlnT)/s.cp,
        std::abs((p.state.delta-m.state.delta)/(2*h)-e.ddelta_dlnT)/s.delta,
        std::abs((p.state.grad_ad-m.state.grad_ad)/(2*h)-e.dgrad_ad_dlnT)/s.grad_ad});
    const auto inner=composition(.54,.06),outer=composition(.56,.06);
    const auto buoyancy=composition_buoyancy(eos,T,s.P,s.delta,-.1,inner,outer,rho);
    stabilizing=stabilizing && buoyancy.B>0;
    constexpr double step=1e-5;
    for(int k=0;k<2;++k) {
      const auto plus=composition_buoyancy(eos,T*std::exp(k==0?step:0),s.P*std::exp(k==1?step:0),s.delta,-.1,inner,outer,rho);
      const auto minus=composition_buoyancy(eos,T*std::exp(k==0?-step:0),s.P*std::exp(k==1?-step:0),s.delta,-.1,inner,outer,rho);
      const double analytic=k==0?buoyancy.dB_dlnT:buoyancy.dB_dlnP;
      buoyancy_response=std::max(buoyancy_response,std::abs((plus.B-minus.B)/(2*step)-analytic));
    }
  }
  check(firstlaw<1e-5,"metal potential preserves first law and Maxwell identities",firstlaw);
  check(derivatives<1e-5,"thermal and composition responses differentiate the same potential",derivatives);
  check(transport<3e-6,"GS98 opacity/conduction composition responses preserve source mapping",transport);
  check(stabilizing && buoyancy_response<3e-6,"metal EOS Ledoux barriers and derivatives agree with independent density inversions",buoyancy_response);
  std::ifstream reference(std::string(EMBER_TEST_DATA_DIR)+"/freeeos300_gs98_reference.dat");
  double x,y,T,rho,P,E,cv,cp,ad,pressure_error=0,energy_error=0,response_error=0;int count=0;
  while(reference>>x>>y>>T>>rho>>P>>E>>cv>>cp>>ad) {
    const auto s=eos.eval(T,rho,composition(x,y));++count;
    pressure_error=std::max(pressure_error,std::abs(s.P/P-1));energy_error=std::max(energy_error,std::abs(s.E/E-1));
    response_error=std::max({response_error,std::abs(s.cv/cv-1),std::abs(s.cp/cp-1),std::abs(s.grad_ad/ad-1)});
  }
  check(count==144 && pressure_error<.001,"144 independent metal-bearing source pressures",pressure_error);
  check(energy_error<.002 && response_error<.003,"independent source energy and thermal responses",response_error);
  auto invalid=c;invalid.basis=AbundanceBasis::atomic_mass;
  check(throws([&]{eos.eval(1e6,1,invalid);}),"metal family rejects wrong abundance basis");
  invalid=c;invalid.metal_inventory=MetalInventory::carried_isotopes;
  check(throws([&]{eos.eval(1e6,1,invalid);}),"metal family requires explicit GS98 inventory");
  invalid=composition(.299,0);
  check(throws([&]{eos.eval(1e6,1,invalid);}),"hydrogen extrapolation rejected");
  invalid=composition(.5,.121);
  check(throws([&]{eos.eval(1e6,1,invalid);}),"helium3 extrapolation rejected");
  Model m;m.M=2;m.m={.1,1,2};m.comp={c,c,c};m.y.assign(3,{0,std::log(200.),std::log(5e6),0});
  const auto burned=burn_and_mix(m,m,nuclear,{{0,3}},1e10);
  check(std::all_of(burned.begin(),burned.end(),[](const Composition& v){return v.metal_inventory==MetalInventory::gs98;}),
        "implicit burning/mixing preserves elemental inventory identity");
  double sum=0,ions=0,electrons=0;
  for(const auto& e:gs98_metals){sum+=e.fraction;ions+=e.fraction/e.mass_number;electrons+=e.fraction*e.charge/e.mass_number;}
  check(std::abs(sum-1)<1e-14 && std::abs(c.mu_ions_inv()-(c.X[0]+c.X[1]/3+c.X[2]/4+c.Z()*ions))<1e-15
      && std::abs(c.mu_elec_inv()-(c.X[0]+2*c.X[1]/3+c.X[2]/2+c.Z()*electrons))<1e-15,
      "GS98 ion/electron inventory closes with conserved baryonic mass");
  return failures?1:0;
}
