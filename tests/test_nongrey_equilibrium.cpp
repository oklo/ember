#include "../examples/stellar_seed.hpp"
#include "ember/eos_helmholtz.hpp"
#include "ember/atmosphere_table.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include "ember/relaxation.hpp"
#include <cstdio>

using namespace ember;
int main() {
  const std::string data=EMBER_DATA_DIR;
  HelmholtzTableEos eos(data+"/eos/freeeos300_hhe_x070_potential.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
  TabulatedAtmosphere atmosphere(eos,data+"/atmosphere/cond_gn93_tau100_solar_proxy.dat",TabulatedAtmosphere::Mixture::allow_documented_proxy);
  AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat"); TopsOpacity high(data+"/opacity");
  BlendedOpacity opacity(low,high,4.4,4.5); PPChains nuclear;
  Physics physics{&eos,&opacity,&nuclear,1.9}; const auto comp=solar_scaled(.7,.02);
  int failures=0;
  auto check=[&](bool ok,const char* label,double value=0.) {
    if(!ok) ++failures;
    std::printf("[%s] %s (%.9g)\n",ok?"PASS":"FAIL",label,value);
  };
  double previous_R=0,previous_L=0,previous_virial=0;
  for(std::size_t n:{1024UL,2048UL,4096UL}) {
    const auto seed=example::stellar_seed(n,.1*constants::Msun,.15*constants::Rsun,comp,nuclear,atmosphere,1.5);
    // Deep matching requires the atmosphere's local T, not Teff.
    const double seed_Teff=std::pow(seed.y.back().L/(4*M_PI*constants::sigma_SB*std::pow(seed.r(n-1),2)),.25);
    check(seed.T(n-1)>1.3*seed_Teff,"deep atmospheric seed uses local matching temperature");
    const auto result=relax(seed,physics,atmosphere);
    check(result.converged && result.residual<1e-9 && result.correction<1e-8,"non-grey potential-EOS equilibrium converges",result.residual);
    if(!result.converged) return 1;
    const auto& m=result.model;
    double L=m.m[0]*nuclear.eval(m.T(0),m.rho(0),comp).eps;
    double pressure=3*m.m[0]*eos.eval(m.T(0),m.rho(0),comp).P/m.rho(0);
    double gravity=.6*constants::G*m.m[0]*m.m[0]/m.r(0);
    bool physical=true; double maxwell=0;
    for(std::size_t i=0;i<n;++i) {
      const auto e=eos.eval(m.T(i),m.rho(i),comp);
      maxwell=std::max(maxwell,std::abs(e.P*e.delta/(m.rho(i)*m.T(i)*e.cp*e.grad_ad)-1));
      physical &= e.cv>0 && e.cp>0 && e.chiRho>0 && std::isfinite(e.E) && m.y[i].L>0;
      if(!i) continue;
      const double dm=m.m[i]-m.m[i-1];
      physical &= dm>0 && m.r(i)>m.r(i-1) && m.T(i)<m.T(i-1) && m.rho(i)<m.rho(i-1);
      L+=.5*dm*(nuclear.eval(m.T(i),m.rho(i),comp).eps+nuclear.eval(m.T(i-1),m.rho(i-1),comp).eps);
      pressure+=1.5*dm*(e.P/m.rho(i)+eos.eval(m.T(i-1),m.rho(i-1),comp).P/m.rho(i-1));
      gravity+=.5*constants::G*dm*(m.m[i]/m.r(i)+m.m[i-1]/m.r(i-1));
    }
    const double surface=4*M_PI*std::pow(m.r(n-1),3)*eos.eval(m.T(n-1),m.rho(n-1),comp).P;
    const double virial=std::abs((pressure-surface)/gravity-1);
    check(physical,"ordered profile with positive heat capacities and finite internal energy");
    check(maxwell<2e-14,"pressure/entropy identity holds throughout the stellar profile",maxwell);
    check(std::abs(L/m.y.back().L-1)<1e-8,"nuclear integral matches surface luminosity",L/m.y.back().L-1);
    if(previous_R>0) {
      check(std::abs(m.r(n-1)/previous_R-1)<2e-4 && std::abs(m.y.back().L/previous_L-1)<.001,
            "mesh doubling changes radius by <0.02% and luminosity by <0.1%");
      check(virial<.3*previous_virial,"independent virial error decreases quadratically",virial/previous_virial);
    }
    previous_R=m.r(n-1);previous_L=m.y.back().L;previous_virial=virial;
  }
  return failures?1:0;
}
