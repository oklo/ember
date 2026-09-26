#include "ember/deuterium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>

using namespace ember;
int main() {
  int failures=0,cases=0;
  auto check=[&](bool okay,const char* label) {
    if(!okay) {++failures;std::printf("FAIL: %s\n",label);}
  };
  std::ifstream in(std::string(EMBER_TEST_DATA_DIR)+"/deuterium_reference.txt");
  check(bool(in),"independent reference opens");
  double T,rate,slope,worst=0;
  while(in>>T>>rate>>slope) {
    auto s=deuterium_bare_rate(T);++cases;
    worst=std::max(worst,std::abs(s.molar_rate/rate-1));
    check(std::abs(s.molar_rate/rate-1)<2e-9,"adaptive energy quadrature reference");
    check(std::abs(s.dlnrate_dlnT-slope)<2e-9*std::max(1.,std::abs(slope)),"temperature moment reference");
    check(pp_bare_rate(T,PPReaction::deuterium_p,PPRates::solar_fusion_iii).molar_rate==s.molar_rate,
          "public pp rate entry point dispatches explicitly");
  }
  check(cases==14,"all independent rates read");
  double worst_derivative=0,worst_mass=0,worst_heat=0;
  for(double xd:{0.,1e-12,2e-5,4e-5,.1}) {
    auto c=solar_scaled(.7-xd,.02);c.basis=AbundanceBasis::baryon_mass;
    c.X[2]-=xd; // Composition excludes the explicitly supplied D inventory.
    for(auto p:{std::pair{5e5,.01},std::pair{1e6,1.},std::pair{3e6,10.}}) {
      const auto [temp,rho]=p;
      const auto s=deuterium_capture(temp,rho,c,xd);
      const double r=s.molar_reactions_per_gram_second;
      check(s.source.eps_neutrino==0,"radiative capture has no neutrino source");
      check(s.source.dXdt[0]<=0 && s.dX_deuterium_dt<=0 && s.source.dXdt[1]>=0,
            "protons and D consumed, He3 produced");
      if(r==0) {check(xd==0 && s.source.eps==0,"zero initial D gives zero added heat");continue;}
      const double total=s.source.dXdt[0]+s.source.dXdt[1]+s.dX_deuterium_dt;
      worst_mass=std::max(worst_mass,std::abs(total)/(6*r));
      check(std::abs(s.source.dXdt[0]+s.dX_deuterium_dt/2+2*s.source.dXdt[1]/3)<1e-14*r,
            "charge conservation");
      const double rest=s.source.dXdt[0]*nuclides[0].A+s.source.dXdt[1]*nuclides[1].A/3
                        +s.dX_deuterium_dt*deuterium_atomic_mass/2;
      worst_heat=std::max(worst_heat,std::abs(-rest*constants::c*constants::c/s.source.eps-1));
      constexpr double h=1e-5;
      const double dt=(std::log(deuterium_capture(temp*std::exp(h),rho,c,xd).source.eps)
                      -std::log(deuterium_capture(temp*std::exp(-h),rho,c,xd).source.eps))/(2*h);
      const double dr=(std::log(deuterium_capture(temp,rho*std::exp(h),c,xd).source.eps)
                      -std::log(deuterium_capture(temp,rho*std::exp(-h),c,xd).source.eps))/(2*h);
      worst_derivative=std::max({worst_derivative,std::abs(dt-s.source.dlneps_dlnT),
                               std::abs(dr-s.source.dlneps_dlnRho)});
      const double q=s.heat_per_mole/constants::NA/(1e6*constants::eV);
      check(std::abs(q-5.4935)<.0002,"capture heat agrees with nuclear mass defect");
    }
  }
  check(worst_mass<1e-15,"baryon conservation");
  check(worst_heat<2e-12,"released heat equals lost rest energy");
  check(worst_derivative<3e-7,"screened thermal and density derivatives");
  check(deuterium_bare_rate(3000).molar_rate==0,"cold surface thermal capture is negligible");
  auto cold=solar_scaled(.69998,.02);cold.basis=AbundanceBasis::baryon_mass;cold.X[2]-=2e-5;
  check(deuterium_capture(3000,1e-6,cold,2e-5).source.eps==0,
        "cold molecular layers do not request ionized screening");
  for(double t:{0.,-1.,2.01e7,std::numeric_limits<double>::quiet_NaN()}) {
    bool rejected=false;try {deuterium_bare_rate(t);}catch(const std::domain_error&) {rejected=true;}
    check(rejected,"unsupported thermal state is rejected");
  }
  bool rejected=false;
  try {pp_bare_rate(1e6,PPReaction::deuterium_p,PPRates::legacy);}
  catch(const std::invalid_argument&) {rejected=true;}
  check(rejected,"no silent legacy D rate");
  auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
  rejected=false;try {deuterium_capture(1e6,1.,c,2e-5);}catch(const std::domain_error&) {rejected=true;}
  check(rejected,"double-counted D mass is rejected");
  std::printf("cases=%d worst_rate=%.4g worst_mass=%.4g worst_heat=%.4g worst_derivative=%.4g failures=%d\n",
              cases,worst,worst_mass,worst_heat,worst_derivative,failures);
  return failures?1:0;
}
