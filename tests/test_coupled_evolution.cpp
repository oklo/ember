#include "../examples/stellar_seed.hpp"
#include "ember/eos_composition.hpp"
#include "ember/evolution.hpp"
#include "ember/evolution_proxies.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include <algorithm>
#include <cstdio>

using namespace ember;
int main() {
  const std::string data=EMBER_DATA_DIR;
  CompositionHelmholtzEos eos(data+"/eos/freeeos300_hhe_composition.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
  AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat");TopsOpacity high(data+"/opacity",TopsOpacity::Grid::composition);
  BlendedOpacity blend(low,high,4.4,4.5);NominalAbundanceOpacity opacity(blend);
  TabulatedAtmosphere cond(eos,data+"/atmosphere/cond_gn93_tau100_solar_proxy.dat",TabulatedAtmosphere::Mixture::allow_documented_proxy);
  FrozenCompositionAtmosphere atmosphere(eos,cond);
  PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
  Physics physics{&eos,&opacity,&nuclear,1.9};
  auto comp=solar_scaled(.7,.02);comp.basis=AbundanceBasis::baryon_mass;
  const auto seed=example::stellar_seed(512,.1*constants::Msun,.15*constants::Rsun,comp,nuclear,atmosphere,1.5);
  const auto equilibrium=relax(seed,physics,atmosphere);
  int failures=0;
  auto check=[&](bool ok,const char* label,double value=0) {
    failures+=!ok;std::printf("[%s] %s (%.9g)\n",ok?"PASS":"FAIL",label,value);
  };
  check(equilibrium.converged,"baryonic composition-dependent starting equilibrium converges");
  if(!equilibrium.converged)return 1;
  const auto& initial=equilibrium.model;
  constexpr double duration=1e9*365.25*86400.;
  double balance=0,restmass=0,abundance=0,nonlinear=0;bool mixed=true,physical=true;
  std::vector<Model> final;
  EvolutionOptions options;options.relaxation.residual_tolerance=1e-10;options.relaxation.correction_tolerance=1e-9;
  auto advance=[&](int steps) {
    auto model=initial;
    for(int k=0;k<steps;++k) {
      const auto result=evolve_step(model,physics,atmosphere,duration/steps,options);
      if(!result.converged) {std::printf("evolution failed: %s\n",result.message.c_str());return false;}
      model=result.model;
      balance=std::max(balance,std::abs(result.luminosity_balance));restmass=std::max(restmass,std::abs(result.nuclear_mass_balance));
      abundance=std::max(abundance,result.abundance_residual);nonlinear=std::max(nonlinear,result.residual);
      mixed&=result.mixed_regions==1 && std::abs(result.convective_mass_fraction-1)<1e-12;
      for(std::size_t i=0;i<model.size();++i) {
        mixed&=model.comp[i].X==model.comp[0].X;
        physical&=std::abs(model.comp[i].sum()-1)<1e-14 && model.comp[i].X[0]<.7 && model.comp[i].X[1]>0
          && model.y[i].L>0 && eos.eval(model.T(i),model.rho(i),model.comp[i]).cv>0;
      }
    }
    final.push_back(model);return true;
  };
  for(int n:{1,2,4}) if(!advance(n))return 1;
  check(balance<1e-8,"coupled surface luminosity equals burning plus caloric/contraction heating",balance);
  check(restmass<1e-7,"integrated abundance change accounts for nuclear heat and neutrinos",restmass);
  check(abundance<1e-12 && nonlinear<1e-10,"composition and structure both converge at the accepted state",abundance);
  check(mixed && physical,"Schwarzschild regions recover a fully mixed, normalized, physical star");
  const double coarse=std::abs(final[0].r(511)-final[1].r(511));
  const double fine=std::abs(final[1].r(511)-final[2].r(511));
  check(coarse>0 && fine/coarse<.8,"halving backward-Euler timestep reduces radius evolution error",fine/coarse);
  check(std::abs(final.back().age/duration-1)<1e-14 && initial.age==0 && initial.comp[0].X[1]==0,
        "accepted steps advance age without mutating the previous model");
  options.max_abundance_change=1e-10;
  const auto rejected=evolve_step(initial,physics,atmosphere,duration,options);
  check(!rejected.converged && rejected.model.age==initial.age && rejected.model.comp[0].X==initial.comp[0].X
    && rejected.model.y[0].lnT==initial.y[0].lnT,"rejected step returns the original composition, structure and age");
  auto edge=initial;
  for(auto& c:edge.comp) {c.X[0]=.695001;c.X[1]=.004999;}
  options.max_abundance_change=.001;
  const auto limited=evolve_step(edge,physics,atmosphere,duration,options);
  check(!limited.converged && limited.message.find("FrozenCompositionAtmosphere")!=std::string::npos
    && limited.model.age==edge.age && limited.model.comp[0].X==edge.comp[0].X
    && limited.model.y[0].lnT==edge.y[0].lnT,"atmosphere composition coverage failure preserves the last supported model");
  return failures?1:0;
}
