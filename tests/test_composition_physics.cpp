#include "ember/eos_composition.hpp"
#include "ember/evolution_proxies.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace ember;
static int failures=0;
static void check(bool ok,const char* label,double value=0) {
  failures+=!ok;std::printf("[%s] %s (%.9g)\n",ok?"PASS":"FAIL",label,value);
}
template<class F> bool throws(F f) {try{f();}catch(const std::exception&){return true;}return false;}
int main() {
  const std::string data=EMBER_DATA_DIR;
  CompositionHelmholtzEos eos(data+"/eos/freeeos300_hhe_composition.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
  HelmholtzTableEos fixed(data+"/eos/freeeos300_hhe_x070_potential.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
  auto comp=solar_scaled(.7,.02);
  check(std::abs(eos.eval(1e6,1,comp).E/fixed.eval(1e6,1,comp).E-1)<1e-13,
        "atomic He3-free family reproduces the existing source plane");
  comp.basis=AbundanceBasis::baryon_mass;
  check(throws([&]{fixed.eval(1e6,1,comp);}),"fixed EOS rejects a silent abundance-basis change");
  comp.X[0]=.675;comp.X[1]=.002;comp.X[2]=1-comp.Z()-comp.X[0]-comp.X[1];
  double derivative=0,firstlaw=0,inversion=0;
  constexpr double step=1e-6;
  for(auto [T,rho]:{std::pair{5500.,1e-4},std::pair{30000.,.03},std::pair{2e6,10.},std::pair{4.6e6,350.}}) {
    const auto e=eos.eval(T,rho,comp);const auto response=eos.composition_response(T,rho,comp);
    inversion=std::max(inversion,std::abs(eos.rho_from_PT(T,e.P,comp)/rho-1));
    for(int j=0;j<2;++j) {
      auto plus=comp,minus=comp;plus.X[j]+=step;plus.X[2]-=step;minus.X[j]-=step;minus.X[2]+=step;
      const auto p=eos.eval(T,rho,plus),m=eos.eval(T,rho,minus);
      derivative=std::max({derivative,std::abs((p.P-m.P)/(2*step)-response.dP[j])/e.P,
        std::abs((p.E-m.E)/(2*step)-response.dE[j])/(T*e.cv)});
    }
    const auto p=eos.eval(T,rho*std::exp(step),comp),m=eos.eval(T,rho*std::exp(-step),comp);
    firstlaw=std::max({firstlaw,std::abs((p.E-m.E)/(2*step)-e.P/rho*(1-e.chiT))/(T*e.cv),
      std::abs((p.S-m.S)/(2*step)+e.P*e.chiT/(rho*T))/e.cv,
      std::abs(e.P*e.delta/(rho*T*e.cp*e.grad_ad)-1)});
  }
  check(derivative<1e-5,"H1/He3 composition derivatives follow actual pressure and energy",derivative);
  check(firstlaw<1e-5,"baryonic isotope EOS preserves first-law and Maxwell identities",firstlaw);
  check(inversion<1e-9,"composition-dependent PT inversion recovers density",inversion);
  auto outside=solar_scaled(.58,.02);outside.basis=AbundanceBasis::baryon_mass;
  check(throws([&]{eos.eval(1e6,1,outside);}),"EOS rejects unsupported hydrogen instead of extrapolating");
  std::ifstream in(std::string(EMBER_TEST_DATA_DIR)+"/freeeos300_composition_reference.dat");
  double X,Y3,T,rho,P,E,cv,cp,ad,values=0,responses=0;int count=0;
  while(in>>X>>Y3>>T>>rho>>P>>E>>cv>>cp>>ad) {
    comp=solar_scaled(X,.02);comp.basis=AbundanceBasis::baryon_mass;comp.X[1]=Y3;comp.X[2]-=Y3;
    const auto e=eos.eval(T,rho,comp);++count;
    values=std::max({values,std::abs(e.P/P-1),std::abs(e.E/E-1)});
    responses=std::max({responses,std::abs(e.cv/cv-1),std::abs(e.cp/cp-1),std::abs(e.grad_ad/ad-1)});
  }
  check(count==28 && values<.001,"independent FreeEOS number-density queries agree in P/E within 0.1%",values);
  check(responses<.02,"off-composition source response audit agrees within 2%",responses);

  AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat");
  TopsOpacity high(data+"/opacity",TopsOpacity::Grid::composition);
  std::ifstream opacity_source(std::string(EMBER_TEST_DATA_DIR)+"/tops_composition_reference.dat");
  double source_kappa=0,source_error=0;int opacity_count=0;
  while(opacity_source>>X>>T>>rho>>source_kappa) {
    ++opacity_count;const auto c=solar_scaled(X,.02);
    source_error=std::max(source_error,std::abs(high.eval(T,rho,c).kappa/source_kappa-1));
  }
  check(opacity_count==60 && source_error<1e-12,"all four opacity planes retain independent original LANL cells",source_error);
  BlendedOpacity blend(low,high,4.4,4.5);NominalAbundanceOpacity opacity(blend);
  comp=solar_scaled(.675,.02);comp.basis=AbundanceBasis::baryon_mass;comp.X[1]=.001;comp.X[2]-=.001;
  double kd=0;
  for(auto [lt,lr]:{std::pair{3.8,-4.},std::pair{4.45,-2.},std::pair{5.65,.3},std::pair{6.6,2.4}}) {
    const auto k=opacity.eval(std::pow(10.,lt),std::pow(10.,lr),comp);
    auto p=comp,m=comp;p.X[0]+=step;p.X[2]-=step;m.X[0]-=step;m.X[2]+=step;
    const double d=(std::log(opacity.eval(std::pow(10.,lt),std::pow(10.,lr),p).kappa)
      -std::log(opacity.eval(std::pow(10.,lt),std::pow(10.,lr),m).kappa))/(2*step);
    kd=std::max(kd,std::abs(d-k.dlnk_dX));
  }
  check(kd<1e-7,"opacity composition derivatives include both temperature blends",kd);
  auto invalid=comp;invalid.X[1]=.006;invalid.X[2]-=.005;
  check(throws([&]{opacity.eval(1e6,1,invalid);}) && throws([&]{high.eval(1e6,1,comp);}),
        "opacity requires explicit bounded nominal-abundance approximation");
  TabulatedAtmosphere cond(eos,data+"/atmosphere/cond_gn93_tau100_solar_proxy.dat",TabulatedAtmosphere::Mixture::allow_documented_proxy);
  FrozenCompositionAtmosphere atmosphere(eos,cond);
  check(throws([&]{atmosphere.eval(2800,1e5,comp);}),"frozen atmosphere rejects significant surface composition evolution");
  comp=solar_scaled(.699,.02);comp.basis=AbundanceBasis::baryon_mass;comp.X[1]=.001;comp.X[2]-=.001;
  const auto atm=atmosphere.eval(2800,1e5,comp);
  check(std::abs(eos.eval(atm.T,atm.rho,comp).P/atm.P-1)<1e-9,
        "frozen atmosphere density uses the actual evolving EOS");
  return failures?1:0;
}
