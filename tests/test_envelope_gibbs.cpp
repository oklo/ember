#include "ember/envelope_gibbs.hpp"
#include "ember/constants.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>

using namespace ember;
namespace {
void table(const std::filesystem::path& path,double lower,double upper,double particles,
           bool hole=false,double shift=0) {
  std::ofstream out(path);out<<std::setprecision(17)<<"EMBER_GIBBS_GAS_V1 3\n";
  const double p0=std::log(1e6),p1=std::log(1e18),t0=std::log(lower),t1=std::log(upper);
  auto axis=[&](double a,double b){out<<"8 ";for(int i=0;i<4;++i)out<<a<<' ';for(int i=0;i<4;++i)out<<b<<' ';out<<'\n';};
  axis(t0,t1);axis(p0,p1);out<<"16 ";
  // Greville coefficients reproduce an ideal monatomic gas exactly.
  for(int i=0;i<4;++i)for(int j=0;j<4;++j)
    out<<particles*((p0+(p1-p0)*j/3.)-2.5*(t0+(t1-t0)*i/3.))+shift<<' ';
  out<<"\n2 "<<t0<<' '<<t1<<"\n5 ";
  for(int j=0;j<5;++j)out<<p0+(p1-p0)*j/4.<<' ';out<<'\n';
  out<<(hole?"1 0 1 1 1 1 1 1 1 1\n":"1 1 1 1 1 1 1 1 1 1\n");
}
struct Top final : Atmosphere {
  AtmosphereState eval(double,double,const Composition&) const override {
    AtmosphereState s{};s.T=3e4;s.P=1e8;s.rho=.01;return s;
  }
  const char* name() const override{return "fixed boundary";}
};
struct Transparent final : Opacity {
  OpacityState eval(double,double,const Composition&) const override {OpacityState s{};s.kappa=1e-20;return s;}
  const char* name() const override{return "transparent layer";}
};
struct NegativeExpansion final : EnvelopeThermodynamics {
  State at_pressure(double,double,const Composition&,double& guess) const override {
    guess=std::log(.01);return {.01,1e8,-.2,-.1,1.};
  }
  double rho_from_PT(double,double,const Composition&,double=0) const override{return .01;}
};
}
int main() {
  int failed=0;
  auto check=[&](bool ok,const char* why) {std::cout<<(ok?"PASS ":"FAIL ")<<why<<'\n';failed+=!ok;};
  auto rejects=[&](auto&& f,const char* why) {bool caught=false;try{f();}catch(const std::domain_error&){caught=true;}check(caught,why);};
  const auto dir=std::filesystem::temp_directory_path()/
      ("ember-gibbs-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(dir);
  try {
    GibbsEnvelope::Files files{(dir/"h").string(),(dir/"he").string(),(dir/"hw").string(),(dir/"hew").string()};
    table(files.hydrogen,1000,1e6,1/1.00782503);
    table(files.helium,1000,1e6,1/4.00260325);
    table(files.hydrogen_warm,9000,2e6,1/1.00782503);
    table(files.helium_warm,1e5,2e6,1/4.00260325);
    GibbsEnvelope eos(files);
    auto c=solar_scaled(.99,0);c.basis=AbundanceBasis::baryon_mass;
    c[Species::He3]=.003;c[Species::He4]-=.003;
    const double gas_constant=constants::R_gas*(c.h1()+c[Species::He3]/3.+c[Species::He4]/4.);
    const double T=20000,P=1e10,rad=constants::a_rad*std::pow(T,4)/3.;
    const auto a=eos.evaluate(T,P,c);const double volume=gas_constant*T/(P-rad);
    check(std::abs(a.pressure.rho*volume-1)<2e-12,"atomic masses convert to the baryon convention");
    check(std::abs(a.pressure.delta-(1+4*rad/(P-rad)))<2e-12,"radiation enters thermal expansion once");
    check(std::abs(a.pressure.chiRho-(P-rad)/P)<2e-12,"radiation enters compressibility once");
    check(std::abs(a.energy-(1.5*gas_constant*T+3*rad*volume))/(gas_constant*T)<2e-12,
          "caloric energy includes gas and radiation");
    constexpr double step=1e-4;
    const auto tp=eos.evaluate(T*std::exp(step),P,c),tm=eos.evaluate(T*std::exp(-step),P,c);
    const auto pp=eos.evaluate(T,P*std::exp(step),c),pm=eos.evaluate(T,P*std::exp(-step),c);
    check(std::abs((tp.entropy-tm.entropy)/(2*step)/a.pressure.cp-1)<1e-7,"Cp differentiates entropy at fixed pressure");
    const double entropy_pressure=(pp.entropy-pm.entropy)/(2*step);
    check(std::abs(entropy_pressure/(P*volume/T)+a.pressure.delta)<1e-7,"Maxwell relation matches density and entropy");
    check(std::abs(-entropy_pressure/a.pressure.cp/a.pressure.grad_ad-1)<1e-7,"adiabatic gradient differentiates the same entropy");
    for(double join:{9000.,1e5}) {
      const auto below=eos.evaluate(join*std::exp(-1e-8),P,c).pressure;
      const auto above=eos.evaluate(join*std::exp(1e-8),P,c).pressure;
      check(std::abs(above.cp/below.cp-1)<1e-6,"warm component join preserves heat capacity");
    }
    auto bad=c;bad.basis=AbundanceBasis::atomic_mass;
    rejects([&]{eos.evaluate(T,P,bad);},"atomic mass input is not mistaken for baryon fractions");
    bad=solar_scaled(.7,0);bad.basis=AbundanceBasis::baryon_mass;
    rejects([&]{eos.evaluate(T,P,bad);},"unassessed composition is refused");
    rejects([&]{eos.evaluate(999,P,c);},"outside temperature coverage is refused");
    table(files.hydrogen_warm,9000,2e6,1/1.00782503,false,.01);
    rejects([&]{GibbsEnvelope discontinuous(files);},"discontinuous component joins are refused");
    table(files.hydrogen_warm,9000,2e6,1/1.00782503);
    table(files.hydrogen,1000,1e6,1/1.00782503,true);
    rejects([&]{GibbsEnvelope masked(files);masked.evaluate(4000,P,c);},"missing source cells cannot be interpolated across");
    Top top;Transparent transparent;NegativeExpansion negative;
    const double mass=.1*constants::Msun,radius=.1*constants::Rsun;
    EnvelopeAtmosphere envelope(top,transparent,negative,1.9,mass,.000001*mass,20);
    envelope.integration_tolerance(1e-7);
    const auto surface=envelope.from_photosphere(3000,radius,c);
    check(std::abs(surface.T_base/3e4-1)<1e-7 && surface.r_base<radius,
          "negative expansion with outward radiation integrates without convective transport");
    const double te=3000*std::sqrt(radius/surface.r_base),g=constants::G*mass/std::pow(surface.r_base,2);
    envelope.evaluation_threads(2);const auto boundary=envelope.eval(te,g,c);
    check(std::isfinite(boundary.dlnT_dlnTeff+boundary.dlnP_dlng),"parallel boundary derivatives use the same implicit integration");
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';++failed;}
  std::filesystem::remove_all(dir);return failed?1:0;
}
