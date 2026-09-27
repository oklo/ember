#include "ember/atmosphere_fixed_metal.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace ember;
namespace {
class Gas final:public Eos {
public:
  mutable Composition last;
  EosState eval(double t,double rho,const Composition& c)const override {
    last=c;EosState s;const double pg=constants::R_gas*c.mu_ions_inv()*rho*t;
    const double pr=constants::a_rad*std::pow(t,4)/3;
    s.P=pg+pr;s.chiRho=pg/s.P;s.chiT=(pg+4*pr)/s.P;return s;
  }
  const char* name()const override{return "test ideal ions plus radiation";}
};
Composition composition(double z) {
  auto c=solar_scaled(.62,z);c.X[1]=.002;c.X[2]-=.002;
  c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.cn_molality=initial_gs98_cn(c);return explicit_cn_material(c);
}
std::string fixture(int missing=-1) {
  std::ostringstream o;o<<std::setprecision(17)
    <<"EMBER_COMPOSITION_ATMOSPHERE 2\nsource \"synthetic test\"\n"
    <<"approximation \"test\"\nbasis baryon_mass\ntau 100\nmetals";
  const auto c=composition(.02);
  for(std::size_t j=METAL_BEGIN;j<METAL_END;++j)o<<' '<<c.X[j];
  o<<"\nhydrogen 2 .6 .7\nhelium3 2 0 .01\nlog_teff 2 3.4 3.6\nlog_g 2 4.5 5.5\ndata\n";
  int i=0;
  for(double h:{.6,.7})for(double y:{0.,.01})for(double t:{3.4,3.6})for(double g:{4.5,5.5}) {
    if(i++==missing)o<<"0\n";
    else o<<"1 "<<t+.15+.2*h-y+.02*g<<' '<<3+.2*t+.5*g-.3*h+y<<'\n';
  }
  return o.str();
}
template<class F>bool rejects(F f){try{f();}catch(const std::exception&){return true;}return false;}
}
int main() {
  int failed=0;
  const auto check=[&](bool ok,const char* label){std::cout<<(ok?"PASS ":"FAIL ")<<label<<'\n';failed+=!ok;};
  try {
    Gas eos;std::istringstream input(fixture());
    const auto proxy=CompositionAtmosphereGrid::Mixture::allow_documented_proxy;
    CompositionAtmosphereGrid source(eos,input,proxy);
    FixedMetalAtmosphere boundary(eos,source,1e-8);
    const double t=2900,g=1.5e5;const auto ref=composition(.02);
    const auto a=source.eval(t,g,ref),b=boundary.eval(t,g,ref);
    check(a.T==b.T && a.P==b.P && a.rho==b.rho,"reference composition remains exact");
    for(double dz:{-9e-9,2.2e-12,9e-9}) {
      auto c=composition(.02+dz);const auto unchanged=c;
      const auto s=boundary.eval(t,g,c);
      check(c==unchanged && eos.last==c,"actual CN ledger and composition reach EOS unchanged");
      check(s.T==a.T && s.P==a.P && s.Pgas==a.Pgas,"only boundary T/P use source metal abundance");
      check(std::abs(eos.eval(s.T,s.rho,c).P/s.P-1)<2e-10,"actual-composition pressure inversion");
      check(!source.covers(t,g,c),"general source reader retains strict metal guard");
      constexpr double h=1e-6;
      const auto p=boundary.eval(t*std::exp(h),g,c),m=boundary.eval(t*std::exp(-h),g,c);
      check(std::abs(std::log(p.T/m.T)/(2*h)-s.dlnT_dlnTeff)<2e-8 &&
          std::abs(std::log(p.P/m.P)/(2*h)-s.dlnP_dlnTeff)<2e-8,"thermal derivatives match selected boundary");
    }
    auto c=composition(.02+5e-9);c.X[0]-=2e-5;c[Species::H2]=2e-5;
    const auto before=c;TraceDeuteriumAtmosphere with_D(eos,boundary);
    const auto d=with_D.eval(t,g,c);
    check(c==before && eos.last==c && std::abs(eos.eval(d.T,d.rho,c).P/d.P-1)<2e-10,
        "trace-D and metal approximations compose without changing isotope inventory");
    check(rejects([&]{boundary.eval(t,g,c);}),"unmapped D rejected");
    for(double dz:{-1.01e-8,1.01e-8})check(rejects([&]{boundary.eval(t,g,composition(.02+dz));}),"approximation bound enforced");
    c=ref;c.X[3]+=1e-10;c.X[4]-=1e-10;
    check(rejects([&]{boundary.eval(t,g,c);}),"changed metal pattern rejected even at unchanged Z");
    c=ref;(*c.cn_molality)[0]=1;
    check(rejects([&]{boundary.eval(t,g,c);}),"invalid nuclear ledger rejected");
    check(rejects([&]{boundary.eval(2000,g,ref);}),"source temperature bounds unchanged");
    c=ref;c.X[0]=.59;c.X[2]+=.03;
    check(rejects([&]{boundary.eval(t,g,c);}),"source hydrogen bounds unchanged");
    std::istringstream incomplete(fixture(0));CompositionAtmosphereGrid masked(eos,incomplete,proxy);
    FixedMetalAtmosphere missing(eos,masked,1e-8);
    check(rejects([&]{missing.eval(t,g,ref);}),"incomplete source stencil rejected");
    FixedMetalAtmosphere measured(eos,source,2e-5);
    c=composition(.020019);const auto measured_state=measured.eval(t,g,c);
    check(eos.last==c && measured_state.T==a.T && measured_state.P==a.P,
        "explicit measured small-metal allowance preserves actual composition");
    check(rejects([&]{measured.eval(t,g,composition(.020021));}),"measured allowance remains bounded");
    check(rejects([&]{FixedMetalAtmosphere invalid(eos,source,2.1e-4);}),"large fixed-metal approximation rejected");
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
  return failed?1:0;
}
