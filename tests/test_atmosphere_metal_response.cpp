#include "ember/atmosphere_metal_response.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace ember;
namespace {
int failures{};
void check(bool value, const char *message) {
  std::cout << (value ? "PASS " : "FAIL ") << message << '\n';
  failures += !value;
}
bool near(double a, double b, double tol=1e-9) {
  return std::abs(a-b) <= tol*std::max({std::abs(a), std::abs(b), 1.});
}
template<class F> bool throws(F f) { try { f(); } catch(const std::exception &) { return true; } return false; }
Composition composition(double z, double y=.001) {
  auto c=solar_scaled(.62,z); c.X[1]=y; c.X[2]-=y;
  c.basis=AbundanceBasis::baryon_mass; c.metal_inventory=MetalInventory::gs98; return c;
}
class Gas final : public Eos {
public:
  EosState eval(double t,double rho,const Composition &c) const override {
    EosState s{};
    const double pg=constants::R_gas*c.mu_ions_inv()*rho*t;
    const double pr=constants::a_rad*std::pow(t,4)/3;
    s.P=pg+pr; s.chiRho=pg/s.P; s.chiT=(pg+4*pr)/s.P; return s;
  }
  const char* name() const override { return "test ideal ions plus radiation"; }
};
class Reference final : public Atmosphere {
  const Eos &eos_;
public:
  explicit Reference(const Eos &eos):eos_(eos){}
  AtmosphereState eval(double t,double g,const Composition &c)const override {
    if(std::abs(c.Z()-.02)>1e-12)throw std::invalid_argument("reference composition changed");
    AtmosphereState s{};
    s.T=t*1.4*std::pow(g/1e5,.04)*std::exp(.1*c.X[0]+2*c.X[1]);
    // Radiation is appreciable in this fixture, exposing total/gas pressure mistakes.
    s.Pgas=20*std::pow(t/3500.,-1.7)*std::pow(g/1e5,.8)*std::exp(-.2*c.X[0]+3*c.X[1]);
    const double pr=constants::a_rad*std::pow(s.T,4)/3;
    s.P=s.Pgas+pr; s.tau=100;
    s.dlnT_dlnTeff=1; s.dlnT_dlng=.04;
    s.dlnP_dlnTeff=(-1.7*s.Pgas+4*pr)/s.P;
    s.dlnP_dlng=(.8*s.Pgas+.16*pr)/s.P;
    s.rho=eos_.rho_from_PT(s.T,s.P,c); return s;
  }
  const char* name()const override{return "analytic reference atmosphere";}
};
double dt(double x,double t,double g){return -.1+.02*x+.01*t-.005*g+.004*x*t*g;}
double dp(double x,double t,double g){return .4-.05*x-.02*t+.015*g-.008*x*t*g;}
std::string fixture(int missing=-1){
  std::ostringstream out;out<<std::setprecision(17)
    <<"EMBER_METAL_ATMOSPHERE_RESPONSE 1\nsource \"synthetic response\"\n"
    <<"approximation \"test separable response\"\nbasis baryon_mass\n"
    <<"reference_Z .02\nsource_Z .01\nmaximum_helium3 .003\ntau 100\n"
    <<"hydrogen 3 .3 .6 .9\nlog_teff 2 3.5 3.7\nlog_g 2 4.9 5.4\ndata\n";
  int i=0;
  for(double x:{.3,.6,.9})for(double t:{3.5,3.7})for(double g:{4.9,5.4}){
    if(i++==missing)out<<"0\n";
    else out<<"1 "<<dt(x,t,g)<<' '<<dp(x,t,g)<<'\n';
  }
  return out.str();
}
}
int main(){
  Gas eos;Reference reference(eos);std::istringstream input(fixture());
  const auto approximation=MetalResponseAtmosphere::Approximation::separable_gs98_response;
  MetalResponseAtmosphere atmosphere(eos,reference,input,approximation,{.00001});
  const double t=3700,g=1.2e5;
  for(double z:{.01,.015,.02,.020006}){
    auto c=composition(z);const auto s=atmosphere.eval(t,g,c);
    const auto base=reference.eval(t,g,composition(.02));
    const double f=(c.Z()-.02)/(-.01);
    check(near(s.T,base.T*std::exp(f*dt(c.X[0],std::log10(t),std::log10(g)))),"source response and bounded metal continuation");
    check(near(s.Pgas,base.Pgas*std::exp(f*dp(c.X[0],std::log10(t),std::log10(g)))),"gas pressure response is applied before adding radiation");
    const double expected=s.Pgas/(constants::R_gas*c.mu_ions_inv()*s.T);
    check(std::abs(s.rho/expected-1)<1e-8,"density uses the actual metal and isotope abundances");
    const double h=1e-5;
    const auto tp=atmosphere.eval(t*std::exp(h),g,c),tm=atmosphere.eval(t*std::exp(-h),g,c);
    const auto gp=atmosphere.eval(t,g*std::exp(h),c),gm=atmosphere.eval(t,g*std::exp(-h),c);
    check(near(std::log(tp.T/tm.T)/(2*h),s.dlnT_dlnTeff,1e-8)&&near(std::log(gp.T/gm.T)/(2*h),s.dlnT_dlng,1e-8),"temperature derivatives match the combined atmosphere");
    check(near(std::log(tp.P/tm.P)/(2*h),s.dlnP_dlnTeff,1e-8)&&near(std::log(gp.P/gm.P)/(2*h),s.dlnP_dlng,1e-8),"total-pressure derivatives include radiation response");
  }
  check(throws([&]{atmosphere.eval(t,g,composition(.009));}),"metal depletion beyond source coverage rejected");
  check(throws([&]{atmosphere.eval(t,g,composition(.02002));}),"unrequested metal continuation rejected");
  check(throws([&]{atmosphere.eval(t,g,composition(.015,.004));}),"large isotope correction rejected");
  check(throws([&]{atmosphere.eval(3000,g,composition(.015));}),"temperature outside coverage rejected");
  std::istringstream incomplete(fixture(8));
  MetalResponseAtmosphere masked(eos,reference,incomplete,approximation,{});
  check(throws([&]{masked.eval(t,g,composition(.015));}),"missing derivative corner rejected");
  auto edge=composition(.015);edge.X[0]=.6;edge.X[2]+=.02;
  check(!throws([&]{masked.eval(t,g,edge);}),"complete one-sided cell retained at an exact knot");
  return failures ? 1 : 0;
}
