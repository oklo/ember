#include "ember/atmosphere_helium_blend.hpp"
#include <algorithm>
#include <iostream>

using namespace ember;
namespace {
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
bool close(double a, double b, double tol=1e-8) {
  return std::abs(a-b) <= tol*std::max({1., std::abs(a), std::abs(b)});
}
class Gas final : public Eos {
public:
  EosState eval(double t, double rho, const Composition& c) const override {
    const double pg=constants::R_gas*c.mu_ions_inv()*rho*t;
    const double pr=constants::a_rad*std::pow(t,4)/3;
    EosState s{};s.P=pg+pr;s.chiT=(pg+4*pr)/s.P;s.chiRho=pg/s.P;return s;
  }
  const char* name() const override { return "analytic test gas"; }
};
Composition comp(double helium, double f3=.3) {
  auto c=solar_scaled(1-helium,0);c.X[1]=helium*f3;c.X[2]=helium-c.X[1];
  c.basis=AbundanceBasis::baryon_mass;return c;
}
class Source final : public Atmosphere {
public:
  const Eos& eos; bool high; mutable int calls=0; bool reject=false; double depth=100;
  Source(const Eos& e,bool h):eos(e),high(h){}
  AtmosphereState eval(double t,double g,const Composition& c) const override {
    ++calls;if(reject)throw std::domain_error("source deliberately missing");
    const double h=c.X[1]+c.X[2], k=high?1.:0.;
    AtmosphereState s{};
    s.T=6000*std::pow(t/4600.,.8+.3*k)*std::pow(g/4e5,.03+.02*k)*std::exp(k*(.1+2*h+c.X[1]));
    s.Pgas=1e7*std::pow(t/4600.,-.2+.1*k)*std::pow(g/4e5,.9-.2*k)*std::exp(k*(-.2-h));
    const double pr=constants::a_rad*std::pow(s.T,4)/3;
    s.P=s.Pgas+pr;s.tau=depth;
    s.dlnT_dlnTeff=.8+.3*k;s.dlnT_dlng=.03+.02*k;
    s.dlnP_dlnTeff=(s.Pgas*(-.2+.1*k)+4*pr*s.dlnT_dlnTeff)/s.P;
    s.dlnP_dlng=(s.Pgas*(.9-.2*k)+4*pr*s.dlnT_dlng)/s.P;
    s.rho=eos.rho_from_PT(s.T,s.P,c);return s;
  }
  const char* name() const override { return "analytic source"; }
};
template<class F>bool rejects(F f){try{f();return false;}catch(const std::exception&){return true;}}
}
int main(){try{
  Gas eos;Source a(eos,false),b(eos,true);
  HeliumBlendAtmosphere grid(eos,a,b,.0045,.005);
  const double t=4600,g=4e5;
  for(double h:{.004,.0045,.005,.006}){
    auto c=comp(h);auto expected=(h<=.0045?a:b).eval(t,g,c);
    a.calls=b.calls=0;auto s=grid.eval(t,g,c);
    require(s.T==expected.T&&s.P==expected.P&&s.rho==expected.rho,"outside blend changed source");
    require((h<=.0045?b.calls:a.calls)==0,"unused source evaluated outside overlap");
  }
  auto c=comp(.00475);auto s=grid.eval(t,g,c);auto sa=a.eval(t,g,c),sb=b.eval(t,g,c);
  require(close(s.T,std::sqrt(sa.T*sb.T),1e-13)&&close(s.Pgas,std::sqrt(sa.Pgas*sb.Pgas),1e-13),"midpoint is not logarithmic mean");
  require(close(eos.eval(s.T,s.rho,c).P,s.P,1e-12),"density not inverted at actual composition");
  for(int axis=0;axis<2;++axis){
    double h=1e-5;auto lo=grid.eval(t*std::exp(axis==0?-h:0),g*std::exp(axis==1?-h:0),c);
    auto hi=grid.eval(t*std::exp(axis==0?h:0),g*std::exp(axis==1?h:0),c);
    require(close(std::log(hi.T/lo.T)/(2*h),axis==0?s.dlnT_dlnTeff:s.dlnT_dlng),"thermal derivative mismatch");
    require(close(std::log(hi.P/lo.P)/(2*h),axis==0?s.dlnP_dlnTeff:s.dlnP_dlng),"pressure derivative mismatch");
  }
  // The composition slope must approach the original source on both sides.
  for(double edge:{.0045,.005}){
    const double h=1e-8;
    auto mid=grid.eval(t,g,comp(edge)),lo=grid.eval(t,g,comp(edge-h)),hi=grid.eval(t,g,comp(edge+h));
    require(std::abs(std::log(mid.T/lo.T)/h-std::log(hi.T/mid.T)/h)<2e-5,"composition slope discontinuity in temperature");
    require(std::abs(std::log(mid.P/lo.P)/h-std::log(hi.P/mid.P)/h)<2e-5,"composition slope discontinuity in pressure");
  }
  b.reject=true;
  require(rejects([&]{grid.eval(t,g,c);}),"missing overlap source silently replaced");
  require(!rejects([&]{grid.eval(t,g,comp(.004));}),"inactive source incorrectly required");
  b.reject=false;b.depth=10;
  require(rejects([&]{grid.eval(t,g,c);}),"inconsistent optical depths accepted");
  require(rejects([&]{HeliumBlendAtmosphere bad(eos,a,b,.005,.0045);}),"reversed interval accepted");
  std::cout<<"Helium blend: source preservation, actual EOS, thermal derivatives, continuous composition slopes and domain rejection passed.\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
