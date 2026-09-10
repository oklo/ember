#include "ember/convection.hpp"
#include "ember/eos_composite.hpp"
#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* label,double error=0) {
  failures+=!ok;std::printf("[%s] %s (%.8g)\n",ok?"PASS":"FAIL",label,error);
}
class Gas final:public Eos {
public:
  EosState eval(double T,double rho,const Composition& c) const override {
    const double R=constants::R_gas*(c.mu_ions_inv()+c.mu_elec_inv());
    EosState e{};e.P=R*rho*T;e.E=1.5*R*T;e.cv=1.5*R;e.cp=2.5*R;
    e.delta=e.chiT=e.chiRho=1;e.grad_ad=.4;return e;
  }
  EosResponse eval_with_derivatives(double T,double rho,const Composition& c) const override {
    EosResponse r{};r.state=eval(T,rho,c);return r;
  }
  const char* name() const override {return "test ideal mixture";}
};
class ConstantOpacity final:public Opacity {
public:
  OpacityState eval(double,double,const Composition&) const override {return {1,0,0};}
  const char* name() const override {return "test opacity";}
};
double jacobian_error(const Model& m,const Physics& p) {
  const auto a=zone_residual(m,0,p,0),b=zone_residual_numerical(m,0,p,0,nullptr,1e-5);
  double worst=0;
  for(std::size_t k=0;k<NVAR;++k) {
    double scale=0,error=0;
    for(int e=0;e<2;++e)for(std::size_t j=0;j<NVAR;++j) {
      const double unit=j==3?std::max(std::abs(m.y[0].L),1.):1.;
      const double v=(e?a.dfdy_hi:a.dfdy_lo)[k][j]*unit,w=(e?b.dfdy_hi:b.dfdy_lo)[k][j]*unit;
      scale=std::max({scale,std::abs(v),std::abs(w)});error=std::max(error,std::abs(v-w));
    }
    worst=std::max(worst,error/std::max(scale,1e-300));
  }
  return worst;
}
}
int main() {
  Gas gas;CompositeEos composite;ConstantOpacity opacity;PPChains nuclear;
  auto inner=solar_scaled(.3,.02),outer=solar_scaled(.7,.02);
  inner.basis=outer.basis=AbundanceBasis::baryon_mass;
  const double T=5e6,P=1e17,contrast=-.25;
  const auto b=composition_buoyancy(gas,T,P,1,contrast,inner,outer);
  const double expected=std::log((inner.mu_ions_inv()+inner.mu_elec_inv())/
      (outer.mu_ions_inv()+outer.mu_elec_inv()))/contrast;
  check(std::abs(b.B-expected)<1e-13 && b.B>0,"helium-rich center gives the ideal-gas Ledoux mu gradient",b.B-expected);
  check(composition_buoyancy(gas,T,P,1,0,inner,inner).B==0,"homogeneous composition gives exactly zero buoyancy");
  bool threw=false;try{composition_buoyancy(gas,T,P,1,0,inner,outer);}catch(const std::domain_error&){threw=true;}
  check(threw,"unresolved pressure contrast is rejected");
  const auto stable=ledoux_mixing_length_gradient(.6,.4,.3,1e-5);
  check(!stable.unstable && stable.grad==.6,"stabilizing composition barrier keeps diffusive transport");
  const auto mixed=ledoux_mixing_length_gradient(.8,.4,.1,1e-7);
  check(mixed.unstable && mixed.grad>=.5 && mixed.grad<.501,"efficient Ledoux convection approaches buoyancy-neutral gradient");
  const auto inverted=ledoux_mixing_length_gradient(.3,.4,-.2,.1);
  check(inverted.unstable && inverted.grad>=.2 && inverted.grad<.3,"inverted composition can destabilize Schwarzschild-stable layer");
  double derivative=0;
  const auto response=composition_buoyancy(composite,T,P,.6,contrast,inner,outer);
  constexpr double h=1e-5;
  for(int k=0;k<4;++k) {
    const auto plus=composition_buoyancy(composite,T*std::exp(k==0?h:0),P*std::exp(k==1?h:0),
        .6+(k==2?h:0),contrast+(k==3?h:0),inner,outer);
    const auto minus=composition_buoyancy(composite,T*std::exp(k==0?-h:0),P*std::exp(k==1?-h:0),
        .6-(k==2?h:0),contrast-(k==3?h:0),inner,outer);
    const std::array<double,4> d{response.dB_dlnT,response.dB_dlnP,response.dB_ddelta,response.dB_dpressure_contrast};
    derivative=std::max(derivative,std::abs((plus.B-minus.B)/(2*h)-d[k])/std::max(1.,std::abs(d[k])));
  }
  check(derivative<2e-7,"nonideal EOS buoyancy derivatives match independent perturbations",derivative);
  Model m;m.M=.1*constants::Msun;m.m={.4*m.M,m.M};m.comp={inner,outer};
  m.y={{std::log(3e9),std::log(gas.rho_from_PT(T,P,inner)),std::log(T),0},
       {std::log(4e9),std::log(gas.rho_from_PT(.98*T,P*std::exp(contrast),outer)),std::log(.98*T),0}};
  Physics physics{&gas,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux};
  auto luminosity=[&](double rad) {
    const double Tb=.5*(m.T(0)+m.T(1)),Pb=.5*(gas.eval(m.T(0),m.rho(0),m.comp[0]).P+gas.eval(m.T(1),m.rho(1),m.comp[1]).P);
    const double L=rad*16*M_PI*constants::a_rad*constants::c*constants::G*.5*(m.m[0]+m.m[1])*std::pow(Tb,4)/(3*Pb);
    m.y[0].L=m.y[1].L=L;
  };
  luminosity(.6);
  check(schwarzschild_mixing_regions(m,physics).size()==1 && convective_mixing_regions(m,physics).size()==2,
        "mixing partition honors Ledoux barrier while Schwarzschild control remains explicit");
  check(jacobian_error(m,physics)<2e-6,"stable stratified structure Jacobian",jacobian_error(m,physics));
  luminosity(.6+b.B);
  check(convective_mixing_regions(m,physics).size()==1,"supercritical flux connects convective mixing region");
  check(jacobian_error(m,physics)<2e-6,"unstable stratified structure Jacobian",jacobian_error(m,physics));
  m.comp[1]=inner;luminosity(.8);
  auto schwarzschild=physics;schwarzschild.criterion=ConvectiveCriterion::schwarzschild;
  check(zone_equations(m,0,physics,0)==zone_equations(m,0,schwarzschild,0),"homogeneous Ledoux exactly preserves previous transport equations");
  physics.eos=&composite;m.comp[1]=outer;
  check(jacobian_error(m,physics)<3e-6,"nonideal stratified structure Jacobian",jacobian_error(m,physics));
  return failures?1:0;
}
