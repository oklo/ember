#include "ember/opacity_conductive_interior.hpp"
#include <iostream>
#include <limits>

using namespace ember;
namespace {
struct Source final : Opacity {
  OpacityState eval(double T,double rho,const Composition& c) const override {
    if(T<9e5 || T>1e7 || rho>1e4 || c.h1()<0)
      throw std::domain_error("test source support");
    const double t=std::log(T/3e6),r=std::log(rho/6000.);
    return {2e4*std::exp(-2*t+.7*r+.02*t*r+.04*r*r+(.3+.02*r)*c.h1()+(.4+.03*r)*c.X[1]),
      -2+.02*r,.7+.02*t+.08*r+.02*c.h1()+.03*c.X[1],.3+.02*r,0,.4+.03*r};
  }
  const char* name() const override {return "analytic curved opacity";}
  std::optional<DensityRange> density_range(double,const Composition&) const override {return DensityRange{1,1e4};}
};
struct Heat final : Conduction {
  double opacity;
  explicit Heat(double k):opacity(k){}
  OpacityState eval(double,double,const Composition&) const override {return {opacity,0,0};}
  const char* name() const override {return "constant conductivity equivalent opacity";}
};
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
}
int main(){try{
  Source source;Heat heat(.01);ConductiveInteriorOpacity continued(source,heat);
  Composition c;c.X[0]=.01;c.X[1]=.005;c.X[2]=.985;
  require(continued.eval(3e6,5000,c).kappa==source.eval(3e6,5000,c).kappa,"ordinary source query changed");
  require(continued.eval(6e6,9000,c).kappa==source.eval(6e6,9000,c).kappa,"warm source query changed");
  for(double T:{2e6,3.4e6,3.5e6,3.6e6})for(double rho:{8499.,8500.,9000.,9500.,9800.}) {
    const auto a=continued.eval(T,rho,c);constexpr double h=1e-6;
    const auto logk=[&](double t,double r,const Composition& cc){return std::log(continued.eval(t,r,cc).kappa);};
    require(std::abs((logk(T*std::exp(h),rho,c)-logk(T*std::exp(-h),rho,c))/(2*h)-a.dlnk_dlnT)<2e-6,"temperature Jacobian mismatch");
    require(std::abs((logk(T,rho*std::exp(h),c)-logk(T,rho*std::exp(-h),c))/(2*h)-a.dlnk_dlnRho)<2e-6,"density Jacobian mismatch");
    auto p=c,m=c;p.X[0]+=h;p.X[2]-=h;m.X[0]-=h;m.X[2]+=h;
    require(std::abs((logk(T,rho,p)-logk(T,rho,m))/(2*h)-a.dlnk_dX)<2e-6,"composition Jacobian mismatch");
    p=c;m=c;p.X[1]+=h;p.X[2]-=h;m.X[1]-=h;m.X[2]+=h;
    require(std::abs((logk(T,rho,p)-logk(T,rho,m))/(2*h)-a.dlnk_dY3)<2e-6,"helium-3 Jacobian mismatch");
  }
  ConductiveInteriorOpacity::Domain cold_domain;
  cold_domain.cold_source_T=1.2e6;cold_domain.cold_full_T=2.2e6;
  ConductiveInteriorOpacity cold(source,heat,.001,1.,cold_domain);
  for(double T:{1.1e6,1.2e6,1.4e6,2e6,2.2e6})for(double rho:{8499.,8500.,9000.,9500.,9800.}) {
    const auto v=cold.eval(T,rho,c);constexpr double h=1e-6;
    const auto logk=[&](double t,double r){return std::log(cold.eval(t,r,c).kappa);};
    require(std::abs((logk(T*std::exp(h),rho)-logk(T*std::exp(-h),rho))/(2*h)-v.dlnk_dlnT)<2e-6,"cold join T derivative");
    require(std::abs((logk(T,rho*std::exp(h))-logk(T,rho*std::exp(-h)))/(2*h)-v.dlnk_dlnRho)<2e-6,"cold join density derivative");
    if(T<=cold_domain.cold_source_T)
      require(v.kappa==source.eval(T,rho,c).kappa,"cold source values changed");
    if(T>=cold_domain.cold_full_T)
      require(v.kappa==continued.eval(T,rho,c).kappa,"warm continuation changed");
  }
  bool cold_rejected=false;try{cold.eval(1.5e6,3e4,c);}catch(const std::domain_error&){cold_rejected=true;}
  require(cold_rejected,"cold overlap ignored source coverage");
  const auto extended=continued.eval(3e6,3e4,c);
  require(extended.kappa>0 && std::isfinite(extended.dlnk_dlnT),"dense continuation unavailable");
  require(continued.maximum_transport_uncertainty()<.001 && continued.continued_evaluations()>0,"contribution diagnostic missing");
  for(double scale:{.1,10.}) {
    ConductiveInteriorOpacity variation(source,heat,.001,scale);
    require(std::abs(variation.eval(3e6,3e4,c).kappa/extended.kappa-scale)<1e-12,"uncertainty arm incorrectly scaled");
  }
  ConductiveInteriorOpacity nominal(source,heat),lower(source,heat,.001,.1);
  const double kn=nominal.eval(3.5e6,9000,c).kappa,kl=lower.eval(3.5e6,9000,c).kappa;
  const double kc=ConductiveInteriorOpacity::conductivity_margin*heat.opacity;
  const double uncertainty=(1/kl+1/kc)/(1/kn+1/kc)-1;
  require(std::abs(uncertainty-nominal.maximum_transport_uncertainty())<1e-12,
      "admission bound differs from the directly evaluated heat-transport change");
  Heat weak(1e4);ConductiveInteriorOpacity unsafe(source,weak);
  bool rejected=false;try{unsafe.eval(3e6,3e4,c);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"radiatively important extension was accepted");
  Heat disabled(std::numeric_limits<double>::infinity());
  ConductiveInteriorOpacity radiative_only(source,disabled);
  const auto direct=source.eval(3e6,9000.,c),fallback=radiative_only.eval(3e6,9000.,c);
  require(fallback.kappa==direct.kappa && fallback.dlnk_dlnT==direct.dlnk_dlnT
      && fallback.dlnk_dlnRho==direct.dlnk_dlnRho && fallback.dlnk_dX==direct.dlnk_dX
      && fallback.dlnk_dZ==direct.dlnk_dZ && fallback.dlnk_dY3==direct.dlnk_dY3,
      "disabled conduction changed the supported radiative source");
  rejected=false;try{radiative_only.eval(3e6,3e4,c);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"disabled conduction permitted radiative extrapolation");
  require(radiative_only.continued_evaluations()==0,"source fallback counted as continuation");
  for(double invalid:{-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN(),0.,-1.}) {
    Heat invalid_heat(invalid);ConductiveInteriorOpacity invalid_continuation(source,invalid_heat);
    rejected=false;try{invalid_continuation.eval(3e6,9000.,c);}catch(const std::domain_error&){rejected=true;}
    require(rejected,"invalid conduction was silently accepted");
  }
  rejected=false;try{continued.eval(8e5,3e4,c);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"missing temperature support was extrapolated");
  // The independent envelope domain uses the same derivative and heat bound.
  struct EnvelopeSource final : Opacity {
    Source source;
    OpacityState eval(double T,double rho,const Composition& c) const override {
      if(T>=std::pow(10.,5.6) && std::log10(rho)-3*std::log10(T/1e6)>3.5)
        throw std::domain_error("test hot hydrogen source corner");
      return source.eval(T*6,rho*20,c);
    }
    const char* name() const override {return "scaled analytic envelope source";}
  } envelope_source;
  Composition hc;hc.X[0]=.99;hc.X[1]=.003;hc.X[2]=.007;
  const ConductiveInteriorOpacity::Domain domain{180.,190.,3.1e5,5.5e5,6e5,1000.,.97,1e-8,100.};
  ConductiveInteriorOpacity ec(envelope_source,heat,.001,1.,domain);
  for(double T:{3.95e5,4e5,5.5e5,5.7e5,6e5})for(double rho:{179.,180.,185.,190.,200.,210.,300.}) {
    auto v=ec.eval(T,rho,hc);constexpr double h=1e-6;
    auto f=[&](double t,double r){return std::log(ec.eval(t,r,hc).kappa);};
    require(std::abs((f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h)-v.dlnk_dlnT)<2e-6,"envelope T derivative");
    require(std::abs((f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h)-v.dlnk_dlnRho)<2e-6,"envelope density derivative");
  }
  require(ec.eval(5e5,160.,hc).kappa==envelope_source.eval(5e5,160.,hc).kappa,"envelope source range changed");
  require(ec.eval(5e5,300.,c).kappa==envelope_source.eval(5e5,300.,c).kappa,"H-poor state continued");
  rejected=false;try{ec.eval(3e5,300.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"envelope temperature bound ignored");
  ConductiveInteriorOpacity ew(envelope_source,weak,.001,1.,domain);
  rejected=false;try{ew.eval(5e5,300.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"radiatively important envelope accepted");
  // A caller may select a larger local transport bound after a stellar
  // sensitivity check. It changes admission, not the opacity prescription.
  const double kr=ec.eval(5e5,300.,hc).kappa;
  Heat intermediate(.002*kr/(2*(domain.uncertainty-1-.002)));
  ConductiveInteriorOpacity strict(envelope_source,intermediate,.001,1.,domain);
  ConductiveInteriorOpacity assessed(envelope_source,intermediate,.003,1.,domain);
  ConductiveInteriorOpacity low_opacity(envelope_source,intermediate,.003,.01,domain);
  rejected=false;try{strict.eval(5e5,300.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"default envelope uncertainty limit changed");
  require(assessed.eval(5e5,300.,hc).kappa==kr,"larger bound changed nominal opacity");
  const double low=low_opacity.eval(5e5,300.,hc).kappa,kcond=2*intermediate.opacity;
  require(std::abs((1/low+1/kcond)/(1/kr+1/kcond)-1-.002)<1e-12,
      "selected envelope bound differs from independent transport change");
  require(std::abs(assessed.maximum_transport_uncertainty()-.002)<1e-12,
      "larger envelope bound was not recorded");
  struct HotBoundedSource final : Opacity {
    OpacityState eval(double T,double rho,const Composition& mixture) const override {
      if(mixture.h1()>.75 && std::log10(rho)-3*std::log10(T/1e6)>3.5)
        throw std::domain_error("hot source log R edge");
      const double r=std::log(rho/180.);
      return {1e3*std::exp(.7*r+.02*r*r)*std::pow(T/7e5,-.4),-.4,.7+.04*r,0,0,0};
    }
    const char* name() const override {return "source with hot density bound";}
  } bounded;
  Heat strong(1e-6);
  ConductiveInteriorOpacity safe_join(bounded,strong,.001,1.,
      ConductiveInteriorOpacity::hydrogen_envelope_domain());
  for(double T:{5.5e5,6e5,7e5,7.5e5,8e5})for(double rho:{190.,500.,900.,1000.}) {
    const auto v=safe_join.eval(T,rho,hc);constexpr double h=1e-6;
    auto f=[&](double t,double r){return std::log(safe_join.eval(t,r,hc).kappa);};
    require(std::abs((f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h)-v.dlnk_dlnT)<2e-6,
        "upper hydrogen overlap temperature derivative");
    if(rho<1000.)
      require(std::abs((f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h)-v.dlnk_dlnRho)<2e-6,
          "upper hydrogen overlap density derivative");
  }
  for(double X:{.699,.70,.71,.725,.74,.745,.75,.85,.97,.99}) {
    auto mixture=hc;mixture.X[0]=X;mixture.X[2]=1-X-mixture.X[1];
    constexpr double T=5.5e5,rho=900.,h=1e-6;
    const auto v=safe_join.eval(T,rho,mixture);
    auto plus=mixture,minus=mixture;
    plus.X[0]+=h;plus.X[2]-=h;minus.X[0]-=h;minus.X[2]+=h;
    const double fd=(std::log(safe_join.eval(T,rho,plus).kappa)
        -std::log(safe_join.eval(T,rho,minus).kappa))/(2*h);
    require(std::abs(fd-v.dlnk_dX)<2e-6,"hydrogen composition join derivative");
    if(X<=.70)require(v.kappa==bounded.eval(T,rho,mixture).kappa,"H-poor source changed");
  }
  // Following the table's log R edge keeps the overlap inside the source
  // when density grows. The uncertainty check remains local and unchanged.
  ConductiveInteriorOpacity coordinate(bounded,strong,.001,1.,
      ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
  for(double T:{4.1e5,6e5,1e6})for(double logR:{3.0999,3.1,3.1001,3.25,3.3999,3.4,3.4001,4.}) {
    const double rho=std::pow(10.,logR)*std::pow(T/1e6,3);
    const auto v=coordinate.eval(T,rho,hc);constexpr double h=1e-6;
    auto f=[&](double t,double r){return std::log(coordinate.eval(t,r,hc).kappa);};
    require(std::abs((f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h)-v.dlnk_dlnT)<2e-6,
        "source-coordinate temperature derivative");
    require(std::abs((f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h)-v.dlnk_dlnRho)<2e-6,
        "source-coordinate density derivative");
    if(logR<=3.1)require(v.kappa==bounded.eval(T,rho,hc).kappa,"supported low log R source changed");
  }
  require(coordinate.eval(6e5,1002.,hc).kappa>0,"density extension unavailable");
  ConductiveInteriorOpacity coordinate_weak(bounded,weak,.001,1.,
      ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
  rejected=false;try{coordinate_weak.eval(6e5,2000.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"source-coordinate extension ignored radiation uncertainty");
  struct CoolBoundedSource final : Opacity {
    HotBoundedSource source;
    OpacityState eval(double T,double rho,const Composition& c) const override {
      if(T<4e5 && rho>250.)throw std::domain_error("cool source density edge");
      return source.eval(T*4,rho,c);
    }
    const char* name() const override{return "source with cool density bound";}
  } cool_bounded;
  ConductiveInteriorOpacity cool_ionized(cool_bounded,strong,.001,1.,
      ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
  for(double T:{2.001e5,2.5e5,2.999e5,3e5,3.001e5,3.1e5,3.2e5,3.5e5}) {
    constexpr double rho=350.,h=1e-6;
    const auto v=cool_ionized.eval(T,rho,hc);
    auto f=[&](double t,double r){return std::log(cool_ionized.eval(t,r,hc).kappa);};
    require(std::abs((f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h)-v.dlnk_dlnT)<2e-6,
        "ionized cool continuation temperature derivative");
    require(std::abs((f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h)-v.dlnk_dlnRho)<2e-6,
        "ionized cool continuation density derivative");
  }
  rejected=false;try{cool_ionized.eval(199999.,350.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"ionized continuation lost its lower temperature bound");
  ConductiveInteriorOpacity cool_weak(cool_bounded,weak,.001,1.,
      ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
  rejected=false;try{cool_weak.eval(2.5e5,350.,hc);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"cool continuation ignored heat-transport uncertainty");
  // A nearly complete blend must still query its source; rounding the C2
  // polynomial above unity must never skip a source-support check.
  struct NarrowSource final : Opacity {
    HotBoundedSource source;
    OpacityState eval(double T,double rho,const Composition& c) const override {
      if(rho>200.)throw std::domain_error("narrow density support");
      return source.eval(T,rho,c);
    }
    const char* name() const override{return "narrow source";}
  } narrow;
  ConductiveInteriorOpacity edge(narrow,strong,.001,1.,
      ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
  auto edge_c=hc;edge_c.X[0]=.745-1e-6;edge_c.X[2]=1-edge_c.X[0]-edge_c.X[1];
  rejected=false;try{edge.eval(6e5,2000.,edge_c);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"nearly complete composition join bypassed source guard");
  std::cout<<"radiative continuation, derivatives, support and contribution checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
