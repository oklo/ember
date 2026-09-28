#include "ember/opacity_conductive_interior.hpp"
#include <iostream>

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
  rejected=false;try{continued.eval(8e5,3e4,c);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"missing temperature support was extrapolated");
  std::cout<<"radiative continuation, derivatives, support and contribution checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
