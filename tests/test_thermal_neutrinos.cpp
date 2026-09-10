#include "ember/structure.hpp"
#include "ember/boundary.hpp"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>

using namespace ember;
class TestOpacity final : public Opacity {
public:
  OpacityState eval(double,double,const Composition&) const override { return {1,0,0}; }
  const char* name() const override { return "test constant"; }
};
class TestNuclear final : public Nuclear {
public:
  NuclearState eval(double,double,const Composition&) const override { NuclearState s;s.eps=.2;return s; }
  const char* name() const override { return "test heating"; }
};
class TestLoss final : public NeutrinoLosses {
public:
  LossState eval(double T,double rho,const Composition&) const override {
    return {.015*std::pow(T/1e7,3.2)*std::pow(rho/1e3,.4),3.2,.4};
  }
  const char* name() const override { return "test analytic cooling"; }
};
int main() {
  int failures=0;
  auto check=[&](bool ok,const char* text,double value=0) {
    failures+=!ok;std::printf("[%s] %s %.12g\n",ok?"PASS":"FAIL",text,value);
  };
  PlasmaNeutrinoLosses plasma;
  std::ifstream input(std::string(EMBER_TEST_DATA_DIR)+"/hrw_plasma_reference.dat");
  std::string line;std::size_t count=0;double worst=0,derivative_error=0;
  while(std::getline(input,line)) {
    if(line.empty() || line[0]=='#')continue;
    double T,rho,ye,reference;std::istringstream row(line);row>>T>>rho>>ye>>reference;
    if(!row)return 1;
    Composition c;c.basis=AbundanceBasis::baryon_mass;c.X[0]=2*ye-1;c.X[2]=2*(1-ye);
    const auto result=plasma.eval(T,rho,c);++count;
    if(reference==0)check(result.eps==0,"exponentially suppressed cold plasma underflows to zero");
    else {
      worst=std::max(worst,std::abs(result.eps/reference-1));
      constexpr double h=1e-5;
      const double t=(std::log(plasma.eval(T*std::exp(h),rho,c).eps)
                     -std::log(plasma.eval(T*std::exp(-h),rho,c).eps))/(2*h);
      const double r=(std::log(plasma.eval(T,rho*std::exp(h),c).eps)
                     -std::log(plasma.eval(T,rho*std::exp(-h),c).eps))/(2*h);
      derivative_error=std::max({derivative_error,std::abs(t-result.dlneps_dlnT)/(1+std::abs(t)),
                                std::abs(r-result.dlneps_dlnRho)/(1+std::abs(r))});
    }
  }
  check(count==21 && worst<2e-12,"published-fit port agrees with independent scalar reference",worst);
  check(derivative_error<1e-7,"plasma temperature and density responses differentiate the fit",derivative_error);
  IdealEos eos;TestOpacity opacity;TestNuclear nuclear;TestLoss loss;NoNeutrinoLosses zero;
  Physics physics{&eos,&opacity,&nuclear};
  Model model;model.M=2e29;model.m={1e29,2e29};model.comp.assign(2,solar_scaled(.7,.02));
  model.y={{std::log(1e8),std::log(1e3),std::log(1e7),1e28},
           {std::log(2e8),std::log(8e2),std::log(9e6),2e28}};
  const auto baseline=zone_residual(model,0,physics,0);
  const auto center=central_residual(model,physics,0);
  physics.neutrino_losses=&zero;
  const auto unchanged=zone_residual(model,0,physics,0);
  check(unchanged.f==baseline.f && unchanged.dfdy_lo==baseline.dfdy_lo
        && unchanged.dfdy_hi==baseline.dfdy_hi,"explicit zero-loss module preserves all values and derivatives exactly");
  physics.neutrino_losses=&loss;
  const auto cooling=zone_residual(model,0,physics,0);
  const double eps0=loss.eval(model.T(0),model.rho(0),model.comp[0]).eps;
  const double eps1=loss.eval(model.T(1),model.rho(1),model.comp[1]).eps;
  check(std::abs((cooling.f[2]-baseline.f[2])/(.5*(eps0+eps1))-1)<1e-13,
        "thermal cooling has the sink sign and trapezoidal zone weights");
  const auto cooling_center=central_residual(model,physics,0);
  check(std::abs((cooling_center.f[1]-center.f[1])/(model.m[0]*eps0)-1)<1e-13,
        "central mass element includes the same thermal sink");
  const auto numerical=zone_residual_numerical(model,0,physics,0,nullptr,1e-5);
  double zone_error=0,center_error=0;
  for(std::size_t j=0;j<NVAR;++j) {
    for(bool high:{false,true}) {
      const auto& a=high?cooling.dfdy_hi:cooling.dfdy_lo;
      const auto& b=high?numerical.dfdy_hi:numerical.dfdy_lo;
      zone_error=std::max(zone_error,std::abs(a[2][j]-b[2][j])/std::max({std::abs(a[2][j]),std::abs(b[2][j]),1e-300}));
    }
    auto plus=model,minus=model;const auto v=static_cast<Var>(j);
    const double step=j==3?1e23:1e-5;plus.y[0][v]+=step;minus.y[0][v]-=step;
    const double finite=(central_residual(plus,physics,0).f[1]-central_residual(minus,physics,0).f[1])/(2*step);
    center_error=std::max(center_error,std::abs(finite-cooling_center.dfdy[1][j])/
                         std::max({std::abs(finite),std::abs(cooling_center.dfdy[1][j]),1e-300}));
  }
  check(zone_error<1e-7 && center_error<1e-7,"central and zone energy Jacobians include loss derivatives",std::max(zone_error,center_error));
  return failures?1:0;
}
