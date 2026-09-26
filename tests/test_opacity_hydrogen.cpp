#include "ember/opacity_hydrogen.hpp"
#include <cstdio>

using namespace ember;
struct AnalyticOpacity : Opacity {
  double scale, low, high;
  explicit AnalyticOpacity(double s, double a, double b):scale(s),low(a),high(b){}
  OpacityState eval(double T, double rho, const Composition& c) const override {
    if (!(rho>=low && rho<=high)) throw std::domain_error("analytic density");
    const double x=c.h1(), z=c.Z();
    return {std::exp(scale*(x*x*std::log(T)+x*std::log(rho)+z*x)),
            scale*x*x,scale*x,scale*(2*x*std::log(T)+std::log(rho)+z),scale*x};
  }
  std::optional<DensityRange> density_range(double, const Composition&) const override {
    return DensityRange{low,high};
  }
  const char* name() const override{return "analytic opacity control";}
};
int main() {
  int failures=0;
  auto check=[&](bool ok,const char* name){std::printf("%s %s\n",ok?"PASS":"FAIL",name);failures+=!ok;};
  AnalyticOpacity retained(.3,.01,100),source(.7,.1,10);
  HydrogenOpacityExtension extension(retained,source,.75,.95);
  auto c=solar_scaled(.7,.02);
  const auto a=retained.eval(3.,2.,c), b=extension.eval(3.,2.,c);
  check(a.kappa==b.kappa && a.dlnk_dX==b.dlnk_dX && a.dlnk_dZ==b.dlnk_dZ,
        "retained composition response is exact");
  c=solar_scaled(.85,.02);
  const auto s=extension.eval(3.,2.,c);
  const double expected=std::exp(.3*(.75*.75*std::log(3.)+.75*std::log(2.)+.02*.75)
      +.7*((.85*.85-.75*.75)*std::log(3.)+(.85-.75)*std::log(2.)+.02*(.85-.75)));
  check(std::abs(s.kappa/expected-1)<2e-15,"independent analytic composition dependence");
  constexpr double h=1e-5;
  auto ln=[&](double t,double r,double x,double z){return std::log(extension.eval(t,r,solar_scaled(x,z)).kappa);};
  check(std::abs((ln(3*std::exp(h),2,.85,.02)-ln(3*std::exp(-h),2,.85,.02))/(2*h)-s.dlnk_dlnT)<1e-9,
        "temperature derivative at fixed density");
  check(std::abs((ln(3,2*std::exp(h),.85,.02)-ln(3,2*std::exp(-h),.85,.02))/(2*h)-s.dlnk_dlnRho)<1e-9,
        "density derivative");
  check(std::abs((ln(3,2,.85+h,.02)-ln(3,2,.85-h,.02))/(2*h)-s.dlnk_dX)<1e-9,
        "hydrogen derivative with helium displaced");
  check(std::abs((ln(3,2,.85,.02+h)-ln(3,2,.85,.02-h))/(2*h)-s.dlnk_dZ)<1e-9,
        "metallicity derivative with helium displaced");
  const auto range=extension.density_range(3.,c).value();
  check(range.min==.1 && range.max==10,"source density intervals intersect");
  check(std::abs(ln(3,2,.75+1e-10,.02)-ln(3,2,.75,.02))<1e-9,"value continuous at anchor");
  bool rejected=false;
  try{extension.eval(3.,2.,solar_scaled(.951,.02));}catch(const std::domain_error&){rejected=true;}
  check(rejected,"hydrogen outside explicit extension rejected");
  struct CombinedSource : AnalyticOpacity {
    CombinedSource():AnalyticOpacity(.3,.01,100){}
    bool includes_conduction() const override{return true;}
  } combined;
  rejected=false;
  try{HydrogenOpacityExtension bad(combined,source,.75,.95);}
  catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"combined conduction cannot enter a radiative source ratio");
  return failures?1:0;
}
