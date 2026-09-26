#pragma once
#include "ember/opacity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
// Bounded composition continuation where retained warm/hot sources lack
// lower Z. Actual source temperature and density bounds remain in force.
// Molecular source tables, when available, should be blended outside this
// wrapper. The alternative log-kappa law measures prescription sensitivity.
class LowerMetalOpacityContinuation final : public Opacity {
public:
  enum class Method { linear_kappa, logarithmic_kappa };
  LowerMetalOpacityContinuation(const Opacity& source, double anchor,
      double minimum, Method method, double comparison_step=.01)
      : source_(source), anchor_(anchor), minimum_(minimum),
        step_(comparison_step), method_(method) {
    if (!std::isfinite(anchor+minimum+comparison_step) || minimum<=0 ||
        anchor<=minimum || comparison_step<=0 || anchor+comparison_step>=1 ||
        source.includes_conduction())
      throw std::invalid_argument("LowerMetalOpacityContinuation: invalid radiative interval");
  }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    check(c);
    if(c.Z()>=anchor_-2e-14) return source_.eval(T,rho,c);
    auto out=source_.eval(T,rho,at(c,anchor_));
    const auto high=source_.eval(T,rho,at(c,anchor_+step_));
    const double f=(anchor_-c.Z())/step_;
    if(method_==Method::linear_kappa) {
      const double a=(1+f)*out.kappa,b=-f*high.kappa,k=a+b;
      if(!(k>0) || !std::isfinite(k))
        throw std::domain_error("LowerMetalOpacityContinuation: nonpositive linear opacity");
      out.dlnk_dlnT=(a*out.dlnk_dlnT+b*high.dlnk_dlnT)/k;
      out.dlnk_dlnRho=(a*out.dlnk_dlnRho+b*high.dlnk_dlnRho)/k;
      out.dlnk_dX=(a*out.dlnk_dX+b*high.dlnk_dX)/k;
      out.dlnk_dZ=(high.kappa-out.kappa)/(step_*k);
      out.kappa=k;
    } else {
      const double response=std::log(high.kappa/out.kappa);
      out.kappa*=std::exp(-f*response);
      out.dlnk_dlnT=(1+f)*out.dlnk_dlnT-f*high.dlnk_dlnT;
      out.dlnk_dlnRho=(1+f)*out.dlnk_dlnRho-f*high.dlnk_dlnRho;
      out.dlnk_dX=(1+f)*out.dlnk_dX-f*high.dlnk_dX;
      out.dlnk_dZ=response/step_;
    }
    out.dlnk_dY3=0.;
    if(!(out.kappa>0) || !std::isfinite(out.kappa) ||
       !std::isfinite(out.dlnk_dlnT) || !std::isfinite(out.dlnk_dlnRho) ||
       !std::isfinite(out.dlnk_dX) || !std::isfinite(out.dlnk_dZ))
      throw std::domain_error("LowerMetalOpacityContinuation: invalid response");
    return out;
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    check(c);
    if(c.Z()>=anchor_-2e-14) return source_.density_range(T,c);
    const auto a=source_.density_range(T,at(c,anchor_));
    const auto b=source_.density_range(T,at(c,anchor_+step_));
    if(!a || !b) throw std::domain_error("LowerMetalOpacityContinuation: missing density bounds");
    DensityRange r{std::max(a->min,b->min),std::min(a->max,b->max)};
    if(!(r.min<r.max)) throw std::domain_error("LowerMetalOpacityContinuation: no density overlap");
    return r;
  }
  const char* name() const override{return "bounded lower-metal radiative opacity continuation";}
private:
  const Opacity& source_;double anchor_,minimum_,step_;Method method_;
  void check(const Composition& c) const {
    if(c.basis!=AbundanceBasis::atomic_mass || c.X[1]!=0 ||
       !std::isfinite(c.Z()) || c.Z()<minimum_-2e-14 || c.Z()>=1 ||
       std::abs(c.sum()-1)>1e-10)
      throw std::domain_error("LowerMetalOpacityContinuation: composition outside interval");
    for(double x:c.X) if(!std::isfinite(x) || x<0)
      throw std::domain_error("LowerMetalOpacityContinuation: invalid composition");
  }
  static Composition at(const Composition& c,double z) {
    auto out=c;const double old=c.Z();
    for(std::size_t j=3;j<METAL_END;++j)out.X[j]*=z/old;
    out.X[2]+=old-z;
    if(out.X[2]<0) throw std::domain_error("LowerMetalOpacityContinuation: comparison helium negative");
    return out;
  }
};
} // namespace ember
