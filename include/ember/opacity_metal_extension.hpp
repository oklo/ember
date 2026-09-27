#pragma once
#include "ember/opacity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {

// Continue a retained radiative opacity above its largest source Z. The
// primary choice transfers a separately tabulated composition ratio. A
// linear-in-kappa alternative measures sensitivity to that choice. Neither
// method extrapolates temperature or density, or changes the metal pattern.
class MetalOpacityExtension final : public Opacity {
public:
  enum class Method { source_ratio, linear_kappa };
  MetalOpacityExtension(const Opacity& source, const Opacity& dependence,
                        double anchor, double maximum, Method method,
                        double comparison_step=.01)
      : source_(source), dependence_(dependence), anchor_(anchor), maximum_(maximum),
        step_(comparison_step), method_(method) {
    if (!std::isfinite(anchor+maximum+comparison_step) || anchor<=comparison_step ||
        comparison_step<=0 || maximum<=anchor || maximum>=1 ||
        source.includes_conduction() || dependence.includes_conduction())
      throw std::invalid_argument("MetalOpacityExtension: invalid radiative extension");
  }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    check(c);
    if(c.Z()<=anchor_+2e-14) return source_.eval(T,rho,c);
    const auto fixed=at(c,anchor_);
    auto out=source_.eval(T,rho,fixed);
    if(method_==Method::source_ratio) {
      const auto high=dependence_.eval(T,rho,c), low=dependence_.eval(T,rho,fixed);
      out.kappa*=high.kappa/low.kappa;
      out.dlnk_dlnT+=high.dlnk_dlnT-low.dlnk_dlnT;
      out.dlnk_dlnRho+=high.dlnk_dlnRho-low.dlnk_dlnRho;
      out.dlnk_dX+=high.dlnk_dX-low.dlnk_dX;
      out.dlnk_dZ=high.dlnk_dZ;
    } else {
      const auto low=source_.eval(T,rho,at(c,anchor_-step_));
      const double f=(c.Z()-anchor_)/step_, a=(1+f)*out.kappa, b=-f*low.kappa;
      const double k=a+b;
      out.dlnk_dlnT=(a*out.dlnk_dlnT+b*low.dlnk_dlnT)/k;
      out.dlnk_dlnRho=(a*out.dlnk_dlnRho+b*low.dlnk_dlnRho)/k;
      out.dlnk_dX=(a*out.dlnk_dX+b*low.dlnk_dX)/k;
      out.dlnk_dZ=(out.kappa-low.kappa)/(step_*k);
      out.kappa=k;
    }
    out.dlnk_dY3=0.;
    if (!(out.kappa>0) || !std::isfinite(out.kappa) ||
        !std::isfinite(out.dlnk_dlnT) || !std::isfinite(out.dlnk_dlnRho) ||
        !std::isfinite(out.dlnk_dX) || !std::isfinite(out.dlnk_dZ))
      throw std::domain_error("MetalOpacityExtension: invalid opacity response");
    return out;
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    check(c);
    if(c.Z()<=anchor_+2e-14) return source_.density_range(T,c);
    auto range=source_.density_range(T,at(c,anchor_));
    auto intersect=[&](std::optional<DensityRange> r) {
      if(!range || !r) throw std::domain_error("MetalOpacityExtension: missing density bounds");
      range=DensityRange{std::max(range->min,r->min),std::min(range->max,r->max)};
      if(!(range->min<range->max)) throw std::domain_error("MetalOpacityExtension: no density overlap");
    };
    if(method_==Method::source_ratio) {
      intersect(dependence_.density_range(T,c));
      intersect(dependence_.density_range(T,at(c,anchor_)));
    } else intersect(source_.density_range(T,at(c,anchor_-step_)));
    return range;
  }
  const char* name() const override{return "bounded upper-metal radiative opacity extension";}
private:
  const Opacity& source_;const Opacity& dependence_;
  double anchor_,maximum_,step_;Method method_;
  void check(const Composition& c) const {
    if(c.basis!=AbundanceBasis::atomic_mass || c.X[1]!=0 ||
       !std::isfinite(c.Z()) || c.Z()<0 || c.Z()>maximum_+2e-14 ||
       std::abs(c.sum()-1)>1e-10)
      throw std::domain_error("MetalOpacityExtension: composition outside interval");
    for(double x:c.X) if(!std::isfinite(x) || x<0)
      throw std::domain_error("MetalOpacityExtension: invalid composition");
  }
  static Composition at(const Composition& c,double z) {
    auto out=c;const double original=c.Z();
    for(std::size_t j=3;j<METAL_END;++j)out.X[j]*=z/original;
    out.X[2]+=original-z;return out;
  }
};
} // namespace ember
