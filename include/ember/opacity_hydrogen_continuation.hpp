#pragma once
#include "ember/opacity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {

// Explicit, bounded linear continuation of log opacity in hydrogen abundance.
// Temperature, density and metallicity are queried only within the original
// source domain. The continuation needs separate physical and stellar checks.
class HydrogenOpacityContinuation final : public Opacity {
public:
  HydrogenOpacityContinuation(const Opacity& source, double anchor, double step, double maximum)
      : source_(source), anchor_(anchor), step_(step), maximum_(maximum) {
    if (!std::isfinite(anchor+step+maximum) || step<=0 || anchor<step ||
        maximum<=anchor || maximum>1 || source.includes_conduction())
      throw std::invalid_argument("HydrogenOpacityContinuation: invalid radiative continuation");
  }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    check(c);
    if(c.h1()<=anchor_) return source_.eval(T,rho,c);
    const auto a=source_.eval(T,rho,at(c,anchor_));
    const auto b=source_.eval(T,rho,at(c,anchor_-step_));
    const double f=(c.h1()-anchor_)/step_;
    auto out=a;
    out.kappa=std::exp(std::log(a.kappa)+f*std::log(a.kappa/b.kappa));
    out.dlnk_dlnT=a.dlnk_dlnT+f*(a.dlnk_dlnT-b.dlnk_dlnT);
    out.dlnk_dlnRho=a.dlnk_dlnRho+f*(a.dlnk_dlnRho-b.dlnk_dlnRho);
    out.dlnk_dX=std::log(a.kappa/b.kappa)/step_;
    out.dlnk_dZ=a.dlnk_dZ+f*(a.dlnk_dZ-b.dlnk_dZ);
    out.dlnk_dY3=0.;
    if(!std::isfinite(out.kappa) || !(out.kappa>0))
      throw std::domain_error("HydrogenOpacityContinuation: invalid opacity");
    return out;
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    check(c);
    if(c.h1()<=anchor_) return source_.density_range(T,c);
    const auto a=source_.density_range(T,at(c,anchor_));
    const auto b=source_.density_range(T,at(c,anchor_-step_));
    if(!a || !b) throw std::domain_error("HydrogenOpacityContinuation: missing source density range");
    DensityRange r{std::max(a->min,b->min),std::min(a->max,b->max)};
    if(!(r.min<r.max)) throw std::domain_error("HydrogenOpacityContinuation: no density overlap");
    return r;
  }
  const char* name() const override{return "bounded log-opacity hydrogen continuation";}
private:
  const Opacity& source_;
  double anchor_,step_,maximum_;
  void check(const Composition& c) const {
    if(c.basis!=AbundanceBasis::atomic_mass || c.X[1]!=0 ||
       !std::isfinite(c.h1()) || c.h1()<0 || c.h1()>maximum_)
      throw std::domain_error("HydrogenOpacityContinuation: composition outside interval");
  }
  static Composition at(const Composition& c,double x) {
    auto out=c;out.X[2]+=out.X[0]-x;out.X[0]=x;return out;
  }
};

} // namespace ember
