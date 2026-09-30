#pragma once
#include "ember/conduction.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <string>

namespace ember {

// Density continuation of radiative opacity, only where electron conduction
// makes the declared opacity uncertainty a small heat-transport term.
// The fixed overlap avoids following the irregular edge of a source stencil.
class ConductiveInteriorOpacity final : public Opacity {
public:
  static constexpr double anchor_density=8500., full_density=9500.;
  static constexpr double minimum_temperature=8e5, full_temperature=3.4e6,
      source_temperature=3.6e6, maximum_density=1e6;
  static constexpr double uncertainty_factor=10.;
  // The microscopic heat law and the tabulated conduction reference differ.
  // Their measured dense-core ratio is 0.7449--0.8475; use half the reference
  // conductivity in this admission bound and review it as the star cools.
  static constexpr double conductivity_margin=2.;

  struct Domain {
    double anchor{anchor_density},full{full_density},minimum_T{minimum_temperature},
        full_T{full_temperature},source_T{source_temperature},maximum_rho{maximum_density};
    double minimum_X{},maximum_Z{1.},uncertainty{uncertainty_factor};
    // Optional lower-temperature overlap with the supported radiative source.
    double cold_source_T{},cold_full_T{};
  };
  static Domain hydrogen_envelope_domain() {
    // At rho <= 1000 the entire 700–800 kK overlap lies below log R=3.5.
    // A 550–600 kK overlap can ask for the old source beyond that edge.
    return {180.,190.,2e5,7e5,8e5,1000.,.97,1e-8,100.,3e5,3.2e5};
  }

  ConductiveInteriorOpacity(const Opacity& source,const Conduction& conduction,
      double maximum_transport_uncertainty=.001,double opacity_scale=1.)
      : ConductiveInteriorOpacity(source,conduction,maximum_transport_uncertainty,opacity_scale,Domain{}) {}

  ConductiveInteriorOpacity(const Opacity& source,const Conduction& conduction,
      double maximum_transport_uncertainty,double opacity_scale,Domain domain)
      : source_(source),conduction_(conduction),limit_(maximum_transport_uncertainty),scale_(opacity_scale),domain_(domain) {
    if(source.includes_conduction() || !std::isfinite(limit_+scale_) || limit_<=0 || limit_>.01
        || !(domain_.anchor>0 && domain_.full>domain_.anchor && domain_.maximum_rho>=domain_.full
             && domain_.minimum_T>0 && domain_.full_T>domain_.minimum_T && domain_.source_T>domain_.full_T
             && domain_.minimum_X>=0 && domain_.minimum_X<=1 && domain_.maximum_Z>=0 && domain_.maximum_Z<=1
             && domain_.uncertainty>=1 && domain_.uncertainty<=100
             && ((domain_.cold_source_T==0 && domain_.cold_full_T==0)
                 || (domain_.cold_source_T>=domain_.minimum_T
                     && domain_.cold_full_T>domain_.cold_source_T
                     && domain_.cold_full_T<domain_.full_T)))
        || scale_<1/domain_.uncertainty || scale_>domain_.uncertainty)
      throw std::invalid_argument("ConductiveInteriorOpacity: invalid radiative continuation");
  }

  OpacityState eval(double T,double rho,const Composition& c) const override {
    if(rho<=domain_.anchor || T>=domain_.source_T || T<=domain_.cold_source_T
        || c.h1()<domain_.minimum_X || c.Z()>domain_.maximum_Z)
      return source_.eval(T,rho,c);
    if(!(T>=domain_.minimum_T && rho<=domain_.maximum_rho))
      throw std::domain_error("ConductiveInteriorOpacity: outside selected temperature/density bounds");
    const double conduction_opacity=conduction_.eval(T,rho,c).kappa;
    // Infinite conductive opacity means that this heat channel is disabled.
    // It cannot justify an extrapolation, but the original radiative source
    // remains usable within its own domain (and still rejects outside it).
    if(std::isinf(conduction_opacity) && conduction_opacity>0)
      return source_.eval(T,rho,c);
    constexpr double h=.05;
    const auto a=source_.eval(T,domain_.anchor,c);
    const auto b=source_.eval(T,domain_.anchor*std::exp(-h),c);
    const double x=std::log(rho/domain_.anchor),slope=std::log(a.kappa/b.kappa)/h;
    const double nominal=a.kappa*std::exp(slope*x);
    const double kc=conductivity_margin*conduction_opacity;
    if(!(nominal>0 && std::isfinite(nominal) && kc>0 && std::isfinite(kc))) {
      std::ostringstream why;why<<std::scientific<<std::setprecision(3)
          <<"ConductiveInteriorOpacity: invalid anchor or conductivity at T="<<T
          <<", rho="<<rho<<", X="<<c.h1()<<", Z="<<c.Z()
          <<"; anchor="<<a.kappa<<", lower_anchor="<<b.kappa
          <<", continued="<<nominal<<", conduction="<<kc;
      throw std::domain_error(why.str());
    }
    double nominal_blend=nominal;

    OpacityState out{nominal*scale_,a.dlnk_dlnT+x*(a.dlnk_dlnT-b.dlnk_dlnT)/h,slope,
      a.dlnk_dX+x*(a.dlnk_dX-b.dlnk_dX)/h,
      a.dlnk_dZ+x*(a.dlnk_dZ-b.dlnk_dZ)/h,
      a.dlnk_dY3+x*(a.dlnk_dY3-b.dlnk_dY3)/h};
    const double width=std::log(domain_.full/domain_.anchor);
    const auto [wr,dwr]=ramp(x/width);
    const double twidth=std::log(domain_.source_T/domain_.full_T);
    const auto [wt,dwt]=ramp(std::log(domain_.source_T/T)/twidth);
    double wc=1.,dwc=0.;
    if(domain_.cold_source_T>0) {
      const double cold_width=std::log(domain_.cold_full_T/domain_.cold_source_T);
      const auto cold=ramp(std::log(T/domain_.cold_source_T)/cold_width);
      wc=cold.first;dwc=cold.second/cold_width;
    }
    const double w=wr*wt*wc;
    const double dwT=wr*(wt*dwc-wc*dwt/twidth),dwR=wt*wc*dwr/width;
    if(w<1.) {
      const auto old=source_.eval(T,rho,c);
      nominal_blend=std::exp((1-w)*std::log(old.kappa)+w*std::log(nominal));
      const double contrast=std::log(out.kappa/old.kappa);
      out.kappa=std::exp(std::log(old.kappa)+w*contrast);
      out.dlnk_dlnT=(1-w)*old.dlnk_dlnT+w*out.dlnk_dlnT+dwT*contrast;
      out.dlnk_dlnRho=(1-w)*old.dlnk_dlnRho+w*out.dlnk_dlnRho+dwR*contrast;
      out.dlnk_dX=(1-w)*old.dlnk_dX+w*out.dlnk_dX;
      out.dlnk_dZ=(1-w)*old.dlnk_dZ+w*out.dlnk_dZ;
      out.dlnk_dY3=(1-w)*old.dlnk_dY3+w*out.dlnk_dY3;
    }
    if(!(out.kappa>0 && std::isfinite(out.kappa+out.dlnk_dlnT+out.dlnk_dlnRho)))
      throw std::domain_error("ConductiveInteriorOpacity: invalid radiative response");
    // The unknown continued opacity enters log(kappa) with weight w. Lowering
    // it by the selected uncertainty factor changes total conductivity by
    // this fraction, exactly. The
    // supported part of the overlap carries no extrapolation uncertainty.
    const double bound=std::expm1(w*std::log(domain_.uncertainty))*kc/(kc+nominal_blend);
    if(!(bound<=limit_))
      throw std::domain_error("ConductiveInteriorOpacity: radiation uncertainty="+std::to_string(bound)
          +" exceeds selected heat fraction at T="+std::to_string(T)+", rho="+std::to_string(rho));
    calls_.fetch_add(1,std::memory_order_relaxed);
    double previous=maximum_.load(std::memory_order_relaxed);
    while(previous<bound && !maximum_.compare_exchange_weak(previous,bound,std::memory_order_relaxed)) {}
    return out;
  }

  // The physical uncertainty test depends on density; retaining the source
  // interval here gives constrained atmosphere solvers a guaranteed range.
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    return source_.density_range(T,c);
  }
  const char* name() const override {return "radiative density continuation in a conductive interior";}
  std::size_t continued_evaluations() const {return calls_.load(std::memory_order_relaxed);}
  double maximum_transport_uncertainty() const {return maximum_.load(std::memory_order_relaxed);}

private:
  const Opacity& source_;const Conduction& conduction_;double limit_,scale_;Domain domain_;
  mutable std::atomic<std::size_t> calls_{};
  mutable std::atomic<double> maximum_{};
  static std::pair<double,double> ramp(double u) {
    if(u<=0)return {0,0};if(u>=1)return {1,0};
    return {u*u*u*(10+u*(-15+6*u)),30*u*u*(1-u)*(1-u)};
  }
};
} // namespace ember
