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
// The overlap stays within the source support; missing source cells still reject.
class ConductiveInteriorOpacity final : public Opacity {
public:
  static constexpr double anchor_density=8500., full_density=9500.;
  static constexpr double minimum_temperature=8e5, full_temperature=3.4e6,
      source_temperature=3.6e6, maximum_density=1e6;
  static constexpr double uncertainty_factor=10.;
  // The microscopic heat law and the tabulated conduction reference differ.
  // Cold-profile comparisons give microscopic/reference conductivity ratios
  // down to 0.256. Use one quarter of the reference in this admission bound;
  // this is a measured range, not a guarantee outside the assessed conditions.
  static constexpr double conductivity_margin=4.;

  struct Domain {
    double anchor{anchor_density},full{full_density},minimum_T{minimum_temperature},
        full_T{full_temperature},source_T{source_temperature},maximum_rho{maximum_density};
    double minimum_X{},maximum_Z{1.},uncertainty{uncertainty_factor};
    // Optional lower-temperature overlap with the supported radiative source.
    double cold_source_T{},cold_full_T{};
    // Optional composition overlap; zero preserves the older selection rule.
    double full_X{};
    // Optional upper join in the radiative source coordinate
    // log R = log10(rho) - 3 log10(T/1e6), instead of fixed temperatures.
    double source_logR{},full_logR{};
  };
  static Domain hydrogen_envelope_domain() {
    // At rho <= 1000 the entire 700–800 kK overlap lies below log R=3.5.
    // A 550–600 kK overlap can ask for the old source beyond that edge.
    // Finish the composition join below the source's atomic X=.75 boundary,
    // including the conversion from baryonic H/He fractions.
    return {180.,190.,2e5,7e5,8e5,1000.,.70,1e-8,100.,3e5,3.2e5,.745};
  }

  static Domain ionized_hydrogen_envelope_domain() {
    auto domain=hydrogen_envelope_domain();
    domain.maximum_rho=1e5;
    domain.minimum_T=1.3e5;
    // Below log T = 5.7 the radiative source is the warm/bridge blend, whose ATOMIC bridge is supported
    // only to log rho = 2.4 there (log rho = log R + 3 log10(T/1e6)). The part of this overlap that still
    // reads the source at the actual state (log R < full_logR) must stay inside it up to log T = 5.7:
    // full_logR <= 3.3, kept at 3.25 for the interpolation stencil.
    domain.source_logR=3.0;domain.full_logR=3.25;
    // Direct ionized conduction remains active below 300 kK. Returning to
    // the unextended radiative table there can leave its density support.
    // Source anchors and H-layer charge/conductivity checks support 130 kK.
    // Keep the contribution bound and explicit lower temperature limit.
    domain.cold_source_T=0;domain.cold_full_T=0;
    // The same fully ionized, low-metal H/He layer becomes helium-rich at
    // depth. Admission depends on radiative heat transport, not hydrogen
    // abundance; both source anchors and the uncertainty bound still apply.
    domain.minimum_X=0;domain.full_X=0;
    // Use the same metal range as the radiative source. Its actual anchor
    // coverage and the local heat-transport bound remain required, including
    // the source conversion between baryonic and atomic mass fractions.
    domain.maximum_Z=.16;
    return domain;
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
             && (domain_.full_X==0 || (domain_.full_X>domain_.minimum_X && domain_.full_X<=1))
             && ((domain_.source_logR==0 && domain_.full_logR==0)
                 || (std::isfinite(domain_.source_logR+domain_.full_logR)
                     && domain_.full_logR>domain_.source_logR))
             && domain_.uncertainty>=1 && domain_.uncertainty<=100
             && ((domain_.cold_source_T==0 && domain_.cold_full_T==0)
                 || (domain_.cold_source_T>=domain_.minimum_T
                     && domain_.cold_full_T>domain_.cold_source_T
                     && domain_.cold_full_T<domain_.full_T)))
        || scale_<1/domain_.uncertainty || scale_>domain_.uncertainty)
      throw std::invalid_argument("ConductiveInteriorOpacity: invalid radiative continuation");
  }

  OpacityState eval(double T,double rho,const Composition& c) const override {
    const bool coordinate_join=domain_.full_logR>domain_.source_logR;
    const double logR=coordinate_join?std::log10(rho)-3*std::log10(T/1e6):0.;
    if(rho<=domain_.anchor || (coordinate_join?logR<=domain_.source_logR:T>=domain_.source_T)
        || T<=domain_.cold_source_T
        || c.h1()<domain_.minimum_X || c.Z()>domain_.maximum_Z)
      return source_.eval(T,rho,c);
    double wx=1.,dwx=0.;
    if(domain_.full_X>0) {
      const double width=domain_.full_X-domain_.minimum_X;
      const auto composition=ramp((c.h1()-domain_.minimum_X)/width);
      wx=composition.first;dwx=composition.second/width;
      if(wx==0)return source_.eval(T,rho,c);
    }
    if(!(T>=domain_.minimum_T && rho<=domain_.maximum_rho)) {
      std::ostringstream why;
      why<<"ConductiveInteriorOpacity: outside selected temperature/density bounds: "
          <<std::scientific<<std::setprecision(3)<<"T="<<T<<", rho="<<rho
          <<", X="<<c.h1()<<", Z="<<c.Z()<<"; minimum T="<<domain_.minimum_T
          <<", maximum rho="<<domain_.maximum_rho;
      throw std::domain_error(why.str());
    }
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
    const auto upper=coordinate_join
        ?ramp((logR-domain_.source_logR)/(domain_.full_logR-domain_.source_logR))
        :ramp(std::log(domain_.source_T/T)/twidth);
    const double wt=upper.first;
    const double upper_derivative=coordinate_join
        ?upper.second/((domain_.full_logR-domain_.source_logR)*std::log(10.)):0.;
    const double dwtT=coordinate_join?-3*upper_derivative:-upper.second/twidth;
    const double dwtR=upper_derivative;
    double wc=1.,dwc=0.;
    if(domain_.cold_source_T>0) {
      const double cold_width=std::log(domain_.cold_full_T/domain_.cold_source_T);
      const auto cold=ramp(std::log(T/domain_.cold_source_T)/cold_width);
      wc=cold.first;dwc=cold.second/cold_width;
    }
    const double w=wr*wt*wc*wx;
    const double dwT=wx*wr*(wt*dwc+wc*dwtT);
    const double dwR=wx*wc*(wt*dwr/width+wr*dwtR);
    const double dwX=wr*wt*wc*dwx;
    if(w<1.) {
      const auto old=source_.eval(T,rho,c);
      nominal_blend=std::exp((1-w)*std::log(old.kappa)+w*std::log(nominal));
      const double contrast=std::log(out.kappa/old.kappa);
      out.kappa=std::exp(std::log(old.kappa)+w*contrast);
      out.dlnk_dlnT=(1-w)*old.dlnk_dlnT+w*out.dlnk_dlnT+dwT*contrast;
      out.dlnk_dlnRho=(1-w)*old.dlnk_dlnRho+w*out.dlnk_dlnRho+dwR*contrast;
      out.dlnk_dX=(1-w)*old.dlnk_dX+w*out.dlnk_dX+dwX*contrast;
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
    const double v=1-u;
    const double value=u<=.5?u*u*u*(10+u*(-15+6*u))
        :1-v*v*v*(10+v*(-15+6*v));
    return {value,30*u*u*v*v};
  }
};
} // namespace ember
