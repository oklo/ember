#pragma once
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {

// Interpolate between matching states at two fixed physical gravities.
// Each source is queried only at its declared endpoint within the interval.
class GravityIntervalAtmosphere final : public Atmosphere {
public:
  GravityIntervalAtmosphere(const Eos& eos, const Atmosphere& lower,
      const Atmosphere& upper, double log_g_low, double log_g_high)
      : eos_(eos), lower_(lower), upper_(upper), low_(log_g_low), high_(log_g_high),
        g_low_(std::pow(10., low_)), g_high_(std::pow(10., high_)) {
    if (!std::isfinite(low_+high_+g_low_+g_high_) || !(g_low_>0) || !(high_>low_))
      throw std::invalid_argument("gravity atmosphere interval: invalid coordinates");
  }

  AtmosphereState eval(double teff, double gravity, const Composition& c) const override {
    if (!(gravity>0) || !std::isfinite(gravity))
      throw std::domain_error("gravity atmosphere interval: invalid gravity");
    const double logg=std::log10(gravity);
    if (logg<=low_) return lower_.eval(teff, gravity, c);
    if (logg>=high_) return upper_.eval(teff, gravity, c);
    const auto a=lower_.eval(teff, g_low_, c), b=upper_.eval(teff, g_high_, c);
    if (a.tau!=b.tau || !(a.T>0 && b.T>0 && a.Pgas>0 && b.Pgas>0))
      throw std::domain_error("gravity atmosphere interval: inconsistent joining states");
    const double f=(logg-low_)/(high_-low_);
    const auto blend=[&](double x,double y){return (1-f)*x+f*y;};
    const double width=(high_-low_)*std::log(10.);
    AtmosphereState s{};
    s.T=std::exp(blend(std::log(a.T),std::log(b.T)));
    s.Pgas=std::exp(blend(std::log(a.Pgas),std::log(b.Pgas)));
    s.tau=a.tau;
    const double ar=constants::a_rad*std::pow(a.T,4)/3;
    const double br=constants::a_rad*std::pow(b.T,4)/3;
    const double pr=constants::a_rad*std::pow(s.T,4)/3;
    s.P=s.Pgas+pr;
    s.dlnT_dlnTeff=blend(a.dlnT_dlnTeff,b.dlnT_dlnTeff);
    s.dlnT_dlng=std::log(b.T/a.T)/width;
    const double pg_t=blend((a.P*a.dlnP_dlnTeff-4*ar*a.dlnT_dlnTeff)/a.Pgas,
                            (b.P*b.dlnP_dlnTeff-4*br*b.dlnT_dlnTeff)/b.Pgas);
    s.dlnP_dlnTeff=(s.Pgas*pg_t+4*pr*s.dlnT_dlnTeff)/s.P;
    s.dlnP_dlng=(s.Pgas*std::log(b.Pgas/a.Pgas)/width+4*pr*s.dlnT_dlng)/s.P;
    s.rho=eos_.rho_from_PT(s.T,s.P,c,a.rho);
    return s;
  }
  const char* name() const override {return "interpolated gravity atmosphere interval";}
private:
  const Eos& eos_;
  const Atmosphere& lower_;
  const Atmosphere& upper_;
  double low_, high_, g_low_, g_high_;
};
} // namespace ember
