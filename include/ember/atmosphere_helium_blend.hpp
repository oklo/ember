#pragma once
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {

// Explicitly join two atmosphere prescriptions over a covered helium interval.
// Both sources must accept the actual composition throughout the overlap.
// The quintic weight has zero first and second derivatives at either endpoint.
// This numerical connection does not establish agreement between the sources.
class HeliumBlendAtmosphere final : public Atmosphere {
public:
  HeliumBlendAtmosphere(const Eos& eos, const Atmosphere& low_helium,
                       const Atmosphere& high_helium, double low, double high)
      : eos_(eos), low_source_(low_helium), high_source_(high_helium),
        low_(low), high_(high) {
    if (!std::isfinite(low) || !std::isfinite(high) || low < 0 || low >= high || high >= 1)
      throw std::invalid_argument("HeliumBlendAtmosphere: invalid helium interval");
  }

  AtmosphereState eval(double teff, double gravity, const Composition& c) const override {
    const double helium = c.X[1]+c.X[2];
    if (!std::isfinite(helium) || helium < 0)
      throw std::domain_error("HeliumBlendAtmosphere: invalid helium abundance");
    // Preserve each original atmosphere exactly outside the explicit overlap.
    if (helium <= low_) return low_source_.eval(teff, gravity, c);
    if (helium >= high_) return high_source_.eval(teff, gravity, c);
    const auto a = low_source_.eval(teff, gravity, c);
    const auto b = high_source_.eval(teff, gravity, c);
    if (a.tau != b.tau || !(a.T > 0 && b.T > 0 && a.Pgas > 0 && b.Pgas > 0))
      throw std::domain_error("HeliumBlendAtmosphere: incompatible source states");
    const double u = (helium-low_)/(high_-low_);
    const double w = u*u*u*(10+u*(-15+6*u));
    const auto mix = [w](double x, double y) { return (1-w)*x+w*y; };
    AtmosphereState s{};
    s.T = std::exp(mix(std::log(a.T), std::log(b.T)));
    s.Pgas = std::exp(mix(std::log(a.Pgas), std::log(b.Pgas)));
    s.tau = a.tau;
    const double ar = constants::a_rad*std::pow(a.T, 4)/3;
    const double br = constants::a_rad*std::pow(b.T, 4)/3;
    const double pr = constants::a_rad*std::pow(s.T, 4)/3;
    s.P = s.Pgas+pr;
    if (!std::isfinite(s.P) || s.P <= pr || !std::isfinite(s.T))
      throw std::domain_error("HeliumBlendAtmosphere: invalid blended pressure or temperature");
    s.dlnT_dlnTeff = mix(a.dlnT_dlnTeff, b.dlnT_dlnTeff);
    s.dlnT_dlng = mix(a.dlnT_dlng, b.dlnT_dlng);
    const auto pressure = [&](double ap, double at, double bp, double bt, double dt) {
      const double dpg = mix((a.P*ap-4*ar*at)/a.Pgas, (b.P*bp-4*br*bt)/b.Pgas);
      return (s.Pgas*dpg+4*pr*dt)/s.P;
    };
    s.dlnP_dlnTeff = pressure(a.dlnP_dlnTeff, a.dlnT_dlnTeff,
                             b.dlnP_dlnTeff, b.dlnT_dlnTeff, s.dlnT_dlnTeff);
    s.dlnP_dlng = pressure(a.dlnP_dlng, a.dlnT_dlng,
                          b.dlnP_dlng, b.dlnT_dlng, s.dlnT_dlng);
    s.rho = eos_.rho_from_PT(s.T, s.P, c, a.rho);
    return s;
  }
  const char* name() const override { return "helium-dependent atmosphere connection"; }

private:
  const Eos& eos_;
  const Atmosphere& low_source_;
  const Atmosphere& high_source_;
  double low_, high_;
};

} // namespace ember
