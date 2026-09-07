#pragma once
#include "ember/atmosphere.hpp"
#include "ember/model.hpp"
#include <array>

namespace ember {

struct SurfaceResidual {
  std::array<double, 2> f{};  // ln T - ln T_atm, ln P - ln P_atm
  std::array<std::array<double, NVAR>, 2> dfdy{};
  AtmosphereState atmosphere{};
  double Teff{}, gravity{};
};

// The outer mesh point is the atmospheric matching surface. In the thin
// atmosphere approximation its radius is the photospheric radius and its
// enclosed mass is the total mass; no fixed interior mass fraction is used.
// L>0 is required here, even though interior points may carry inward flux.
SurfaceResidual surface_residual(const Point& surface, double mass,
                                 const Composition&, const Eos&, const Atmosphere&);

} // namespace ember
