#pragma once
#include "ember/atmosphere.hpp"
#include "ember/model.hpp"
#include "ember/structure.hpp"
#include <array>

namespace ember {

struct CentralResidual {
  std::array<double, 2> f{};  // central sphere radius, L - enclosed luminosity mass * eps_total
  std::array<std::array<double, NVAR>, 2> dfdy{};
};

// Leading regular-center expansion over the unresolved sphere 0 < m < m[0].
// The innermost point has finite radius and strictly positive enclosed mass.
// Density and heating are approximated by their values at that point; check
// convergence as the central sphere shrinks. Positive dt requires prev on
// the same mesh, with exactly the energy convention used by the zone rows.
// With face luminosities, the energy volume extends to (m[0]+m[1])/2;
// the radius condition still describes the sphere inside m[0].
CentralResidual central_residual(const Model&, const Physics&, double dt = 0.0,
                                 const Model* prev = nullptr);

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
                                 const Composition&, const Eos&, const Atmosphere&, bool derivatives = true);

} // namespace ember
