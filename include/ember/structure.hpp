#pragma once
#include "ember/eos.hpp"
#include "ember/model.hpp"
#include "ember/nuclear.hpp"
#include "ember/opacity.hpp"
#include <array>

namespace ember {

// The four difference equations of a zone, and their derivatives with respect
// to the eight variables that bound it.
//
//   (1) mass:      d(ln r)/dm  = 1/(4 pi r^3 rho)
//   (2) hydrostatic: d(ln P)/dm = -G m /(4 pi r^4 P)
//   (3) energy:    dL/dm       = eps_deposited - du/dt - P d(1/rho)/dt
//   (4) transport: d(ln T)/dm  = grad * d(ln P)/dm
//
// Nuclear neutrinos are already removed from eps_deposited. Internal energy
// includes its composition dependence but excludes nuclear rest mass.
// Time differences are backward Euler; nodal sources use trapezoidal mass
// weights, including both endpoints of the zone between points i and i+1. The
// gradient in (4) is radiative in Schwarzschild-stable zones, and the
// Bohm-Vitense mixing-length result in unstable zones. The choice of form
// matters as much as the value: scaling the radiative equation by a convective
// efficiency is only well conditioned while that efficiency is of order unity, a lesson
// this code inherits at the cost of several days.
struct ZoneResidual {
  std::array<double, NVAR> f{};                              // residuals
  std::array<std::array<double, NVAR>, NVAR> dfdy_lo{};      // d f / d y_i
  std::array<std::array<double, NVAR>, NVAR> dfdy_hi{};      // d f / d y_{i+1}
};

enum class ConvectiveCriterion { schwarzschild, ledoux };
struct Physics {
  const Eos* eos{};
  const Opacity* opacity{};
  const Nuclear* nuclear{};
  double alpha_mlt{1.9};             // l/H_P, strictly positive; not calibrated
  ConvectiveCriterion criterion{ConvectiveCriterion::schwarzschild};
  double alpha_semiconvection{};     // Langer mixing-only closure; zero disables
  double alpha_thermohaline{};       // Kippenhahn closure; zero disables
};

// Evaluate one zone.  `dt` <= 0 means a static model: the time-dependent term
// in (3) is dropped rather than divided by zero. Positive dt requires prev on
// the same mass mesh. The Jacobian uses analytic module derivatives and the
// chain rule; no state perturbations or whole-model copies are made.
ZoneResidual zone_residual(const Model& mdl, std::size_t i,
                           const Physics& phys, double dt,
                           const Model* prev = nullptr);

// Values only, and the independent finite-difference reference Jacobian.
std::array<double, NVAR> zone_equations(const Model&, std::size_t i,
                                      const Physics&, double dt, const Model* prev = nullptr);
ZoneResidual zone_residual_numerical(const Model&, std::size_t i,
                                    const Physics&, double dt, const Model* prev = nullptr,
                                    double relative_step = 1e-5);

} // namespace ember
