#pragma once
#include <array>
#include <functional>
#include <cstddef>

namespace ember {
// Coordinates: ln T, ln rho, X_H, X_He3, Z; He4 is the remainder.
using PhaseCoordinates=std::array<double,5>;
// Dimensionless F/(R_gas T), and derivatives in PhaseCoordinates.
struct PhasePotential {
  double value{};
  PhaseCoordinates gradient{};
  std::array<PhaseCoordinates,5> hessian{};
};
using PhaseEvaluator=std::function<PhasePotential(const PhaseCoordinates&)>;
struct PhaseSplit {
  double log_density_first{},log_density_second{};
  std::array<double,3> composition_first{},composition_second{};
  double second_mass_fraction{};
};
struct PhaseEquilibrium {
  PhasePotential potential;
  PhaseSplit split;
  double residual{};
  std::size_t iterations{};
};
// Find a stationary pair of the specified phases, conserving volume and every
// component. Evaluators must reject unsupported states with domain_error.
// The caller selects the globally stable phase pair and supplies a nearby
// split. This does not search the phase diagram or model separation kinetics.
// Requires positive species fractions and nonzero amounts of both phases.
// The returned Hessian differentiates the equilibrium constraints; phase
// fractions and compositions are not frozen in the thermal response.
PhaseEquilibrium equilibrate_phases(const PhaseCoordinates&,const PhaseSplit&,
                                   const PhaseEvaluator&,const PhaseEvaluator&);
} // namespace ember
