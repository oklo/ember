#pragma once
#include "ember/eos_helmholtz.hpp"
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

struct PhaseMaterial {
  PhaseEquilibrium equilibrium;
  // Material F/T, including the supplied phase mixing terms, without radiation.
  // Channels: value, X, He3, Z, XX, XHe3, XZ, He3He3, He3Z, ZZ.
  // Thermal/density derivatives are populated through total order three.
  std::array<HelmholtzJet,10> jets{};
};
// Connect coexistence thermodynamics to the material EOS representation.
// Four nearby coexistence solves differentiate the implicit Hessian for its
// thermal/density third derivatives. The same distinct phases must remain
// supported throughout this stencil; a failed solve is not silently replaced
// by a single phase. Phase selection and its boundaries remain the caller's job.
// Intended for material evaluation/tabulation, not unseeded searches in each
// stellar Newton iteration. Do not restore mixing or add latent heat again.
PhaseMaterial phase_equilibrium_material(const PhaseCoordinates&,const PhaseSplit&,
    const PhaseEvaluator&,const PhaseEvaluator&,double log_step=1e-4);
} // namespace ember
