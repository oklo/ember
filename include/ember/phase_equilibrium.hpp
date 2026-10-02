#pragma once
#include "ember/eos_helmholtz.hpp"
#include <array>
#include <functional>
#include <cstddef>
#include <span>
#include <vector>

namespace ember {
// Coordinates: ln T, ln rho, X_H, X_He3, Z; He4 is the remainder.
using PhaseCoordinates=std::array<double,5>;
// Dimensionless F/(R_gas T), and derivatives in PhaseCoordinates.
struct PhasePotential {
  double value{};
  PhaseCoordinates gradient{};
  std::array<PhaseCoordinates,5> hessian{};
  // Optional third derivatives. When supplied by every phase, the material
  // interface differentiates coexistence analytically instead of differencing.
  std::array<std::array<PhaseCoordinates,5>,5> third{};
  bool has_third{};
};
using PhaseEvaluator=std::function<PhasePotential(const PhaseCoordinates&)>;
struct PhaseState {
  double log_density{};
  std::array<double,3> composition{};
  double mass_fraction{};
};
struct PhaseMixture {
  PhasePotential potential;
  std::vector<PhaseState> phases;
  // Changes of each phase's mass fraction with the five bulk coordinates.
  std::vector<PhaseCoordinates> fraction_response;
  // Per phase: ln rho, ln X, ln He3, ln Z responses. These also supply a
  // starting guess for nearby material derivatives without extra phase solves.
  std::vector<std::array<PhaseCoordinates,4>> state_response;
  double residual{};
  std::size_t iterations{};
};
// One homogeneous phase, or two or three distinct coexisting phases.
// A phase may use the same evaluator as another (two solid compositions, for
// example). Coexisting phase fractions and all bulk abundances must be positive.
// This finds a stationary coexistence state; global selection is the caller's
// responsibility. No phase is discarded or merged during this solve. For one
// phase the bulk state is used directly; its supplied seed is not a constraint.
PhaseMixture equilibrate_phases(const PhaseCoordinates&,std::span<const PhaseState>,
                               std::span<const PhaseEvaluator>);
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
// Supplied third derivatives differentiate the constraints analytically;
// otherwise four nearby solves differentiate the implicit Hessian. The same
// distinct phases must remain supported throughout that stencil; a failed solve is not silently replaced
// by a single phase. Phase selection and its boundaries remain the caller's job.
// Intended for material evaluation/tabulation, not unseeded searches in each
// stellar Newton iteration. Do not restore mixing or add latent heat again.
PhaseMaterial phase_equilibrium_material(const PhaseCoordinates&,const PhaseSplit&,
    const PhaseEvaluator&,const PhaseEvaluator&,double log_step=1e-4);

struct PhaseMixtureMaterial {
  PhaseMixture equilibrium;
  std::array<HelmholtzJet,10> jets{};
  std::array<double,2> log_steps{};
  bool analytic_third{};
};
// Without supplied third derivatives, the difference step is an upper bound.
// Fraction and state responses shorten it near a phase boundary.
// Failed/vanishing phases are reported, not replaced
// by a different phase set. Radiation and phase heat must not be added again.
PhaseMixtureMaterial phase_equilibrium_material(const PhaseCoordinates&,
    std::span<const PhaseState>,std::span<const PhaseEvaluator>,double log_step=1e-4);
} // namespace ember
