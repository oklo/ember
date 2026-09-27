#pragma once
#include "ember/relaxation.hpp"
#include <span>
#include <utility>

namespace ember {
using MixingRegions = std::vector<std::pair<std::size_t,std::size_t>>; // [first,last), includes isolated points
std::vector<double> nodal_mass_weights(const Model&);
MixingRegions schwarzschild_mixing_regions(const Model&,const Physics&);
MixingRegions convective_mixing_regions(const Model&,const Physics&,std::size_t threads=1);
// Implicit pp burning with instantaneous homogeneous mixing within each
// supplied connected region. Isolated points burn locally. Metals are inert.
std::vector<Composition> burn_and_mix(const Model& thermal,const Model& previous,
    const Nuclear&,const MixingRegions&,double dt,double tolerance=1e-13);
// Slow mixing coefficients at mesh faces, cm^2/s. Diffusive heat transport
// remains radiative+conductive in these Ledoux-stable layers.
std::vector<double> secular_mixing_diffusivities(const Model&,const Physics&);
// Conventional Ledoux/Schwarzschild MLT, using the same thermal luminosity
// and microscopic heat split as structure. No temperature cutoff is applied.
// The composition-to-buoyancy ratio diagnoses the closure approximation;
// it is not a proof of a consistent parcel treatment at large composition
// gradients (see docs/CONVECTION_COMPOSITION_CHECK.md).
struct ConvectiveMixingFace {
  double diffusivity{},velocity{},length{},buoyancy_contrast{},composition_term{};
};
std::vector<ConvectiveMixingFace> convective_mixing_faces(const Model&,const Physics&,std::size_t threads=1);
// (4 pi r^2 rho)^2 D / Delta m, in g/s, including selected slow mixing.
std::vector<double> finite_mixing_conductances(const Model&,const Physics&);
// Conservative backward-Euler diffusion and burning, with instantaneous
// convection collapsed into mass-weighted regions. Zero exterior flux.
std::vector<Composition> burn_and_transport(const Model& thermal,const Model& previous,
    const Nuclear&,const MixingRegions&,const std::vector<double>& diffusivity,
    double dt,double tolerance=1e-13);

// Finite choices currently require the physical-metal CN solver, matching
// total-species heat, and volume-face luminosities. Initial-D microscopic
// transport still needs its isotope-aware provider; it is not omitted here.
enum class ConvectiveMixing { instantaneous, finite_implicit, finite_lagged };
struct EvolutionOptions {
  RelaxationOptions relaxation{};
  // Optional starting structure on the same mass mesh. Composition, age and
  // thermal history still come from the accepted previous model.
  const Model* initial_structure_guess{};
  std::size_t max_coupling_iterations{30};
  double abundance_tolerance{1e-12};
  // Optional tighter correction bound when instantaneous convection joins
  // the entire star into one abundance unknown. Zero uses the bound above.
  // This is a solver setting, independent of the global inventory audit.
  double homogeneous_abundance_tolerance{};
  // Optional outer composition bound, independent of the inner species and
  // inventory solves. Zero retains the abundance tolerance above.
  double coupling_stop_tolerance{};
  // Final structure verification after the last composition update. Zero
  // retains the corresponding relaxation tolerance.
  double verification_residual_tolerance{},verification_correction_tolerance{};
  double max_abundance_change{.001};
  // Convergence of total-species enthalpy transport outside the local thermal
  // Jacobian, normalized at each face to max(|L_lo|,|L_hi|,1e-12*max|L|).
  double material_heat_tolerance{1e-10};
  ConvectiveMixing convective_mixing{ConvectiveMixing::instantaneous};
  // In a finite mode, collapse only convective faces with either endpoint
  // cooler than this temperature (K). Zero leaves all faces finite.
  // This selects a mixing approximation, not a microscopic validity limit;
  // exposed interfaces must still satisfy the chosen transport provider.
  double instantaneous_mixing_below_T{};
  // Accepted total face rates from the preceding interval, matching previous.
  // A lagged finite coefficient needs these to recover the previous thermal
  // luminosity. Only an exactly homogeneous initial model can omit them:
  // its macroscopic species flux is zero. Save these with production restarts.
  std::span<const SpeciesVector> previous_species_heat_rates;
  std::span<const std::array<double,3>> previous_metal_heat_rates;
};
// Numerical composition blocks, including single cells. Radiative faces
// always separate blocks, independently of the temperature selection.
MixingRegions instantaneous_mixing_regions(const Model&,const Physics&,const EvolutionOptions&);
struct EvolutionStep {
  Model model; // previous model on failure; never a partially accepted step
  bool converged{};
  std::size_t coupling_iterations{};
  std::string message;
  double residual{},correction{},abundance_residual{};
  double material_heat_residual{};
  // Rates held by the final converged thermal solve; empty for providers that
  // do not request total-species enthalpy. Together with model these reproduce
  // its thermal residuals. The outer residual above bounds their update.
  std::vector<SpeciesVector> total_species_rates;
  // H1, He3 and total Z for an explicitly selected moving-metal provider.
  std::vector<std::array<double,3>> total_metal_species_rates;
  // Relative discrete first-law and nuclear rest-mass audits for this step.
  // Zero nuclear release reports zero mass imbalance (cold no-burning limit).
  double luminosity_balance{},nuclear_mass_balance{},convective_mass_fraction{};
  double nuclear_luminosity{},gravitational_luminosity{},neutrino_luminosity{};
  double thermal_neutrino_luminosity{}; // positive sink; separate from nuclear neutrinos
  std::size_t mixed_regions{};
  // Largest |B|/q^2 at an active convective face; diagnostic, not an error bound.
  double convection_composition_ratio{};
  std::size_t convection_composition_face{};
};
// Converged block iteration of burning/mixing and implicit thermal structure.
// Fixed baryonic mass mesh; the selected convection criterion and optional
// slow/microscopic transport participate in the same update. No remeshing.
EvolutionStep evolve_step(const Model&,const Physics&,const Atmosphere&,double dt,
                          const EvolutionOptions& = {});
} // namespace ember
