#pragma once
#include "ember/relaxation.hpp"
#include <utility>

namespace ember {
using MixingRegions = std::vector<std::pair<std::size_t,std::size_t>>; // [first,last), includes isolated points
std::vector<double> nodal_mass_weights(const Model&);
MixingRegions schwarzschild_mixing_regions(const Model&,const Physics&);
MixingRegions convective_mixing_regions(const Model&,const Physics&);
// Implicit pp burning with instantaneous homogeneous mixing within each
// supplied connected region. Isolated points burn locally. Metals are inert.
std::vector<Composition> burn_and_mix(const Model& thermal,const Model& previous,
    const Nuclear&,const MixingRegions&,double dt,double tolerance=1e-13);
// Slow mixing coefficients at mesh faces, cm^2/s. Diffusive heat transport
// remains radiative+conductive in these Ledoux-stable layers.
std::vector<double> secular_mixing_diffusivities(const Model&,const Physics&);
// Conservative backward-Euler diffusion and burning, with instantaneous
// convection collapsed into mass-weighted regions. Zero exterior flux.
std::vector<Composition> burn_and_transport(const Model& thermal,const Model& previous,
    const Nuclear&,const MixingRegions&,const std::vector<double>& diffusivity,
    double dt,double tolerance=1e-13);

struct EvolutionOptions {
  RelaxationOptions relaxation{};
  std::size_t max_coupling_iterations{30};
  double abundance_tolerance{1e-12};
  double max_abundance_change{.001};
};
struct EvolutionStep {
  Model model; // previous model on failure; never a partially accepted step
  bool converged{};
  std::size_t coupling_iterations{};
  std::string message;
  double residual{},correction{},abundance_residual{};
  // Relative discrete first-law and nuclear rest-mass audits for this step.
  // Zero nuclear release reports zero mass imbalance (cold no-burning limit).
  double luminosity_balance{},nuclear_mass_balance{},convective_mass_fraction{};
  double nuclear_luminosity{},gravitational_luminosity{},neutrino_luminosity{};
  std::size_t mixed_regions{};
};
// Converged block iteration of burning/mixing and implicit thermal structure.
// Fixed baryonic mass mesh, Schwarzschild boundaries, no remeshing/diffusion.
EvolutionStep evolve_step(const Model&,const Physics&,const Atmosphere&,double dt,
                          const EvolutionOptions& = {});
} // namespace ember
