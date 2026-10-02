#pragma once
#include "ember/evolution.hpp"
#include "ember/species_flux.hpp"
#include <functional>
#include <span>

namespace ember {
using SpeciesFlux = std::function<SpeciesFaceResponse(
    std::size_t face, const Composition& left, const Composition& right, bool derivatives)>;
struct SpeciesTransportOptions {
  // Bound the Newton abundance correction; also the default balance bound.
  // A local rate residual alone is not an abundance error in stiff diffusion.
  double abundance_tolerance{1e-13};
  std::size_t max_iterations{60}, max_backtracks{50};
  // Optional cell-wise Newton guess, mass-averaged within each mixed region.
  // This never replaces previous abundances or changes the conservative RHS.
  // It permits an interior guess when old cells contain an absent species.
  std::span<const Composition> initial_guess{};
  // For a connected transport region requiring positive new-time fractions,
  // seed a locally absent but globally present species by a conservative
  // mixture with the region-wide mean. Globally absent species remain zero.
  bool seed_present_species{false};
  // CN source regions and microscopic faces may be evaluated independently.
  // Opt-in parallel callers must supply a thread-safe const flux callback.
  // Matrix assembly and conservation sums retain their serial order.
  std::size_t evaluation_threads{1};
  // Independent bound on the integrated species balance. Zero uses
  // abundance_tolerance. Local correction accuracy need not equal the
  // accuracy of a mass-weighted conservation sum.
  double integrated_balance_tolerance{};
  // Optional bound on the binding energy of the integrated species residual,
  // in erg/g over the step. Representable abundance roundoff is added to it.
  double integrated_binding_tolerance{};
};
struct SpeciesBoundaryFlux {std::size_t face{};SpeciesVector rate{};};
struct SpeciesTransportResult {
  std::vector<Composition> composition;
  // Only region boundaries: internal microscopic fluxes cancel in a mixed
  // region's species balance and must be evaluated separately for heat flow.
  std::vector<SpeciesBoundaryFlux> boundary_fluxes;
  SpeciesVector integrated_balance{}; // delta mass fraction minus integrated nuclear source
  double residual{};
  double abundance_correction{};
  std::size_t iterations{};
  std::vector<double> residual_history;
  bool interior_guess_adjusted{};
};

struct SpeciesFluxReconstruction {
  // Total outward H1/He3 rates, g/s, on every interior face. Within a mixed
  // region these include the redistribution required by homogeneous mixing.
  // At region boundaries the supplied rates are retained exactly.
  std::vector<SpeciesVector> face_rates;
  // Closed exterior fluxes are used. Balances are normalized to stellar mass:
  // w/M * (Xnew-Xold-dt*source) + dt/M * (rate_out-rate_in).
  // An inconsistent region is reported, never repaired by changing its source
  // or boundary flux. Its accumulated mismatch occurs in its final cell.
  std::vector<SpeciesVector> cell_balances, region_balances;
  SpeciesVector integrated_balance{};
};

// Reconstruct the total material redistribution implied by a converged
// backward-Euler abundance update. New abundances must be exactly homogeneous
// in each supplied mixed region; previous abundances need not be homogeneous.
// Nuclear sources are evaluated at the new thermal/composition state. The
// boundaries must match the partition and be ordered from center to surface.
// Subtract a separately evaluated microscopic rate to obtain the convective
// redistribution. This does not prescribe an energy flux or change a model.
SpeciesFluxReconstruction reconstruct_species_fluxes(const Model& current,
    const Model& previous,const Nuclear&,const MixingRegions&,
    std::span<const SpeciesBoundaryFlux>,double dt);

// Solve E_i x_i + f_i-f_{i-1}=b_i, f_i=A_i x_i-B_i x_{i+1}+d_i.
// Closed ends. General two-component Jacobians are allowed. The elimination
// retains E explicitly, without subtracting large diffusive diagonals.
std::vector<SpeciesVector> solve_species_flux_chain(std::span<const SpeciesMatrix> E,
    std::span<const SpeciesMatrix> A,std::span<const SpeciesMatrix> B,
    std::span<const SpeciesVector> b,std::span<const SpeciesVector> d);

// Simultaneous backward-Euler burning and microscopic flux at fixed thermal
// structure. Convective regions are constrained to one homogeneous composition
// with exact nodal mass weights. All metal carriers remain fixed. No abundance
// clipping or renormalization. An implicit burn/mix predictor supplies a guess,
// not the conservative right side. The flux callback must be defined at every
// accepted/initial candidate, and reevaluates the complete trial flux.
// Convergence requires both a small Newton correction and a small integrated
// balance. The unpreconditioned local residual is retained as a diagnostic;
// stiff face terms can amplify roundoff in nearly equal chemical potentials.
// This does not update thermal energy or age: couple to the structure solve
// with the full material-energy flux and original previous model.
SpeciesTransportResult burn_and_diffuse(const Model& thermal,const Model& previous,
    const Nuclear&,const MixingRegions&,const SpeciesFlux&,double dt,
    const SpeciesTransportOptions& = {});
} // namespace ember
