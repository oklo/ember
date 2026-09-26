#pragma once
#include "ember/cn_burning.hpp"
#include "ember/species_transport.hpp"

namespace ember {
using CNSpeciesFlux=std::function<CNSpeciesFaceResponse(
    std::size_t,const Composition&,const Composition&,bool)>;
struct CNBoundaryFlux {std::size_t face{};CNSpeciesVector rate{};};
struct CNTransportResult {
  std::vector<Composition> composition;
  std::vector<CNBoundaryFlux> boundary_fluxes;
  CNSpeciesVector integrated_balance{};
  double residual{},abundance_correction{};
  std::size_t iterations{};
  std::vector<double> residual_history;
};
// General five-species conservative implicit solve. Local catalyst number
// may vary; fixed-GS98 material slots remain lookup coordinates. The physical
// ledger supplies actual He4/Z and nuclear binding energy. Flux derivatives
// use independent H1, He3, C12, C13, N14 mass fractions at fixed lookup Z.
CNTransportResult burn_cn_and_diffuse(const Model&,const Model&,const PPCNNetwork&,
    const MixingRegions&,const CNSpeciesFlux&,double dt,const SpeciesTransportOptions& = {},
    std::span<const double> common_mixing_rates = {});
// Optional common_mixing_rates contains the linear coefficient (g/s per
// abundance difference) at every mesh face. Keeping this term separate
// preserves the finite storage/reaction equation when mixing is very stiff;
// the callback supplies only the remaining, possibly nonlinear flux.
CNSpeciesVector cn_transport_abundances(const Composition&);
Composition cn_transport_composition(const Composition& lookup,const CNSpeciesVector&);
// Explicit trace approximations for coupling an H/He microscopic provider.
// Helium velocity includes CN in the helium-group mass flux, preserving zero
// total baryonic flux. It omits separate C/N settling and its feedback on the
// provider's collision coefficients. Secular mixing is added separately.
CNSpeciesFaceResponse trace_cn_flux(const SpeciesFaceResponse&,
    const Composition&,const Composition&,CNMicroscopicApproximation,bool derivatives);
} // namespace ember
