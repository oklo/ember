#pragma once
#include "ember/evolution.hpp"
#include "ember/nuclear_cn.hpp"

namespace ember {
struct CNBurnResult {
  std::vector<Composition> lookup;
  std::vector<CNAbundances> catalysts;
  std::size_t maximum_iterations{};
  double maximum_equation_residual{};
  double nuclear_luminosity{},neutrino_luminosity{};
  double nuclear_mass_balance{};
};
// Backward-Euler pp/CN burning with simultaneous instantaneous convection.
// The independent variables are H1, He3, C12 and C13; N14 follows from
// conserved catalyst number. Lookup He4 closes the fixed GS98 composition,
// while cn_physical_ledger supplies actual He4, metals and binding energy.
// No secular or microscopic diffusion, nor thermal/structural update here.
CNBurnResult burn_cn_and_mix(const Model& thermal,const Model& previous,
    const std::vector<CNAbundances>& old_catalysts,const PPChains&,const CNNetwork&,
    const MixingRegions&,double dt,double tolerance=1e-12);
// Simultaneous burning, convection and finite secular mixing. D is the
// common diffusivity of every isotope at each internal face; fixed uniform
// total catalyst number is preserved. Microscopic drift is not included.
std::vector<Composition> burn_cn_and_transport(const Model& thermal,const Model& previous,
    const PPCNNetwork&,const MixingRegions&,const std::vector<double>& D,
    double dt,double tolerance=1e-12);
} // namespace ember
