#pragma once
#include "ember/composition.hpp"

namespace ember {

// Yakovlev et al. (2006), uniform-mixture phenomenological rates. The high
// choice is a sensitivity scenario, not a rigorous bound for light impurities.
enum class DenseNuclearModel { none, uniform_optimal, uniform_high };
void set_dense_nuclear_model(DenseNuclearModel);
DenseNuclearModel dense_nuclear_model();

struct NuclearPair {
  double charge1{},charge2{},mass1{},mass2{}; // masses in atomic mass units
  double s0{}; // MeV barn; constant low-energy S factor
};
struct NuclearRateResponse {
  double molar_rate{}; // N_A <sigma v>, no symmetry or abundance factors
  double dlnrate_dlnT{},dlnrate_dlnRho{};
  std::array<double,NSPEC> dlnrate_dX{};
};
// Fully ionized background, normalized abundances for physical evaluation.
// Unconstrained abundance partials extend number densities linearly in X.
// This supplies a coefficient only; existing networks retain their exact
// stoichiometry, nuclear mass differences and neutrino accounting.
NuclearRateResponse dense_nuclear_rate(double T,double rho,const Composition&,
                                      const NuclearPair&,DenseNuclearModel);

} // namespace ember
