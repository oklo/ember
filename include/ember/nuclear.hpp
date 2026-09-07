#pragma once
#include "ember/composition.hpp"
#include <array>

namespace ember {

// Nuclear energy generation and the composition change that produces it.
struct NuclearState {
  double eps{};                       // erg/g/s
  double dlneps_dlnT{};
  double dlneps_dlnRho{};
  std::array<double, NSPEC> dXdt{};   // 1/s, mass fractions
};

class Nuclear {
public:
  virtual ~Nuclear() = default;
  virtual NuclearState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
};

// The proton-proton chains, which are the whole of the energy budget for a
// star of a tenth of a solar mass.
//
// He3 is followed explicitly rather than assumed to be in equilibrium.  In a
// fully convective star of this mass the He3 abundance rises for a trillion
// years before it is consumed, and the resulting change in mean molecular
// weight expands the core and drives a good part of the evolution: assuming
// equilibrium removes the phenomenon rather than approximating it.
//
// Rates follow the Adelberger et al. (2011) compilation of S-factors in the
// standard non-resonant form.  Screening is the weak-limit Salpeter factor,
// which is what these densities call for; the strong-screening and
// pycnonuclear regimes a cold dense remnant would need are a separate
// implementation of this same interface.
class PPChains final : public Nuclear {
public:
  NuclearState eval(double T, double rho, const Composition&) const override;
  const char* name() const override { return "pp chains (Adelberger+2011)"; }
};

} // namespace ember
