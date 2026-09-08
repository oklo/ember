#pragma once
#include "ember/composition.hpp"
#include <array>
#include <stdexcept>

namespace ember {

// Nuclear energy generation and the composition change that produces it.
struct NuclearState {
  double eps{};                       // erg/g/s
  double dlneps_dlnT{};
  double dlneps_dlnRho{};
  std::array<double, NSPEC> dXdt{};   // 1/s, in the input composition's abundance basis
  double eps_neutrino{};             // escaping nuclear neutrinos, erg/g/s
};

struct NuclearResponse {
  NuclearState state;
  std::array<std::array<double,NSPEC>,NSPEC> d_dXdt_dX{};
  std::array<double,NSPEC> deps_dX{};
};

class Nuclear {
public:
  virtual ~Nuclear() = default;
  virtual NuclearState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  virtual NuclearResponse composition_response(double,double,const Composition&) const {
    throw std::logic_error("Nuclear: composition derivatives unavailable");
  }
};

// The proton-proton chains, which are the whole of the energy budget for a
// star of a tenth of a solar mass.
//
// He3 is followed explicitly rather than assumed to be in equilibrium.  In a
// fully convective star of this mass the He3 abundance rises for a trillion
// years before it reaches equilibrium, and the resulting change in mean molecular
// weight expands the core and drives a good part of the evolution: assuming
// equilibrium removes the phenomenon rather than approximating it.
//
// Rates follow the Adelberger et al. (2011) compilation of S-factors in the
// standard non-resonant form. Screening is the classical weak-limit Salpeter
// factor with a retained exp(2) cap, an approximation requiring further
// qualification in partially degenerate/intermediate-coupling matter. The strong-screening and
// pycnonuclear regimes a cold dense remnant would need are a separate
// implementation of this same interface.
class PPChains final : public Nuclear {
public:
  NuclearState eval(double T, double rho, const Composition&) const override;
  NuclearResponse composition_response(double T,double rho,const Composition&) const override;
  const char* name() const override { return "pp chains (Adelberger+2011)"; }
};

} // namespace ember
