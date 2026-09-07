#pragma once
#include "ember/composition.hpp"

namespace ember {

// Volumetric energy sinks: thermal neutrino emission.
//
// Not needed by a 0.1 Msun star, whose core never gets hot or dense enough for
// it to matter, and so not implemented yet.  It is declared here because it is
// the term that governs the *early* cooling of a massive white dwarf: plasmon
// neutrinos carry away more than photons do until the star has faded well
// below a thousandth of solar, and a cooling age computed without them is
// wrong by a factor of several.  A solar remnant or a Chandrasekhar-mass
// remnant will need it on the first day it is asked for.
struct LossState {
  double eps{};          // erg/g/s, positive is a sink
  double dlneps_dlnT{};
  double dlneps_dlnRho{};
};

class NeutrinoLosses {
public:
  virtual ~NeutrinoLosses() = default;
  virtual LossState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
};

// The explicit "no losses" choice, so that omitting them is a decision in the
// configuration rather than an absence in the code.
class NoNeutrinoLosses final : public NeutrinoLosses {
public:
  LossState eval(double, double, const Composition&) const override { return {}; }
  const char* name() const override { return "none"; }
};

} // namespace ember
