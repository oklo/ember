#pragma once
#include "ember/envelope_atmosphere.hpp"
#include <memory>

namespace ember {

// Additive Gibbs potentials at common gas pressure, with radiation included once.
// Tables contain G/(R T) on the atomic-mass basis; evaluation converts to baryon
// mass. The He3 term uses the He4 reference at equal number density.
// This implementation is assessed for X>=0.985 and trace metals only. It is not
// a phase-separation model or a substitute for full-composition PMS coverage.
class GibbsEnvelope final : public EnvelopeThermodynamics {
public:
  struct Files { std::string hydrogen, helium, hydrogen_warm, helium_warm; };
  explicit GibbsEnvelope(const Files&);
  ~GibbsEnvelope();
  struct Thermodynamics { State pressure; double entropy, energy; };
  Thermodynamics evaluate(double T, double P, const Composition&) const;
  State at_pressure(double lnT, double lnP, const Composition&, double& guess) const override;
  double rho_from_PT(double T, double P, const Composition&, double rho_guess=0) const override;
private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace ember
