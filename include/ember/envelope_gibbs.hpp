#pragma once
#include "ember/envelope_atmosphere.hpp"
#include <memory>

namespace ember {

// Additive Gibbs potentials at common gas pressure, with radiation included once.
// Tables contain G/(R T) on the atomic-mass basis; evaluation converts to baryon
// mass. He3 and trace D use the classical equal-number isotope approximation.
// Optional ideal GS98 metals are neutral or fully ionized; these are declared
// approximations, not an ionization or phase-separation calculation. Component
// additivity omits nonideal H/He mixing and needs an uncertainty assessment.
class GibbsEnvelope final : public EnvelopeThermodynamics {
public:
  struct Files { std::string hydrogen, helium, hydrogen_warm, helium_warm; };
  explicit GibbsEnvelope(const Files&,EnvelopeMetals = EnvelopeMetals::reject);
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
