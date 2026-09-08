#pragma once
#include "ember/composition.hpp"
#include <stdexcept>
#include <optional>
#include <array>

namespace ember {

// Thermodynamic state at (T, rho).  Everything a Henyey solver needs, with
// analytic derivatives: the F77 ancestor formed these by centred differences,
// five equation-of-state calls per zone per Newton iteration, which cost both
// accuracy and a factor of five in time.
struct EosState {
  double P{};        // total pressure, gas + radiation           [dyn/cm^2]
  double E{};        // specific internal energy, if available     [erg/g]
  double S{};        // specific entropy                          [erg/g/K]
  double chiT{};     // dlnP/dlnT at constant rho
  double chiRho{};   // dlnP/dlnrho at constant T
  double cv{};       // T*dS/dT at constant rho (=dE/dT if E supplied) [erg/g/K]
  double cp{};       // specific heat at constant P               [erg/g/K]
  double grad_ad{};  // dlnT/dlnP at constant S
  double Gamma1{};   // dlnP/dlnrho at constant S
  double delta{};    // -dlnrho/dlnT at constant P
  double mu{};       // mean molecular weight
  double free_e{};   // free electrons per nucleon

  // Derived, so callers never re-derive them inconsistently.
  double dPdT_rho(double T)     const { return P * chiT / T; }
  double dPdRho_T(double rho)   const { return P * chiRho / rho; }
};

// Additional state derivatives needed by an analytic structure Jacobian.
// Ordinary eval() remains available without computing second derivatives.
struct EosResponse {
  EosState state{};
  double dE_dlnRho{};
  double dcp_dlnT{}, dcp_dlnRho{};
  double ddelta_dlnT{}, ddelta_dlnRho{};
  double dgrad_ad_dlnT{}, dgrad_ad_dlnRho{};
};

// Composition directions replace He4 by H1 or He3, at fixed T, rho and metals.
struct EosCompositionResponse {
  std::array<double,2> dP{}, dE{};
};

class Eos {
public:
  virtual ~Eos() = default;
  virtual EosState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  // A pressure/entropy table can support static structure without a reliable
  // caloric energy. Such implementations return false and E=NaN; time-
  // dependent structure must reject them before using E or its derivatives.
  virtual bool has_internal_energy() const { return true; }
  struct DensityRange { double min, max; };
  virtual std::optional<DensityRange> density_range(double, const Composition&) const { return {}; }
  virtual EosResponse eval_with_derivatives(double, double, const Composition&) const {
    throw std::logic_error("Eos: analytic transport derivatives are not implemented");
  }
  virtual EosCompositionResponse composition_response(double, double, const Composition&) const {
    throw std::logic_error("Eos: composition responses are not implemented");
  }

  // Invert to density at given (T, P).  Needed by the atmosphere, which
  // integrates in pressure.  Newton on ln rho using the analytic chiRho, so it
  // converges quadratically instead of the fixed-point crawl the F77 used.
  // Invalid states and failure to converge throw; no approximate density is
  // returned as if the pressure constraint had been satisfied.
  virtual double rho_from_PT(double T, double P, const Composition&,
                     double rho_guess = 0.0) const;
};

// Ideal gas + radiation + analytic Fermi-Dirac electron degeneracy.
//
// This is the fallback and the sanity check, not the production equation of
// state: it has no partial ionisation, no H2, and no Coulomb corrections.  Its
// job is to be exactly right in the limits (ideal, fully degenerate,
// radiation-dominated) so the table-based equations of state can be tested
// against it where they should agree.
class IdealEos final : public Eos {
public:
  EosState eval(double T, double rho, const Composition&) const override;
  EosResponse eval_with_derivatives(double T, double rho, const Composition&) const override;
  const char* name() const override { return "ideal+rad+degeneracy"; }
};

} // namespace ember
