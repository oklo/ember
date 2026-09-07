#pragma once
#include "ember/composition.hpp"
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

namespace ember {

// One additive contribution to the equation of state.
//
// The thermodynamics of a star is a sum: ions, electrons, radiation, and - as
// the material gets cold and dense - the Coulomb energy of the ion lattice,
// the latent heat of crystallisation, and phase separation of a carbon-oxygen
// mixture.  Those last three are what govern how a white dwarf cools over the
// very long run, and they are the reason this is a sum of components rather
// than one closed formula: adding them must be an addition, not a rewrite.
//
// Each term reports its own analytic derivatives.  Pressures and energies add;
// so do their derivatives; chi_T, Gamma_1 and the rest are then formed once,
// from the totals, where they cannot be made inconsistent term by term.
struct EosTerm {
  double P{};          // dyn/cm^2
  double E{};          // erg/g
  double S{};          // erg/g/K   (optional; used for cooling budgets)
  double dP_dlnT{};    // at fixed rho
  double dP_dlnRho{};  // at fixed T
  double dE_dlnT{};    // at fixed rho
  double dE_dlnRho{};  // at fixed T
  // Populated by eval_with_derivatives(). These add before cp, delta, and
  // grad_ad are differentiated; derivatives of those ratios do not add.
  double d2P_dlnT2{}, d2P_dlnTdlnRho{}, d2P_dlnRho2{};
  double d2E_dlnT2{}, d2E_dlnTdlnRho{};

  EosTerm& operator+=(const EosTerm& o) {
    P += o.P; E += o.E; S += o.S;
    dP_dlnT += o.dP_dlnT; dP_dlnRho += o.dP_dlnRho;
    dE_dlnT += o.dE_dlnT; dE_dlnRho += o.dE_dlnRho;
    d2P_dlnT2 += o.d2P_dlnT2; d2P_dlnTdlnRho += o.d2P_dlnTdlnRho;
    d2P_dlnRho2 += o.d2P_dlnRho2;
    d2E_dlnT2 += o.d2E_dlnT2; d2E_dlnTdlnRho += o.d2E_dlnTdlnRho;
    return *this;
  }
};

class EosComponent {
public:
  virtual ~EosComponent() = default;
  virtual EosTerm eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  virtual EosTerm eval_with_derivatives(double, double, const Composition&) const {
    throw std::logic_error("EosComponent: second derivatives are not implemented");
  }
};

// Ideal non-degenerate ions.  Replaced by a Coulomb-corrected version once the
// one-component-plasma terms go in; the interface does not change.
class IonGas final : public EosComponent {
public:
  EosTerm eval(double T, double rho, const Composition&) const override;
  EosTerm eval_with_derivatives(double T, double rho, const Composition& c) const override {
    return eval(T, rho, c);
  }
  const char* name() const override { return "ions (ideal)"; }
};

class Radiation final : public EosComponent {
public:
  EosTerm eval(double T, double rho, const Composition&) const override;
  EosTerm eval_with_derivatives(double T, double rho, const Composition& c) const override {
    return eval(T, rho, c);
  }
  const char* name() const override { return "radiation"; }
};

// Electrons of arbitrary degeneracy and arbitrary relativity.  This one is
// already good enough for a Chandrasekhar-mass white dwarf; what such a star
// additionally needs is the ion physics above, not better electrons.
class ElectronGas final : public EosComponent {
public:
  EosTerm eval(double T, double rho, const Composition&) const override;
  EosTerm eval_with_derivatives(double T, double rho, const Composition&) const override;
  const char* name() const override { return "electrons (relativistic FD)"; }
};

} // namespace ember
