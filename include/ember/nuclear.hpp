#pragma once
#include "ember/composition.hpp"
#include "ember/dense_nuclear.hpp"
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
  // Physical nuclear rest energy minus that represented by the X slots,
  // erg/g. Its change belongs in the rest-mass audit, not as an additional
  // heat source: eval already supplies the physical nuclear energy release.
  virtual double rest_energy_correction(const Composition&) const { return 0; }
  virtual NuclearResponse composition_response(double,double,const Composition&) const {
    throw std::logic_error("Nuclear: composition derivatives unavailable");
  }
};

// The reduced proton-proton network. Missing branches and CNO burning must
// be assessed separately before interpreting a complete stellar trajectory.
//
// He3 is followed explicitly rather than assumed to be in equilibrium.  In a
// fully convective star of this mass the He3 abundance rises for a trillion
// years before it reaches equilibrium, and the resulting change in mean molecular
// weight expands the core and drives a good part of the evolution: assuming
// equilibrium removes the phenomenon rather than approximating it.
//
// The legacy default preserves historical static benchmarks. Evolution selects
// Solar Fusion II S-factor quadrature and finite-degeneracy Salpeter--Van Horn
// screening explicitly. These remain a reduced pp network, without pep, hep,
// ppIII or CNO, and are not a prescription for pycnonuclear burning.
enum class PPRates { legacy, solar_fusion_ii, solar_fusion_iii };
enum class PPScreening { legacy_weak, debye_fermi, salpeter_van_horn };
enum class PPReaction { pp, he3_he3, he3_he4, deuterium_p };
struct ThermonuclearRate {
  double molar_rate{}; // N_A <sigma v>, cm^3 mol^-1 s^-1; no symmetry factor
  double dlnrate_dlnT{};
};
struct ScreeningState {
  double log_factor{}, dlog_dlnT{}, dlog_dlnRho{};
  std::array<double,NSPEC> dlog_dX{}; // unconstrained abundance partials
  double electron_eta{}, electron_susceptibility{}; // kT/ne * dne/dmu
  double gamma_e{}, zeta{}; // electron-sphere coupling; 3 Gamma_12/tau
  bool reaction_omitted{}; // explicit trace-fuel approximation; not a screening fit
};
ThermonuclearRate pp_bare_rate(double T, PPReaction, PPRates);
// Combined coefficients avoid multiplying an underflowed thermal rate by a
// divergent screening factor in the quantum limit. Dense model is opt-in.
NuclearRateResponse pp_rate_response(double T,double rho,const Composition&,PPReaction,PPRates,PPScreening);
ScreeningState pp_screening(double T,double rho,const Composition&,PPReaction,PPScreening);
// Optional process-wide first-order reuse of the electron screening
// susceptibility on a grid of spacing h in (ln T, ln n_e); h=0 is exact.
void set_screening_reuse(double h);
// Optional Chugunov & DeWitt (2009, PRC 80, 014611; eqs. 23-25, A4) quantum-tunnelling
// correction added to the Salpeter--Van Horn exponent. zeta_max=0 disables it (default;
// results unchanged). When enabled, the classical zeta<=0.2 stop is replaced by the
// mean-field domain zeta<=zeta_max (<=1.6) and Gamma12<=200; beyond it, an explicit
// thermo-pycnonuclear domain error. The correction vanishes as zeta->0.
void set_quantum_screening(double zeta_max);
// Optional omission of reactions outside the quantum screening domain, only
// when X(H1)+X(D)+X(He3) <= max_fraction. Zero (default) retains the domain error.
// This preserves fuel and omits heat and all rate derivatives together. The
// caller must assess the missing power; a small inventory alone is not a rate.
void set_quantum_burning_fuel_limit(double max_fraction);

// N14(p,gamma) bottleneck of a closed CN cycle. The Solar Fusion III choice
// uses its S(0) with Solar Fusion II's first and second derivatives.
ThermonuclearRate cn_bare_rate(double T, PPRates);
ScreeningState cn_screening(double T,double rho,const Composition&,PPScreening);
// Resolved CN captures; the two-argument overload above retains N14(p,gamma).
enum class CNReaction { c12_p, c13_p, n14_p };
ThermonuclearRate cn_bare_rate(double T,CNReaction,PPRates);
NuclearRateResponse cn_rate_response(double T,double rho,const Composition&,CNReaction,PPRates,PPScreening);
ScreeningState cn_screening(double T,double rho,const Composition&,CNReaction,PPScreening);

// Fixed GS98 catalyst-number approximation; does not convert the inert metal
// slots into literal isotopes or account for prior C-to-N fuel consumption.
// converted_carbon=0 supplies the nitrogen-only comparison, =1 the CN limit.
class CNCycle final : public Nuclear {
public:
  explicit CNCycle(PPRates rates=PPRates::solar_fusion_iii,
                   PPScreening screening=PPScreening::salpeter_van_horn,
                   double converted_carbon=1.);
  NuclearState eval(double,double,const Composition&) const override;
  NuclearResponse composition_response(double,double,const Composition&) const override;
  const char* name() const override { return "closed CN cycle with fixed GS98 catalyst number"; }
private:
  PPRates rates_;
  PPScreening screening_;
  double converted_carbon_;
};

class PPChains final : public Nuclear {
public:
  explicit PPChains(PPRates rates=PPRates::legacy,
                    PPScreening screening=PPScreening::legacy_weak)
      : rates_(rates), screening_(screening) {}
  NuclearState eval(double T, double rho, const Composition&) const override;
  NuclearResponse composition_response(double T,double rho,const Composition&) const override;
  const char* name() const override;
private:
  PPRates rates_;
  PPScreening screening_;
};

class PPCNO final : public Nuclear {
public:
  explicit PPCNO(PPRates pp_rates=PPRates::solar_fusion_ii,
                 PPScreening screening=PPScreening::salpeter_van_horn,
                 PPRates cn_rates=PPRates::solar_fusion_iii,
                 double converted_carbon=1.)
      : pp_(pp_rates,screening), cn_(cn_rates,screening,converted_carbon) {}
  NuclearState eval(double,double,const Composition&) const override;
  NuclearResponse composition_response(double,double,const Composition&) const override;
  const char* name() const override { return "reduced pp network plus closed CN cycle"; }
private:
  PPChains pp_;
  CNCycle cn_;
};

} // namespace ember
