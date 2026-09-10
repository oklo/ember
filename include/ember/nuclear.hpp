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
// The legacy default preserves historical static benchmarks. Evolution selects
// Solar Fusion II S-factor quadrature and finite-degeneracy Salpeter--Van Horn
// screening explicitly. These remain a reduced pp network, without pep, hep,
// ppIII or CNO, and are not a prescription for pycnonuclear burning.
enum class PPRates { legacy, solar_fusion_ii };
enum class PPScreening { legacy_weak, debye_fermi, salpeter_van_horn };
enum class PPReaction { pp, he3_he3, he3_he4 };
struct ThermonuclearRate {
  double molar_rate{}; // N_A <sigma v>, cm^3 mol^-1 s^-1; no symmetry factor
  double dlnrate_dlnT{};
};
struct ScreeningState {
  double log_factor{}, dlog_dlnT{}, dlog_dlnRho{};
  std::array<double,NSPEC> dlog_dX{}; // unconstrained abundance partials
  double electron_eta{}, electron_susceptibility{}; // kT/ne * dne/dmu
  double gamma_e{}, zeta{}; // electron-sphere coupling; 3 Gamma_12/tau
};
ThermonuclearRate pp_bare_rate(double T, PPReaction, PPRates);
ScreeningState pp_screening(double T,double rho,const Composition&,PPReaction,PPScreening);

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

} // namespace ember
