#pragma once
#include "ember/nuclear.hpp"
#include "ember/deuterium.hpp"

namespace ember {
// Actual catalyst molalities, in mol/g of baryon mass. These are separate
// from Composition's fixed GS98 metal slots, which are not isotope abundances.
using CNAbundances = std::array<double,3>; // C12, C13, N14
CNAbundances initial_gs98_cn(const Composition& lookup);

struct CNNetworkResponse {
  // Physical H1/He4/C12/C13/N14 mass-fraction sources. d_dXdt_dX and
  // deps_dX are partials in the independent LOOKUP composition at fixed Y.
  // Do not apply the CNO entries directly to a fixed-mixture EOS composition.
  NuclearResponse physical;
  std::array<double,3> frequency{}; // proton captures per catalyst per second
  std::array<double,3> reaction_rate{}; // mol/g/s
  std::array<double,3> deps_dY{};
  std::array<std::array<double,3>,NSPEC> d_dXdt_dY{};
};

// C12(p,gamma)N13(beta+)C13, C13(p,gamma)N14, and
// N14(p,gamma)O15(beta+)N15(p,alpha)C12. Short beta decays and N15 capture
// are eliminated; C12, C13 and N14 remain time dependent. Oxygen branches
// and catalyst diffusion are not included. Screening uses the supplied bulk
// composition; the error of a fixed-mixture screening lookup is separate.
class CNNetwork {
public:
  explicit CNNetwork(PPRates=PPRates::solar_fusion_iii,
                     PPScreening=PPScreening::salpeter_van_horn);
  CNNetworkResponse response(double T,double rho,const Composition& lookup,
                             const CNAbundances&) const;
private:
  PPRates rates_;
  PPScreening screening_;
};

// Initial deuterium, pp chains and the local, time-dependent CN inventory.
// Only the primordial D inventory is explicit; pp-created D remains folded
// into the reduced pp reactions. The returned
// dXdt/derivatives describe the lookup X slots at fixed cn_molality. CN isotope
// derivatives and updates are supplied by CNNetwork and burn_cn_and_mix.
class PPCNNetwork final : public Nuclear {
public:
  explicit PPCNNetwork(PPRates pp_rates=PPRates::solar_fusion_ii,
      PPScreening screening=PPScreening::salpeter_van_horn,
      PPRates cn_rates=PPRates::solar_fusion_iii):pp_(pp_rates,screening),cn_(cn_rates,screening) {}
  NuclearState eval(double,double,const Composition&) const override;
  NuclearResponse composition_response(double,double,const Composition&) const override;
  double rest_energy_correction(const Composition&) const override;
  const char* name() const override {return "initial deuterium, reduced pp and time-dependent CN network";}
  const PPChains& pp() const {return pp_.pp();}
  const PPDeuterium& light() const {return pp_;}
  const CNNetwork& cn() const {return cn_;}
private:
  NuclearResponse response(double,double,const Composition&,bool initial_deuterium) const;
  PPDeuterium pp_;
  CNNetwork cn_;
};

struct CNPhysicalLedger {
  CNAbundances molality{};
  double hydrogen{},helium3{},helium4{},metal_fraction{};
  double extra_metal_mass{}; // physical metal mass minus lookup Z
  double ion_molality_difference{},electron_molality_difference{};
  // Physical nuclear rest energy minus the fixed-mixture He4 proxy, erg/g.
  // Its time difference must enter any nuclear mass-defect audit using lookup X.
  double rest_energy_difference{};
};
CNPhysicalLedger cn_physical_ledger(const Composition& lookup,const CNAbundances&);
// Remaining metal mass after CN; allow only relative floating-point roundoff
// when mass fractions and catalyst molalities describe the same zero remainder.
double cn_inert_metal_fraction(double total_metals,const CNAbundances&);
// Re-express an existing isotope inventory with actual total metal and He4
// mass fractions. Preserves H, He3, CN numbers and inert metal mass. This is
// an explicit material-model conversion, not an evolutionary timestep.
Composition explicit_cn_material(const Composition&);

// Conservative backward-Euler catalyst update at specified end-of-step
// capture frequencies. Can be eliminated inside a hydrogen/composition solve;
// this function alone does not update fuel, heat, screening or stellar structure.
CNAbundances cn_backward_euler(const CNAbundances& old,
                             const std::array<double,3>& frequency,double dt);
} // namespace ember
