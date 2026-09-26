#pragma once
#include "ember/nuclear.hpp"

namespace ember {

// Initial deuterium is a separate fuel inventory. The reduced pp chain already
// includes capture of its freshly produced, short-lived deuterium; never apply
// this source to that same material a second time.
inline constexpr double deuterium_atomic_mass = nuclides[static_cast<std::size_t>(Species::H2)].A;

// Solar Fusion III, Table IV: linear interpolation of the recommended S12(E).
// The integral is restricted to the published 0--120 keV interval, and the
// thermal domain to at most 20 MK. Below 0.01 MK the negligible thermal
// capture source is explicitly zero; this is not a pycnonuclear model.
// Returns the bare molar capture rate.
ThermonuclearRate deuterium_bare_rate(double T);

struct DeuteriumCapture {
  NuclearState source; // dXdt contains the physical H1 sink and He3 source
  double dX_deuterium_dt{};
  double molar_reactions_per_gram_second{};
  double heat_per_mole{};
};

// `other` contains true proton, helium and metal baryon fractions EXCLUDING
// X_deuterium, so other.sum()+X_deuterium=1. It is not an EOS composition.
// This local source is not yet attached to the stellar composition solver.
DeuteriumCapture deuterium_capture(double T, double rho,
    const Composition& other, double X_deuterium,
    PPScreening screening=PPScreening::salpeter_van_horn);

// Explicit initial D plus the existing reduced pp chain. H1 is always the
// true proton mass fraction, H2 the separate initial-D inventory. Deuterium
// formed by pp remains the already folded short-lived intermediate.
class PPDeuterium final : public Nuclear {
 public:
  explicit PPDeuterium(PPRates rates=PPRates::solar_fusion_iii,
      PPScreening screening=PPScreening::salpeter_van_horn):pp_(rates,screening),screening_(screening) {}
  NuclearState eval(double,double,const Composition&) const override;
  NuclearResponse composition_response(double,double,const Composition&) const override;
  const char* name() const override {return "explicit initial deuterium plus reduced pp chains";}
  const PPChains& pp() const {return pp_;}
 private:
  PPChains pp_;
  PPScreening screening_;
};

} // namespace ember
