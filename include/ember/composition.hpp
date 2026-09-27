#pragma once
#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string_view>
#include "ember/gs98_mixture.hpp"

namespace ember {

// Species slots for isotope mixtures. With MetalInventory::gs98, the metal
// slots instead carry a fixed elemental table mixture; actual CN isotopes
// are kept separately below. CN equilibration is time dependent.
enum class Species : std::size_t {
  H1 = 0, He3, He4, C12, C13, N14, O16, Zrest, H2, COUNT
};
inline constexpr std::size_t NSPEC = static_cast<std::size_t>(Species::COUNT);
// Historical material tables and checkpoint versions 1--3 have eight slots.
// Appending D preserves every existing isotope index and the five metal slots.
inline constexpr std::size_t TABLE_NSPEC = 8;
inline constexpr std::size_t METAL_BEGIN = 3, METAL_END = 8;
inline constexpr std::size_t NMETALS = METAL_END-METAL_BEGIN;
constexpr bool is_metal_species(std::size_t i) { return i>=METAL_BEGIN && i<METAL_END; }
inline constexpr std::array<double, NSPEC> mass_numbers{1,3,4,12,13,14,16,20,2};
enum class AbundanceBasis { atomic_mass, baryon_mass };
enum class MetalInventory { carried_isotopes, gs98 };
enum class CNMassConvention { fixed_metal_proxy, explicit_metal_mass };

struct Nuclide { double A; double Z; std::string_view name; };

inline constexpr std::array<Nuclide, NSPEC> nuclides{{
  {1.00782503,  1.0, "h1"},   {3.01602932,  2.0, "he3"},
  {4.00260325,  2.0, "he4"},  {12.0000000,  6.0, "c12"},
  {13.00335484, 6.0, "c13"},  {14.0030740,  7.0, "n14"},
  {15.9949146,  8.0, "o16"},  {20.0,       10.0, "zrest"},
  {2.01410177812,1.0,"h2"}
}};

// Mass fractions.  Kept as a plain aggregate so it lives in registers and
// copies freely; a stellar model carries one per zone.
struct Composition {
  std::array<double, NSPEC> X{};
  // Legacy static tables use atomic mass fractions. Evolution uses X_i=A_i Y_i
  // with integer mass numbers, so sum(X)=1 conserves baryons; nuclear rest
  // mass changes are accounted for separately in the released energy.
  AbundanceBasis basis{AbundanceBasis::atomic_mass};
  // GS98 selects the fixed elemental distribution used by atmosphere/EOS
  // sources. The five inert metal slots continue to carry conserved mass;
  // individual CNO reactions must not interpret those proxy labels literally.
  MetalInventory metal_inventory{MetalInventory::carried_isotopes};
  // Optional actual C12/C13/N14 molalities (mol/g of baryonic mass).
  // X remains the fixed-GS98 EOS/opacity lookup. cn_physical_ledger supplies
  // the physical helium/metal fractions and the nuclear binding correction.
  // Inactive by default; activation requires a compatible evolution network.
  std::optional<std::array<double,3>> cn_molality;
  // Explicit mode: X[He4] and Z() are the actual total mass fractions; the
  // GS98 pattern is only the material approximation within that metal mass.
  CNMassConvention cn_mass_convention{CNMassConvention::fixed_metal_proxy};
  constexpr bool operator==(const Composition&) const = default;
  constexpr double abundance_weight(std::size_t i) const {
    return basis == AbundanceBasis::baryon_mass ? mass_numbers[i] : nuclides[i].A;
  }

  constexpr double operator[](Species s) const { return X[static_cast<std::size_t>(s)]; }
  constexpr double& operator[](Species s)      { return X[static_cast<std::size_t>(s)]; }

  constexpr double sum() const {
    double s = 0.0; for (double v : X) s += v; return s;
  }
  constexpr double h1()  const { return (*this)[Species::H1]; }
  constexpr double he4() const { return (*this)[Species::He4]; }
  // Metallicity: everything heavier than helium.
  constexpr double Z() const {
    double z = 0.0;
    for (std::size_t i = METAL_BEGIN; i < METAL_END; ++i) z += X[i];
    return z;
  }
  constexpr double Y() const { return (*this)[Species::He3] + (*this)[Species::He4]; }
  constexpr double metal_ion_moment(int power) const {
    return gs98_ion_moment(power)/(basis==AbundanceBasis::baryon_mass?1.:gs98_atomic_mass_scale());
  }

  // Moles of ions and of electrons per gram, for the fully ionised mixture.
  constexpr double mu_ions_inv() const {
    double s = 0.0;
    for (std::size_t i = 0; i < NSPEC; ++i)
      if(metal_inventory!=MetalInventory::gs98 || !is_metal_species(i)) s += X[i] / abundance_weight(i);
    if(metal_inventory==MetalInventory::gs98)s+=Z()*metal_ion_moment(0);
    return s;
  }
  constexpr double mu_elec_inv() const {
    double s = 0.0;
    for (std::size_t i = 0; i < NSPEC; ++i)
      if(metal_inventory!=MetalInventory::gs98 || !is_metal_species(i)) s += X[i] * nuclides[i].Z / abundance_weight(i);
    if(metal_inventory==MetalInventory::gs98)s+=Z()*metal_ion_moment(1);
    return s;
  }
};

// Asplund, Amarsi & Grevesse (2021) photospheric mixture, the current
// consensus solar composition, scaled to a requested (X, Z).  The metal
// pattern within Z is AAG21; only the four CNO isotopes are resolved.
Composition solar_scaled(double X, double Z);

// A face composition must average both the material lookup and actual CN
// inventory. Mixing conventions here would give a different physical state.
inline Composition mean_composition(const Composition& a,const Composition& b) {
  if(a.basis!=b.basis || a.metal_inventory!=b.metal_inventory ||
      a.cn_mass_convention!=b.cn_mass_convention ||
      a.cn_molality.has_value()!=b.cn_molality.has_value())
    throw std::domain_error("face composition: incompatible abundance conventions");
  auto c=a;
  for(std::size_t k=0;k<NSPEC;++k)c.X[k]=.5*(a.X[k]+b.X[k]);
  if(c.cn_molality)for(std::size_t k=0;k<3;++k)
    (*c.cn_molality)[k]=.5*((*a.cn_molality)[k]+(*b.cn_molality)[k]);
  return c;
}

} // namespace ember
