#pragma once
#include <array>
#include <cstddef>
#include <string_view>
#include "ember/gs98_mixture.hpp"

namespace ember {

// Isotopes carried explicitly.  The pp chain needs H1/He3/He4; the CNO cycle
// needs C12/C13/N14/O16 out of equilibrium, because in a low-mass star the
// approach to CN equilibrium happens on the pre-main-sequence and never
// repeats.  Everything else in Z is inert ballast that contributes to the mean
// molecular weight and to the opacity's metallicity, and is carried as one
// lumped species.
enum class Species : std::size_t {
  H1 = 0, He3, He4, C12, C13, N14, O16, Zrest, COUNT
};
inline constexpr std::size_t NSPEC = static_cast<std::size_t>(Species::COUNT);
inline constexpr std::array<double, NSPEC> mass_numbers{1,3,4,12,13,14,16,20};
enum class AbundanceBasis { atomic_mass, baryon_mass };
enum class MetalInventory { carried_isotopes, gs98 };

struct Nuclide { double A; double Z; std::string_view name; };

inline constexpr std::array<Nuclide, NSPEC> nuclides{{
  {1.00782503,  1.0, "h1"},   {3.01602932,  2.0, "he3"},
  {4.00260325,  2.0, "he4"},  {12.0000000,  6.0, "c12"},
  {13.00335484, 6.0, "c13"},  {14.0030740,  7.0, "n14"},
  {15.9949146,  8.0, "o16"},  {20.0,       10.0, "zrest"}
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
    for (std::size_t i = static_cast<std::size_t>(Species::C12); i < NSPEC; ++i) z += X[i];
    return z;
  }
  constexpr double Y() const { return (*this)[Species::He3] + (*this)[Species::He4]; }
  constexpr double metal_ion_moment(int power) const {
    return gs98_ion_moment(power)/(basis==AbundanceBasis::baryon_mass?1.:gs98_atomic_mass_scale());
  }

  // Moles of ions and of electrons per gram, for the fully ionised mixture.
  constexpr double mu_ions_inv() const {
    double s = 0.0;
    for (std::size_t i = 0; i < (metal_inventory==MetalInventory::gs98?3:NSPEC); ++i) s += X[i] / abundance_weight(i);
    if(metal_inventory==MetalInventory::gs98)s+=Z()*metal_ion_moment(0);
    return s;
  }
  constexpr double mu_elec_inv() const {
    double s = 0.0;
    for (std::size_t i = 0; i < (metal_inventory==MetalInventory::gs98?3:NSPEC); ++i) s += X[i] * nuclides[i].Z / abundance_weight(i);
    if(metal_inventory==MetalInventory::gs98)s+=Z()*metal_ion_moment(1);
    return s;
  }
};

// Asplund, Amarsi & Grevesse (2021) photospheric mixture, the current
// consensus solar composition, scaled to a requested (X, Z).  The metal
// pattern within Z is AAG21; only the four CNO isotopes are resolved.
Composition solar_scaled(double X, double Z);

} // namespace ember
