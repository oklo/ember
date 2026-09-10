#include "../examples/stellar_seed.hpp"
#include "ember/atmosphere_composition.hpp"
#include "ember/conduction_table.hpp"
#include "ember/eos_composition.hpp"
#include "ember/evolution.hpp"
#include "ember/opacity_mixture.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace ember;
int main() {
  try {
    const std::string data = EMBER_DATA_DIR;
    CompositionHelmholtzEos eos(data + "/eos/freeeos300_hhe_extended.dat",
                                HelmholtzTableEos::Mixture::allow_documented_proxy);
    auto rad = std::make_shared<StellarMixtureOpacity>(data + "/opacity");
    TabulatedConduction cond(data + "/conduction/condtab21wd.dat");
    CombinedOpacity opacity(rad, std::make_shared<HotConduction>(cond));
    TabulatedAtmosphere table(eos, data + "/atmosphere/cond_gn93_tau100_solar_proxy.dat",
                              TabulatedAtmosphere::Mixture::allow_documented_proxy);
    ConvectiveAtmosphere column(eos, *rad);
    auto ref = solar_scaled(.7, .02);
    ref.basis = AbundanceBasis::baryon_mass;
    CompositionCorrectedAtmosphere atmosphere(eos, table, column, ref);
    PPChains nuclear(PPRates::solar_fusion_ii, PPScreening::salpeter_van_horn);
    Physics physics{&eos, &opacity, &nuclear, 1.9};
    auto c = solar_scaled(.58, .02);
    c.basis = AbundanceBasis::baryon_mass;
    c.X[1] = .085;
    c.X[2] -= .085;
    auto seed = example::stellar_seed(256, .1 * constants::Msun, .15 * constants::Rsun, ref, nuclear,
                                      atmosphere, 1.5);
    auto solar = relax(seed, physics, atmosphere);
    if (!solar.converged)
      throw std::runtime_error(solar.message);
    seed = solar.model;
    for (auto& mixture : seed.comp)
      mixture = c;
    auto equilibrium = relax(seed, physics, atmosphere);
    if (!equilibrium.converged)
      throw std::runtime_error(equilibrium.message);
    constexpr double dt = 1e8 * 365.25 * 86400;
    const auto full = evolve_step(equilibrium.model, physics, atmosphere, dt);
    const auto half = evolve_step(equilibrium.model, physics, atmosphere, dt / 2);
    const auto fine = half.converged ? evolve_step(half.model, physics, atmosphere, dt / 2) : half;
    if (!full.converged || !fine.converged)
      throw std::runtime_error(full.converged ? fine.message : full.message);
    double spread = 0, difference = 0;
    for (std::size_t i = 0; i < full.model.size(); ++i)
      for (std::size_t j = 0; j < NSPEC; ++j) {
        spread = std::max(spread, std::abs(full.model.comp[i].X[j] - full.model.comp[0].X[j]));
        difference = std::max(difference, std::abs(full.model.comp[i].X[j] - fine.model.comp[i].X[j]));
      }
    bool pass = full.model.age == dt && full.model.comp[0].X[0] < c.X[0] &&
                full.convective_mass_fraction > .99999 && spread < 1e-12 && difference < 1e-7 &&
                std::abs(full.luminosity_balance) < 2e-7 && std::abs(full.nuclear_mass_balance) < 2e-6;
    std::printf("[%s] helium-rich coupled evolution: convective mass %.9g, abundance difference %.9g, energy "
                "%.9g, rest mass %.9g\n",
                pass ? "PASS" : "FAIL", full.convective_mass_fraction, difference, full.luminosity_balance,
                full.nuclear_mass_balance);
    return pass ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "extended evolution: %s\n", e.what());
    return 1;
  }
}
