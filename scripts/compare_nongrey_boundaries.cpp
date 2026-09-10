// Static and coupled-step checks at compositions on the existing 0.1 Msun
// track.
#include "../examples/stellar_seed.hpp"
#include "ember/atmosphere_composition.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/eos_composition.hpp"
#include "ember/evolution.hpp"
#include "ember/opacity_mixture.hpp"
#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv) {
  using namespace ember;
  try {
    if (argc < 2 || argc > 3)
      throw std::invalid_argument("grid path and optional mesh size required");
    const std::size_t points = argc == 3 ? std::stoul(argv[2]) : 512;
    CompositionHelmholtzEos eos(
        "data/eos/freeeos300_hhe_extended.dat",
        HelmholtzTableEos::Mixture::allow_documented_proxy);
    auto rad = std::make_shared<StellarMixtureOpacity>("data/opacity");
    TabulatedConduction wd("data/conduction/condtab21wd.dat");
    auto hot = std::make_shared<HotConduction>(wd, 3e5, 1e6);
    CombinedOpacity opacity(rad, hot);
    TabulatedAtmosphere cond(
        eos, "data/atmosphere/cond_gn93_tau100_solar_proxy.dat",
        TabulatedAtmosphere::Mixture::allow_documented_proxy);
    auto ref = solar_scaled(.7, .02);
    ref.basis = AbundanceBasis::baryon_mass;
    ConvectiveAtmosphere column(eos, *rad);
    CompositionCorrectedAtmosphere corrected(eos, cond, column, ref);
    CompositionAtmosphereGrid grid(
        eos, argv[1],
        CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    PPChains nuclear(PPRates::solar_fusion_ii, PPScreening::salpeter_van_horn);
    Physics physics{&eos, &opacity, &nuclear, 1.9};
    Model continuation;
    std::printf("{\"points\":%zu,\"columns\":[\"X\",\"Y3\",\"mode\",\"R_Rsun\","
                "\"L_Lsun\","
                "\"Teff_K\",\"Tc_K\",\"rhoc\",\"residual\"],\"rows\":[\n",
                points);
    bool first = true;
    for (auto [x, y] :
         {std::pair{.7, 0.}, std::pair{.6, .08},
          std::pair{.53254996944256, .09343622410422}, std::pair{.46, .08}}) {
      auto c = solar_scaled(x, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.X[1] = y;
      c.X[2] -= y;
      for (int mode = 0; mode < 2; ++mode) {
        const Atmosphere &atmosphere =
            mode ? static_cast<const Atmosphere &>(grid)
                 : static_cast<const Atmosphere &>(corrected);
        auto seed = continuation.size()
                        ? continuation
                        : example::stellar_seed(points, .1 * constants::Msun,
                                                .15 * constants::Rsun, ref,
                                                nuclear, atmosphere, 1.5);
        for (auto &composition : seed.comp)
          composition = c;
        auto result = relax(seed, physics, atmosphere);
        if (!result.converged)
          throw std::runtime_error(result.message);
        if (!mode)
          continuation = result.model;
        const auto &m = result.model;
        const double radius = m.r(m.size() - 1), lum = m.y.back().L;
        std::printf(
            "%s[%.17g,%.17g,\"%s\",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",
            first ? "" : ",\n", x, y, mode ? "nongrey" : "cond-corrected",
            radius / constants::Rsun, lum / constants::Lsun,
            std::pow(lum / (4 * M_PI * constants::sigma_SB * radius * radius),
                     .25),
            m.T(0), m.rho(0), result.residual);
        std::fflush(stdout);
        first = false;
        if (mode && y > 0) {
          EvolutionOptions options;
          const auto step = evolve_step(m, physics, atmosphere,
                                        1e7 * 365.25 * 86400, options);
          if (!step.converged)
            throw std::runtime_error("helium-rich coupled step: " +
                                     step.message);
          std::fprintf(stderr,
                       "X=%g X3=%g: non-grey static and coupled step passed\n",
                       x, y);
        }
      }
    }
    std::puts("\n]}");
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
