// Static sensitivity experiments at fixed homogeneous compositions.
#include "../examples/stellar_seed.hpp"
#include "ember/atmosphere_composition.hpp"
#include "ember/conduction_table.hpp"
#include "ember/eos_composition.hpp"
#include "ember/evolution.hpp"
#include "ember/opacity_mixture.hpp"
#include <cstdio>
int main() {
  using namespace ember;
  try {
    CompositionHelmholtzEos eos("data/eos/freeeos300_hhe_extended.dat",
                                HelmholtzTableEos::Mixture::allow_documented_proxy);
    auto rad = std::make_shared<StellarMixtureOpacity>("data/opacity");
    TabulatedConduction wd("data/conduction/condtab21wd.dat"), classic("data/conduction/condtab21_I.dat"),
        undamped("data/conduction/condtab21nd.dat");
    TabulatedAtmosphere cond(eos, "data/atmosphere/cond_gn93_tau100_solar_proxy.dat",
                             TabulatedAtmosphere::Mixture::allow_documented_proxy);
    auto ref = solar_scaled(.7, .02);
    ref.basis = AbundanceBasis::baryon_mass;
    PPChains nuclear(PPRates::solar_fusion_ii, PPScreening::salpeter_van_horn);
    std::puts("{\"description\":\"static boundary/transport sensitivity at fixed homogeneous compositions; "
              "not age-matched "
              "tracks\",\"columns\":[\"X\",\"Y3\",\"mode\",\"R_Rsun\",\"L_Lsun\",\"Teff_K\",\"Tc_K\","
              "\"rhoc\",\"residual\"],\"rows\":[");
    bool first = true;
    Model continuation;
    for (auto [x, y3] : {std::pair{.7, 0.}, std::pair{.60, .08}, std::pair{.5, .08}, std::pair{.45, .08}}) {
      auto c = solar_scaled(x, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.X[1] = y3;
      c.X[2] -= y3;
      Model baseline;
      for (const std::string mode : {"reference", "raw-grey", "y076", "alpha15", "top002", "no-conduction",
                                     "classic", "undamped", "hot-join"}) {
        ConvectiveAtmosphereOptions options;
        if (mode == "y076")
          options.henyey_y = .076;
        if (mode == "alpha15")
          options.alpha = 1.5;
        if (mode == "top002")
          options.tau_top = .002;
        ConvectiveAtmosphere column(eos, *rad, options);
        CompositionCorrectedAtmosphere corrected(eos, cond, column, ref);
        const Atmosphere& atmosphere = mode == "raw-grey" ? static_cast<const Atmosphere&>(column)
                                                          : static_cast<const Atmosphere&>(corrected);
        const Conduction& source = mode == "classic" ? classic : mode == "undamped" ? undamped : wd;
        auto hot = std::make_shared<HotConduction>(source, mode == "hot-join" ? 6e5 : 3e5,
                                                   mode == "hot-join" ? 2e6 : 1e6);
        CombinedOpacity opacity(rad, mode == "no-conduction" ? nullptr : hot);
        Physics physics{&eos, &opacity, &nuclear, 1.9};
        auto seed = baseline.size() ? baseline
                    : continuation.size()
                        ? continuation
                        : example::stellar_seed(512, .1 * constants::Msun, .15 * constants::Rsun, ref,
                                                nuclear, atmosphere, 1.5);
        for (auto& mixture : seed.comp)
          mixture = c;
        const auto result = relax(seed, physics, atmosphere);
        if (!result.converged)
          throw std::runtime_error(mode + ": " + result.message);
        if (mode == "reference") {
          baseline = result.model;
          continuation = result.model;
        }
        const auto& m = result.model;
        const double radius = m.r(m.size() - 1), lum = m.y.back().L;
        std::printf("%s[%.17g,%.17g,\"%s\",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]", first ? "" : ",\n", x, y3,
                    mode.c_str(), radius / constants::Rsun, lum / constants::Lsun,
                    std::pow(lum / (4 * M_PI * constants::sigma_SB * radius * radius), .25), m.T(0), m.rho(0),
                    result.residual);
        std::fflush(stdout);
        first = false;
        std::fprintf(stderr, "X=%g Y3=%g %s R=%g L=%g\n", x, y3, mode.c_str(), radius / constants::Rsun,
                     lum / constants::Lsun);
      }
    }
    std::puts("\n]}");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
