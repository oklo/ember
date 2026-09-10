#include "ember/atmosphere_composition.hpp"
#include "ember/constants.hpp"
#include "ember/eos_composite.hpp"
#include "ember/eos_composition.hpp"
#include "ember/opacity_mixture.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace ember;
static int failures = 0;
void check(bool pass, const char* label, double value = 0) {
  failures += !pass;
  std::printf("[%s] %s %.10g\n", pass ? "PASS" : "FAIL", label, value);
}
template <class F> bool throws(F f) {
  try {
    f();
  } catch (const std::exception&) {
    return true;
  }
  return false;
}
class ConstantOpacity : public Opacity {
  OpacityState eval(double, double, const Composition&) const override { return {.4, 0, 0}; }
  const char* name() const override { return "constant"; }
};
int main() {
  try {
    // A constant-opacity, ideal ion gas is radiatively stable: recover the
    // Eddington law and hydrostatic pressure including the top radiation term.
    CompositeEos gas(std::vector<std::unique_ptr<EosComponent>>{});
    gas.add<IonGas>();
    gas.add<Radiation>();
    ConstantOpacity k;
    auto ref = solar_scaled(.7, .02);
    ref.basis = AbundanceBasis::baryon_mass;
    ConvectiveAtmosphereOptions opt;
    opt.tolerance = 1e-9;
    opt.sensitivity_tolerance = 1e-9;
    ConvectiveAtmosphere analytic(gas, k, opt);
    const auto exact = analytic.eval(4000, 1e5, ref);
    const double T = 4000 * std::pow(.75 * (100 + 2. / 3), .25),
                 P = 1e5 * 100 / .4 + 2. / 3 * constants::sigma_SB * std::pow(4000., 4) / constants::c;
    check(std::abs(exact.T / T - 1) < 2e-8 && std::abs(exact.P / P - 1) < 2e-8,
          "radiatively stable column recovers analytic Eddington solution", exact.T / T - 1);
    check(std::abs(exact.dlnT_dlnTeff - 1) < 2e-8 && std::abs(exact.dlnT_dlng) < 2e-8 &&
              std::abs(exact.dlnP_dlng - (1e5 * 100 / .4) / P) < 1e-8,
          "constant-opacity column sensitivities recover analytic limits");
    const std::string data = EMBER_DATA_DIR;
    CompositionHelmholtzEos eos(data + "/eos/freeeos300_hhe_extended.dat",
                                HelmholtzTableEos::Mixture::allow_documented_proxy);
    MixtureOpacity source(data + "/opacity/aesopus21_gs98_mixture.dat");
    ElementalOpacity opacity(source);
    TabulatedAtmosphere cond(eos, data + "/atmosphere/cond_gn93_tau100_solar_proxy.dat",
                             TabulatedAtmosphere::Mixture::allow_documented_proxy);
    ConvectiveAtmosphere column(eos, opacity, opt);
    CompositionCorrectedAtmosphere corrected(eos, cond, column, ref);
    const double teff = 2817, g = 1.73e5;
    const auto baseline = cond.eval(teff, g, cond.composition()), b = corrected.eval(teff, g, ref);
    check(std::abs(b.T / baseline.T - 1) < 1e-13 && std::abs(b.P / baseline.P - 1) < 1e-13,
          "reference mixture exactly recovers COND boundary");
    check(std::abs(b.dlnT_dlnTeff - baseline.dlnT_dlnTeff) < 1e-13 &&
              std::abs(b.dlnP_dlng - baseline.dlnP_dlng) < 1e-13,
          "reference composition retains COND derivatives");
    constexpr double h = 3e-5;
    double derivative = 0, change = 0;
    for (auto [x, y3] : {std::pair{.675, .02}, std::pair{.52, .04}, std::pair{.4, .01}, std::pair{.55, .105}}) {
      auto c = solar_scaled(x, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.X[1] = y3;
      c.X[2] -= y3;
      const auto a = corrected.eval(teff, g, c);
      change = std::max(change, std::abs(a.T / b.T - 1));
      for (int j = 0; j < 2; ++j) {
        const auto p = corrected.eval(teff * std::exp(j == 0 ? h : 0), g * std::exp(j == 1 ? h : 0), c),
                   m = corrected.eval(teff * std::exp(j == 0 ? -h : 0), g * std::exp(j == 1 ? -h : 0), c);
        derivative = std::max(
            {derivative, std::abs(std::log(p.T / m.T) / (2 * h) - (j == 0 ? a.dlnT_dlnTeff : a.dlnT_dlng)),
             std::abs(std::log(p.P / m.P) / (2 * h) - (j == 0 ? a.dlnP_dlnTeff : a.dlnP_dlng))});
      }
      const auto radiative = GreyAtmosphere(eos, opacity, {100, .001, 1e-9}).eval(teff, g, c),
                 convective = column.eval(teff, g, c);
      check(convective.T < radiative.T,
            "convective transport lowers deep-column temperature relative to radiative-only closure",
            convective.T / radiative.T);
    }
    check(derivative < 2e-4, "corrected-atmosphere responses follow finite differences through convection",
          derivative);
    check(change > .02, "helium enrichment produces a resolved atmosphere temperature change", change);
    auto c = ref;
    c.X[0] = .55;
    c.X[1] = .03;
    c.X[2] = .40;
    ConvectiveAtmosphereOptions raised = opt;
    raised.tau_top = .002;
    ConvectiveAtmosphere other(eos, opacity, raised);
    CompositionCorrectedAtmosphere topcheck(eos, cond, other, ref);
    const auto s = corrected.eval(teff, g, c), t = topcheck.eval(teff, g, c);
    check(std::abs(t.T / s.T - 1) < 1e-4 && std::abs(t.P / s.P - 1) < 1e-4,
          "composition response converges as starting optical depth is changed", t.P / s.P - 1);
    check(throws([&] { corrected.eval(1500, g, c); }),
          "corrected atmosphere retains non-grey anchor source bounds");
  } catch (const std::exception& e) {
    std::printf("EXCEPTION %s\n", e.what());
    return 1;
  }
  return failures ? 1 : 0;
}
