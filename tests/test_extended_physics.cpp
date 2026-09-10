#include "ember/conduction_table.hpp"
#include "ember/constants.hpp"
#include "ember/eos_composition.hpp"
#include "ember/opacity_mixture.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
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
class UnboundedOpacity final : public Opacity {
public:
  OpacityState eval(double, double, const Composition&) const override { return {1, 0, 0}; }
  const char* name() const override { return "constant opacity"; }
};
int main() {
  try {
    const std::string data = EMBER_DATA_DIR;
    CompositionHelmholtzEos eos(data + "/eos/freeeos300_hhe_extended.dat",
                                HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::ifstream in(std::string(EMBER_TEST_DATA_DIR) + "/freeeos300_extended_reference.dat");
    double X, Y3, T, rho, P, E, cv, cp, ad, values = 0, responses = 0;
    int count = 0;
    while (in >> X >> Y3 >> T >> rho >> P >> E >> cv >> cp >> ad) {
      auto c = solar_scaled(X, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.X[1] = Y3;
      c.X[2] -= Y3;
      const auto e = eos.eval(T, rho, c);
      ++count;
      values = std::max({values, std::abs(e.P / P - 1), std::abs(e.E / E - 1)});
      responses = std::max(
          {responses, std::abs(e.cv / cv - 1), std::abs(e.cp / cp - 1), std::abs(e.grad_ad / ad - 1)});
    }
    check(count == 63 && values < .001, "extended EOS P/E agree with independent FreeEOS queries", values);
    check(responses < .01, "extended EOS thermal responses agree with direct source", responses);
    auto c = solar_scaled(.57, .02);
    c.basis = AbundanceBasis::baryon_mass;
    c.X[1] = .035;
    c.X[2] -= .035;
    const auto bounds = eos.density_range(2355, c).value();
    const double test_rho = bounds.min * 1.1;
    check(std::abs(eos.rho_from_PT(2355, eos.eval(2355, test_rho, c).P, c) / test_rho - 1) < 1e-9,
          "PT inversion respects the table floor when the ideal guess lies below it");
    check(throws([&] { eos.rho_from_PT(2355, .5 * eos.eval(2355, bounds.min, c).P, c); }),
          "PT inversion rejects pressure below table support");
    MixtureOpacity low(data + "/opacity/aesopus21_gs98_mixture.dat");
    ElementalOpacity mapped(low);
    StellarMixtureOpacity stellar(data + "/opacity");
    constexpr double h = 1e-6;
    double derivatives = 0, invariance = 0;
    for (auto [t, r] : {std::pair{3800., 1e-5}, std::pair{6200., 3e-4}, std::pair{24000., .01}}) {
      const auto k = mapped.eval(t, r, c);
      for (int j = 0; j < 2; ++j) {
        auto p = c, m = c;
        p.X[j] += h;
        p.X[2] -= h;
        m.X[j] -= h;
        m.X[2] += h;
        const double fd = std::log(mapped.eval(t, r, p).kappa / mapped.eval(t, r, m).kappa) / (2 * h);
        derivatives = std::max(derivatives, std::abs(fd - (j == 0 ? k.dlnk_dX : k.dlnk_dY3)));
      }
      for (int j = 0; j < 2; ++j) {
        const double fd =
            std::log(mapped.eval(t * std::exp(j == 0 ? h : 0), r * std::exp(j == 1 ? h : 0), c).kappa /
                     mapped.eval(t * std::exp(j == 0 ? -h : 0), r * std::exp(j == 1 ? -h : 0), c).kappa) /
            (2 * h);
        derivatives = std::max(derivatives, std::abs(fd - (j == 0 ? k.dlnk_dlnT : k.dlnk_dlnRho)));
      }
      // Independently construct an atomic, isotope-free mixture with identical atom densities.
      auto src = c;
      src.basis = AbundanceBasis::atomic_mass;
      for (std::size_t j = 0; j < NSPEC; ++j)
        src.X[j] = nuclides[j].A * c.X[j] / mass_numbers[j];
      src.X[2] += nuclides[2].A * c.X[1] / 3;
      src.X[1] = 0;
      const double scale = src.sum();
      for (double& x : src.X)
        x /= scale;
      invariance = std::max(invariance, std::abs(k.kappa / (scale * low.eval(t, r * scale, src).kappa) - 1));
    }
    check(derivatives < 2e-6, "opacity thermal/H/He3 derivatives include density and metallicity mapping",
          derivatives);
    check(invariance < 1e-13, "elemental mapping preserves extinction at equal atom densities", invariance);
    double hot_derivatives = 0;
    for (auto [t, r] : {std::pair{28000., .03}, std::pair{450000., 1.3}, std::pair{4.7e6, 420.}}) {
      const auto k = stellar.eval(t, r, c);
      for (int j = 0; j < 2; ++j) {
        auto p = c, m = c;
        p.X[j] += h;
        p.X[2] -= h;
        m.X[j] -= h;
        m.X[2] += h;
        const double fd = std::log(stellar.eval(t, r, p).kappa / stellar.eval(t, r, m).kappa) / (2 * h);
        hot_derivatives = std::max(hot_derivatives, std::abs(fd - (j == 0 ? k.dlnk_dX : k.dlnk_dY3)));
      }
      for (int j = 0; j < 2; ++j) {
        const double fd =
            std::log(stellar.eval(t * std::exp(j == 0 ? h : 0), r * std::exp(j == 1 ? h : 0), c).kappa /
                     stellar.eval(t * std::exp(j == 0 ? -h : 0), r * std::exp(j == 1 ? -h : 0), c).kappa) /
            (2 * h);
        hot_derivatives = std::max(hot_derivatives, std::abs(fd - (j == 0 ? k.dlnk_dlnT : k.dlnk_dlnRho)));
      }
    }
    check(hot_derivatives < 2e-6,
          "full opacity family includes both temperature joins and isotope derivatives", hot_derivatives);
    auto outside = solar_scaled(.7, .005);
    check(throws([&] { mapped.eval(4000, 1e-5, outside); }), "opacity rejects metallicity extrapolation");
    TabulatedConduction conduction(data + "/conduction/condtab21wd.dat");
    double cd = 0;
    for (auto [t, r] : {std::pair{7.1e5, .2}, std::pair{4.65e6, 383.}, std::pair{1.2e7, 1.2e4}}) {
      const auto k = conduction.eval(t, r, c);
      for (int j = 0; j < 2; ++j) {
        auto p = c, m = c;
        p.X[j] += h;
        p.X[2] -= h;
        m.X[j] -= h;
        m.X[2] += h;
        const double fd = std::log(conduction.eval(t, r, p).kappa / conduction.eval(t, r, m).kappa) / (2 * h);
        cd = std::max(cd, std::abs(fd - (j == 0 ? k.dlnk_dX : k.dlnk_dY3)));
      }
      for (int j = 0; j < 2; ++j) {
        const double fd =
            std::log(conduction.eval(t * std::exp(j == 0 ? h : 0), r * std::exp(j == 1 ? h : 0), c).kappa /
                     conduction.eval(t * std::exp(j == 0 ? -h : 0), r * std::exp(j == 1 ? -h : 0), c).kappa) /
            (2 * h);
        cd = std::max(cd, std::abs(fd - (j == 0 ? k.dlnk_dlnT : k.dlnk_dlnRho)));
      }
    }
    check(cd < 2e-6, "mixture conduction derivatives differentiate the collision sum", cd);
    // Pure H at an original source node: logT=6, log rho=2, log K=11.919.
    Composition hydrogen{};
    hydrogen.basis = AbundanceBasis::baryon_mass;
    hydrogen.X[0] = 1;
    std::ifstream source(data + "/conduction/condtab21wd.dat");
    std::string line;
    std::getline(source, line);
    std::getline(source, line);
    std::getline(source, line);
    double z, mass;
    source >> z >> mass;
    double original = 0;
    for (int ir = 0; ir < 64; ++ir)
      for (int it = 0; it < 19; ++it) {
        double v;
        source >> v;
        if (ir == 32 && it == 9)
          original = v;
      }
    const double K = 16 * constants::sigma_SB * 1e18 / (3 * 100 * conduction.eval(1e6, 100, hydrogen).kappa);
    check(std::abs(std::log10(K) - original) < 1e-12, "pure hydrogen recovers the original conductivity cell",
          K);
    TabulatedConduction classic(data + "/conduction/condtab21_I.dat"),
        undamped(data + "/conduction/condtab21nd.dat");
    std::ifstream original_code(std::string(EMBER_TEST_DATA_DIR) + "/conduction_reference.dat");
    int mode;
    double charge, amass, ct, cr, ck, conductivity_error = 0;
    int samples = 0;
    while (original_code >> mode >> charge >> amass >> ct >> cr >> ck) {
      Composition pure{};
      pure.basis = AbundanceBasis::baryon_mass;
      const std::size_t species = charge == 1   ? 0
                                  : charge == 2 ? 2
                                  : charge == 6 ? 3
                                  : charge == 7 ? 5
                                  : charge == 8 ? 6
                                                : 7;
      pure.X[species] = 1;
      const auto& variant = mode < 0 ? classic : mode == 0 ? undamped : conduction;
      const double interpolated_K =
          16 * constants::sigma_SB * ct * ct * ct / (3 * cr * variant.eval(ct, cr, pure).kappa);
      conductivity_error = std::max(conductivity_error, std::abs(interpolated_K / ck - 1));
      ++samples;
    }
    check(samples == 30 && conductivity_error < .015,
          "interpolated conductivity agrees with independent source-code evaluations", conductivity_error);
    HotConduction hot(conduction);
    check(std::isinf(hot.eval(1e5, 1e-4, c).kappa), "neutral envelope has zero adopted electron conductance");
    const auto conductive_bounds = conduction.density_range(5e6, c).value();
    check(std::isfinite(conduction.eval(5e6, conductive_bounds.min * 1.000001, c).kappa) &&
              std::isfinite(conduction.eval(5e6, conductive_bounds.max / 1.000001, c).kappa) &&
              throws([&] { conduction.eval(5e6, conductive_bounds.min / 1.000001, c); }) &&
              throws([&] { conduction.eval(5e6, conductive_bounds.max * 1.000001, c); }),
          "conductivity density support follows all mapped source ions");
    CombinedOpacity combined(std::make_shared<UnboundedOpacity>(),
                             std::make_shared<HotConduction>(conduction));
    const auto combined_bounds = combined.density_range(5e6, c).value();
    check(combined_bounds.min == conductive_bounds.min && combined_bounds.max == conductive_bounds.max &&
              !combined.density_range(1e5, c) && throws([&] { combined.density_range(1e10, c); }),
          "combined transport retains conductive source bounds only where active");
    const auto k = hot.eval(6e5, .1, c);
    const double fd =
        std::log(hot.eval(6e5 * std::exp(h), .1, c).kappa / hot.eval(6e5 * std::exp(-h), .1, c).kappa) /
        (2 * h);
    check(std::abs(fd - k.dlnk_dlnT) < 2e-6, "conduction activation derivative is included",
          fd - k.dlnk_dlnT);
    check(throws([&] { conduction.eval(1e10, 1, c); }),
          "conduction rejects source temperature extrapolation");
  } catch (const std::exception& e) {
    std::printf("EXCEPTION %s\n", e.what());
    return 1;
  }
  return failures ? 1 : 0;
}
