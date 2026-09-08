// Opacity tests.  Two kinds: values the table itself asserts (a node must
// return its own number), and physics the table must express (the grain
// opacity below 1700 K, whose absence removed the Hayashi limit from cool
// giants in the Fortran line this code replaces).
#include "ember/opacity_ferguson.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_blend.hpp"
#include <fstream>
#include <filesystem>
#include <chrono>
#include <limits>
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>

using namespace ember;
static int failures = 0;
static void check(bool ok, const std::string& what, double got, double want) {
  if (!ok) ++failures;
  std::printf("  [%s] %-52s got %-12.5g want %-12.5g\n",
              ok ? "PASS" : "FAIL", what.c_str(), got, want);
}
static void near(double got, double want, double tol, const std::string& what) {
  check(std::abs(got - want) <= tol * std::abs(want), what, got, want);
}

int main() {
  const std::string data = std::string(EMBER_DATA_DIR) + "/opacity/ferguson_gs98_z020.dat";
  FergusonOpacity op(data);
  const auto r = op.range();
  std::printf("ember opacity - %s\n  log T %.2f..%.2f, log R %.1f..%.1f, X %.2f..%.2f\n\n",
              op.name(), r.logT_min, r.logT_max, r.logD_min, r.logD_max, r.X_min, r.X_max);

  auto comp_with_X = [](double X) {
    Composition c = solar_scaled(X, 0.020);
    return c;
  };

  // 1. Grain opacity below 1700 K.  This is the whole point of the table: the
  //    1983 grid had to be extrapolated here and returned ~1e-9, so the
  //    envelope went transparent and cool giants had no Hayashi limit left.
  {
    const double T = 1500.0, lt = std::log10(T);
    const double rho = std::pow(10.0, -2.0 + 3.0 * (lt - 6.0)); // log R = -2
    const auto s = op.eval(T, rho, comp_with_X(0.70));
    check(s.kappa > 1e-3, "grains at 1500 K give a real opacity", s.kappa, 1e-3);
  }

  // 2. Monotone rise into the grain regime: cooling from 3000 K to 1200 K at
  //    fixed log R must increase the opacity by orders of magnitude.
  {
    auto kap_at = [&](double T) {
      const double lt = std::log10(T);
      const double rho = std::pow(10.0, -3.0 + 3.0 * (lt - 6.0));
      return op.eval(T, rho, comp_with_X(0.70)).kappa;
    };
    const double hot = kap_at(3000.0), cold = kap_at(1200.0);
    check(cold > 10.0 * hot, "opacity climbs steeply into grain condensation",
          cold / hot, 10.0);
  }

  // 3. Analytic derivatives must match centred differences of the code.
  {
    const double T = 4000.0, rho = 1e-8, d = 1e-4;
    const auto c = comp_with_X(0.70);
    const auto s = op.eval(T, rho, c);
    const double kp = op.eval(T * (1 + d), rho, c).kappa;
    const double km = op.eval(T * (1 - d), rho, c).kappa;
    near(s.dlnk_dlnT, (std::log(kp) - std::log(km)) / (2 * d), 1e-4,
         "dlnk/dlnT vs numerical (log R term included)");
    const double rp = op.eval(T, rho * (1 + d), c).kappa;
    const double rm = op.eval(T, rho * (1 - d), c).kappa;
    near(s.dlnk_dlnRho, (std::log(rp) - std::log(rm)) / (2 * d), 1e-4,
         "dlnk/dlnrho vs numerical");
  }

  // 3b. The density derivative must follow the temperature slope limiter.
  {
    double worst = 0.0;
    constexpr double step = 1e-6;
    for (double T : {2450.0, 2850.0, 3430.0, 3700.0, 4300.0}) {
      for (double logR : {-6.3, -4.2, -2.6, -0.3}) {
        const double rho = std::pow(10.0, logR + 3.0 * (std::log10(T) - 6.0));
        const auto c = comp_with_X(0.7);
        const auto s = op.eval(T, rho, c);
        const double dr = (std::log(op.eval(T, rho * std::exp(step), c).kappa)
                          - std::log(op.eval(T, rho * std::exp(-step), c).kappa)) / (2 * step);
        const double dT = (std::log(op.eval(T * std::exp(step), rho, c).kappa)
                          - std::log(op.eval(T * std::exp(-step), rho, c).kappa)) / (2 * step);
        worst = std::max({worst, std::abs(dr - s.dlnk_dlnRho) / std::max(1.0, std::abs(dr)),
                         std::abs(dT - s.dlnk_dlnT) / std::max(1.0, std::abs(dT))});
      }
    }
    check(worst < 2e-5, "density derivative follows the temperature slope limiter", worst, 0.0);
  }

  // 4. Composition dependence is real, not a scaling law.
  {
    const double T = 3000.0, rho = 1e-10;
    const double k1 = op.eval(T, rho, comp_with_X(0.10)).kappa;
    const double k2 = op.eval(T, rho, comp_with_X(0.80)).kappa;
    check(k2 > k1, "opacity increases with hydrogen fraction", k2 / k1, 1.0);
  }

  // 5. Outside the table is an error, never a silent extrapolation.
  {
    bool threw = false;
    try { op.eval(200.0, 1e-10, comp_with_X(0.70)); } catch (const std::exception&) { threw = true; }
    check(threw, "below the table floor throws rather than extrapolating", threw, 1.0);
    threw = false;
    try { op.eval(1e6, 1.0, comp_with_X(0.70)); } catch (const std::exception&) { threw = true; }
    check(threw, "above the table ceiling throws", threw, 1.0);
  }

  OpalOpacity high(std::string(EMBER_DATA_DIR) + "/opacity/opal_gs98_z020.dat");
  BlendedOpacity blend(op, high);
  const auto comp = comp_with_X(0.7);
  // Original GS98hz table 73: logT=4.00, logR=-3.0, log(kappa)=1.264;
  // table 8: X=0, logT=6.50, logR=0, log(kappa)=1.638.
  near(std::log10(high.eval(1e4, 1e-9, comp).kappa), 1.264, 1e-12, "OPAL hydrogen source node (GS98hz table 73)");
  near(std::log10(high.eval(std::pow(10.0, 6.5), std::pow(10.0, 1.5), comp_with_X(0)).kappa),
       1.638, 1e-12, "OPAL helium source node (GS98hz table 8)");
  // Test every stored cell as well, including all boundaries and X planes.
  {
    std::ifstream table(std::string(EMBER_DATA_DIR) + "/opacity/opal_gs98_z020.dat");
    std::size_t nx, nt, nr; table >> nx >> nt >> nr;
    std::string title; std::getline(table, title);
    std::vector<double> lr(nr), lt(nt);
    for (double& v : lr) table >> v;
    for (double& v : lt) table >> v;
    double worst = 0.0;
    for (std::size_t ix = 0; ix < nx; ++ix) {
      double X, Z; table >> X >> Z;
      for (double T : lt) for (double R : lr) {
        double value; table >> value;
        const double rho = std::pow(10.0, R + 3.0 * (T - 6.0));
        const double got = std::log10(high.eval(std::pow(10.0, T), rho, solar_scaled(X, Z)).kappa);
        worst = std::max(worst, std::abs(got - value));
      }
    }
    check(worst < 1e-12, "OPAL reproduces every imported source cell", worst, 0.0);
  }
  {
    const double k = high.eval(1e7, 2e-5, comp).kappa;
    near(k, 0.2 * (1.0 + comp.h1()), 0.04, "hot dilute gas approaches electron scattering opacity");
    double worst = 0.0;
    constexpr double h = 1e-6;
    for (double lt : {4.001, 4.123, 4.271, 4.499, 5.371, 6.573, 7.013})
      for (double lr : {-6.13, -3.37, -0.19, 0.77}) {
        const double T = std::pow(10.0, lt), rho = std::pow(10.0, lr + 3.0 * (lt - 6.0));
        const auto c = comp_with_X(0.63);
        const auto k0 = blend.eval(T, rho, c);
        const double dt = std::log(blend.eval(T * std::exp(h), rho, c).kappa
                                / blend.eval(T * std::exp(-h), rho, c).kappa) / (2 * h);
        const double dr = std::log(blend.eval(T, rho * std::exp(h), c).kappa
                                / blend.eval(T, rho * std::exp(-h), c).kappa) / (2 * h);
        worst = std::max({worst, std::abs(dt - k0.dlnk_dlnT) / std::max(1.0, std::abs(dt)),
                                std::abs(dr - k0.dlnk_dlnRho) / std::max(1.0, std::abs(dr))});
      }
    check(worst < 2e-5, "OPAL and blend derivatives follow actual opacity", worst, 0.0);
    for (double lt : {4.0, 4.5}) {
      const double T = std::pow(10.0, lt), rho = std::pow(10.0, -2.3 + 3 * (lt - 6));
      const auto a = blend.eval(T * std::exp(-h), rho, comp);
      const auto b = blend.eval(T * std::exp(h), rho, comp);
      check(std::abs(std::log(a.kappa / b.kappa)) < 1e-4
            && std::abs(a.dlnk_dlnT - b.dlnk_dlnT) < 1e-3
            && std::abs(a.dlnk_dlnRho - b.dlnk_dlnRho) < 1e-3,
            "blend value and partials are continuous at endpoint", a.kappa, b.kappa);
    }
  }
  {
    int rejected = 0;
    for (const Opacity* source : std::array<const Opacity*, 3>{&op, &high, &blend}) {
      try { source->eval(15000, 1e-8, solar_scaled(.7, .014)); }
      catch (const std::domain_error&) { ++rejected; }
      try { source->density_range(15000, solar_scaled(.7, .014)); }
      catch (const std::domain_error&) { ++rejected; }
    }
    check(rejected == 6, "fixed-Z opacity rejects mismatched metallicity", rejected, 6);
    rejected = 0;
    for (const auto& state : std::array{std::array{1e8, 1.0}, std::array{1e6, 11.0},
                                      std::array{1e6, 1e-10}, std::array{200.0, 1e-10},
                                      std::array{std::numeric_limits<double>::quiet_NaN(), 1.0}}) {
      try { blend.eval(state[0], state[1], comp); }
      catch (const std::domain_error&) { ++rejected; }
    }
    check(rejected == 5, "blend never fills missing temperature/density coverage", rejected, 5);
  }
  {
    const auto path = std::filesystem::temp_directory_path()
        / ("ember-opacity-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".dat");
    auto invalid = [&](std::string_view axis, double second_Z, double cell, bool truncate, bool extra) {
      std::ofstream file(path);
      file << "2 4 4 test\n" << axis << "\n4 5 6 7\n";
      for (int i = 0; i < 2; ++i) {
        file << (i == 0 ? 0.0 : 0.7) << ' ' << (i == 0 ? .02 : second_Z) << '\n';
        for (int j = 0; j < (truncate ? 15 : 16); ++j) file << cell << ' ';
        file << '\n';
      }
      if (extra) file << "extra\n";
      file.close();
      try { OpalOpacity bad(path); return false; }
      catch (const std::runtime_error&) { return true; }
    };
    const bool rejected = invalid("-8 -8 -2 1", .02, 0, false, false)
                       && invalid("-8 -5 -2 1", .03, 0, false, false)
                       && invalid("-8 -5 -2 1", .02, 9.999, false, false)
                       && invalid("-8 -5 -2 1", .02, 0, true, false)
                       && invalid("-8 -5 -2 1", .02, 0, false, true);
    std::filesystem::remove(path);
    check(rejected, "invalid axes, mixed Z, missing cells and truncation throw", rejected, 1);
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
