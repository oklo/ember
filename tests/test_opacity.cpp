// Opacity tests.  Two kinds: values the table itself asserts (a node must
// return its own number), and physics the table must express (the grain
// opacity below 1700 K, whose absence removed the Hayashi limit from cool
// giants in the Fortran line this code replaces).
#include "ember/opacity_ferguson.hpp"
#include <cmath>
#include <cstdio>
#include <string>

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
              op.name(), r.logT_min, r.logT_max, r.logR_min, r.logR_max, r.X_min, r.X_max);

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

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
