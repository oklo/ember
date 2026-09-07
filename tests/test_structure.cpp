// Structure-equation tests.  A stellar model is hard to check as a whole, but
// each equation is a statement that can be checked on its own, and the
// Jacobian can be checked against the residual it claims to differentiate.
#include "ember/structure.hpp"
#include "ember/eos_composite.hpp"
#include "ember/opacity.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include <algorithm>
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

// A smooth analytic opacity, so the structure test is about the equations and
// not about a table's edges.
class PowerLawOpacity final : public Opacity {
public:
  OpacityState eval(double T, double rho, const Composition&) const override {
    OpacityState s{};
    s.kappa = 0.2 * std::pow(rho, 0.5) * std::pow(T / 1e6, -2.0);
    s.dlnk_dlnRho = 0.5; s.dlnk_dlnT = -2.0;
    return s;
  }
  const char* name() const override { return "power law (test)"; }
};

// An exactly ideal, monatomic gas isolates transport-row conditioning from
// changes in the adiabatic gradient of a partially degenerate EOS.
class TestIdealEos final : public Eos {
public:
  EosState eval(double T, double rho, const Composition&) const override {
    const double gas = constants::R_gas / 0.6;
    EosState e{};
    e.P = gas * rho * T; e.E = 1.5 * gas * T;
    e.chiT = e.chiRho = e.delta = 1.0;
    e.cv = 1.5 * gas; e.cp = 2.5 * gas; e.grad_ad = 0.4;
    return e;
  }
  const char* name() const override { return "ideal gas (structure test)"; }
};
class ConstantOpacity final : public Opacity {
public:
  double kappa{1.0};
  OpacityState eval(double, double, const Composition&) const override {
    return {kappa, 0.0, 0.0};
  }
  const char* name() const override { return "constant (structure test)"; }
};

int main() {
  CompositeEos eos; PowerLawOpacity op; PPChains nuc;
  Physics phys{&eos, &op, &nuc, 1.9};
  std::printf("ember structure equations\n\n");

  // Build two points of a plausible low-mass interior.
  Model m;
  m.M = 0.1 * constants::Msun;
  m.m = {0.30 * m.M, 0.32 * m.M};
  m.comp = {solar_scaled(0.70, 0.014), solar_scaled(0.70, 0.014)};
  m.y.resize(2);
  m.y[0] = {std::log(0.05 * constants::Rsun), std::log(80.0), std::log(6.0e6), 1.0e30};
  m.y[1] = {std::log(0.055 * constants::Rsun), std::log(70.0), std::log(5.6e6), 1.05e30};

  const auto R = zone_residual(m, 0, phys, -1.0);

  // 1. The residuals must be finite and of a sane size - a first sanity gate.
  {
    bool finite = true;
    for (double v : R.f) finite = finite && std::isfinite(v);
    check(finite, "residuals are finite", finite ? 1.0 : 0.0, 1.0);
  }

  // 2. The Jacobian must differentiate the residual it is paired with.  Recompute
  //    by a coarser difference and compare; agreement to a few parts in 1e4 is
  //    all a difference-of-differences can offer, and is enough to catch a
  //    wrong sign, a wrong variable, or a missing term.
  {
    double worst = 0.0; std::size_t wk = 0, wv = 0;
    for (std::size_t v = 0; v < NVAR; ++v) {
      Model mp = m, mm = m;
      const Var vv = static_cast<Var>(v);
      const double s = (vv == Var::L) ? 1e-4 * std::abs(m.y[0].L) : 1e-4;
      mp.y[0][vv] += s; mm.y[0][vv] -= s;
      const auto fp = zone_residual(mp, 0, phys, -1.0).f;
      const auto fm = zone_residual(mm, 0, phys, -1.0).f;
      for (std::size_t k = 0; k < NVAR; ++k) {
        const double num = (fp[k] - fm[k]) / (2 * s);
        const double ana = R.dfdy_lo[k][v];
        const double scale = std::max({std::abs(num), std::abs(ana), 1e-300});
        const double rel = std::abs(num - ana) / scale;
        if (rel > worst) { worst = rel; wk = k; wv = v; }
      }
    }
    std::printf("       (worst at equation %zu, variable %zu)\n", wk, wv);
    check(worst < 1e-3, "Jacobian matches its own residual", worst, 0.0);
  }

  // 3. Mass conservation must be exact for a shell built to satisfy it: place
  //    the outer point at the radius the equation demands and the residual
  //    must vanish.
  {
    // The zone has to be thin for this to be a test of the equation rather
    // than of the integration used to build the shell: the equation takes a
    // centred difference, the construction below a one-point one, and across a
    // zone where the radius changes ten percent those differ at the ten
    // percent level - which is discretisation, not error.
    Model t = m;
    t.m[1] = t.m[0] + 1e-4 * t.M;
    const double r0 = std::exp(t.y[0].lnr), rho = std::exp(t.y[0].lnrho);
    const double dm = t.m[1] - t.m[0];
    // d(ln r)/dm = 1/(4 pi r^3 rho), integrated crudely over the zone
    const double dlnr = dm / (4.0 * M_PI * r0 * r0 * r0 * rho);
    t.y[1].lnr = t.y[0].lnr + dlnr;
    t.y[1].lnrho = t.y[0].lnrho;   // uniform density shell
    const auto Rt = zone_residual(t, 0, phys, -1.0);
    check(std::abs(Rt.f[0]) < 1e-2 * std::abs(dlnr / dm),
          "mass equation vanishes for a thin shell built to it",
          std::abs(Rt.f[0]) / std::abs(dlnr / dm), 0.0);
  }

  // 4. Exercise the actual transport row across the convection transition.
  //    Opacity and luminosity set independent U and grad_rad in this fixture.
  {
    TestIdealEos ideal;
    ConstantOpacity constant;
    Physics p{&ideal, &constant, &nuc, 1.9};
    for (const auto [rad, targetU] : {std::pair{0.2, 1.0}, std::pair{5.0, 1.0},
                                    std::pair{1e6, 1e-12}}) {
      Model t = m;
      const auto a = ideal.eval(t.T(0), t.rho(0), t.comp[0]);
      const auto b = ideal.eval(t.T(1), t.rho(1), t.comp[1]);
      const double Tb = 0.5 * (t.T(0) + t.T(1));
      const double rhob = 0.5 * (t.rho(0) + t.rho(1));
      const double rb = 0.5 * (t.r(0) + t.r(1));
      const double mb = 0.5 * (t.m[0] + t.m[1]);
      EosState eb = a; eb.P = 0.5 * (a.P + b.P);
      constant.kappa = mixing_length_U(Tb, rhob, 1.0,
          constants::G * mb / (rb * rb), eb, p.alpha_mlt) / targetU;
      const double L = rad * 16.0 * M_PI * constants::a_rad * constants::c
          * constants::G * mb * std::pow(Tb, 4) / (3.0 * constant.kappa * eb.P);
      t.y[0].L = t.y[1].L = L;
      const double dm = t.m[1] - t.m[0];
      const double dlnP = std::log(b.P) - std::log(a.P);
      auto implied_grad = [&](const ZoneResidual& residual) {
        return (t.y[1].lnT - t.y[0].lnT - dm * residual.f[3]) / dlnP;
      };
      const auto r = zone_residual(t, 0, p, -1.0);
      const double grad = implied_grad(r);
      if (rad < 0.4) {
        check(std::abs(grad / rad - 1.0) < 1e-12,
              "stable transport row uses the radiative gradient", grad, rad);
      } else if (rad == 5.0) {
        check(grad > 0.4 && grad < rad,
              "unstable row retains finite superadiabaticity", grad, 0.4);
        Physics longer = p; longer.alpha_mlt *= 2.0;
        const double grad_longer = implied_grad(zone_residual(t, 0, longer, -1.0));
        check(grad_longer < grad && grad_longer > 0.4,
              "alpha_mlt changes the structure's temperature gradient", grad_longer, grad);
      } else {
        check(grad / rad < 1e-6 && grad > 0.4,
              "giant-efficiency fixture reaches grad/grad_rad < 1e-6", grad / rad, 1e-6);
        // Change outer T at fixed P by compensating rho. The coefficient of
        // dlnT in dm*f_transport stays one even when grad/grad_rad is tiny.
        const double coefficient = dm * (r.dfdy_hi[3][static_cast<std::size_t>(Var::lnT)]
                                        - r.dfdy_hi[3][static_cast<std::size_t>(Var::lnrho)]);
        check(std::abs(coefficient - 1.0) < 1e-4,
              "efficient convection leaves temperature row unscaled", coefficient, 1.0);
      }

      // Both endpoints and all variables, with a coarser step than the
      // implementation. A row-scaled comparison respects small physical
      // derivatives without turning roundoff in them into a false failure.
      double worst = 0.0;
      for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
        double scale = 0.0, error = 0.0;
        for (std::size_t v = 0; v < NVAR; ++v) {
          const Var vv = static_cast<Var>(v);
          const double unit = vv == Var::L ? std::abs(L) : 1.0;
          const double step = 1e-4 * unit;
          Model plus = t, minus = t;
          plus.y[endpoint][vv] += step; minus.y[endpoint][vv] -= step;
          const double num = dm * (zone_residual(plus, 0, p, -1.0).f[3]
                                - zone_residual(minus, 0, p, -1.0).f[3]) / (2e-4);
          const auto& jac = endpoint == 0 ? r.dfdy_lo : r.dfdy_hi;
          const double actual = dm * jac[3][v] * unit;
          error = std::max(error, std::abs(num - actual));
          scale = std::max({scale, std::abs(num), std::abs(actual)});
        }
        worst = std::max(worst, error / scale);
      }
      check(worst < 2e-6, "transport Jacobian matches at both zone endpoints", worst, 0.0);
    }
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
