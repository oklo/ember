// Structure-equation tests.  A stellar model is hard to check as a whole, but
// each equation is a statement that can be checked on its own, and the
// Jacobian can be checked against the residual it claims to differentiate.
#include "ember/structure.hpp"
#include "ember/eos_composite.hpp"
#include "ember/opacity.hpp"
#include "ember/opacity_ferguson.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
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
  EosResponse eval_with_derivatives(double T, double rho, const Composition& comp) const override {
    EosResponse r{}; r.state = eval(T, rho, comp);
    return r;  // E_rho=0; cp, delta and grad_ad are constants in this EOS.
  }
};
class ConstantOpacity final : public Opacity {
public:
  double kappa{1.0};
  OpacityState eval(double, double, const Composition&) const override {
    return {kappa, 0.0, 0.0};
  }
  const char* name() const override { return "constant (structure test)"; }
};

// Compare every row and both endpoints after converting the luminosity
// column to a local luminosity unit. Row normalization prevents tiny physical
// terms from turning subtraction roundoff into a spurious relative failure.
static void check_jacobian(const Model& model, const Physics& phys, double dt,
                           const Model* prev, const std::string& label) {
  const auto analytic = zone_residual(model, 0, phys, dt, prev);
  const auto numerical = zone_residual_numerical(model, 0, phys, dt, prev, 3e-5);
  const double dm = model.m[1] - model.m[0];
  const double Lunit = std::max({std::abs(model.y[0].L), std::abs(model.y[1].L),
      std::abs(model.y[1].L - model.y[0].L - dm * analytic.f[2]), 1.0});
  double worst = 0.0, value_error = 0.0;
  std::size_t worst_row = 0;
  for (std::size_t k = 0; k < NVAR; ++k) {
    double scale = 0.0, error = 0.0;
    for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
      const auto& a = endpoint == 0 ? analytic.dfdy_lo : analytic.dfdy_hi;
      const auto& n = endpoint == 0 ? numerical.dfdy_lo : numerical.dfdy_hi;
      for (std::size_t v = 0; v < NVAR; ++v) {
        const double unit = v == static_cast<std::size_t>(Var::L) ? Lunit : 1.0;
        scale = std::max({scale, std::abs(a[k][v] * unit), std::abs(n[k][v] * unit)});
        error = std::max(error, std::abs((a[k][v] - n[k][v]) * unit));
      }
    }
    const double relative = error / std::max(scale, 1e-300);
    if (relative > worst) { worst = relative; worst_row = k; }
    value_error = std::max(value_error, std::abs(analytic.f[k] - numerical.f[k])
        / std::max({std::abs(analytic.f[k]), scale, 1e-300}));
  }
  std::printf("       %s: worst Jacobian row %zu\n", label.c_str(), worst_row);
  check(worst < 2e-6, label + ": all Jacobian entries", worst, 0.0);
  check(value_error < 1e-12, label + ": values unchanged", value_error, 0.0);
}

template<class F> static void rejects(F&& f, const std::string& label) {
  bool threw = false;
  try { f(); } catch (const std::exception&) { threw = true; }
  check(threw, label, threw ? 1.0 : 0.0, 1.0);
}

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

  // 2. The production Jacobian must differentiate the value-only equations.
  check_jacobian(m, phys, -1.0, nullptr, "low-mass interior");
  {
    FergusonOpacity low(std::string(EMBER_DATA_DIR) + "/opacity/ferguson_gs98_z020.dat");
    OpalOpacity high(std::string(EMBER_DATA_DIR) + "/opacity/opal_gs98_z020.dat");
    BlendedOpacity tables(low, high);
    Physics real{&eos, &tables, &nuc, 1.9};
    Model t = m;
    t.comp.assign(2, solar_scaled(.7, .02));
    check_jacobian(t, real, -1.0, nullptr, "OPAL + FD + pp interior");
    for (double logT : {4.0, 4.23, 4.5}) {
      for (std::size_t i = 0; i < 2; ++i) {
        const double lt = logT + (i == 0 ? .001 : -.001);
        t.y[i].lnT = std::log(10.0) * lt;
        t.y[i].lnrho = std::log(10.0) * (-2.37 + 3 * (lt - 6));
      }
      check_jacobian(t, real, -1.0, nullptr, "opacity transition logT=" + std::to_string(logT));
    }
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

      // Both endpoints and all variables, with independently perturbed
      // value-only equations. A row-scaled comparison respects small physical
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
          const double num = dm * (zone_equations(plus, 0, p, -1.0)[3]
                                - zone_equations(minus, 0, p, -1.0)[3]) / (2e-4);
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

  // 5. Check the complete matrix under contraction and in degenerate matter.
  //    The previous model stays fixed when the new model is perturbed.
  {
    Model previous = m;
    for (auto& point : previous.y) {
      point.lnT -= 0.01; point.lnrho -= 0.02; point.lnr += 0.01;
    }
    check_jacobian(m, phys, 1e11, &previous, "contracting interior");

    Model wd = m;
    Composition helium{}; helium[Species::He4] = 1.0;
    wd.comp = {helium, helium};
    wd.y[0] = {std::log(1e8), std::log(1e6), std::log(1e6), 1e30};
    wd.y[1] = {std::log(1.02e8), std::log(0.95e6), std::log(0.98e6), 1.01e30};
    check_jacobian(wd, phys, 0.0, nullptr, "degenerate helium");
    previous = wd;
    previous.y[0].lnrho -= 0.01; previous.y[0].lnT -= 0.02;
    check_jacobian(wd, phys, 1e12, &previous, "contracting helium");

    // Force finite-efficiency convection with the real EOS. Unlike the
    // ideal-gas fixtures, cp, delta and grad_ad all respond to the state.
    for (const auto& model : {m, wd}) {
      Model t = model;
      ConstantOpacity constant;
      Physics p{&eos, &constant, &nuc, 1.9};
      const auto a = eos.eval(t.T(0), t.rho(0), t.comp[0]);
      const auto b = eos.eval(t.T(1), t.rho(1), t.comp[1]);
      EosState mid{};
      mid.P = 0.5 * (a.P + b.P); mid.cp = 0.5 * (a.cp + b.cp);
      mid.delta = 0.5 * (a.delta + b.delta);
      const double T = 0.5 * (t.T(0) + t.T(1)), rho = 0.5 * (t.rho(0) + t.rho(1));
      const double r = 0.5 * (t.r(0) + t.r(1)), mass = 0.5 * (t.m[0] + t.m[1]);
      constant.kappa = mixing_length_U(T, rho, 1.0, constants::G * mass / (r * r), mid, p.alpha_mlt);
      t.y[0].L = t.y[1].L = 5.0 * 16.0 * M_PI * constants::a_rad * constants::c
          * constants::G * mass * std::pow(T, 4) / (3.0 * constant.kappa * mid.P);
      check_jacobian(t, p, 0.0, nullptr, t.comp[0].h1() > 0.0
          ? "finite-efficiency interior" : "finite-efficiency helium");
    }

    wd.y[0].L = wd.y[1].L = 0.0;
    check_jacobian(wd, phys, 0.0, nullptr, "zero luminosity");
    wd.y[0].L = wd.y[1].L = -1e28;
    check_jacobian(wd, phys, 0.0, nullptr, "inward luminosity");

    // In an ideal gas E has no density response. Compress at fixed T to
    // isolate the sign and magnitude of the P d(1/rho)/dt contribution.
    TestIdealEos ideal;
    Physics p{&ideal, &op, &nuc, 1.9};
    previous = m; previous.y[0].lnrho -= std::log(2.0);
    const double dt = 1e11;
    const auto stat = zone_residual(m, 0, p, 0.0);
    const auto dynamic = zone_residual(m, 0, p, dt, &previous);
    const double expected = -0.5 * constants::R_gas / 0.6 * m.T(0) / dt;
    check(std::abs((dynamic.f[2] - stat.f[2]) / expected - 1.0) < 1e-12,
          "isothermal compression releases gravitational heat", dynamic.f[2] - stat.f[2], expected);
    const auto rho_var = static_cast<std::size_t>(Var::lnrho);
    const double derivative = dynamic.dfdy_lo[2][rho_var] - stat.dfdy_lo[2][rho_var];
    // The exact finite volume change gives R*T*(rho/rho_old-1).
    // At rho/rho_old=2 its log-density derivative is twice its value.
    check(std::abs(derivative / (2*expected) - 1.0) < 1e-12,
          "compression derivative holds the old density fixed", derivative, 2*expected);
  }

  // 6. Reject undefined zones and time-dependent states that cannot be
  //    compared on the same Lagrangian mesh.
  {
    rejects([&] { (void)zone_residual(m, 1, phys, 0.0); }, "out-of-range zone rejected");
    Model invalid = m; invalid.m[1] = invalid.m[0];
    rejects([&] { (void)zone_residual(invalid, 0, phys, 0.0); }, "zero mass interval rejected");
    rejects([&] { (void)zone_residual(m, 0, phys, 1.0); }, "positive dt requires previous model");
    invalid = m; invalid.m[1] *= 1.01;
    rejects([&] { (void)zone_residual(m, 0, phys, 1.0, &invalid); }, "changed previous mesh rejected");
    rejects([&] { (void)zone_residual(m, 0, phys, std::numeric_limits<double>::quiet_NaN()); },
            "non-finite time step rejected");
    class ValueOnlyEos final : public Eos {
      TestIdealEos ideal;
    public:
      EosState eval(double T, double rho, const Composition& comp) const override {
        return ideal.eval(T, rho, comp);
      }
      const char* name() const override { return "value-only test EOS"; }
    } value_only;
    Physics p{&value_only, &op, &nuc, 1.9};
    const auto reference = zone_residual_numerical(m, 0, p, 0.0);
    check(std::isfinite(reference.f[3]), "numerical reference supports a value-only EOS", reference.f[3], 0.0);
    rejects([&] { (void)zone_residual(m, 0, p, 0.0); }, "missing EOS derivatives are never silently zero");
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
