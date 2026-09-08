#include "ember/atmosphere.hpp"
#include "ember/atmosphere_table.hpp"
#include "ember/boundary.hpp"
#include "ember/constants.hpp"
#include "ember/eos_composite.hpp"
#include "ember/opacity_ferguson.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_blend.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

using namespace ember;
static int failures = 0;
static void check(bool ok, const std::string& what, double got, double want) {
  if (!ok) ++failures;
  std::printf("  [%s] %-61s got %-12.5g want %-12.5g\n",
              ok ? "PASS" : "FAIL", what.c_str(), got, want);
}
static void near(double got, double want, double tol, const std::string& what) {
  check(std::abs(got - want) <= tol * std::max(std::abs(want), 1e-30), what, got, want);
}
template<class F> static bool throws(F&& f) {
  try { f(); } catch (const std::exception&) { return true; }
  return false;
}

class GasRadiationEos final : public Eos {
public:
  EosState eval(double T, double rho, const Composition& comp) const override {
    if (!(T > 0.0) || !(rho > 0.0)) throw std::domain_error("test EOS: invalid state");
    EosState e{};
    const double Pg = gas(comp) * rho * T;
    const double Pr = constants::a_rad * std::pow(T, 4) / 3.0;
    e.P = Pg + Pr; e.chiRho = Pg / e.P; e.chiT = (Pg + 4.0 * Pr) / e.P;
    return e;
  }
  static double gas(const Composition& comp) {
    return constants::R_gas * (comp.mu_ions_inv() + comp.mu_elec_inv());
  }
  const char* name() const override { return "ideal gas + radiation (atmosphere test)"; }
};
class ConstantOpacity final : public Opacity {
public:
  OpacityState eval(double, double, const Composition&) const override { return {0.4, 0.0, 0.0}; }
  const char* name() const override { return "constant (atmosphere test)"; }
};
class PressureOpacity final : public Opacity {
public:
  OpacityState eval(double T, double rho, const Composition& comp) const override {
    const auto e = eos.eval(T, rho, comp);
    return {1e-7 * e.P, e.chiT, e.chiRho};
  }
  const char* name() const override { return "kappa proportional to total P (test)"; }
  GasRadiationEos eos;
};
class VariableOpacity final : public Opacity {
public:
  OpacityState eval(double T, double rho, const Composition&) const override {
    return {0.01 * std::sqrt(rho / 1e-8) * std::pow(T / 4000.0, 3), 3.0, 0.5};
  }
  const char* name() const override { return "density and temperature dependent (test)"; }
};

static void check_surface(const Eos& eos, const Atmosphere& atm, const Composition& comp) {
  const double mass = 0.1 * constants::Msun, gravity = std::pow(10.0, 4.3), Teff = 3000.0;
  const double radius = std::sqrt(constants::G * mass / gravity);
  const auto a = atm.eval(Teff, gravity, comp);
  Point p{std::log(radius), std::log(a.rho), std::log(a.T),
          4.0 * M_PI * constants::sigma_SB * radius * radius * std::pow(Teff, 4)};
  const auto r = surface_residual(p, mass, comp, eos, atm);
  const double error = std::max(std::abs(r.f[0]), std::abs(r.f[1]));
  check(error < 1e-10, "surface equations vanish on an atmosphere-matched point", error, 0.0);

  // Check off the solution as well: the EOS derivatives belong to the mesh
  // point, while the atmosphere derivatives belong to Teff and g.
  p.lnT += 0.03; p.lnrho -= 0.1;
  const auto off = surface_residual(p, mass, comp, eos, atm);
  double worst = 0.0;
  for (std::size_t v = 0; v < NVAR; ++v) {
    const Var var = static_cast<Var>(v);
    const double unit = var == Var::L ? p.L : 1.0, step = 1e-4 * unit;
    Point plus = p, minus = p; plus[var] += step; minus[var] -= step;
    const auto fp = surface_residual(plus, mass, comp, eos, atm).f;
    const auto fm = surface_residual(minus, mass, comp, eos, atm).f;
    for (std::size_t k = 0; k < 2; ++k) {
      const double num = (fp[k] - fm[k]) / (2e-4);
      worst = std::max(worst, std::abs(num - off.dfdy[k][v] * unit));
    }
  }
  check(worst < 2e-7, "analytic surface Jacobian includes Teff and gravity chains", worst, 0.0);
}

int main() {
  std::printf("ember atmosphere boundaries\n\n");
  GasRadiationEos eos; ConstantOpacity constant;
  Composition comp{}; comp[Species::H1] = 0.7; comp[Species::He4] = 0.3;

  // Constant opacity is exactly soluble, including the finite surface
  // radiation pressure and the radiation force on the material.
  for (double tau : {0.1, 2.0 / 3.0, 10.0}) {
    const GreyAtmosphere grey(eos, constant, {.tau_match = tau});
    const double Teff = 4000.0, gravity = 1e4;
    const auto a = grey.eval(Teff, gravity, comp);
    const double Q = constants::sigma_SB * std::pow(Teff, 4) / constants::c;
    const double P0 = 2.0 * Q / 3.0, P = P0 + gravity * tau / 0.4;
    near(a.T, Teff * std::pow(0.75 * (tau + 2.0 / 3.0), 0.25), 1e-14,
         "grey temperature obeys the Eddington relation");
    near(a.P, P, 2e-10, "hydrostatic pressure includes radiation at tau=0");
    near(a.Pgas, tau * (gravity / 0.4 - Q), 2e-10,
         "gas pressure accounts for outward radiation force");
    near(a.rho, a.Pgas / (GasRadiationEos::gas(comp) * a.T), 2e-10,
         "boundary density comes from gas plus radiation EOS");
    near(a.dlnP_dlnTeff, 4.0 * P0 / P, 2e-8, "Teff pressure sensitivity includes the top boundary");
    near(a.dlnP_dlng, (P - P0) / P, 2e-10, "gravity pressure sensitivity matches hydrostatics");
  }

  // Radiation must also work where it is a substantial fraction of pressure.
  {
    const double Teff = 4000.0, gravity = 1.0;
    const double Q = constants::sigma_SB * std::pow(Teff, 4) / constants::c;
    const auto a = GreyAtmosphere(eos, constant).eval(Teff, gravity, comp);
    const double P = (2.0 / 3.0) * (Q + gravity / 0.4);
    near(a.P, P, 2e-9, "constant-opacity solution remains correct with substantial radiation");
    near(a.dlnP_dlnTeff, (8.0 / 3.0) * Q / P, 2e-8,
         "radiation pressure sensitivity is correct beyond the gas limit");
  }
  {
    PressureOpacity pressure;
    // With kappa=k0*P, P(tau)^2=P(0)^2+2*g*tau/k0. Refining the finite top
    // checks its truncation error independently of the ODE error estimator.
    const double Teff = 4000, gravity = 1e4, tau = 2.0 / 3.0;
    const double P0 = (2.0 / 3.0) * constants::sigma_SB * std::pow(Teff, 4) / constants::c;
    const double expected = std::sqrt(P0 * P0 + 2.0 * gravity * tau / 1e-7);
    const auto coarse = GreyAtmosphere(eos, pressure, {.tau_top = 1e-3}).eval(Teff, gravity, comp);
    const auto fine = GreyAtmosphere(eos, pressure, {.tau_top = 1e-8, .tolerance = 1e-10}).eval(Teff, gravity, comp);
    near(fine.P, expected, 2e-8, "varying-opacity integration matches an analytic atmosphere");
    check(std::abs(fine.P - expected) < 1e-3 * std::abs(coarse.P - expected),
          "reducing tau_top removes the finite-top truncation error",
          std::abs(fine.P - expected) / std::abs(coarse.P - expected), 1e-3);
  }
  {
    VariableOpacity variable;
    GreyAtmosphere grey(eos, variable);
    const double Teff = 3500, gravity = 1e4, step = 1e-4;
    const auto a = grey.eval(Teff, gravity, comp);
    const auto tp = grey.eval(Teff * std::exp(step), gravity, comp);
    const auto tm = grey.eval(Teff * std::exp(-step), gravity, comp);
    const auto gp = grey.eval(Teff, gravity * std::exp(step), comp);
    const auto gm = grey.eval(Teff, gravity * std::exp(-step), comp);
    near(a.dlnP_dlnTeff, (std::log(tp.P) - std::log(tm.P)) / (2 * step), 2e-7,
         "integrated Teff sensitivity includes EOS and opacity response");
    near(a.dlnP_dlng, (std::log(gp.P) - std::log(gm.P)) / (2 * step), 2e-7,
         "integrated gravity sensitivity includes opacity response");
    check_surface(eos, grey, comp);
  }

  const std::string path = std::string(EMBER_TEST_DATA_DIR) + "/synthetic_atmosphere.dat";
  const TabulatedAtmosphere table(eos, path);
  {
    const std::string cond_path=std::string(EMBER_DATA_DIR)+"/atmosphere/cond_gn93_tau100_solar_proxy.dat";
    const auto solar=solar_scaled(.7,.02);
    check(throws([&]{TabulatedAtmosphere exact(eos,cond_path);}),
          "native solar-mixture mismatch requires an explicit proxy",1,1);
    const TabulatedAtmosphere cond(eos,cond_path,TabulatedAtmosphere::Mixture::allow_documented_proxy);
    // Original source cells, independent of the converted log table.
    const auto a=cond.eval(2800,1e5,solar);
    near(a.T,4081.407,2e-14,"COND source temperature is the tau=100 layer, not Teff");
    near(a.Pgas,1.555799e7,2e-14,"COND source material pressure includes electrons, excludes radiation");
    near(a.tau,100,1e-15,"physical atmosphere matches at Rosseland tau=100");
    const double Teff=2845,g=1.8e5,h=1e-5;
    const auto b=cond.eval(Teff,g,solar),tp=cond.eval(Teff*std::exp(h),g,solar),
               tm=cond.eval(Teff*std::exp(-h),g,solar),gp=cond.eval(Teff,g*std::exp(h),solar),
               gm=cond.eval(Teff,g*std::exp(-h),solar);
    near(b.dlnP_dlnTeff,std::log(tp.P/tm.P)/(2*h),2e-8,"non-grey pressure derivative follows the actual interpolation");
    near(b.dlnT_dlng,std::log(gp.T/gm.T)/(2*h),2e-8,"non-grey temperature includes its gravity derivative");
    check(throws([&]{cond.eval(3400,g,solar);}) && throws([&]{cond.eval(2800,1e3,solar);}),
          "grey-filled and source-transition regions are excluded",1,1);
  }
  {
    double worst = 0.0;
    for (double lt : {3.0, 3.2, 3.5, 3.8, 4.0}) {
      for (double lg : {3.0, 3.3, 4.0, 4.4, 5.0}) {
        const auto a = table.eval(std::pow(10.0, lt), std::pow(10.0, lg), comp);
        const double T = std::pow(10.0, 0.9 * lt + 0.1 * lg + 0.1);
        const double Pg = std::pow(10.0, 2.0 * lg - 0.5 * lt - 0.2);
        const double Pr = constants::a_rad * std::pow(T, 4) / 3.0;
        worst = std::max({worst, std::abs(a.T / T - 1.0), std::abs(a.Pgas / Pg - 1.0),
                         std::abs(a.P / (Pg + Pr) - 1.0), std::abs(a.dlnT_dlnTeff - 0.9),
                         std::abs(a.dlnT_dlng - 0.1),
                         std::abs(a.dlnP_dlnTeff - (-0.5 * Pg + 3.6 * Pr) / (Pg + Pr)),
                         std::abs(a.dlnP_dlng - (2.0 * Pg + 0.4 * Pr) / (Pg + Pr))});
      }
    }
    check(worst < 2e-13, "table nodes, cell interiors, and derivatives reproduce power laws", worst, 0.0);
    check_surface(eos, table, comp);
  }
  {
    Composition changed = comp; changed[Species::H1] -= 0.01; changed[Species::He4] += 0.01;
    bool rejected = throws([&] { table.eval(999.0, 1e4, comp); })
        && throws([&] { table.eval(10001.0, 1e4, comp); })
        && throws([&] { table.eval(3000.0, 1e2, comp); })
        && throws([&] { table.eval(3000.0, 1e6, comp); })
        && throws([&] { table.eval(3000.0, 1e4, changed); });
    check(rejected, "table edges and composition mismatches throw", rejected, 1.0);

    std::ifstream in(path); std::ostringstream buffer; buffer << in.rdbuf();
    const std::string good = buffer.str();
    auto bad_table = [&](const std::string& content) {
      return throws([&] { std::istringstream stream(content); TabulatedAtmosphere bad(eos, stream); });
    };
    std::string duplicate = good;
    duplicate.replace(duplicate.find("3.0 3.5 4.0"), 11, "3.0 3.0 4.0");
    std::string bad_comp = good;
    bad_comp.replace(bad_comp.find("0.7 0 0.3"), 9, "0.7 0 0.2");
    rejected = bad_table(duplicate) && bad_table(bad_comp) && bad_table(good + "extra\n")
        && bad_table(good.substr(0, good.find("data") + 4)) && bad_table("EMBER_ATMOSPHERE 2\n");
    check(rejected, "malformed and incomplete tables cannot be loaded", rejected, 1.0);
  }

  // Real modules must interoperate as well as the analytic test physics.
  {
    CompositeEos real_eos;
    FergusonOpacity ferguson(std::string(EMBER_DATA_DIR) + "/opacity/ferguson_gs98_z020.dat");
    GreyAtmosphere grey(real_eos, ferguson);
    const auto mixture = solar_scaled(0.70, 0.020);
    const auto a = grey.eval(4000.0, 100.0, mixture);
    near(real_eos.eval(a.T, a.rho, mixture).P, a.P, 1e-10,
         "grey boundary works with CompositeEos and Ferguson opacity");
    check(std::isfinite(a.dlnP_dlnTeff) && std::isfinite(a.dlnP_dlng),
          "production-module boundary derivatives are finite", a.dlnP_dlng, 1.0);
    constexpr double step = 1e-4;
    const double dT = (std::log(grey.eval(4000.0 * std::exp(step), 100.0, mixture).P)
                     - std::log(grey.eval(4000.0 * std::exp(-step), 100.0, mixture).P)) / (2 * step);
    const double dg = (std::log(grey.eval(4000.0, 100.0 * std::exp(step), mixture).P)
                     - std::log(grey.eval(4000.0, 100.0 * std::exp(-step), mixture).P)) / (2 * step);
    near(a.dlnP_dlnTeff, dT, 2e-4, "Ferguson boundary Teff derivative matches an independent perturbation");
    near(a.dlnP_dlng, dg, 2e-4, "Ferguson boundary gravity derivative matches an independent perturbation");
    check(throws([&] { grey.eval(100.0, 100.0, mixture); }),
          "opacity table errors propagate through atmosphere", 1.0, 1.0);
  }
  {
    CompositeEos real_eos;
    FergusonOpacity low(std::string(EMBER_DATA_DIR) + "/opacity/ferguson_gs98_z020.dat");
    OpalOpacity high(std::string(EMBER_DATA_DIR) + "/opacity/opal_gs98_z020.dat");
    BlendedOpacity blend(low, high);
    GreyAtmosphere grey(real_eos, blend);
    const auto mixture = solar_scaled(.7, .02);
    constexpr double Teff = 16000, gravity = 2e5, step = 1e-5;
    const auto a = grey.eval(Teff, gravity, mixture);
    near(real_eos.eval(a.T, a.rho, mixture).P, a.P, 1e-10, "blended opacity atmosphere closes its EOS pressure");
    const double dt = std::log(grey.eval(Teff * std::exp(step), gravity, mixture).P
                            / grey.eval(Teff * std::exp(-step), gravity, mixture).P) / (2 * step);
    const double dg = std::log(grey.eval(Teff, gravity * std::exp(step), mixture).P
                            / grey.eval(Teff, gravity * std::exp(-step), mixture).P) / (2 * step);
    near(a.dlnP_dlnTeff, dt, 3e-4, "grey atmosphere carries opacity blend's temperature response");
    near(a.dlnP_dlng, dg, 3e-4, "grey atmosphere carries opacity blend's density response");
    check(throws([&] { grey.eval(3000, gravity, mixture); }),
          "cool high-gravity grey atmosphere reports missing density coverage", 1, 1);
  }
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    GreyAtmosphere grey(eos, constant);
    bool rejected = throws([&] { grey.eval(0, 1e4, comp); })
        && throws([&] { grey.eval(4000, nan, comp); })
        && throws([&] { grey.eval(4000, 0.1, comp); })
        && throws([&] { GreyAtmosphere bad(eos, constant, {.tau_top = 1.0}); })
        && throws([&] { GreyAtmosphere(eos, constant, {.max_steps = 1}).eval(4000, 1e4, comp); })
        && throws([&] { surface_residual(Point{}, constants::Msun, comp, eos, grey); });
    check(rejected, "invalid states and exhausted integration throw", rejected, 1.0);
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
