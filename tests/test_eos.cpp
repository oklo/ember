// Limit tests for the equation of state.  Each one is a place where the answer
// is known independently of the code, so a regression shows up as a physics
// statement rather than a changed number.
#include "ember/eos.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <cstdio>
#include <string>

using namespace ember;
using namespace ember::constants;

static int failures = 0;
static void check(bool ok, const std::string& what, double got, double want, double tol) {
  const double rel = std::abs(got - want) / (std::abs(want) > 0 ? std::abs(want) : 1.0);
  if (!ok) ++failures;
  std::printf("  [%s] %-46s got %-13.6g want %-13.6g (rel %.2e, tol %.0e)\n",
              ok ? "PASS" : "FAIL", what.c_str(), got, want, rel, tol);
}
static void near(double got, double want, double tol, const std::string& what) {
  check(std::abs(got - want) <= tol * std::abs(want), what, got, want, tol);
}

int main() {
  IdealEos eos;
  const Composition h_pure = [] { Composition cc{}; cc[Species::H1] = 1.0; return cc; }();
  const Composition solar = solar_scaled(0.70, 0.014);

  std::printf("ember equation of state - limit tests (%s)\n\n", eos.name());

  // 1. Non-degenerate, gas-dominated: P -> (n_i + n_e) R rho T exactly.
  //    The point has to earn the name.  At T=1e6, rho=1e-4 radiation is 13%
  //    of the pressure and chi_T is legitimately 1.40, which looks like a bug
  //    and is not one; here P_rad/P_gas ~ 1e-4 and e^eta ~ 4e-4, so both
  //    corrections sit below the tolerance.
  {
    const double T = 1.0e5, rho = 1.0e-4;
    const auto s = eos.eval(T, rho, h_pure);
    const double Pid = (h_pure.mu_ions_inv() + h_pure.mu_elec_inv()) * R_gas * rho * T;
    near(s.P, Pid, 2e-3, "ideal limit: P = (ni+ne) R rho T");
    near(s.chiT, 1.0, 2e-3, "ideal limit: chi_T = 1");
    near(s.chiRho, 1.0, 2e-3, "ideal limit: chi_rho = 1");
    near(s.Gamma1, 5.0 / 3.0, 3e-3, "ideal limit: Gamma_1 = 5/3");
    near(s.grad_ad, 0.4, 3e-3, "ideal limit: grad_ad = 2/5");
  }

  // 2. Radiation-dominated: Gamma_1 -> 4/3, grad_ad -> 1/4.
  {
    const double T = 1.0e8, rho = 1.0e-6;
    const auto s = eos.eval(T, rho, h_pure);
    near(s.Gamma1, 4.0 / 3.0, 2e-2, "radiation limit: Gamma_1 -> 4/3");
    near(s.grad_ad, 0.25, 3e-2, "radiation limit: grad_ad -> 1/4");
  }

  // 3. Strong non-relativistic degeneracy: P_e -> K rho^{5/3}, so chi_rho -> 5/3
  //    and the pressure becomes nearly independent of temperature.
  {
    const double T = 1.0e6, rho = 1.0e4;
    const auto s = eos.eval(T, rho, h_pure);
    near(s.chiRho, 5.0 / 3.0, 3e-2, "degenerate limit: chi_rho -> 5/3");
    check(s.chiT < 0.1, "degenerate limit: chi_T << 1", s.chiT, 0.0, 0.1);
    // Absolute check against the zero-temperature Chandrasekhar constant.
    const double K = 1.0036e13; // dyn cm^-2 (g/cm^3)^{-5/3}, mu_e = 1
    const double mue = 1.0 / h_pure.mu_elec_inv();
    const double Pdeg = K * std::pow(rho / mue, 5.0 / 3.0);
    near(s.P, Pdeg, 5e-2, "degenerate limit: P = K (rho/mu_e)^5/3");
  }

  // 4. Thermodynamic consistency: the Maxwell relation behind grad_ad.
  //    cp - cv = P delta^2 /(rho T chi_rho) must hold identically.
  {
    const double T = 3.0e6, rho = 1.0e-2;
    const auto s = eos.eval(T, rho, solar);
    const double lhs = s.cp - s.cv;
    // cp - cv = P chi_T^2 /(rho T chi_rho).  (Writing delta^2 here instead of
    // chi_T^2 is only right when chi_rho = 1, which is exactly the regime that
    // would hide the error.)
    const double rhs = s.P * s.chiT * s.chiT / (rho * T * s.chiRho);
    near(lhs, rhs, 1e-10, "consistency: cp - cv = P chi_T^2/(rho T chi_rho)");
  }

  // 5. Analytic derivatives must match centred differences of the code itself.
  {
    const double T = 2.0e6, rho = 1.0e-3, d = 1e-5;
    const auto s = eos.eval(T, rho, solar);
    const double Pp = eos.eval(T * (1 + d), rho, solar).P;
    const double Pm = eos.eval(T * (1 - d), rho, solar).P;
    near(s.chiT, (std::log(Pp) - std::log(Pm)) / (2 * d), 1e-6, "chi_T vs numerical");
    const double Rp = eos.eval(T, rho * (1 + d), solar).P;
    const double Rm = eos.eval(T, rho * (1 - d), solar).P;
    near(s.chiRho, (std::log(Rp) - std::log(Rm)) / (2 * d), 1e-6, "chi_rho vs numerical");
    const double Ep = eos.eval(T * (1 + d), rho, solar).E;
    const double Em = eos.eval(T * (1 - d), rho, solar).E;
    near(s.cv, (Ep - Em) / (2 * d * T), 1e-6, "c_v vs numerical");
  }

  // 6. The (T,P) inversion must return the density it was given.
  {
    const double T = 5.0e5, rho = 3.7e-3;
    const auto s = eos.eval(T, rho, solar);
    near(eos.rho_from_PT(T, s.P, solar), rho, 1e-10, "inversion: rho(T, P(T,rho)) = rho");
  }

  // 7. Composition bookkeeping.
  {
    const auto cs = solar_scaled(0.70, 0.014);
    near(cs.sum(), 1.0, 1e-12, "composition: mass fractions sum to one");
    near(cs.Z(), 0.014, 1e-12, "composition: Z as requested");
    near(cs.h1(), 0.70, 1e-12, "composition: X as requested");
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
