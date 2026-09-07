#include "ember/convection.hpp"
#include "ember/constants.hpp"
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
  std::printf("  [%s] %-58s got %-12.5g want %-12.5g\n",
              ok ? "PASS" : "FAIL", what.c_str(), got, want);
}
static void near(double got, double want, double tol, const std::string& what) {
  check(std::abs(got - want) <= tol * std::abs(want), what, got, want);
}
template<class F> static bool domain_error(F&& f) {
  try { f(); } catch (const std::domain_error&) { return true; }
  return false;
}

int main() {
  std::printf("ember mixing-length convection\n\n");

  // Schwarzschild stability includes zero and inward flux. Alpha/U must
  // never change a stable layer's radiative gradient or luminosity response.
  for (double rad : {-1.0, 0.0, 0.2, 0.4}) {
    const auto s = mixing_length_gradient(rad, 0.4, 1e-10);
    check(!s.unstable && s.grad == rad && s.convective_excess == 0.0
          && s.dgrad_dgrad_rad == 1.0 && s.dgrad_dgrad_ad == 0.0
          && s.dgrad_dlnU == 0.0, "stable layer is exactly radiative", s.grad, rad);
  }

  // A hand-solvable instance: U=1 and q=2 imply element cooling=4,
  // element contrast=4, and convective flux (in gradient units)=9.
  {
    const auto s = mixing_length_gradient(17.4, 0.4, 1.0);
    near(s.grad, 8.4, 2e-14, "finite-efficiency background gradient");
    near(s.grad_element, 4.4, 2e-14, "finite-efficiency element gradient");
    near(s.convective_excess, 9.0, 2e-14, "finite-efficiency convective flux");
  }

  // Over 500 decades in U, the physical gradient is bounded, the partition
  // of the luminosity closes, and the response to increasing U is monotone.
  {
    bool bounded = true, monotone = true, derivatives = true;
    double worst = 0.0;
    for (double W : {1e-12, 1.0, 1e6, 1e12}) {
      const double rad = 0.4 + W, actualW = rad - 0.4;
      double previous = 0.0;
      for (int exponent = -250; exponent <= 250; ++exponent) {
        const auto s = mixing_length_gradient(rad, 0.4, std::pow(10.0, exponent));
        bounded = bounded && s.unstable && std::isfinite(s.grad)
            && s.grad >= 0.4 && s.grad <= rad
            && s.grad_element >= 0.4 && s.grad_element <= s.grad
            && s.element_contrast >= 0.0 && s.superadiabaticity >= 0.0
            && s.convective_excess >= 0.0;
        monotone = monotone && s.superadiabaticity >= previous * (1.0 - 1e-13);
        previous = s.superadiabaticity;
        derivatives = derivatives && s.dgrad_dgrad_rad >= 0.0
            && s.dgrad_dgrad_rad <= 1.0 && s.dgrad_dgrad_ad >= 0.0
            && s.dgrad_dgrad_ad <= 1.0 && std::isfinite(s.dgrad_dlnU)
            && s.dgrad_dlnU >= 0.0;
        worst = std::max(worst,
            std::abs((s.superadiabaticity + s.convective_excess) / actualW - 1.0));
      }
    }
    check(bounded, "gradient stays between adiabatic and radiative", bounded, 1.0);
    check(monotone, "larger radiative loss raises superadiabaticity", monotone, 1.0);
    check(derivatives, "analytic response has the physical signs", derivatives, 1.0);
    check(worst < 2e-15, "radiation and convection share the imposed flux", worst, 0.0);
  }

  // Efficient limit: q^3 = 8*U*W/9 and grad-grad_ad ~ q^2. Retain this
  // small excess even when adding it to grad_ad rounds back to grad_ad.
  {
    const double U = 1e-30, W = 1.0;
    const auto s = mixing_length_gradient(1.4, 0.4, U);
    near(s.superadiabaticity, std::pow(8.0 * U * W / 9.0, 2.0 / 3.0), 2e-12,
         "efficient-limit superadiabatic excess survives rounding");
    near(s.dgrad_dlnU / s.superadiabaticity, 2.0 / 3.0, 2e-12,
         "efficient-limit excess scales as U^(2/3)");
  }
  // Inefficient limit: q=W/(2U) and radiation carries all but an O(U^-4)
  // fraction. Testing the stored flux catches cancellation hidden by grad.
  {
    const double U = 1e8;
    const auto s = mixing_length_gradient(1.4, 0.4, U);
    near(s.element_contrast, 1.0 / (4.0 * U * U), 2e-12,
         "inefficient-limit element contrast survives rounding");
    near(s.convective_excess, 9.0 / (64.0 * std::pow(U, 4)), 2e-12,
         "inefficient-limit convective flux scales as U^-4");
    near(s.grad, 1.4, 2e-15, "inefficient convection approaches radiation");
  }

  // Check the physical, dimensional equations independently of the cubic:
  // buoyancy velocity, enthalpy flux, and diffusive cooling of the element.
  // BV58 coefficients (1/8, 1/2, 24): Salaris & Cassisi 2008, eqs. 1-3.
  {
    const double T = 1e5, rho = 1e-6, kappa = 10.0, gravity = 1e4, alpha = 1.9;
    EosState e{};
    e.P = 1e7; e.cp = 2.5e8; e.delta = 1.0; e.grad_ad = 0.4;
    const double Hp = e.P / (rho * gravity), ell = alpha * Hp;
    const double U = mixing_length_U(T, rho, kappa, gravity, e, alpha);
    const auto s = mixing_length_gradient(5.0, e.grad_ad, U);
    const double v = std::sqrt(ell * ell * gravity * e.delta * s.element_contrast
                              / (8.0 * Hp));
    const double Fc = rho * v * e.cp * T * ell * s.element_contrast / (2.0 * Hp);
    const double K = 16.0 * constants::sigma_SB * std::pow(T, 4)
                   / (3.0 * kappa * rho * Hp);
    near(K * s.grad + Fc, K * 5.0, 2e-12, "dimensional radiation + enthalpy flux = total flux");
    const double gamma = e.cp * rho * rho * ell * v * kappa
                        / (24.0 * constants::sigma_SB * std::pow(T, 3));
    near(s.element_contrast / (s.grad_element - e.grad_ad), gamma, 2e-12,
         "element cooling agrees with radiative diffusion");
    near(mixing_length_U(T, rho, kappa, gravity, e, 2.0 * alpha), U / 4.0, 2e-14,
         "doubling mixing length lowers U by four");
  }

  // Independent differences across the transition between the limits.
  {
    double worst = 0.0;
    constexpr double h = 1e-4;
    for (double U : {1e-5, 0.01, 0.3, 1.0, 10.0}) {
      const auto s = mixing_length_gradient(1.4, 0.4, U);
      const double dr = (mixing_length_gradient(1.4 + h, 0.4, U).grad
                       - mixing_length_gradient(1.4 - h, 0.4, U).grad) / (2.0 * h);
      const double da = (mixing_length_gradient(1.4, 0.4 + h, U).grad
                       - mixing_length_gradient(1.4, 0.4 - h, U).grad) / (2.0 * h);
      const double du = (mixing_length_gradient(1.4, 0.4, U * std::exp(h)).grad
                       - mixing_length_gradient(1.4, 0.4, U * std::exp(-h)).grad) / (2.0 * h);
      worst = std::max({worst, std::abs(dr - s.dgrad_dgrad_rad) / s.dgrad_dgrad_rad,
                       std::abs(da - s.dgrad_dgrad_ad) / s.dgrad_dgrad_ad,
                       std::abs(du - s.dgrad_dlnU) / s.dgrad_dlnU});
    }
    check(worst < 2e-6, "analytic gradient derivatives match independent differences", worst, 0.0);
    const auto edge = mixing_length_gradient(0.4 + 1e-10, 0.4, 0.1);
    near(edge.dgrad_dgrad_rad, 1.0, 1e-12, "radiative derivative is continuous at Schwarzschild edge");
  }

  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    bool rejected = true;
    for (double U : {0.0, -1.0, nan, inf})
      rejected = rejected && domain_error([&] { mixing_length_gradient(0.2, 0.4, U); });
    rejected = rejected && domain_error([&] { mixing_length_gradient(nan, 0.4, 1.0); })
        && domain_error([&] { mixing_length_gradient(1.0, -0.4, 1.0); });
    EosState e{}; e.P = 1e7; e.cp = 2.5e8; e.delta = 1.0;
    for (double bad : {0.0, -1.0, nan, inf}) {
      rejected = rejected && domain_error([&] { mixing_length_U(1e5, 1e-6, 10, bad, e, 1.9); })
          && domain_error([&] { mixing_length_U(1e5, 1e-6, 10, 1e4, e, bad); });
    }
    check(rejected, "invalid states throw, including in stable layers", rejected, 1.0);
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
