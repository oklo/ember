#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ember {
namespace {
bool positive(double x) { return std::isfinite(x) && x > 0.0; }
using State = std::array<double, 3>;  // ln Pgas, dlnPgas/dlnTeff, dlnPgas/dlng

template<class F> State rk4(const F& rhs, double x, const State& y, double step) {
  auto shifted = [&](const State& slope, double factor) {
    State out{};
    for (std::size_t i = 0; i < y.size(); ++i) out[i] = y[i] + factor * slope[i];
    return out;
  };
  const auto k1 = rhs(x, y);
  const auto k2 = rhs(x + 0.5 * step, shifted(k1, 0.5 * step));
  const auto k3 = rhs(x + 0.5 * step, shifted(k2, 0.5 * step));
  const auto k4 = rhs(x + step, shifted(k3, step));
  State out{};
  for (std::size_t i = 0; i < y.size(); ++i)
    out[i] = y[i] + step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
  return out;
}
} // namespace

GreyAtmosphere::GreyAtmosphere(const Eos& eos, const Opacity& opacity,
                               GreyAtmosphereOptions options)
    : eos_(eos), opacity_(opacity), options_(options) {
  if (!positive(options.tau_match) || !positive(options.tau_top)
      || options.tau_top >= options.tau_match || !positive(options.tolerance)
      || options.tolerance < 1e-12 || options.tolerance > 1e-2 || options.max_steps == 0)
    throw std::invalid_argument("GreyAtmosphere: invalid integration options");
}

AtmosphereState GreyAtmosphere::eval(double Teff, double gravity, const Composition& comp) const {
  if (!positive(Teff) || !positive(gravity))
    throw std::domain_error("GreyAtmosphere: Teff and gravity must be positive and finite");
  const double Q = constants::sigma_SB * std::pow(Teff, 4) / constants::c;
  if (!positive(Q)) throw std::domain_error("GreyAtmosphere: flux outside representable range");
  auto temperature = [&](double tau) { return Teff * std::pow(0.75 * (tau + 2.0 / 3.0), 0.25); };

  struct Local {
    double Pg, P, Pr, T, rho, kappa, k_u, k_T;
  };
  auto local = [&](double tau, double logPg) {
    Local s{};
    s.Pg = std::exp(logPg);
    s.T = temperature(tau);
    s.Pr = Q * (tau + 2.0 / 3.0);
    s.P = s.Pg + s.Pr;
    if (!positive(s.Pg) || !positive(s.P) || s.P == s.Pr)
      throw std::domain_error("GreyAtmosphere: gas pressure is not representable");
    const double guess = s.Pg / (constants::R_gas * s.T
        * (comp.mu_ions_inv() + comp.mu_elec_inv()));
    s.rho = eos_.rho_from_PT(s.T, s.P, comp, guess);
    const auto e = eos_.eval(s.T, s.rho, comp);
    OpacityState k;
    try { k = opacity_.eval(s.T, s.rho, comp); }
    catch (const std::domain_error& e) {
      throw std::domain_error(std::string("GreyAtmosphere: opacity lookup: ") + e.what());
    }
    if (!positive(e.chiRho) || !std::isfinite(e.chiT) || !positive(k.kappa)
        || !std::isfinite(k.dlnk_dlnRho) || !std::isfinite(k.dlnk_dlnT))
      throw std::domain_error("GreyAtmosphere: invalid EOS or opacity derivatives");
    s.kappa = k.kappa;
    // Opacity partials at fixed tau, with u=ln Pgas. The Teff partial holds
    // Pgas fixed, not total P: radiation contributes 4*Pr/P to dlnP/dlnTeff.
    s.k_u = k.dlnk_dlnRho * (s.Pg / s.P) / e.chiRho;
    s.k_T = k.dlnk_dlnT + k.dlnk_dlnRho * (4.0 * s.Pr / s.P - e.chiT) / e.chiRho;
    return s;
  };

  const double tau0 = options_.tau_top;
  State y{std::log(gravity) + std::log(tau0), 0.0, 0.0};  // kappa=1 initial guess
  auto top_residual = [&](const Local& s) {
    return std::log(s.Pg + tau0 * Q) + std::log(s.kappa) - std::log(tau0) - std::log(gravity);
  };
  double lo = -std::numeric_limits<double>::infinity();
  double hi = std::numeric_limits<double>::infinity();
  auto bounds = opacity_.density_range(temperature(tau0), comp);
  if (const auto eos_bounds = eos_.density_range(temperature(tau0), comp)) {
    if (bounds) {
      bounds->min = std::max(bounds->min, eos_bounds->min);
      bounds->max = std::min(bounds->max, eos_bounds->max);
    } else bounds = Opacity::DensityRange{eos_bounds->min, eos_bounds->max};
  }
  if (bounds) {
    if (!positive(bounds->min) || !positive(bounds->max) || bounds->min >= bounds->max)
      throw std::domain_error("GreyAtmosphere: EOS/opacity density domains do not overlap");
    const double Pr = Q * (tau0 + 2.0 / 3.0);
    // Move slightly inside the declared edges so the EOS inversion's rounding
    // cannot turn an endpoint into an out-of-table density.
    const double Pg_lo = eos_.eval(temperature(tau0), bounds->min * (1.0 + 1e-6), comp).P - Pr;
    const double Pg_hi = eos_.eval(temperature(tau0), bounds->max * (1.0 - 1e-6), comp).P - Pr;
    if (!positive(Pg_lo) || !positive(Pg_hi) || Pg_lo >= Pg_hi)
      throw std::domain_error("GreyAtmosphere: opacity bounds cannot bracket gas pressure");
    lo = std::log(Pg_lo); hi = std::log(Pg_hi);
    if (top_residual(local(tau0, lo)) > 0.0 || top_residual(local(tau0, hi)) < 0.0)
      throw std::domain_error("GreyAtmosphere: top pressure outside EOS/opacity support; check tau_top");
    y[0] = std::clamp(y[0], lo, hi);
  }
  bool found_top = false;
  for (int it = 0; it < 100; ++it) {
    const auto s = local(tau0, y[0]);
    // Positive-term residual for Pg = tau0*(g/kappa-Q), avoiding a negative
    // pressure during Newton's trial evaluations.
    const double f = top_residual(s);
    const double df = s.Pg / (s.Pg + tau0 * Q) + s.k_u;
    if (std::abs(f) < 1e-11) { found_top = true; break; }
    if (!(df > 0.0) || !std::isfinite(df))
      throw std::domain_error("GreyAtmosphere: top pressure has no monotone local root");
    if (f > 0.0) hi = y[0]; else lo = y[0];
    const double next = y[0] - std::clamp(f / df, -2.0, 2.0);
    y[0] = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
  }
  if (!found_top) throw std::runtime_error("GreyAtmosphere: top pressure did not converge");

  // The ODE and its independent partials in x=ln tau, u=ln Pgas.
  auto coefficients = [&](double tau, const Local& s) {
    const double radiative = gravity / s.kappa;
    if (!(radiative > Q))
      throw std::domain_error("GreyAtmosphere: radiative acceleration reaches gravity");
    const double f = tau * (radiative - Q) / s.Pg;
    const double f_u = -f - tau * radiative * s.k_u / s.Pg;
    const double f_T = -tau * (radiative * s.k_T + 4.0 * Q) / s.Pg;
    const double f_g = tau * radiative / s.Pg;
    return std::array{f, f_u, f_T, f_g};
  };
  {
    const auto a = coefficients(tau0, local(tau0, y[0]));
    // Differentiate the top condition f=1, rather than treating the starting
    // pressure as independent of Teff and gravity.
    y[1] = -a[2] / a[1];
    y[2] = -a[3] / a[1];
  }
  auto rhs = [&](double x, const State& state) {
    const double tau = std::exp(x);
    const auto a = coefficients(tau, local(tau, state[0]));
    return State{a[0], a[1] * state[1] + a[2], a[1] * state[2] + a[3]};
  };

  double x = std::log(tau0), end = std::log(options_.tau_match);
  double step = std::min(0.1, end - x);
  std::size_t steps = 0;
  while (x < end) {
    if (++steps > options_.max_steps)
      throw std::runtime_error("GreyAtmosphere: integration step limit exceeded");
    step = std::min(step, end - x);
    if (x + step == x)
      throw std::runtime_error("GreyAtmosphere: integration step underflow");
    const auto coarse = rk4(rhs, x, y, step);
    const auto half = rk4(rhs, x, y, 0.5 * step);
    const auto fine = rk4(rhs, x + 0.5 * step, half, 0.5 * step);
    double error = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i) {
      if (!std::isfinite(coarse[i]) || !std::isfinite(fine[i]))
        throw std::runtime_error("GreyAtmosphere: non-finite integration state");
      const double scale = i == 0 ? 1.0 : std::max(1.0, std::abs(fine[i]));
      error = std::max(error, std::abs(fine[i] - coarse[i]) / (15.0 * scale));
    }
    if (error <= options_.tolerance) {
      for (std::size_t i = 0; i < y.size(); ++i)
        y[i] = fine[i] + (fine[i] - coarse[i]) / 15.0;
      x += step;
    }
    const double factor = error == 0.0 ? 2.0
        : std::clamp(0.9 * std::pow(options_.tolerance / error, 0.2), 0.2, 2.0);
    step = std::min(0.5, step * factor);
  }

  const auto s = local(options_.tau_match, y[0]);
  AtmosphereState out{};
  out.T = s.T; out.P = s.P; out.Pgas = s.Pg; out.rho = s.rho; out.tau = options_.tau_match;
  out.dlnT_dlnTeff = 1.0;
  out.dlnP_dlnTeff = (s.Pg / s.P) * y[1] + 4.0 * s.Pr / s.P;
  out.dlnP_dlng = (s.Pg / s.P) * y[2];
  return out;
}

} // namespace ember
