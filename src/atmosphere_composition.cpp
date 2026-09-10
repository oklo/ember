#include "ember/atmosphere_composition.hpp"
#include "differential.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace ember {
namespace {
using D = detail::Differential<2>;
using State = std::array<D, 2>; // ln(total pressure), ln T, with Teff/g sensitivities
bool positive(double x) {
  return std::isfinite(x) && x > 0;
}
template <class F> State rk4(const F& rhs, double x, const State& y, double h) {
  const auto shift = [&](const State& k, double s) { return State{y[0] + s * k[0], y[1] + s * k[1]}; };
  const auto a = rhs(x, y), b = rhs(x + h / 2, shift(a, h / 2)), c = rhs(x + h / 2, shift(b, h / 2)),
             d = rhs(x + h, shift(c, h));
  return {y[0] + h / 6 * (a[0] + 2 * b[0] + 2 * c[0] + d[0]),
          y[1] + h / 6 * (a[1] + 2 * b[1] + 2 * c[1] + d[1])};
}

// q=sqrt(grad-grad_element). Flux conservation and element cooling give
// W=grad_rad-grad_ad=q^2+B*q+C*q^3. All three terms are positive.
D gradient(const D& rad, const D& ad, const D& logB, const D& logC) {
  if (rad.value <= ad.value)
    return rad;
  const double W = rad.value - ad.value, logW = std::log(W);
  const double l2 = .5 * logW, l1 = logW - logB.value, l3 = (logW - logC.value) / 3;
  const double lq = std::min({l2, l1, l3});
  const double A = std::exp(2 * (lq - l2)), B = std::exp(lq - l1), C = std::exp(3 * (lq - l3));
  double lo = 0, hi = 1, t = .5;
  bool converged = false;
  for (int it = 0; it < 80; ++it) {
    const double f = ((C * t + A) * t + B) * t - 1;
    if (std::abs(f) < 2e-15) {
      converged = true;
      break;
    }
    if (f > 0)
      hi = t;
    else
      lo = t;
    const double next = t - f / ((3 * C * t + 2 * A) * t + B);
    t = next > lo && next < hi ? next : .5 * (lo + hi);
  }
  if (!converged)
    throw std::runtime_error("ConvectiveAtmosphere: MLT root failed");
  const double total = ((C * t + A) * t + B) * t, a = A * t * t / total, b = B * t / total,
               c = C * t * t * t / total;
  const double denom = 2 * a + b + 3 * c, r = (2 * a + b) / denom;
  D out(ad.value + W * (a + b));
  for (std::size_t j = 0; j < 2; ++j)
    out.d[j] =
        r * rad.d[j] + (1 - r) * ad.d[j] + W * c / denom * (3 * b * logB.d[j] - (2 * a + b) * logC.d[j]);
  return out;
}
std::array<double, 2> gas_derivatives(const AtmosphereState& s) {
  const double Pr = s.P - s.Pgas;
  return {(s.P * s.dlnP_dlnTeff - 4 * Pr * s.dlnT_dlnTeff) / s.Pgas,
          (s.P * s.dlnP_dlng - 4 * Pr * s.dlnT_dlng) / s.Pgas};
}
} // namespace
ConvectiveAtmosphere::ConvectiveAtmosphere(const Eos& e, const Opacity& k, ConvectiveAtmosphereOptions o)
    : eos_(e), opacity_(k), options_(o) {
  for (double x :
       {o.tau_match, o.tau_start, o.tau_top, o.alpha, o.henyey_y, o.tolerance, o.sensitivity_tolerance})
    if (!positive(x))
      throw std::invalid_argument("ConvectiveAtmosphere: invalid options");
  if (o.tau_top >= o.tau_start || o.tau_start >= o.tau_match || o.tau_start > .1 || o.tolerance < 1e-12 ||
      o.tolerance > 1e-3 || o.sensitivity_tolerance < 1e-12 || o.sensitivity_tolerance > 1e-3 || !o.max_steps)
    throw std::invalid_argument("ConvectiveAtmosphere: invalid integration domain");
}
AtmosphereState ConvectiveAtmosphere::eval(double Teff, double gravity, const Composition& c) const {
  if (!positive(Teff) || !positive(gravity))
    throw std::domain_error("ConvectiveAtmosphere: invalid Teff/g");
  GreyAtmosphere top(eos_, opacity_,
                     {options_.tau_start, options_.tau_top, options_.tolerance, options_.max_steps});
  const auto initial = top.eval(Teff, gravity, c);
  State y{D(std::log(initial.P)), D(std::log(initial.T))};
  y[0].d = {initial.dlnP_dlnTeff, initial.dlnP_dlng};
  y[1].d = {initial.dlnT_dlnTeff, initial.dlnT_dlng};
  const auto lng = D::variable(std::log(gravity), 1), lnTeff = D::variable(std::log(Teff), 0);
  const D flux = constants::sigma_SB * exp(4 * lnTeff);
  const auto rhs = [&](double x, const State& s) {
    const D P = exp(s[0]), T = exp(s[1]);
    const double rho = eos_.rho_from_PT(T.value, P.value, c);
    const auto e = eos_.eval_with_derivatives(T.value, rho, c);
    const auto& v = e.state;
    const auto k = opacity_.eval(T.value, rho, c);
    if (!positive(v.chiRho) || !positive(v.cp) || !positive(v.delta) || !positive(k.kappa))
      throw std::domain_error("ConvectiveAtmosphere: invalid material response");
    D lnrho(std::log(rho));
    for (std::size_t j = 0; j < 2; ++j)
      lnrho.d[j] = (s[0].d[j] - v.chiT * s[1].d[j]) / v.chiRho;
    const auto material = [&](double value, double dT, double drho) {
      D a(value);
      for (std::size_t j = 0; j < 2; ++j)
        a.d[j] = dT * s[1].d[j] + drho * lnrho.d[j];
      return a;
    };
    const D lnk = material(std::log(k.kappa), k.dlnk_dlnT, k.dlnk_dlnRho);
    const D lncp = material(std::log(v.cp), e.dcp_dlnT / v.cp, e.dcp_dlnRho / v.cp);
    const D lndelta = material(std::log(v.delta), e.ddelta_dlnT / v.delta, e.ddelta_dlnRho / v.delta);
    const D ad = material(v.grad_ad, e.dgrad_ad_dlnT, e.dgrad_ad_dlnRho);
    const D lnHp = s[0] - lnrho - lng;
    const D lnell = std::log(options_.alpha) + lnHp;
    const D lnv0 = std::log(options_.alpha) + .5 * (lng + lnHp + lndelta - std::log(8.));
    const D lntauell = lnrho + lnk + lnell;
    // Henyey (nu=8), Gustafsson et al. 2008 eqs. 13--18.
    const D logB = std::log(8 * constants::sigma_SB) + 3 * s[1] + lntauell - lnv0 - lnrho - lncp -
                   log(1 + options_.henyey_y * exp(2 * lntauell));
    const D logC = std::log(3 * options_.alpha / (32 * constants::sigma_SB)) + lnk + s[0] + lnrho + lncp +
                   lnv0 - 3 * s[1] - lng;
    const D rad = 3 * exp(lnk) * P * flux / (16 * constants::sigma_SB * exp(4 * s[1] + lng));
    const D grad = gradient(rad, ad, logB, logC);
    const D hydro = exp(x + lng - lnk - s[0]);
    return State{hydro, hydro * grad};
  };
  double x = std::log(options_.tau_start), end = std::log(options_.tau_match), h = .05;
  std::size_t steps = 0;
  while (x < end) {
    if (++steps > options_.max_steps)
      throw std::runtime_error("ConvectiveAtmosphere: step limit");
    h = std::min(h, end - x);
    if (x + h == x)
      throw std::runtime_error("ConvectiveAtmosphere: step underflow");
    const auto coarse = rk4(rhs, x, y, h), half = rk4(rhs, x, y, h / 2),
               fine = rk4(rhs, x + h / 2, half, h / 2);
    double error = 0;
    for (std::size_t i = 0; i < 2; ++i) {
      error = std::max(error, std::abs(fine[i].value - coarse[i].value) / (15 * options_.tolerance));
      for (std::size_t j = 0; j < 2; ++j)
        error =
            std::max(error, std::abs(fine[i].d[j] - coarse[i].d[j]) /
                                (15 * options_.sensitivity_tolerance * std::max(1., std::abs(fine[i].d[j]))));
    }
    if (!std::isfinite(error))
      throw std::runtime_error("ConvectiveAtmosphere: nonfinite integration");
    if (error <= 1) {
      for (std::size_t i = 0; i < 2; ++i)
        y[i] = fine[i] + (fine[i] - coarse[i]) / 15;
      x += h;
    }
    h = std::min(.25, h * (error == 0 ? 2 : std::clamp(.9 * std::pow(1 / error, .2), .2, 2.)));
  }
  AtmosphereState out{};
  out.P = std::exp(y[0].value);
  out.T = std::exp(y[1].value);
  out.Pgas = out.P - constants::a_rad * std::pow(out.T, 4) / 3;
  out.tau = options_.tau_match;
  if (!positive(out.Pgas))
    throw std::domain_error("ConvectiveAtmosphere: radiation-dominated boundary");
  out.rho = eos_.rho_from_PT(out.T, out.P, c);
  out.dlnP_dlnTeff = y[0].d[0];
  out.dlnP_dlng = y[0].d[1];
  out.dlnT_dlnTeff = y[1].d[0];
  out.dlnT_dlng = y[1].d[1];
  return out;
}
CompositionCorrectedAtmosphere::CompositionCorrectedAtmosphere(const Eos& e, const TabulatedAtmosphere& t,
                                                               const Atmosphere& a, Composition c)
    : eos_(e), table_(t), column_(a), reference_(c) {
  if (std::abs(c.sum() - 1) > 1e-10)
    throw std::invalid_argument("CompositionCorrectedAtmosphere: reference not normalized");
  for (double v : c.X)
    if (!std::isfinite(v) || v < 0)
      throw std::invalid_argument("CompositionCorrectedAtmosphere: invalid reference");
}
AtmosphereState CompositionCorrectedAtmosphere::eval(double Teff, double gravity,
                                                     const Composition& c) const {
  if (c.basis != reference_.basis)
    throw std::domain_error("CompositionCorrectedAtmosphere: changed abundance basis");
  for (std::size_t j = 3; j < NSPEC; ++j)
    if (std::abs(c.X[j] - reference_.X[j]) > 1e-10)
      throw std::domain_error("CompositionCorrectedAtmosphere: correction requires fixed metal pattern");
  const auto anchor = table_.eval(Teff, gravity, table_.composition());
  const auto ref = column_.eval(Teff, gravity, reference_), actual = column_.eval(Teff, gravity, c);
  if (ref.tau != actual.tau || ref.tau != anchor.tau)
    throw std::domain_error("CompositionCorrectedAtmosphere: optical-depth mismatch");
  const auto a = gas_derivatives(anchor), b = gas_derivatives(ref), d = gas_derivatives(actual);
  AtmosphereState out{};
  out.tau = anchor.tau;
  out.T = anchor.T * actual.T / ref.T;
  out.Pgas = anchor.Pgas * actual.Pgas / ref.Pgas;
  out.dlnT_dlnTeff = anchor.dlnT_dlnTeff + actual.dlnT_dlnTeff - ref.dlnT_dlnTeff;
  out.dlnT_dlng = anchor.dlnT_dlng + actual.dlnT_dlng - ref.dlnT_dlng;
  const double Pr = constants::a_rad * std::pow(out.T, 4) / 3;
  out.P = out.Pgas + Pr;
  out.dlnP_dlnTeff = (out.Pgas * (a[0] + d[0] - b[0]) + 4 * Pr * out.dlnT_dlnTeff) / out.P;
  out.dlnP_dlng = (out.Pgas * (a[1] + d[1] - b[1]) + 4 * Pr * out.dlnT_dlng) / out.P;
  out.rho = eos_.rho_from_PT(out.T, out.P, c);
  return out;
}
} // namespace ember
