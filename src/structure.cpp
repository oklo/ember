#include "ember/structure.hpp"
#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include "differential.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace ember {
namespace {
using constants::G;
using detail::Differential;
using detail::exp;
using detail::log;

struct Previous { double E, rho; };
std::optional<Previous> prepare(const Model& m, std::size_t i, const Physics& phys,
                                double dt, const Model* prev) {
  if (!phys.eos || !phys.opacity || !phys.nuclear)
    throw std::invalid_argument("zone_residual: missing physics module");
  if (m.size() < 2 || i >= m.size() - 1 || m.m.size() != m.size() || m.comp.size() != m.size())
    throw std::invalid_argument("zone_residual: invalid zone index or model arrays");
  if (!std::isfinite(m.m[i]) || !std::isfinite(m.m[i + 1]) || m.m[i] < 0.0
      || !(m.m[i + 1] > m.m[i]) || !std::isfinite(dt))
    throw std::domain_error("zone_residual: invalid mass interval or time step");
  if (dt <= 0.0) return std::nullopt;
  if (!prev || prev->size() != m.size() || prev->comp.size() != m.size() || prev->m != m.m)
    throw std::invalid_argument("zone_residual: time dependence requires a previous model on the same mesh");
  const double rho = prev->rho(i);
  return Previous{phys.eos->eval(prev->T(i), rho, prev->comp[i]).E, rho};
}

template<std::size_t N> struct Local {
  Differential<N> lnr, lnT, r, rho, T, L, P, E, cp, delta, grad_ad, kappa, eps;
};

template<std::size_t N>
Local<N> gather(const Point& point, const Composition& comp, const Physics& phys, std::size_t offset) {
  using D = Differential<N>;
  for (double x : {point.lnr, point.lnrho, point.lnT, point.L})
    if (!std::isfinite(x)) throw std::domain_error("zone_residual: non-finite mesh point");
  Local<N> q{};
  q.lnr = D::variable(point.lnr, offset);
  q.r = exp(q.lnr);
  q.rho = exp(D::variable(point.lnrho, offset + 1));
  q.lnT = D::variable(point.lnT, offset + 2);
  q.T = exp(q.lnT);
  q.L = D::variable(point.L, offset + 3);
  EosResponse response{};
  if constexpr (N > 0) response = phys.eos->eval_with_derivatives(q.T.value, q.rho.value, comp);
  else response.state = phys.eos->eval(q.T.value, q.rho.value, comp);
  const auto& e = response.state;
  const auto k = phys.opacity->eval(q.T.value, q.rho.value, comp);
  const auto n = phys.nuclear->eval(q.T.value, q.rho.value, comp);
  if (!(e.P > 0.0) || !(e.cp > 0.0) || !(e.delta > 0.0) || !(k.kappa > 0.0))
    throw std::domain_error("zone_residual: invalid thermodynamic state or opacity");
  auto material = [&](double value, double dT, double drho) {
    D out(value);
    if (!std::isfinite(value)) throw std::domain_error("zone_residual: non-finite physics value");
    if constexpr (N > 0) {
      if (!std::isfinite(dT) || !std::isfinite(drho))
        throw std::domain_error("zone_residual: non-finite physics derivative");
      out.d[offset + 2] = dT; out.d[offset + 1] = drho;
    }
    return out;
  };
  q.P = material(e.P, e.P * e.chiT, e.P * e.chiRho);
  q.E = material(e.E, e.cv * q.T.value, response.dE_dlnRho);
  q.cp = material(e.cp, response.dcp_dlnT, response.dcp_dlnRho);
  q.delta = material(e.delta, response.ddelta_dlnT, response.ddelta_dlnRho);
  q.grad_ad = material(e.grad_ad, response.dgrad_ad_dlnT, response.dgrad_ad_dlnRho);
  q.kappa = material(k.kappa, k.kappa * k.dlnk_dlnT, k.kappa * k.dlnk_dlnRho);
  q.eps = material(n.eps, n.eps * n.dlneps_dlnT, n.eps * n.dlneps_dlnRho);
  return q;
}

template<std::size_t N>
std::array<Differential<N>, NVAR> equations(const Model& model, std::size_t i,
    const Point& lo, const Point& hi, const Physics& phys, double dt, const std::optional<Previous>& prev) {
  using D = Differential<N>;
  const auto a = gather<N>(lo, model.comp[i], phys, 0);
  const auto b = gather<N>(hi, model.comp[i + 1], phys, NVAR);
  const double dm = model.m[i + 1] - model.m[i];
  const double mb = 0.5 * (model.m[i] + model.m[i + 1]);
  const D rb = 0.5 * (a.r + b.r), rhob = 0.5 * (a.rho + b.rho);
  const D Tb = 0.5 * (a.T + b.T), Pb = 0.5 * (a.P + b.P);
  const D kb = 0.5 * (a.kappa + b.kappa), Lb = 0.5 * (a.L + b.L);
  const D cp = 0.5 * (a.cp + b.cp), delta = 0.5 * (a.delta + b.delta);
  const D grad_ad = 0.5 * (a.grad_ad + b.grad_ad);
  const D dlnP = (log(b.P) - log(a.P)) / dm;

  std::array<D, NVAR> f{};
  f[0] = (b.lnr - a.lnr) / dm - 1.0 / (4.0 * M_PI * rb * rb * rb * rhob);
  f[1] = dlnP + G * mb / (4.0 * M_PI * rb * rb * rb * rb * Pb);
  D eps_grav{};
  if (prev) eps_grav = -(a.E - prev->E - a.P / (a.rho * a.rho) * (a.rho - prev->rho)) / dt;
  // Retains the existing left-endpoint backward energy difference. A future
  // time integrator must address its order separately from Jacobian assembly.
  f[2] = (b.L - a.L) / dm - (0.5 * (a.eps + b.eps) + eps_grav);

  const D gravity = G * mb / (rb * rb);
  const D grad_rad = 3.0 * kb * Lb * Pb
      / (16.0 * M_PI * constants::a_rad * constants::c * G * mb * Tb * Tb * Tb * Tb);
  EosState midpoint{};
  midpoint.P = Pb.value; midpoint.cp = cp.value; midpoint.delta = delta.value;
  const double U = mixing_length_U(Tb.value, rhob.value, kb.value, gravity.value, midpoint, phys.alpha_mlt);
  const auto convection = mixing_length_gradient(grad_rad.value, grad_ad.value, U);
  D grad(convection.grad);
  if constexpr (N > 0) {
    // Constant terms in ln U do not contribute. Carry the state dependence of
    // H_P=P/(rho*g), cp, delta, and opacity, rather than freezing efficiency.
    const D logHp = log(Pb) - log(rhob) - log(gravity);
    const D logU = 3.0 * log(Tb) - log(cp) - 2.0 * log(rhob) - log(kb)
                 - 1.5 * logHp - 0.5 * log(gravity) - 0.5 * log(delta);
    for (std::size_t v = 0; v < N; ++v)
      grad.d[v] = convection.dgrad_dgrad_rad * grad_rad.d[v]
                + convection.dgrad_dgrad_ad * grad_ad.d[v] + convection.dgrad_dlnU * logU.d[v];
  }
  // The temperature row is never multiplied by grad/grad_rad.
  f[3] = (b.lnT - a.lnT) / dm - grad * dlnP;
  for (const auto& row : f) {
    if (!std::isfinite(row.value)) throw std::domain_error("zone_residual: non-finite residual");
    for (double d : row.d)
      if (!std::isfinite(d)) throw std::domain_error("zone_residual: non-finite Jacobian");
  }
  return f;
}
} // namespace

std::array<double, NVAR> zone_equations(const Model& model, std::size_t i,
    const Physics& phys, double dt, const Model* prev) {
  const auto old = prepare(model, i, phys, dt, prev);
  const auto f = equations<0>(model, i, model.y[i], model.y[i + 1], phys, dt, old);
  std::array<double, NVAR> out{};
  for (std::size_t k = 0; k < NVAR; ++k) out[k] = f[k].value;
  return out;
}

ZoneResidual zone_residual(const Model& model, std::size_t i,
    const Physics& phys, double dt, const Model* prev) {
  const auto old = prepare(model, i, phys, dt, prev);
  const auto f = equations<2 * NVAR>(model, i, model.y[i], model.y[i + 1], phys, dt, old);
  ZoneResidual out{};
  for (std::size_t k = 0; k < NVAR; ++k) {
    out.f[k] = f[k].value;
    for (std::size_t v = 0; v < NVAR; ++v) {
      out.dfdy_lo[k][v] = f[k].d[v]; out.dfdy_hi[k][v] = f[k].d[NVAR + v];
    }
  }
  return out;
}

ZoneResidual zone_residual_numerical(const Model& model, std::size_t i,
    const Physics& phys, double dt, const Model* prev, double step) {
  if (!std::isfinite(step) || !(step > 0.0))
    throw std::invalid_argument("zone_residual_numerical: invalid difference step");
  const auto old = prepare(model, i, phys, dt, prev);
  const auto values = equations<0>(model, i, model.y[i], model.y[i + 1], phys, dt, old);
  ZoneResidual out{};
  for (std::size_t k = 0; k < NVAR; ++k) out.f[k] = values[k].value;
  const double dm = model.m[i + 1] - model.m[i];
  // Scale luminosity from this zone's flux and heating, with a 1 erg/s floor;
  // a fixed fraction of solar luminosity is inappropriate for a cold remnant.
  const double Lscale = std::max({std::abs(model.y[i].L), std::abs(model.y[i + 1].L),
      std::abs(model.y[i + 1].L - model.y[i].L - dm * out.f[2]), 1.0});
  for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
    auto& jac = endpoint == 0 ? out.dfdy_lo : out.dfdy_hi;
    for (std::size_t v = 0; v < NVAR; ++v) {
      const auto variable = static_cast<Var>(v);
      const double h = step * (variable == Var::L ? Lscale : 1.0);
      Point lo = model.y[i], hi = model.y[i + 1];
      Point& point = endpoint == 0 ? lo : hi;
      point[variable] += h;
      const auto plus = equations<0>(model, i, lo, hi, phys, dt, old);
      point[variable] -= 2.0 * h;
      const auto minus = equations<0>(model, i, lo, hi, phys, dt, old);
      for (std::size_t k = 0; k < NVAR; ++k) jac[k][v] = (plus[k].value - minus[k].value) / (2.0 * h);
    }
  }
  return out;
}

} // namespace ember
