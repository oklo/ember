#include "ember/boundary.hpp"
#include "ember/constants.hpp"
#include "differential.hpp"
#include "energy.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {

CentralResidual central_residual(const Model& model, const Physics& phys, double dt, const Model* prev) {
  if (!phys.eos || !phys.nuclear || model.size() == 0 || model.m.size() != model.size()
      || model.comp.size() != model.size())
    throw std::invalid_argument("central_residual: missing physics or invalid model arrays");
  const double mass = model.m.front();
  const auto& point = model.y.front();
  if (!std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(dt))
    throw std::domain_error("central_residual: inner mass must be positive and dt finite");
  for (double v : {point.lnr, point.lnrho, point.lnT, point.L})
    if (!std::isfinite(v)) throw std::domain_error("central_residual: non-finite state");
  using D = detail::Differential<NVAR>;
  const D rho = detail::exp(D::variable(point.lnrho, 1));
  const double T = std::exp(point.lnT);
  if (!(rho.value > 0.0) || !std::isfinite(rho.value) || !(T > 0.0) || !std::isfinite(T))
    throw std::domain_error("central_residual: invalid density or temperature");
  const auto n = phys.nuclear->eval(T, rho.value, model.comp.front());
  D heating(n.eps);
  heating.d[1] = n.eps * n.dlneps_dlnRho;
  heating.d[2] = n.eps * n.dlneps_dlnT;
  if (dt > 0.0) {
    if (!phys.eos->has_internal_energy())
      throw std::logic_error("central_residual: EOS has no validated internal energy for time dependence");
    if (!prev || prev->size() != model.size() || prev->comp.size() != model.size() || prev->m != model.m)
      throw std::invalid_argument("central_residual: time dependence requires previous model on the same mesh");
    const auto e = phys.eos->eval_with_derivatives(T, rho.value, model.comp.front());
    const double old_rho = prev->rho(0);
    const double old_E = phys.eos->eval(prev->T(0), old_rho, prev->comp.front()).E;
    D P(e.state.P), E(e.state.E);
    P.d[1] = P.value * e.state.chiRho; P.d[2] = P.value * e.state.chiT;
    E.d[1] = e.dE_dlnRho; E.d[2] = e.state.cv * T;
    heating = heating + detail::gravitational_heating(E, P, rho, old_E, old_rho, dt);
  }
  const std::array<D, 2> f{
      D::variable(point.lnr, 0) - (std::log(3.0 / (4.0 * M_PI)) + std::log(mass)
                                  - D::variable(point.lnrho, 1)) / 3.0,
      D::variable(point.L, 3) - mass * heating};
  CentralResidual out{};
  for (std::size_t k = 0; k < 2; ++k) {
    out.f[k] = f[k].value; out.dfdy[k] = f[k].d;
    if (!std::isfinite(out.f[k])) throw std::domain_error("central_residual: non-finite residual");
    for (double v : out.dfdy[k])
      if (!std::isfinite(v)) throw std::domain_error("central_residual: non-finite derivative");
  }
  return out;
}

SurfaceResidual surface_residual(const Point& surface, double mass, const Composition& comp,
                                 const Eos& eos, const Atmosphere& atmosphere) {
  if (!std::isfinite(surface.lnr) || !std::isfinite(surface.lnrho) || !std::isfinite(surface.lnT)
      || !std::isfinite(surface.L) || !(surface.L > 0.0) || !std::isfinite(mass) || !(mass > 0.0))
    throw std::domain_error("surface_residual: invalid surface state or mass");
  SurfaceResidual r{};
  r.Teff = std::exp(0.25 * (std::log(surface.L) - std::log(4.0 * M_PI * constants::sigma_SB)
                          - 2.0 * surface.lnr));
  r.gravity = std::exp(std::log(constants::G) + std::log(mass) - 2.0 * surface.lnr);
  r.atmosphere = atmosphere.eval(r.Teff, r.gravity, comp);
  const auto e = eos.eval(std::exp(surface.lnT), std::exp(surface.lnrho), comp);
  const auto& a = r.atmosphere;
  if (!(e.P > 0.0) || !std::isfinite(e.P) || !std::isfinite(e.chiT) || !(e.chiRho > 0.0)
      || !std::isfinite(e.chiRho) || !(a.T > 0.0) || !(a.P > 0.0)
      || !std::isfinite(a.T) || !std::isfinite(a.P) || !std::isfinite(a.dlnT_dlnTeff)
      || !std::isfinite(a.dlnT_dlng) || !std::isfinite(a.dlnP_dlnTeff) || !std::isfinite(a.dlnP_dlng))
    throw std::domain_error("surface_residual: invalid EOS or atmosphere result");
  r.f = {surface.lnT - std::log(a.T), std::log(e.P) - std::log(a.P)};
  // ln Teff = (ln L - 2 ln r - ln(4*pi*sigma))/4; ln g = ln(G*M)-2 ln r.
  r.dfdy[0] = {0.5 * a.dlnT_dlnTeff + 2.0 * a.dlnT_dlng, 0.0, 1.0,
               -a.dlnT_dlnTeff / (4.0 * surface.L)};
  r.dfdy[1] = {0.5 * a.dlnP_dlnTeff + 2.0 * a.dlnP_dlng, e.chiRho, e.chiT,
               -a.dlnP_dlnTeff / (4.0 * surface.L)};
  return r;
}

} // namespace ember
