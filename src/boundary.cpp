#include "ember/boundary.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {

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
