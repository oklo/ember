#include "ember/structure.hpp"
#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include <cmath>

namespace ember {
using namespace constants;

namespace {

// Everything the zone equations need at a point, gathered once.
struct Local {
  double lnr, lnT, r, rho, T, L, P, lnP;
  double chiT, chiRho, grad_ad, cp, delta, kappa, dlnk_dlnT, dlnk_dlnRho;
  double eps, dlneps_dlnT, dlneps_dlnRho;
  double E;
};

Local gather(const Model& m, std::size_t i, const Physics& p) {
  Local q{};
  q.lnr = m.y[i].lnr;
  q.r   = std::exp(q.lnr);
  q.rho = std::exp(m.y[i].lnrho);
  q.lnT = m.y[i].lnT;
  q.T   = std::exp(q.lnT);
  q.L   = m.y[i].L;
  const auto e = p.eos->eval(q.T, q.rho, m.comp[i]);
  q.P = e.P; q.lnP = std::log(e.P);
  q.chiT = e.chiT; q.chiRho = e.chiRho; q.grad_ad = e.grad_ad; q.cp = e.cp;
  q.delta = e.delta;
  q.E = e.E;
  const auto k = p.opacity->eval(q.T, q.rho, m.comp[i]);
  q.kappa = k.kappa; q.dlnk_dlnT = k.dlnk_dlnT; q.dlnk_dlnRho = k.dlnk_dlnRho;
  const auto n = p.nuclear->eval(q.T, q.rho, m.comp[i]);
  q.eps = n.eps; q.dlneps_dlnT = n.dlneps_dlnT; q.dlneps_dlnRho = n.dlneps_dlnRho;
  return q;
}

} // namespace

ZoneResidual zone_residual(const Model& mdl, std::size_t i,
                           const Physics& phys, double dt,
                           const Model* prev) {
  ZoneResidual R{};
  const std::size_t j = i + 1;
  const double dm = mdl.m[j] - mdl.m[i];
  if (!(dm > 0.0)) return R;

  // Numerical Jacobian for now, taken by perturbing the eight variables that
  // bound the zone.  It is written this way deliberately: the analytic form
  // will replace it once the equations themselves are settled, and having both
  // lets the analytic version be checked against something. The MLT gradient
  // supplies its partials; full assembly will also need state derivatives of
  // the EOS's cp, delta, and grad_ad, beyond their current returned values.
  auto residual = [&](const Point& yl, const Point& yh) {
    Model tmp = mdl;               // cheap enough at this stage; the solver
    tmp.y[i] = yl; tmp.y[j] = yh;  // will not do this per iteration
    const Local a = gather(tmp, i, phys);
    const Local b = gather(tmp, j, phys);

    const double rb   = 0.5 * (a.r + b.r);
    const double rhob = 0.5 * (a.rho + b.rho);
    const double Tb   = 0.5 * (a.T + b.T);
    const double Pb   = 0.5 * (a.P + b.P);
    const double mb   = 0.5 * (mdl.m[i] + mdl.m[j]);

    std::array<double, NVAR> f{};
    // (1) mass conservation
    f[0] = (b.lnr - a.lnr) / dm - 1.0 / (4.0 * M_PI * rb * rb * rb * rhob);
    // (2) hydrostatic equilibrium
    f[1] = (b.lnP - a.lnP) / dm + G * mb / (4.0 * M_PI * std::pow(rb, 4) * Pb);
    // (3) energy: nuclear minus the heat taken to warm the material
    double eps_grav = 0.0;
    if (dt > 0.0 && prev) {
      // -T ds/dt, written as -(dE/dt - (P/rho^2) drho/dt) so it needs no
      // entropy from the equation of state.
      const auto pa = phys.eos->eval(std::exp(prev->y[i].lnT),
                                     std::exp(prev->y[i].lnrho), prev->comp[i]);
      const double dE   = a.E - pa.E;
      const double drho = a.rho - std::exp(prev->y[i].lnrho);
      eps_grav = -(dE - (a.P / (a.rho * a.rho)) * drho) / dt;
    }
    f[2] = (b.L - a.L) / dm - (0.5 * (a.eps + b.eps) + eps_grav);
    // (4) Schwarzschild criterion and optically thick mixing-length transport.
    const double kb = 0.5 * (a.kappa + b.kappa);
    const double Lb = 0.5 * (a.L + b.L);
    const double grad_rad = 3.0 * kb * Lb * Pb
                          / (16.0 * M_PI * a_rad * c * G * mb * std::pow(Tb, 4));
    const double grad_ad  = 0.5 * (a.grad_ad + b.grad_ad);
    EosState eb{};
    eb.P = Pb;
    eb.cp = 0.5 * (a.cp + b.cp);
    eb.delta = 0.5 * (a.delta + b.delta);
    const double gravity = G * mb / (rb * rb);
    const double U = mixing_length_U(Tb, rhob, kb, gravity, eb, phys.alpha_mlt);
    const double grad = mixing_length_gradient(grad_rad, grad_ad, U).grad;
    // Plain gradient form at every efficiency. Never multiply this row by
    // grad/grad_rad: in a giant that factor can be 1e-6 or smaller.
    f[3] = (b.lnT - a.lnT) / dm - grad * ((b.lnP - a.lnP) / dm);
    return f;
  };

  R.f = residual(mdl.y[i], mdl.y[j]);
  constexpr double h = 1e-6;
  for (std::size_t v = 0; v < NVAR; ++v) {
    const Var vv = static_cast<Var>(v);
    {
      Point yl = mdl.y[i];
      const double s = (vv == Var::L) ? h * (std::abs(yl.L) + 1e-3 * constants::Lsun)
                                      : h;
      yl[vv] += s;
      const auto fp = residual(yl, mdl.y[j]);
      yl[vv] -= 2 * s;
      const auto fm = residual(yl, mdl.y[j]);
      for (std::size_t k = 0; k < NVAR; ++k) R.dfdy_lo[k][v] = (fp[k] - fm[k]) / (2 * s);
    }
    {
      Point yh = mdl.y[j];
      const double s = (vv == Var::L) ? h * (std::abs(yh.L) + 1e-3 * constants::Lsun)
                                      : h;
      yh[vv] += s;
      const auto fp = residual(mdl.y[i], yh);
      yh[vv] -= 2 * s;
      const auto fm = residual(mdl.y[i], yh);
      for (std::size_t k = 0; k < NVAR; ++k) R.dfdy_hi[k][v] = (fp[k] - fm[k]) / (2 * s);
    }
  }
  return R;
}

} // namespace ember
