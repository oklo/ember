#include "ember/eos_component.hpp"
#include "ember/eos_composite.hpp"
#include "ember/conduction.hpp"
#include "ember/constants.hpp"
#include "fermi.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
using namespace constants;

EosTerm IonGas::eval(double T, double rho, const Composition& c) const {
  const double nI = c.mu_ions_inv();
  EosTerm t{};
  t.P = nI * R_gas * rho * T;
  t.E = 1.5 * nI * R_gas * T;
  t.dP_dlnT = t.P;          // P ~ T
  t.dP_dlnRho = t.P;        // P ~ rho
  t.dE_dlnT = t.E;          // E ~ T
  t.dE_dlnRho = 0.0;
  return t;
}

EosTerm Radiation::eval(double T, double rho, const Composition&) const {
  EosTerm t{};
  t.P = a_rad * T * T * T * T / 3.0;
  t.E = 3.0 * t.P / rho;
  t.dP_dlnT = 4.0 * t.P;
  t.dP_dlnRho = 0.0;
  t.dE_dlnT = 4.0 * t.E;
  t.dE_dlnRho = -t.E;       // E = 3P/rho at fixed T
  return t;
}

EosTerm ElectronGas::eval(double T, double rho, const Composition& c) const {
  const double mc2  = me * c_light * c_light;
  const double beta = kB * T / mc2;
  const double A    = 8.0 * M_PI * std::pow(me * c_light / h, 3.0);
  const double ne   = c.mu_elec_inv() * NA * rho;

  double eta;
  {
    const double lam = h / std::sqrt(2.0 * M_PI * me * kB * T);
    const double nd  = std::log(std::max(ne * lam * lam * lam / 2.0, 1e-300));
    const double xF  = std::cbrt(std::max(3.0 * ne / (8.0 * M_PI), 1e-300)) * (h / (me * c_light));
    const double dg  = (std::sqrt(1.0 + xF * xF) - 1.0) / beta;
    eta = (nd < 0.0) ? nd : dg;
    for (int it = 0; it < 200; ++it) {
      const auto f = fermi::evaluate(eta, beta);
      const double r = A * f.In - ne, dr = A * f.dIn_deta;
      if (!(dr > 0.0)) break;
      double step = r / dr;
      const double cap = 4.0 * (1.0 + std::abs(eta));
      step = std::clamp(step, -cap, cap);
      eta -= step;
      if (std::abs(step) < 1e-13 * (1.0 + std::abs(eta))) break;
    }
  }
  const auto f = fermi::evaluate(eta, beta);
  const double dEta_dlnT   = -f.dIn_dlnb / f.dIn_deta;
  const double dEta_dlnRho =  f.In       / f.dIn_deta;

  EosTerm t{};
  t.P = A * mc2 * f.Ip / 3.0;
  const double u = A * mc2 * f.Iu;
  t.E = u / rho;
  t.dP_dlnT   = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnT + f.dIp_dlnb);
  t.dP_dlnRho = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnRho);
  t.dE_dlnT   = A * mc2 * (f.dIu_deta * dEta_dlnT + f.dIu_dlnb) / rho;
  // E = u(eta,beta)/rho, so at fixed T the explicit 1/rho contributes -E.
  t.dE_dlnRho = A * mc2 * (f.dIu_deta * dEta_dlnRho) / rho - t.E;
  return t;
}

CompositeEos::CompositeEos() {
  parts_.push_back(std::make_unique<IonGas>());
  parts_.push_back(std::make_unique<Radiation>());
  parts_.push_back(std::make_unique<ElectronGas>());
  name_ = "composite(ions + radiation + electrons)";
}

EosState CompositeEos::eval(double T, double rho, const Composition& c) const {
  if (!(T > 0.0) || !(rho > 0.0))
    throw std::domain_error("CompositeEos: non-positive T or rho");
  EosTerm sum{};
  for (const auto& p : parts_) sum += p->eval(T, rho, c);

  EosState s{};
  s.P = sum.P;
  s.E = sum.E;
  s.S = sum.S;
  s.chiT   = sum.dP_dlnT / sum.P;
  s.chiRho = sum.dP_dlnRho / sum.P;
  s.cv     = sum.dE_dlnT / T;
  s.Gamma1 = s.chiRho + s.chiT * s.chiT * s.P / (rho * T * s.cv);
  s.grad_ad = s.chiT * s.P / (rho * T * s.cv * s.Gamma1);
  s.cp     = s.cv * s.Gamma1 / s.chiRho;
  s.delta  = s.chiT / s.chiRho;
  s.mu     = 1.0 / (c.mu_ions_inv() + c.mu_elec_inv());
  s.free_e = c.mu_elec_inv() / c.mu_ions_inv();
  return s;
}

OpacityState CombinedOpacity::eval(double T, double rho, const Composition& c) const {
  const OpacityState r = rad_->eval(T, rho, c);
  if (!cond_) return r;
  const OpacityState k = cond_->eval(T, rho, c);
  // 1/kappa = 1/kappa_rad + 1/kappa_cond; differentiate the reciprocal sum so
  // the combined derivatives stay exact rather than being re-differenced.
  const double ir = 1.0 / r.kappa, ic = 1.0 / k.kappa, it = ir + ic;
  OpacityState s{};
  s.kappa = 1.0 / it;
  const double wr = ir / it, wc = ic / it;   // weights sum to one
  s.dlnk_dlnT   = wr * r.dlnk_dlnT   + wc * k.dlnk_dlnT;
  s.dlnk_dlnRho = wr * r.dlnk_dlnRho + wc * k.dlnk_dlnRho;
  return s;
}

} // namespace ember
