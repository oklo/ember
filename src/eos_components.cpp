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
  t.d2P_dlnT2 = t.d2P_dlnTdlnRho = t.d2P_dlnRho2 = t.P;
  t.d2E_dlnT2 = t.E;
  // Sackur--Tetrode translational entropy, with unit nuclear internal
  // partition functions. Distinct species have their own number densities,
  // hence ideal mixing is included. Kinematic masses follow the declared
  // representative-isotope convention; matching another EOS's entropy
  // zero/internal partition functions is a separate caloric-reference task.
  const auto add_species=[&](double number_per_baryon_mass,double mass) {
    if(number_per_baryon_mass==0.)return;
    const double log_lambda=std::log(h)-.5*std::log(2*M_PI*mass*amu*kB*T);
    t.S+=R_gas*number_per_baryon_mass*
      (2.5-std::log(number_per_baryon_mass*NA*rho)-3*log_lambda);
  };
  for(std::size_t i=0;i<(c.metal_inventory==MetalInventory::gs98?3:NSPEC);++i)
    add_species(c.X[i]/c.abundance_weight(i),c.abundance_weight(i));
  if(c.metal_inventory==MetalInventory::gs98) {
    const double scale=c.basis==AbundanceBasis::baryon_mass?1.:gs98_atomic_mass_scale();
    for(const auto& element:gs98_metals)
      add_species(c.Z()*element.fraction/(element.mass_number*scale),
        c.basis==AbundanceBasis::baryon_mass?element.mass_number:element.atomic_weight);
  }
  return t;
}

EosTerm Radiation::eval(double T, double rho, const Composition&) const {
  EosTerm t{};
  t.P = a_rad * T * T * T * T / 3.0;
  t.E = 3.0 * t.P / rho;
  t.S = 4.0 * t.P / (rho*T);
  t.dP_dlnT = 4.0 * t.P;
  t.dP_dlnRho = 0.0;
  t.dE_dlnT = 4.0 * t.E;
  t.dE_dlnRho = -t.E;       // E = 3P/rho at fixed T
  t.d2P_dlnT2 = 16.0 * t.P;
  t.d2E_dlnT2 = 16.0 * t.E;
  t.d2E_dlnTdlnRho = -4.0 * t.E;
  return t;
}

static EosTerm electron_term(double T, double rho, const Composition& comp, bool second) {
  const double mc2  = me * c_light * c_light;
  const double beta = kB * T / mc2;
  const double A    = 8.0 * M_PI * std::pow(me * c_light / h, 3.0);
  const double ne   = comp.mu_elec_inv() * NA * rho;

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
  t.S = A*kB*fermi::entropy(eta,beta,f)/rho;
  t.dP_dlnT   = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnT + f.dIp_dlnb);
  t.dP_dlnRho = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnRho);
  t.dE_dlnT   = A * mc2 * (f.dIu_deta * dEta_dlnT + f.dIu_dlnb) / rho;
  // E = u(eta,beta)/rho, so at fixed T the explicit 1/rho contributes -E.
  t.dE_dlnRho = A * mc2 * (f.dIu_deta * dEta_dlnRho) / rho - t.E;
  if (second || eta>100.) {
    const auto d = fermi::density_response(eta, beta, f, second);
    // Equivalent quadrature at strong degeneracy, avoiding loss of the
    // small thermal response. Warm-state arithmetic is retained.
    if(eta>100.) {
      t.dP_dlnT = (A * mc2 / 3.0) * d.dIp_dlnT;
      t.dE_dlnT = A * mc2 * d.dIu_dlnT / rho;
    }
    if(second) {
      t.d2P_dlnT2 = (A * mc2 / 3.0) * d.d2Ip_dlnT2;
      t.d2P_dlnTdlnRho = (A * mc2 / 3.0) * d.d2Ip_dlnTdlnRho;
      t.d2P_dlnRho2 = (A * mc2 / 3.0) * d.d2Ip_dlnRho2;
      t.d2E_dlnT2 = (A * mc2 / rho) * d.d2Iu_dlnT2;
      // Specific energy contains an explicit rho^-1, unlike the integrals.
      t.d2E_dlnTdlnRho = (A * mc2 / rho) * d.d2Iu_dlnTdlnRho - t.dE_dlnT;
    }
  }
  return t;
}

EosTerm ElectronGas::eval(double T, double rho, const Composition& comp) const {
  return electron_term(T, rho, comp, false);
}
EosTerm ElectronGas::eval_with_derivatives(double T, double rho, const Composition& comp) const {
  return electron_term(T, rho, comp, true);
}

CompositeEos::CompositeEos() {
  parts_.push_back(std::make_unique<IonGas>());
  parts_.push_back(std::make_unique<Radiation>());
  parts_.push_back(std::make_unique<ElectronGas>());
  name_ = "composite(ions + radiation + electrons)";
}

static EosState assemble(const EosTerm& sum, double T, double rho, const Composition& comp) {
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
  s.mu     = 1.0 / (comp.mu_ions_inv() + comp.mu_elec_inv());
  s.free_e = comp.mu_elec_inv() / comp.mu_ions_inv();
  return s;
}

EosState CompositeEos::eval(double T, double rho, const Composition& comp) const {
  if (!std::isfinite(T) || !std::isfinite(rho) || !(T > 0.0) || !(rho > 0.0))
    throw std::domain_error("CompositeEos: invalid T or rho");
  EosTerm sum{};
  for (const auto& p : parts_) sum += p->eval(T, rho, comp);
  return assemble(sum, T, rho, comp);
}

EosResponse CompositeEos::eval_with_derivatives(double T, double rho, const Composition& comp) const {
  if (!std::isfinite(T) || !std::isfinite(rho) || !(T > 0.0) || !(rho > 0.0))
    throw std::domain_error("CompositeEos: invalid T or rho");
  EosTerm sum{};
  for (const auto& p : parts_) sum += p->eval_with_derivatives(T, rho, comp);
  EosResponse out{};
  out.state = assemble(sum, T, rho, comp);
  out.dE_dlnRho = sum.dE_dlnRho;
  const auto& s = out.state;
  // delta=P_T/P_rho, cp=E_T/T + P_T*delta/(rho*T),
  // grad_ad=P*delta/(rho*T*cp); subscripts here denote logarithmic partials.
  for (int axis = 0; axis < 2; ++axis) {
    const double Px = axis == 0 ? sum.dP_dlnT : sum.dP_dlnRho;
    const double PTx = axis == 0 ? sum.d2P_dlnT2 : sum.d2P_dlnTdlnRho;
    const double PRx = axis == 0 ? sum.d2P_dlnTdlnRho : sum.d2P_dlnRho2;
    const double ETx = axis == 0 ? sum.d2E_dlnT2 : sum.d2E_dlnTdlnRho;
    const double dx = (PTx - s.delta * PRx) / sum.dP_dlnRho;
    const double cvx = (ETx - (axis == 0 ? sum.dE_dlnT : 0.0)) / T;
    const double cpx = cvx + (PTx * s.delta + sum.dP_dlnT * dx - sum.dP_dlnT * s.delta) / (rho * T);
    const double gx = s.grad_ad * (Px / s.P - 1.0 - cpx / s.cp)
                    + s.P * dx / (rho * T * s.cp);
    if (axis == 0) {
      out.dcp_dlnT = cpx; out.ddelta_dlnT = dx; out.dgrad_ad_dlnT = gx;
    } else {
      out.dcp_dlnRho = cpx; out.ddelta_dlnRho = dx; out.dgrad_ad_dlnRho = gx;
    }
  }
  return out;
}

std::optional<Opacity::DensityRange> CombinedOpacity::density_range(double T, const Composition& comp) const {
  const auto r = rad_->density_range(T, comp);
  if (!cond_) return r;
  const auto k = cond_->density_range(T, comp);
  if (!r) return k;
  if (!k) return r;
  const DensityRange range{std::max(r->min, k->min), std::min(r->max, k->max)};
  if (range.min >= range.max)
    throw std::domain_error("CombinedOpacity: no common source density support");
  return range;
}

OpacityState CombinedOpacity::eval(double T, double rho, const Composition& comp) const {
  const OpacityState r = rad_->eval(T, rho, comp);
  if (!cond_) return r;
  const OpacityState k = cond_->eval(T, rho, comp);
  // 1/kappa = 1/kappa_rad + 1/kappa_cond; differentiate the reciprocal sum so
  // the combined derivatives stay exact rather than being re-differenced.
  const double ir = 1.0 / r.kappa, ic = 1.0 / k.kappa, it = ir + ic;
  OpacityState s{};
  s.kappa = 1.0 / it;
  const double wr = ir / it, wc = ic / it;   // weights sum to one
  if(wc==0) return r;
  s.dlnk_dlnT   = wr * r.dlnk_dlnT   + wc * k.dlnk_dlnT;
  s.dlnk_dlnRho = wr * r.dlnk_dlnRho + wc * k.dlnk_dlnRho;
  s.dlnk_dX = wr*r.dlnk_dX+wc*k.dlnk_dX;
  s.dlnk_dZ = wr*r.dlnk_dZ+wc*k.dlnk_dZ;
  s.dlnk_dY3 = wr*r.dlnk_dY3+wc*k.dlnk_dY3;
  return s;
}

} // namespace ember
