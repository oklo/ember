#include "ember/eos.hpp"
#include "ember/eos_composite.hpp"
#include "ember/constants.hpp"
#include "fermi.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
using namespace constants;

EosState IdealEos::eval(double T, double rho, const Composition& comp) const {
  if (!(T > 0.0) || !(rho > 0.0)) throw std::domain_error("IdealEos: non-positive T or rho");

  const double nI = comp.mu_ions_inv();          // mol/g of ions
  const double nE = comp.mu_elec_inv();          // mol/g of electrons

  // --- ions: ideal, non-degenerate -----------------------------------------
  const double P_ion = nI * R_gas * rho * T;
  const double E_ion = 1.5 * nI * R_gas * T;

  // --- radiation ------------------------------------------------------------
  const double P_rad = a_rad * T * T * T * T / 3.0;
  const double E_rad = 3.0 * P_rad / rho;

  // --- electrons: relativistic Fermi-Dirac ---------------------------------
  const double mc2  = me * c * c;
  const double beta = kB * T / mc2;
  const double A    = 8.0 * M_PI * std::pow(me * c / h, 3.0);
  const double ne   = nE * NA * rho;

  // Solve A * In(eta, beta) = ne for eta.  First guess from the two limits:
  // non-degenerate Maxwell-Boltzmann, or the zero-temperature Fermi momentum.
  double eta;
  {
    const double lam = h / std::sqrt(2.0 * M_PI * me * kB * T);
    const double nd  = std::log(std::max(ne * lam * lam * lam / 2.0, 1e-300));
    const double xF  = std::cbrt(std::max(3.0 * ne / (8.0 * M_PI), 1e-300)) * (h / (me * c));
    const double dg  = (std::sqrt(1.0 + xF * xF) - 1.0) / beta;
    eta = (nd < 0.0) ? nd : dg;
    for (int it = 0; it < 200; ++it) {
      const auto f = fermi::evaluate(eta, beta);
      const double r  = A * f.In - ne;
      const double dr = A * f.dIn_deta;
      if (!(dr > 0.0)) break;
      double step = r / dr;
      const double cap = 4.0 * (1.0 + std::abs(eta));
      if (step >  cap) step =  cap;
      if (step < -cap) step = -cap;
      eta -= step;
      if (std::abs(step) < 1e-13 * (1.0 + std::abs(eta))) break;
    }
  }
  const auto f = fermi::evaluate(eta, beta);
  const double P_e = A * mc2 * f.Ip / 3.0;
  const double u_e = A * mc2 * f.Iu;            // energy density
  const double E_e = u_e / rho;

  // eta responds to T and rho through the constraint A In = ne:
  //   (deta/dlnT)_rho  = -(dIn/dlnbeta)/(dIn/deta)      [beta ~ T]
  //   (deta/dlnrho)_T  =  In/(dIn/deta)
  const double dEta_dlnT   = -f.dIn_dlnb / f.dIn_deta;
  const double dEta_dlnRho =  f.In       / f.dIn_deta;

  double dPe_dlnT   = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnT + f.dIp_dlnb);
  const double dPe_dlnRho = (A * mc2 / 3.0) * (f.dIp_deta * dEta_dlnRho);
  double due_dlnT   = A * mc2 * (f.dIu_deta * dEta_dlnT + f.dIu_dlnb);
  if(eta>100.) {
    const auto response=fermi::density_response(eta,beta,f,false);
    dPe_dlnT=(A*mc2/3.)*response.dIp_dlnT;
    due_dlnT=A*mc2*response.dIu_dlnT;
  }

  // --- assemble -------------------------------------------------------------
  EosState s{};
  s.P = P_ion + P_rad + P_e;
  s.E = E_ion + E_rad + E_e;

  const double dPdlnT   = P_ion + 4.0 * P_rad + dPe_dlnT;
  const double dPdlnRho = P_ion + dPe_dlnRho;
  s.chiT   = dPdlnT / s.P;
  s.chiRho = dPdlnRho / s.P;

  // E_e = u_e/rho, so (dE_e/dlnT)_rho = (du_e/dlnT)/rho.
  const double dEdlnT = 1.5 * nI * R_gas * T + 4.0 * E_rad + due_dlnT / rho;
  s.cv = dEdlnT / T;

  s.Gamma1  = s.chiRho + s.chiT * s.chiT * s.P / (rho * T * s.cv);
  s.grad_ad = s.chiT * s.P / (rho * T * s.cv * s.Gamma1);
  s.cp      = s.cv * s.Gamma1 / s.chiRho;
  s.delta   = s.chiT / s.chiRho;
  s.mu      = 1.0 / (nI + nE);
  s.free_e  = nE / nI;
  s.S = IonGas{}.eval(T,rho,comp).S+4*P_rad/(rho*T)
        +A*kB*fermi::entropy(eta,beta,f)/rho;
  return s;
}

EosResponse IdealEos::eval_with_derivatives(double T, double rho, const Composition& comp) const {
  // The monolithic value implementation remains an independent assembly
  // check. Its second derivatives are those of the same three components.
  static const CompositeEos components;
  return components.eval_with_derivatives(T, rho, comp);
}

double Eos::rho_from_PT(double T, double P, const Composition& comp,
                        double rho_guess) const {
  if (!std::isfinite(T) || !(T > 0.0) || !std::isfinite(P) || !(P > 0.0)
      || !std::isfinite(rho_guess) || rho_guess < 0.0)
    throw std::domain_error("Eos::rho_from_PT: invalid T, P, or density guess");
  double rho = rho_guess > 0.0 ? rho_guess
             : P / (comp.mu_ions_inv() + comp.mu_elec_inv()) / (R_gas * T);
  double lo=-std::numeric_limits<double>::infinity(),hi=std::numeric_limits<double>::infinity();
  if(const auto bounds=density_range(T,comp)) {
    if(!(bounds->min>0) || !std::isfinite(bounds->max) || bounds->min>=bounds->max)
      throw std::domain_error("Eos::rho_from_PT: invalid density bounds");
    lo=std::log(bounds->min);hi=std::log(bounds->max);
    if(eval(T,bounds->min,comp).P>P || eval(T,bounds->max,comp).P<P)
      throw std::domain_error("Eos::rho_from_PT: pressure outside supported density interval");
    rho=std::clamp(rho,bounds->min,bounds->max);
  }
  for (int it = 0; it < 200; ++it) {
    if (!std::isfinite(rho) || !(rho > 0.0))
      throw std::domain_error("Eos::rho_from_PT: density outside representable range");
    const EosState s = eval(T, rho, comp);
    if (!std::isfinite(s.P) || !(s.P > 0.0) || !std::isfinite(s.chiRho) || !(s.chiRho > 0.0))
      throw std::domain_error("Eos::rho_from_PT: invalid EOS pressure or density derivative");
    if (s.chiRho < 32.0 * std::numeric_limits<double>::epsilon())
      throw std::domain_error("Eos::rho_from_PT: density is unresolved by total pressure");
    const double f = std::abs(s.P - P) < 0.5 * P
        ? std::log1p((s.P - P) / P) : std::log(s.P) - std::log(P);
    double step = f / s.chiRho;
    // A small relative pressure residual alone is insufficient when radiation
    // dominates: require convergence in the inferred density as well.
    if (std::abs(f) < 1e-12 && std::abs(step) < 1e-10) return rho;
    if (step >  0.7) step =  0.7;
    if (step < -0.7) step = -0.7;
    const double r=std::log(rho);
    if(f<0)lo=r;else hi=r;
    const double next=r-step;
    rho=std::exp(next>lo && next<hi?next:.5*(lo+hi));
  }
  throw std::runtime_error("Eos::rho_from_PT: density inversion did not converge");
}

} // namespace ember
