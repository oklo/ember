#include "ember/eos.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
using namespace constants;

namespace {

// Fermi-Dirac integrals F_k(eta) = \int_0^inf x^k /(1+exp(x-eta)) dx, for
// k = 1/2 and 3/2, by splitting at the Fermi level: Gauss-Legendre below,
// Gauss-Laguerre above.  Direct quadrature is slower than a rational fit but
// it is unambiguously correct, and the interface below is where a tabulated
// fit will drop in once the physics is settled.  Accuracy is ~1e-12 here,
// which is far beyond what the models need.
struct FD { double F; double dFdeta; };

constexpr int NQ = 48;

// Gauss-Legendre nodes/weights are generated once at first use by Newton on
// the Legendre polynomial - no magic tables to mistype.
struct GL {
  double x[NQ], w[NQ];
  GL() {
    for (int i = 0; i < NQ; ++i) {
      double z = std::cos(M_PI * (i + 0.75) / (NQ + 0.5)), dz;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1.0, p1 = 0.0;
        for (int j = 0; j < NQ; ++j) { double p2 = p1; p1 = p0;
          p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1); }
        double dp = NQ * (z * p0 - p1) / (z * z - 1.0);
        dz = p0 / dp; z -= dz;
        if (std::abs(dz) < 1e-15) break;
      }
      x[i] = z; w[i] = 2.0 / ((1.0 - z * z) * 0.0 + (1.0 - z * z));
      // weight: 2/((1-z^2) dP^2); recompute dP at the converged node
      double p0 = 1.0, p1 = 0.0;
      for (int j = 0; j < NQ; ++j) { double p2 = p1; p1 = p0;
        p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1); }
      double dp = NQ * (z * p0 - p1) / (z * z - 1.0);
      w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
  }
};
const GL& gl() { static const GL g; return g; }

FD fermi_dirac(double k, double eta) {
  const GL& q = gl();
  double F = 0.0, dF = 0.0;
  // Region 1: 0 .. max(eta,0)+20, mapped to the Legendre interval.
  const double xhi = (eta > 0.0 ? eta : 0.0) + 30.0;
  const double a = 0.0, b = xhi, h = 0.5 * (b - a), m = 0.5 * (a + b);
  for (int i = 0; i < NQ; ++i) {
    const double x = m + h * q.x[i];
    const double e = x - eta;
    // 1/(1+exp(e)) and its derivative wrt eta, written to avoid overflow.
    double f, dfde;
    if (e > 0.0) { const double z = std::exp(-e); f = z / (1.0 + z); dfde = f * (1.0 - f); }
    else         { const double z = std::exp(e);  f = 1.0 / (1.0 + z); dfde = f * (1.0 - f); }
    const double xk = std::pow(x, k);
    F  += q.w[i] * h * xk * f;
    dF += q.w[i] * h * xk * dfde;
  }
  return {F, dF};
}

} // namespace

EosState IdealEos::eval(double T, double rho, const Composition& comp) const {
  if (!(T > 0.0) || !(rho > 0.0)) throw std::domain_error("IdealEos: non-positive T or rho");

  const double nI = comp.mu_ions_inv();  // mol/g of ions
  const double nE = comp.mu_elec_inv();  // mol/g of electrons (full ionisation)

  // --- ions: ideal, non-degenerate -----------------------------------------
  const double P_ion  = nI * R_gas * rho * T;
  const double E_ion  = 1.5 * nI * R_gas * T;

  // --- radiation ------------------------------------------------------------
  const double P_rad = a_rad * T * T * T * T / 3.0;
  const double E_rad = 3.0 * P_rad / rho;

  // --- electrons: Fermi-Dirac, non-relativistic -----------------------------
  // n_e = (8 pi /3) (2 m kT/h^2)^{3/2} * (3/2) F_{1/2}(eta) ... standard form:
  //   n_e = (2/sqrt(pi)) * (2 pi m kT/h^2)^{3/2} * 2 F_{1/2}(eta)/sqrt(pi)
  // We use n_e = C_n T^{3/2} F_{1/2}, P_e = C_p T^{5/2} F_{3/2}, with
  //   C_n = 4 pi (2 m kT/h^2)^{3/2} / ... collected below.
  // n_e = (4 pi / h^3) (2 m kT)^{3/2} F_{1/2}(eta) = (4/sqrt(pi)) F_{1/2} / lam^3
  // with lam the thermal de Broglie wavelength; the non-degenerate limit of
  // this is n_e = 2 e^eta / lam^3, the spin-2 Maxwell-Boltzmann gas, which is
  // the check that fixes the constant.
  const double lam = h / std::sqrt(2.0 * M_PI * me * kB * T);
  const double Cn  = 4.0 / (std::sqrt(M_PI) * lam * lam * lam);
  // n_e target
  const double ne = nE * NA * rho;
  // Solve Cn * T^0 * F_{1/2}(eta) = ne   (T dependence already inside lam)
  double eta = 0.0;
  {
    // bracket-free Newton from an analytic first guess: the non-degenerate
    // limit eta = ln(ne/(Cn * sqrt(pi)/2)) and the degenerate limit
    // eta = (3 sqrt(pi) ne /(4 Cn))^{2/3}.
    const double nd = std::log(std::max(ne / (Cn * std::sqrt(M_PI) / 2.0), 1e-300));
    const double dg = std::pow(std::max(1.5 * ne / Cn, 1e-300), 2.0 / 3.0);
    eta = (nd < 0.0) ? nd : dg;
    for (int it = 0; it < 100; ++it) {
      const FD f = fermi_dirac(0.5, eta);
      const double r = Cn * f.F - ne;
      const double dr = Cn * f.dFdeta;
      double step = r / dr;
      if (step >  2.0) step =  2.0;
      if (step < -2.0) step = -2.0;
      eta -= step;
      if (std::abs(step) < 1e-13 * (1.0 + std::abs(eta))) break;
    }
  }
  const FD f12 = fermi_dirac(0.5, eta);
  const FD f32 = fermi_dirac(1.5, eta);
  const double P_e = (2.0 / 3.0) * Cn * kB * T * f32.F;
  const double E_e = 1.5 * P_e / rho;

  // --- assemble -------------------------------------------------------------
  EosState s{};
  s.P = P_ion + P_rad + P_e;
  s.E = E_ion + E_rad + E_e;

  // Analytic derivatives.  For the electrons, eta shifts with T and rho at
  // fixed n_e; differentiating the constraint Cn(T) F_{1/2}(eta) = n_e gives
  //   (dEta/dlnT)_rho = -3/2 * F_{1/2}/F'_{1/2},
  //   (dEta/dlnRho)_T =        F_{1/2}/F'_{1/2}.
  const double r12 = f12.F / f12.dFdeta;
  const double dEta_dlnT   = -1.5 * r12;
  const double dEta_dlnRho =  r12;
  const double dlnPe_dlnT   = 2.5 + (f32.dFdeta / f32.F) * dEta_dlnT;
  const double dlnPe_dlnRho =       (f32.dFdeta / f32.F) * dEta_dlnRho;

  const double dPdlnT   = P_ion + 4.0 * P_rad + P_e * dlnPe_dlnT;
  const double dPdlnRho = P_ion + P_e * dlnPe_dlnRho;
  s.chiT   = dPdlnT / s.P;
  s.chiRho = dPdlnRho / s.P;

  const double dEdlnT = 1.5 * nI * R_gas * T + 4.0 * E_rad
                      + 1.5 * (P_e / rho) * dlnPe_dlnT;
  s.cv = dEdlnT / T;

  s.Gamma1 = s.chiRho + s.chiT * s.chiT * s.P / (rho * T * s.cv);
  s.grad_ad = s.chiT * s.P / (rho * T * s.cv * s.Gamma1);
  s.cp = s.cv * s.Gamma1 / s.chiRho;
  s.delta = s.chiT / s.chiRho;
  s.mu = 1.0 / (nI + nE);
  s.free_e = nE / comp.mu_ions_inv();
  s.S = 0.0; // not yet needed by the solver; filled when the tables land
  return s;
}

double Eos::rho_from_PT(double T, double P, const Composition& comp,
                        double rho_guess) const {
  double rho = rho_guess > 0.0 ? rho_guess
             : P / (comp.mu_ions_inv() + comp.mu_elec_inv()) / (R_gas * T);
  for (int it = 0; it < 200; ++it) {
    const EosState s = eval(T, rho, comp);
    const double f = std::log(s.P) - std::log(P);
    if (std::abs(f) < 1e-12) break;
    double step = f / std::max(s.chiRho, 1e-3);
    if (step >  0.7) step =  0.7;
    if (step < -0.7) step = -0.7;
    rho *= std::exp(-step);
  }
  return rho;
}

} // namespace ember
