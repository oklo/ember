#pragma once
#include <cmath>

namespace ember::fermi {

// Relativistic Fermi-Dirac integrals for an electron gas of arbitrary
// degeneracy and arbitrary relativity, in the dimensionless momentum
// x = p/(m_e c) with beta = kT/(m_e c^2) and eta = mu_kin/kT:
//
//   n = A * I_n,   I_n = \int x^2 f dx
//   P = (A m c^2 /3) * I_P,   I_P = \int x^4/sqrt(1+x^2) f dx
//   u = A m c^2 * I_u,        I_u = \int x^2 (sqrt(1+x^2)-1) f dx
//   f = 1/(1 + exp((sqrt(1+x^2)-1)/beta - eta))
//
// with A = 8 pi (m c/h)^3.  Evaluating these directly, rather than through a
// fit valid in patches, means the non-degenerate, degenerate, non-relativistic
// and relativistic corners are all the same code - and the 0.1 Msun white
// dwarf core, at p_F/mc ~ 0.8, sits in none of them.
//
// Returned alongside each integral is its derivative with respect to eta and
// to ln beta, which is all the chain rule needs to produce chi_T and chi_rho
// analytically.
struct Integrals {
  double In{}, Ip{}, Iu{};
  double dIn_deta{}, dIp_deta{}, dIu_deta{};
  double dIn_dlnb{}, dIp_dlnb{}, dIu_dlnb{};
  // Susceptibility response, needed for finite-degeneracy nuclear screening.
  double d2In_deta2{}, d2In_detadlnb{};
};

Integrals evaluate(double eta, double beta);

// Dimensionless entropy integral -integral x^2 [f ln f+(1-f)ln(1-f)] dx.
// Specific electron entropy is A*kB*Is/rho. The degenerate branch evaluates
// the narrow Fermi shell directly, without subtracting chemical energy.
double entropy(double eta, double beta, const Integrals&);

// Second derivatives after enforcing n proportional to rho. Temperature
// variations keep n fixed; density variations keep T fixed. Centered kernels
// avoid subtracting O(eta^2) integrals to recover a small thermal response.
struct DensityResponse {
  double dIp_dlnT{}, dIu_dlnT{};
  double d2Ip_dlnT2{}, d2Ip_dlnTdlnRho{}, d2Ip_dlnRho2{};
  double d2Iu_dlnT2{}, d2Iu_dlnTdlnRho{};
};
DensityResponse density_response(double eta, double beta, const Integrals&, bool second=true);

} // namespace ember::fermi
