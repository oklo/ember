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
};

Integrals evaluate(double eta, double beta);

} // namespace ember::fermi
