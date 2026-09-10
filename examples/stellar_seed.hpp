#pragma once
#include "lane_emden.hpp"
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include "ember/model.hpp"
#include "ember/nuclear.hpp"

namespace ember::example {
// An n=3 or n=1.5 initial guess only. The solve replaces its transport and
// energy profile with the EOS, tabulated opacity, pp heating and atmosphere.
inline Model stellar_seed(std::size_t points, double mass, double radius,
                          const Composition& comp, const Nuclear& nuclear, const Atmosphere& atmosphere,
                          double index = 3.0, double seed_Teff = 0.0) {
  if (points < 32 || points > 8192 || !std::isfinite(mass) || !std::isfinite(radius)
      || !(mass > 0.0) || !(radius > 0.0))
    throw std::invalid_argument("stellar_seed: invalid mesh, mass or radius");
  if (index != 3 && index != 1.5) throw std::invalid_argument("stellar_seed: supported indices are 3 and 1.5");
  if (!std::isfinite(seed_Teff) || seed_Teff < 0)
    throw std::invalid_argument("stellar_seed: invalid trial effective temperature");
  const double outer_xi = index == 3 ? 6.875 : 3.65;
  std::vector<double> targets(points);
  for (std::size_t i = 0; i < points; ++i) {
    const double f = static_cast<double>(i) / static_cast<double>(points - 1);
    targets[i] = 1e-3 + (outer_xi - 1e-3) * f * f * (3.0 - 2.0 * f);
  }
  const auto le = lane_emden(targets, index);
  const double a = radius / outer_xi;
  const double rhoc = mass / (4.0 * M_PI * a * a * a * le.back().mass);
  const double Pc = 4 * M_PI / (index + 1) * constants::G * a * a * rhoc * rhoc;
  // Ideal-gas seed only; relaxation includes the actual degeneracy pressure.
  const double Tc = Pc / (constants::R_gas * rhoc * (comp.mu_ions_inv() + comp.mu_elec_inv()));
  Model m; m.M = mass;
  for (const auto& s : le) {
    m.m.push_back(mass * s.mass / le.back().mass);
    m.y.push_back({std::log(a * s.xi), std::log(rhoc * std::pow(s.theta, index)),
                   std::log(Tc * s.theta), 0.0});
    m.comp.push_back(comp);
  }
  m.m.back() = mass;
  double eps = nuclear.eval(m.T(0), m.rho(0), comp).eps;
  m.y.front().L = m.m.front() * eps;
  for (std::size_t i = 1; i < points; ++i) {
    const double next = nuclear.eval(m.T(i), m.rho(i), comp).eps;
    m.y[i].L = m.y[i - 1].L + 0.5 * (eps + next) * (m.m[i] - m.m[i - 1]);
    eps = next;
  }
  // A bounded atmosphere may need a trial luminosity inside its support.
  // This only sets the initial guess; relaxation still solves nuclear and
  // thermal balance, with no alteration of the selected atmosphere.
  if (seed_Teff > 0) {
    const double target = 4 * M_PI * constants::sigma_SB * radius * radius * std::pow(seed_Teff,4);
    const double scale_L = target / m.y.back().L;
    for (auto& state : m.y) state.L *= scale_L;
  }
  const double Ts = m.T(points - 1);
  const double Teff = std::pow(m.y.back().L / (4.0 * M_PI * constants::sigma_SB * radius * radius), 0.25);
  const auto surface = atmosphere.eval(Teff, constants::G * mass / (radius * radius), comp);
  const double drho = std::log(surface.rho / m.rho(points - 1));
  // Smoothly adjust the trial envelope toward its atmosphere. This is a
  // starting guess, not hydrostatic or thermal equilibrium. Its envelope
  // density change MUST also change the enclosed-mass mesh below.
  const double join_T = std::sqrt(Tc * Ts);
  for (std::size_t i = 0; i < points; ++i) {
    const double f = std::clamp(std::log(join_T / m.T(i)) / std::log(join_T / Ts), 0.0, 1.0);
    m.y[i].lnrho += f * drho;
    m.y[i].lnT += f * std::log(surface.T / Ts);
  }
  // Match the discrete mass equation and central sphere exactly. Scaling
  // rho and m together sets total mass without a surface mass discontinuity.
  m.m[0] = 4.0 * M_PI / 3.0 * m.rho(0) * std::pow(m.r(0), 3);
  for (std::size_t i = 1; i < points; ++i) {
    const double rb = 0.5 * (m.r(i) + m.r(i - 1));
    const double rhob = 0.5 * (m.rho(i) + m.rho(i - 1));
    m.m[i] = m.m[i - 1] + 4.0 * M_PI * std::pow(rb, 3) * rhob * (m.y[i].lnr - m.y[i - 1].lnr);
  }
  const double scale = mass / m.m.back();
  for (std::size_t i = 0; i < points; ++i) { m.m[i] *= scale; m.y[i].lnrho += std::log(scale); }
  m.m.back() = mass;
  for (std::size_t i = 1; i < points; ++i)
    if (!(m.m[i] > m.m[i - 1]))
      throw std::domain_error("stellar_seed: outer mass intervals are unresolved; reduce mesh points");
  return m;
}
} // namespace ember::example
