#pragma once
// Controlled solver benchmark, not production stellar physics. A monatomic
// ideal gas, constant opacity/heating, and an explicitly artificial atmosphere
// admit an independent n=3 Lane-Emden solution with grad_rad=1/4 everywhere.
#include "ember/atmosphere.hpp"
#include "lane_emden.hpp"
#include "ember/constants.hpp"
#include "ember/model.hpp"
#include "ember/nuclear.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace ember::example {
inline constexpr double gas_constant = constants::R_gas / 0.6;

class PolytropeEos final : public Eos {
public:
  EosState eval(double T, double rho, const Composition&) const override {
    if (!(T > 0.0) || !(rho > 0.0) || !std::isfinite(T) || !std::isfinite(rho))
      throw std::domain_error("polytrope EOS: invalid state");
    EosState e{};
    e.P = gas_constant * rho * T; e.E = 1.5 * gas_constant * T;
    e.chiT = e.chiRho = e.delta = 1.0;
    e.cv = 1.5 * gas_constant; e.cp = 2.5 * gas_constant;
    e.grad_ad = 0.4; e.Gamma1 = 5.0 / 3.0;
    return e;
  }
  EosResponse eval_with_derivatives(double T, double rho, const Composition& comp) const override {
    EosResponse response{}; response.state = eval(T, rho, comp); return response;
  }
  const char* name() const override { return "benchmark ideal gas"; }
};

class PolytropeHeating final : public Nuclear {
public:
  double epsilon{1.0};  // erg/g/s; deliberately independent of T, rho and composition
  NuclearState eval(double, double, const Composition&) const override {
    NuclearState n{}; n.eps = epsilon; return n;
  }
  const char* name() const override { return "benchmark constant heating"; }
};

class PolytropeOpacity final : public Opacity {
public:
  double kappa{};
  OpacityState eval(double, double, const Composition&) const override { return {kappa, 0.0, 0.0}; }
  const char* name() const override { return "benchmark constant opacity"; }
};

class PolytropeAtmosphere final : public Atmosphere {
public:
  double coefficient{};  // P/T^4, chosen to enforce the polytropic entropy
  AtmosphereState eval(double Teff, double gravity, const Composition&) const override {
    if (!(Teff > 0.0) || !(gravity > 0.0) || !std::isfinite(Teff) || !std::isfinite(gravity))
      throw std::domain_error("polytrope boundary: invalid state");
    AtmosphereState a{};
    a.T = Teff; a.P = coefficient * std::pow(Teff, 4); a.Pgas = a.P;
    a.rho = a.P / (gas_constant * Teff);
    a.dlnT_dlnTeff = 1.0; a.dlnP_dlnTeff = 4.0;
    return a;
  }
  const char* name() const override { return "artificial polytrope boundary (not an atmosphere grid)"; }
};

struct RadiativePolytrope {
  PolytropeEos eos;
  PolytropeHeating heating;
  PolytropeOpacity opacity;
  PolytropeAtmosphere atmosphere;
  Model reference;
  double K{}, central_density{}, central_temperature{};

  explicit RadiativePolytrope(std::size_t points, double inner_xi = 1e-3) {
    if (points < 8 || points > 100000 || !(inner_xi >= 1e-6) || !(inner_xi < 4.0))
      throw std::invalid_argument("polytrope benchmark: invalid mesh");
    std::vector<double> xi(points);
    for (std::size_t i = 0; i < points; ++i) {
      const double fraction = static_cast<double>(i) / static_cast<double>(points - 1);
      xi[i] = inner_xi + (4.0 - inner_xi) * fraction * fraction;
    }
    xi.back() = 4.0;
    const auto solution = lane_emden(xi);
    const auto& surface = solution.back();
    const double mass = 0.1 * constants::Msun;
    K = M_PI * constants::G * std::pow(mass / (4.0 * M_PI * surface.mass), 2.0 / 3.0);
    const double length_coefficient = std::sqrt(K / (M_PI * constants::G));
    // Stefan-Boltzmann matching fixes rho_c for the chosen truncation xi=4,
    // mass and constant heating. The atmosphere enforces P/T^4=R_gas^4/K^3.
    const double flux_coefficient = 4.0 * M_PI * constants::sigma_SB
        * length_coefficient * length_coefficient * surface.xi * surface.xi
        * std::pow(K / gas_constant * surface.theta, 4);
    central_density = std::pow(heating.epsilon * mass / flux_coefficient, 1.5);
    central_temperature = K / gas_constant * std::cbrt(central_density);
    opacity.kappa = 4.0 * M_PI * constants::a_rad * constants::c * constants::G * std::pow(K, 3)
                   / (3.0 * heating.epsilon * std::pow(gas_constant, 4));
    atmosphere.coefficient = std::pow(gas_constant, 4) / std::pow(K, 3);
    const double length = length_coefficient / std::cbrt(central_density);
    reference.M = mass;
    for (const auto& s : solution) {
      const double enclosed = mass * s.mass / surface.mass;
      reference.m.push_back(enclosed);
      reference.y.push_back({std::log(length * s.xi), std::log(central_density * std::pow(s.theta, 3)),
                             std::log(central_temperature * s.theta), heating.epsilon * enclosed});
      reference.comp.push_back(solar_scaled(0.7, 0.014));
    }
    reference.m.back() = mass;
  }

  Model initial(double amplitude = 0.08) const {
    Model model = reference;
    for (std::size_t i = 0; i < model.size(); ++i) {
      const double phase = M_PI * static_cast<double>(i) / static_cast<double>(model.size() - 1);
      model.y[i].lnr += amplitude * (1.0 + 0.3 * std::sin(phase));
      model.y[i].lnrho -= 0.8 * amplitude * (1.0 + 0.2 * std::cos(phase));
      model.y[i].lnT += 0.7 * amplitude;
      model.y[i].L *= 1.0 + amplitude * (0.4 + std::sin(phase));
    }
    return model;
  }
};

} // namespace ember::example
