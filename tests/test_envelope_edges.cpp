#include "ember/envelope_atmosphere.hpp"
#include "ember/constants.hpp"
#include "ember/eos.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace ember;
namespace {
struct Gas final : Eos {
  EosState eval(double T, double rho, const Composition& c) const override {
    const double R = constants::R_gas * (c.mu_ions_inv() + c.mu_elec_inv());
    EosState s{};
    s.P = rho * R * T; s.cp = 2.5 * R; s.cv = 1.5 * R; s.E = s.cv * T;
    s.delta = s.chiRho = s.chiT = 1; s.grad_ad = .4;
    return s;
  }
  const char* name() const override { return "ideal gas"; }
};
struct Transparent final : Opacity {
  double maximum_T = std::numeric_limits<double>::infinity();
  OpacityState eval(double T, double, const Composition&) const override {
    if (T > maximum_T) throw std::domain_error("missing interior opacity");
    OpacityState s{}; s.kappa = 1e-20; return s;
  }
  const char* name() const override { return "transparent test layer"; }
};
struct Boundary final : Atmosphere {
  double min_T = 0, max_T = 1e9, min_g = 0, max_g = 1e12;
  Gas gas;
  AtmosphereState eval(double T, double g, const Composition& c) const override {
    if (T < min_T || T > max_T || g < min_g || g > max_g)
      throw std::domain_error("outside surface support");
    AtmosphereState s{};
    s.T = 2 * T; s.P = 1e8; s.rho = gas.rho_from_PT(s.T, s.P, c);
    s.dlnT_dlnTeff = 1;
    return s;
  }
  const char* name() const override { return "bounded surface"; }
};
}

int main() try {
  Gas gas; Transparent opacity; Boundary broad;
  auto c = solar_scaled(.65, .08);
  const double M = .1 * constants::Msun, R = .1 * constants::Rsun;
  const double surface_g = constants::G * M / (R * R);
  int failures = 0;
  auto check = [&](bool ok, const char* name, double error = 0) {
    std::cout << (ok ? "PASS " : "FAIL ") << name << " " << error << '\n';
    failures += !ok;
  };
  for (int test = 0; test < 5; ++test) {
    Boundary limited;
    double T = 2850;
    if (test < 3) {
      limited.min_T = 2500; limited.max_T = 3200;
      T = test == 0 ? 2500.01 : test == 1 ? 2850 : 3199.99;
    } else if (test == 3) limited.min_g = surface_g * (1 - 2e-5);
    else limited.max_g = surface_g * (1 + 2e-5);
    EnvelopeAtmosphere edge(limited, opacity, gas, 1.9, M, 1e-5*M, 40);
    EnvelopeAtmosphere control(broad, opacity, gas, 1.9, M, 1e-5*M, 40);
    edge.integration_tolerance(1e-9); control.integration_tolerance(1e-9);
    edge.evaluation_threads(2);
    const auto surface = edge.from_photosphere(T, R, c);
    const double Tb = T * std::sqrt(R / surface.r_base);
    const double gb = constants::G * M / (surface.r_base * surface.r_base);
    const auto value = edge.eval_value(Tb, gb, c), ref = control.eval(Tb, gb, c);
    try {
      const auto got = edge.eval(Tb, gb, c);
      const double error = std::max({std::abs(got.dlnT_dlnTeff-ref.dlnT_dlnTeff),
          std::abs(got.dlnT_dlng-ref.dlnT_dlng),
          std::abs(got.dlnP_dlnTeff-ref.dlnP_dlnTeff),
          std::abs(got.dlnP_dlng-ref.dlnP_dlng)});
      check(std::abs(got.T/value.T-1) < 1e-7 && std::abs(got.P/value.P-1) < 1e-7
            && error < 2e-4, "supported edge matches unrestricted boundary", error);
    } catch (const std::exception& e) {
      check(false, e.what());
    }
  }
  Boundary limited; limited.min_T = 2500; limited.max_T = 3200;
  EnvelopeAtmosphere edge(limited, opacity, gas, 1.9, M, 1e-5*M, 40);
  bool rejected = false;
  try { edge.from_photosphere(2400, R, c); }
  catch (const std::domain_error&) { rejected = true; }
  check(rejected, "unsupported physical surface is still rejected");

  // The central value is supported; only the hotter derivative query fails
  // inside the envelope. It must not be mistaken for a surface-table edge.
  Transparent incomplete; incomplete.maximum_T = 2*2900*std::exp(4e-5);
  EnvelopeAtmosphere interior(broad, incomplete, gas, 1.9, M, 1e-5*M, 40);
  const auto surface = interior.from_photosphere(2900, R, c);
  const double Tb = 2900 * std::sqrt(R / surface.r_base);
  const double gb = constants::G * M / (surface.r_base * surface.r_base);
  interior.eval_value(Tb, gb, c);
  rejected = false;
  try { interior.eval(Tb, gb, c); }
  catch (const std::domain_error& e) {
    rejected = std::string(e.what()).find("missing interior opacity") != std::string::npos;
  }
  check(rejected, "interior failure during derivative evaluation propagates");
  return failures ? 1 : 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n'; return 2;
}
