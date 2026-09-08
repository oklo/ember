#include "../examples/stellar_seed.hpp"
#include "ember/eos_cms19.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include "ember/relaxation.hpp"
#include <cstdio>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0) {
  if (!ok) ++failures;
  std::printf("[%s] %s (%.9g)\n", ok ? "PASS" : "FAIL", label, value);
}

// Numerical integration and spatial checks of the experimental 0.1 Msun
// equations. These do not calibrate the grey surface or repair the known
// CMS19 pressure/entropy consistency defect. No evolution is attempted.
int main() {
  const std::string data = EMBER_DATA_DIR;
  Cms19Eos eos(data + "/eos/cms19_h_tp.dat", data + "/eos/cms19_he_tp.dat", Cms19Eos::Metals::helium_proxy);
  AesopusOpacity low(data + "/opacity/aesopus21_gs98_z020.dat");
  TopsOpacity high(data + "/opacity");
  BlendedOpacity opacity(low, high, 4.4, 4.5);
  PPChains nuclear;
  GreyAtmosphere atmosphere(eos, opacity, {.tau_top = .001});
  const Physics physics{&eos, &opacity, &nuclear, 1.9};
  const auto comp = solar_scaled(.7, .02);
  std::vector<double> radii, luminosities, virials;
  Model finest;
  for (std::size_t n : {512UL, 1024UL, 2048UL}) {
    const auto seed = example::stellar_seed(n, .1 * constants::Msun, .15 * constants::Rsun,
                                           comp, nuclear, atmosphere, 1.5);
    const auto result = relax(seed, physics, atmosphere);
    check(result.converged && result.residual < 1e-9 && result.correction < 1e-8,
          "0.1 Msun satisfies residual and undamped correction tolerances", result.residual);
    if (!result.converged) { std::printf("%s\n", result.message.c_str()); return 1; }
    const auto& m = result.model;
    double L = m.m[0] * nuclear.eval(m.T(0), m.rho(0), comp).eps;
    double pressure = 3 * m.m[0] * eos.eval(m.T(0), m.rho(0), comp).P / m.rho(0);
    double gravity = .6 * constants::G * m.m[0] * m.m[0] / m.r(0);
    bool ordered = m.y[0].L > 0;
    for (std::size_t i = 1; i < n; ++i) {
      const double dm = m.m[i] - m.m[i - 1];
      L += .5 * dm * (nuclear.eval(m.T(i), m.rho(i), comp).eps
                    + nuclear.eval(m.T(i - 1), m.rho(i - 1), comp).eps);
      pressure += 1.5 * dm * (eos.eval(m.T(i), m.rho(i), comp).P / m.rho(i)
                            + eos.eval(m.T(i - 1), m.rho(i - 1), comp).P / m.rho(i - 1));
      gravity += .5 * constants::G * dm * (m.m[i] / m.r(i) + m.m[i - 1] / m.r(i - 1));
      ordered &= dm > 0 && m.r(i) > m.r(i - 1) && m.rho(i) < m.rho(i - 1)
                 && m.T(i) < m.T(i - 1) && m.y[i].L > 0;
    }
    const double surface = 4 * M_PI * std::pow(m.r(n - 1), 3)
                           * eos.eval(m.T(n - 1), m.rho(n - 1), comp).P;
    const double virial = std::abs((pressure - surface) / gravity - 1);
    check(ordered, "positive mass shells and outward luminosity; ordered radius, density and temperature");
    check(m.m[0] / m.M < 2e-10, "unresolved center is a negligible fraction of stellar mass", m.m[0] / m.M);
    check(std::abs(L / m.y.back().L - 1) < 1e-8,
          "surface luminosity balances integrated nuclear heating", L / m.y.back().L - 1);
    check(virial < 5e-5, "independent virial integral closes within spatial error", virial);
    radii.push_back(m.r(n - 1)); luminosities.push_back(m.y.back().L); virials.push_back(virial);
    std::printf("n=%zu R/Rsun=%.10g L/Lsun=%.10g\n", n, radii.back() / constants::Rsun,
                luminosities.back() / constants::Lsun);
    finest = m;
  }
  check(std::abs(radii[2] / radii[1] - 1) < .0002
         && std::abs(luminosities[2] / luminosities[1] - 1) < .001,
        "1024 and 2048 points agree within 0.02% radius and 0.1% luminosity");
  check(virials[2] < .3 * virials[1] && virials[1] < .3 * virials[0],
        "independent virial error decreases approximately quadratically with spacing");
  // Start from the converged model: this isolates the atmospheric truncation
  // from the seed's small, atmosphere-dependent change in mass mesh.
  GreyAtmosphere deeper(eos, opacity, {.tau_top = .002});
  const auto shifted = relax(finest, physics, deeper);
  check(shifted.converged, "doubling atmosphere starting depth converges on the identical mass mesh");
  if (shifted.converged) {
    check(std::abs(shifted.model.r(finest.size() - 1) / finest.r(finest.size() - 1) - 1) < .0001
           && std::abs(shifted.model.y.back().L / finest.y.back().L - 1) < .0002,
          "finite atmosphere top column has less than 0.02% effect on R and L");
  }
  return failures ? 1 : 0;
}
