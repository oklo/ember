#include "../examples/stellar_seed.hpp"
#include "ember/eos_composite.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/relaxation.hpp"
#include <cstdio>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0) {
  if (!ok) ++failures;
  std::printf("[%s] %s (%.9g)\n", ok ? "PASS" : "FAIL", label, value);
}

// An integration test of the actual modules at fixed composition, NOT an
// observational stellar calibration. Ionization/H2 EOS and non-grey surface
// physics are absent. Refining this case formerly failed the linear check.
int main() {
  const std::string data = EMBER_DATA_DIR;
  CompositeEos eos;
  AesopusOpacity low(data + "/opacity/aesopus21_gs98_z020.dat");
  OpalOpacity high(data + "/opacity/opal_gs98_z020.dat");
  BlendedOpacity opacity(low, high, 4.4, 4.5);
  PPChains nuclear;
  GreyAtmosphere atmosphere(eos, opacity);
  const Physics physics{&eos, &opacity, &nuclear, 1.9};
  const auto composition = solar_scaled(.7, .02);
  std::vector<double> radii, luminosities, virials;
  for (std::size_t n : {128UL, 256UL, 512UL}) {
    std::printf("%zu mesh points\n", n);
    const auto seed = example::stellar_seed(n, .5 * constants::Msun, .6 * constants::Rsun,
                                           composition, nuclear, atmosphere);
    const auto result = relax(seed, physics, atmosphere);
    check(result.converged && result.residual < 1e-9 && result.correction < 1e-8,
          "fixed-composition star satisfies residual and undamped-correction criteria", result.residual);
    if (!result.converged) {
      std::printf("%s\n", result.message.c_str());
      return 1;
    }
    const auto& m = result.model;
    double L = m.m[0] * nuclear.eval(m.T(0), m.rho(0), composition).eps;
    double pressure = 3 * m.m[0] * eos.eval(m.T(0), m.rho(0), composition).P / m.rho(0);
    double gravity = .6 * constants::G * m.m[0] * m.m[0] / m.r(0);
    bool monotone = true;
    for (std::size_t i = 1; i < n; ++i) {
      const double dm = m.m[i] - m.m[i - 1];
      L += .5 * dm * (nuclear.eval(m.T(i), m.rho(i), composition).eps
                    + nuclear.eval(m.T(i - 1), m.rho(i - 1), composition).eps);
      pressure += 1.5 * dm * (eos.eval(m.T(i), m.rho(i), composition).P / m.rho(i)
                            + eos.eval(m.T(i - 1), m.rho(i - 1), composition).P / m.rho(i - 1));
      gravity += .5 * constants::G * dm * (m.m[i] / m.r(i) + m.m[i - 1] / m.r(i - 1));
      monotone &= m.r(i) > m.r(i - 1) && m.T(i) < m.T(i - 1)
               && m.rho(i) < m.rho(i - 1) && m.y[i].L > 0;
    }
    const double boundary = 4 * M_PI * std::pow(m.r(n - 1), 3)
                            * eos.eval(m.T(n - 1), m.rho(n - 1), composition).P;
    const double virial = std::abs((pressure - boundary) / gravity - 1);
    check(monotone, "radius, density and temperature remain ordered with positive outward luminosity");
    check(std::abs(L / m.y.back().L - 1) < 1e-8,
          "surface luminosity equals the integrated nuclear heating", L / m.y.back().L - 1);
    check(virial < 0.01, "independent integrated virial balance closes within spatial error", virial);
    radii.push_back(m.r(n - 1)); luminosities.push_back(m.y.back().L); virials.push_back(virial);
  }
  check(std::abs(radii[2] / radii[1] - 1) < .002
         && std::abs(luminosities[2] / luminosities[1] - 1) < .003,
        "finest two meshes agree within 0.2% in radius and 0.3% in luminosity");
  check(std::abs(radii[2] - radii[1]) < std::abs(radii[1] - radii[0])
         && std::abs(luminosities[2] - luminosities[1]) < std::abs(luminosities[1] - luminosities[0])
         && virials[2] < virials[1] && virials[1] < virials[0],
        "radius, luminosity and independent virial errors decrease on refinement");
  return failures ? 1 : 0;
}
