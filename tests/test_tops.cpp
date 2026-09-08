#include "ember/opacity_tops.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0) {
  if (!ok) ++failures;
  std::printf("[%s] %s (%.9g)\n", ok ? "PASS" : "FAIL", label, value);
}
int main() {
  const std::string data = EMBER_DATA_DIR;
  TopsOpacity high(data + "/opacity");
  AesopusOpacity low(data + "/opacity/aesopus21_gs98_z020.dat");
  BlendedOpacity opacity(low, high, 4.4, 4.5);
  const auto c = solar_scaled(.7, .02);
  const double keV = 1000 * constants::eV / constants::kB;
  for (const auto [energy, rho, kappa] : {std::array{.002, 1e-6, 1957.9},
                                        std::array{.5, 1000., 39.926}, std::array{2., 1000., 1.1659}})
    check(std::abs(high.eval(energy * keV, rho, c).kappa / kappa - 1) < 1e-12,
          "native TOPS opacity reproduces an independent printed source cell", kappa);
  for (const std::string which : {"low", "high"}) {
    const auto path = data + "/opacity/tops_gs98_x070_z020_" + which + ".dat";
    TabulatedOpacity table(path, "TOPS test", TabulatedOpacity::DensityAxis::logRho);
    std::ifstream input(path); std::size_t nx, nt, nr; input >> nx >> nt >> nr;
    std::string title; std::getline(input, title);
    std::vector<double> r(nr), t(nt);
    for (double& v : r) input >> v;
    for (double& v : t) input >> v;
    double X, Z; input >> X >> Z;
    double worst = 0;
    for (double lt : t) for (double lr : r) {
      double want; input >> want;
      const double got = std::log10(table.eval(std::pow(10., lt), std::pow(10., lr), c).kappa);
      worst = std::max(worst, std::abs(got - want));
    }
    check(worst < 1e-12, "every imported native-density cell is reproduced", worst);
  }
  {
    constexpr double h = 1e-6;
    double worst = 0;
    for (double lt : {4.443, 4.489, 4.781, 5.631, 5.687, 6.343, 7.017})
      for (double lr : {-6.13, -2.117, .817, 3.173}) {
        if (lt < 5.7 && lr > 2) continue;
        const double T = std::pow(10., lt), rho = std::pow(10., lr);
        const auto a = opacity.eval(T, rho, c);
        const double dt = std::log(opacity.eval(T * std::exp(h), rho, c).kappa
                                / opacity.eval(T * std::exp(-h), rho, c).kappa) / (2 * h);
        const double dr = std::log(opacity.eval(T, rho * std::exp(h), c).kappa
                                / opacity.eval(T, rho * std::exp(-h), c).kappa) / (2 * h);
        worst = std::max({worst, std::abs(dt - a.dlnk_dlnT), std::abs(dr - a.dlnk_dlnRho)});
      }
    check(worst < 3e-5, "density and temperature derivatives include both blends and the native density coordinate", worst);
  }
  {
    int rejected = 0;
    for (const auto [energy, rho] : {std::pair{.002, 1000.}, std::pair{.01, 10000.},
                                    std::pair{2.01, 1.}, std::pair{.5, 1e-11}})
      try { high.eval(energy * keV, rho, c); } catch (const std::domain_error&) { ++rejected; }
    try { high.eval(.5 * keV, 10, solar_scaled(.71, .02)); } catch (const std::domain_error&) { ++rejected; }
    try { high.eval(.5 * keV, 10, solar_scaled(.7, .021)); } catch (const std::domain_error&) { ++rejected; }
    check(rejected == 6, "server-clamped cells, missing coverage and other compositions are rejected", rejected);
    const auto a = high.eval(.5 * keV, 1e-6, c);
    check(std::abs(a.kappa / .34 - 1) < .05, "dilute hot gas approaches electron-scattering opacity", a.kappa);
  }
  return failures ? 1 : 0;
}
