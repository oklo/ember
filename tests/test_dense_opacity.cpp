#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/eos_composite.hpp"
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0) {
  if (!ok) ++failures;
  std::printf("[%s] %s (%.8g)\n", ok ? "PASS" : "FAIL", label, value);
}

int main() {
  const std::string data = EMBER_DATA_DIR;
  AesopusOpacity low(data + "/opacity/aesopus21_gs98_z020.dat");
  OpalOpacity high(data + "/opacity/opal_gs98_z020.dat");
  BlendedOpacity opacity(low, high, 4.4, 4.5);
  const auto comp = solar_scaled(.7, .02);
  // Independent cells in the source's highR/X=.7 file: logR=3.0.
  for (const auto [lt, want] : {std::pair{3.48, -.5916}, std::pair{4.0, 3.1936}}) {
    const double T = std::pow(10.0, lt), rho = std::pow(10.0, 3.0 + 3 * (lt - 6));
    check(std::abs(std::log10(low.eval(T, rho, comp).kappa) - want) < 1e-12,
          "dense AESOPUS cell reproduces published value");
  }
  {
    std::ifstream file(data + "/opacity/aesopus21_gs98_z020.dat");
    std::size_t nx, nt, nr; file >> nx >> nt >> nr;
    std::string title; std::getline(file, title);
    std::vector<double> rr(nr), tt(nt);
    for (double& v : rr) file >> v;
    for (double& v : tt) file >> v;
    double worst = 0;
    for (std::size_t ix = 0; ix < nx; ++ix) {
      double X, Z; file >> X >> Z;
      for (double lt : tt) for (double lr : rr) {
        double want; file >> want;
        const double T = std::pow(10.0, lt), rho = std::pow(10.0, lr + 3 * (lt - 6));
        const double got = std::log10(low.eval(T, rho, solar_scaled(X, Z)).kappa);
        worst = std::max(worst, std::abs(got - want));
      }
    }
    check(worst < 1e-12, "all 149810 imported cells, including seams and edges, are reproduced", worst);
  }
  {
    constexpr double h = 1e-6;
    double worst = 0;
    for (double lt : {2.893, 3.147, 3.483, 3.697, 4.137, 4.387, 4.443, 4.497}) {
      for (double lr : {-.31, .99, 1.11, 3.17, 5.83}) {
        if (lt > 4.4 && lr > 1) continue;
        const double T = std::pow(10.0, lt), rho = std::pow(10.0, lr + 3 * (lt - 6));
        const auto c = solar_scaled(.63, .02);
        const auto k = opacity.eval(T, rho, c);
        const double dt = std::log(opacity.eval(T * std::exp(h), rho, c).kappa
                                / opacity.eval(T * std::exp(-h), rho, c).kappa) / (2 * h);
        const double dr = std::log(opacity.eval(T, rho * std::exp(h), c).kappa
                                / opacity.eval(T, rho * std::exp(-h), c).kappa) / (2 * h);
        worst = std::max({worst, std::abs(dt - k.dlnk_dlnT) / std::max(1., std::abs(dt)),
                                std::abs(dr - k.dlnk_dlnRho) / std::max(1., std::abs(dr))});
      }
    }
    check(worst < 3e-5, "dense opacity and upper-temperature blend derivatives match their values", worst);
  }
  {
    CompositeEos eos;
    GreyAtmosphere grey(eos, opacity);
    constexpr double g = 2e5, h = 1e-5;
    for (double Teff : {2500., 3000., 4000.}) {
      const auto a = grey.eval(Teff, g, comp);
      const auto finer_top = GreyAtmosphere(eos, opacity, {.tau_top = 1e-7}).eval(Teff, g, comp);
      const double dT = std::log(grey.eval(Teff * std::exp(h), g, comp).P
                              / grey.eval(Teff * std::exp(-h), g, comp).P) / (2 * h);
      const double dg = std::log(grey.eval(Teff, g * std::exp(h), comp).P
                              / grey.eval(Teff, g * std::exp(-h), comp).P) / (2 * h);
      check(std::abs(std::log(eos.eval(a.T, a.rho, comp).P / a.P)) < 1e-10,
            "cool high-gravity atmosphere closes EOS pressure", Teff);
      check(std::abs(std::log(a.P / finer_top.P)) < 1e-5,
            "cool atmosphere pressure is insensitive to a smaller top column", Teff);
      check(std::abs(dT - a.dlnP_dlnTeff) < 3e-4 && std::abs(dg - a.dlnP_dlng) < 3e-4,
            "cool high-gravity atmosphere sensitivities match independent integrations", Teff);
    }
  }
  {
    int rejected = 0;
    for (const auto [lt, lr] : {std::pair{1.9, 0.0}, std::pair{3.5, 6.1},
                              std::pair{4.45, 1.1}, std::pair{5.0, 1.1},
                              std::pair{3.5, 6.0 + 1e-10}, std::pair{2.0 - 1e-10, 0.0}}) {
      try { opacity.eval(std::pow(10., lt), std::pow(10., lr + 3 * (lt - 6)), comp); }
      catch (const std::domain_error&) { ++rejected; }
    }
    check(rejected == 6, "edge roundoff tolerance never authorizes extrapolation or hot gaps", rejected);
  }
  return failures ? 1 : 0;
}
