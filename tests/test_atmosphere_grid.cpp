#include "../examples/stellar_seed.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace ember;
namespace {
int failures{};
void check(bool v, const char *message) {
  std::printf("[%s] %s\n", v ? "PASS" : "FAIL", message);
  if (!v)
    ++failures;
}
bool near(double a, double b, double tol = 1e-10) {
  return std::abs(a - b) <= tol * std::max({1., std::abs(a), std::abs(b)});
}
template <class F> bool throws(F &&f) {
  try {
    f();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}
class Gas final : public Eos {
public:
  EosState eval(double t, double rho, const Composition &c) const override {
    EosState e{};
    const double pg = constants::R_gas * c.mu_ions_inv() * rho * t;
    const double pr = constants::a_rad * std::pow(t, 4) / 3;
    e.P = pg + pr;
    e.chiT = (pg + 4 * pr) / e.P;
    e.chiRho = pg / e.P;
    return e;
  }
  const char *name() const override { return "synthetic atmosphere test gas"; }
};
double temperature(double x, double y, double t, double g) {
  return .25 + .95 * t + .02 * g + .1 * x - .2 * y + .04 * x * y * t * g;
}
double pressure(double x, double y, double t, double g) {
  return 5. + .2 * t + .3 * g - .4 * x + .15 * y - .03 * x * y * t * g;
}
std::string fixture() {
  const auto c = solar_scaled(.7, .02);
  std::ostringstream s;
  s << std::setprecision(17);
  s << "EMBER_COMPOSITION_ATMOSPHERE 1\nsource \"synthetic multiaffine "
       "fixture, not physical data\"\n"
       "approximation \"test source\"\nbasis baryon_mass\ntau 100\nmetals";
  for (std::size_t j = 3; j < NSPEC; ++j)
    s << ' ' << c.X[j];
  s << "\nhydrogen 3 .4 .55 .7\nhelium3 2 0 .12\nlog_teff 2 3.4 3.6\nlog_g 2 "
       "4.5 5.5\ndata\n";
  for (double x : {.4, .55, .7})
    for (double y : {0., .12})
      for (double t : {3.4, 3.6})
        for (double g : {4.5, 5.5})
          s << temperature(x, y, t, g) << ' ' << pressure(x, y, t, g) << '\n';
  return s.str();
}
} // namespace
int main() {
  Gas eos;
  const auto proxy = CompositionAtmosphereGrid::Mixture::allow_documented_proxy;
  std::istringstream in(fixture());
  CompositionAtmosphereGrid grid(eos, in, proxy);
  const auto bounds = grid.support();
  check(bounds.hydrogen == std::array{.4, .7} &&
            bounds.helium3 == std::array{0., .12} &&
            near(bounds.teff[0], std::pow(10., 3.4)) &&
            near(bounds.gravity[1], std::pow(10., 5.5)),
        "physical support exposed for initialization without extrapolation");
  auto c = solar_scaled(.51, .02);
  c.basis = AbundanceBasis::baryon_mass;
  c.X[1] = .08;
  c.X[2] -= .08;
  const double t = std::pow(10., 3.47), g = std::pow(10., 5.14);
  const auto s = grid.eval(t, g, c);
  const double radius = .15 * constants::Rsun;
  PPChains nuclear;
  const auto seed = example::stellar_seed(64, .1 * constants::Msun, radius, c,
                                          nuclear, grid, 1.5, 2800);
  const double seed_teff = std::pow(
      seed.y.back().L / (4 * M_PI * constants::sigma_SB * radius * radius),
      .25);
  check(near(seed_teff, 2800) && seed.m.back() == seed.M,
        "bounded-grid trial luminosity and baryonic mass initialization");
  check(near(std::log10(s.T), temperature(.51, .08, 3.47, 5.14)) &&
            near(std::log10(s.Pgas), pressure(.51, .08, 3.47, 5.14)),
        "independent 4D multiaffine states");
  check(near(eos.eval(s.T, s.rho, c).P, s.P) && s.tau == 100,
        "actual-composition EOS inversion and matching depth");
  constexpr double h = 1e-5;
  const auto tp = grid.eval(t * std::exp(h), g, c),
             tm = grid.eval(t * std::exp(-h), g, c);
  const auto gp = grid.eval(t, g * std::exp(h), c),
             gm = grid.eval(t, g * std::exp(-h), c);
  check(near(s.dlnT_dlnTeff, std::log(tp.T / tm.T) / (2 * h), 1e-8) &&
            near(s.dlnP_dlnTeff, std::log(tp.P / tm.P) / (2 * h), 1e-8) &&
            near(s.dlnT_dlng, std::log(gp.T / gm.T) / (2 * h), 1e-8) &&
            near(s.dlnP_dlng, std::log(gp.P / gm.P) / (2 * h), 1e-8),
        "four thermal/gravity derivatives including radiation");
  const auto response = grid.composition_response(t, g, c);
  for (std::size_t k = 0; k < 2; ++k) {
    auto a = c, b = c;
    a.X[k] += h;
    a.X[2] -= h;
    b.X[k] -= h;
    b.X[2] += h;
    const auto sa = grid.eval(t, g, a), sb = grid.eval(t, g, b);
    check(near(k == 0 ? response.dlnT_dXH : response.dlnT_dX3,
               std::log(sa.T / sb.T) / (2 * h), 1e-8) &&
              near(k == 0 ? response.dlnP_dXH : response.dlnP_dX3,
                   std::log(sa.P / sb.P) / (2 * h), 1e-8),
          "composition derivative along conserved baryonic substitution");
  }
  for (double x : {.4, .55, .7})
    for (double y : {0., .12}) {
      auto a = solar_scaled(x, .02);
      a.basis = c.basis;
      a.X[1] = y;
      a.X[2] -= y;
      for (double lt : {3.4, 3.6})
        for (double lg : {4.5, 5.5}) {
          const auto v = grid.eval(std::pow(10., lt), std::pow(10., lg), a);
          check(near(std::log10(v.T), temperature(x, y, lt, lg)),
                "source knot and closed edges reproduced");
        }
    }
  check(throws([&] { grid.eval(1e3, g, c); }) &&
            throws([&] { grid.eval(t, 1e6, c); }),
        "no Teff/gravity extrapolation");
  auto bad = c;
  bad.basis = AbundanceBasis::atomic_mass;
  check(!grid.covers(t, g, bad), "abundance basis enforced");
  bad = c;
  bad.X[0] = .3;
  bad.X[2] += .21;
  check(!grid.covers(t, g, bad), "no composition extrapolation");
  bad = c;
  bad.X[3] += 1e-7;
  bad.X[2] -= 1e-7;
  check(!grid.covers(t, g, bad), "metal pattern enforced");
  bad = c;
  bad.X[2] -= .01;
  check(!grid.covers(t, g, bad), "normalization enforced");
  bad = c;
  bad.X[0] = std::numeric_limits<double>::quiet_NaN();
  check(!grid.covers(t, g, bad), "nonfinite composition rejected");
  check(throws([&] {
          std::istringstream f(fixture());
          CompositionAtmosphereGrid rejected(eos, f);
        }),
        "source approximations require explicit selection");
  for (const auto &malformed :
       {fixture() + "extra", fixture().substr(0, fixture().size() - 35)})
    check(throws([&] {
            std::istringstream f(malformed);
            CompositionAtmosphereGrid rejected(eos, f, proxy);
          }),
          "truncated or trailing source states rejected");
  return failures ? 1 : 0;
}
