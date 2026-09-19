#include "ember/opacity_table.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* text) {
  if (!ok) ++failures;
  std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", text);
}

int main() {
  const auto root = std::filesystem::temp_directory_path() /
      ("ember-opacity-support-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(root);
  auto write = [&](const char* name, bool version2, bool incomplete, int bad_count = 0) {
    const auto path = root / name;
    std::ofstream file(path);
    file << std::setprecision(17);
    if (version2) file << "EMBER_OPACITY_TABLE 2 ";
    file << "2 8 7 analytic support control\n-3 -2 -1 0 1 2 3\n4 4.5 5 5.5 6 6.5 7 7.5\n";
    for (int ix = 0; ix < 2; ++ix) {
      const double X = .7 * ix;
      file << X << " 0.02\n";
      for (int it = 0; it < 8; ++it) {
        int count = incomplete ? std::min(4 + it, 7) : 7;
        if (incomplete && ix == 1 && it == 6) count = 6;
        if (bad_count && ix == 0 && it == 0) count = bad_count;
        if (version2) file << count << ' ';
        for (int ir = 0; ir < std::min(count, 7); ++ir)
          file << .2 * (ir - 3.) - .3 * (4. + .5 * it) + .1 * X << ' ';
        file << '\n';
      }
    }
    return path;
  };
  using Axis = TabulatedOpacity::DensityAxis;
  TabulatedOpacity old(write("old.dat", false, false), "rectangle", Axis::logRho);
  TabulatedOpacity complete(write("complete.dat", true, false), "complete v2", Axis::logRho);
  TabulatedOpacity limited(write("limited.dat", true, true), "limited v2", Axis::logRho);
  bool exact = true, analytic = true;
  for (double X : {0., .23, .7})
    for (double lt : {4., 4.2, 5.37, 6.3, 7.5})
      for (double lr : {-3., -.31, 0.}) {
        const auto c = solar_scaled(X, .02);
        const double T = std::pow(10., lt), rho = std::pow(10., lr);
        const auto a = old.eval(T, rho, c), b = complete.eval(T, rho, c), d = limited.eval(T, rho, c);
        exact &= a.kappa == b.kappa && a.dlnk_dlnT == b.dlnk_dlnT
            && a.dlnk_dlnRho == b.dlnk_dlnRho && a.dlnk_dX == b.dlnk_dX;
        analytic &= std::abs(std::log10(d.kappa) - (.2 * lr - .3 * lt + .1 * X)) < 1e-13
            && std::abs(d.dlnk_dlnT + .3) < 1e-12 && std::abs(d.dlnk_dlnRho - .2) < 1e-12
            && std::abs(d.dlnk_dX - .1 * std::log(10.)) < 1e-12;
      }
  check(exact, "complete version-2 table preserves all values and derivatives exactly");
  check(analytic, "unequal density support retains analytic opacity and derivatives");
  const auto comp = solar_scaled(.2, .02);
  const auto cool = limited.density_range(std::pow(10., 4.2), comp);
  const auto hot = limited.density_range(std::pow(10., 5.8), comp);
  check(cool && hot && cool->min == 1e-3 && cool->max == 1. && hot->max == 100.,
        "density range includes all temperature and composition derivative inputs");
  const auto warm = limited.eval(std::pow(10., 5.8), 50., comp);
  check(std::abs(warm.dlnk_dlnRho - .2) < 1e-12, "supported hotter rows extend beyond the cool density limit");
  int rejected = 0;
  for (const auto point : {std::array{4.2, 50., .2}, std::array{6.8, 1000., 0.},
                           std::array{5.8, 100.*(1.+1e-10), .2}}) {
    const double T = std::pow(10., point[0]);
    check(!limited.covers(T, point[1], point[2]), "coverage excludes a missing derivative stencil");
    try { limited.eval(T, point[1], solar_scaled(point[2], .02)); }
    catch (const std::domain_error&) { ++rejected; }
  }
  check(rejected == 3, "evaluation rejects unavailable rows and extrapolation");
  const auto boundary = limited.eval(std::pow(10., 5.8), hot->max, comp);
  check(std::isfinite(boundary.kappa), "closed supported density boundary remains usable");
  rejected = 0;
  for (int count : {3, 8}) {
    try { TabulatedOpacity invalid(write("invalid.dat", true, true, count), "invalid", Axis::logRho); }
    catch (const std::runtime_error&) { ++rejected; }
  }
  check(rejected == 2, "truncated and excessive density prefixes are rejected");
  std::filesystem::remove_all(root);
  return failures ? 1 : 0;
}
