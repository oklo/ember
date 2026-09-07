#include "ember/henyey.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0.0) {
  if (!ok) ++failures;
  std::printf("  [%s] %-60s %.6g\n", ok ? "PASS" : "FAIL", label, value);
}

// Independent dense Gauss-Jordan reference, with no use of the Henyey blocks
// after the complete matrix is assembled. Used only for small systems.
static std::vector<double> dense_solve(std::vector<std::vector<double>> a) {
  const std::size_t n = a.size();
  for (auto& row : a) {
    double norm = 0.0;
    for (std::size_t j = 0; j < n; ++j) norm = std::max(norm, std::abs(row[j]));
    for (double& v : row) v /= norm;
  }
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    std::swap(a[k], a[pivot]);
    const double value = a[k][k];
    if (value == 0.0) throw std::runtime_error("singular test reference");
    for (double& v : a[k]) v /= value;
    for (std::size_t i = 0; i < n; ++i) {
      if (i == k) continue;
      const double factor = a[i][k];
      for (std::size_t j = 0; j <= n; ++j) a[i][j] -= factor * a[k][j];
    }
  }
  std::vector<double> result(n);
  for (std::size_t i = 0; i < n; ++i) result[i] = a[i][n];
  return result;
}

int main() {
  std::printf("ember Henyey block elimination\n");
  for (std::size_t points : {1UL, 2UL, 7UL, 31UL, 2048UL}) {
    BoundaryBlock inner{}, outer{};
    // These inner rows have zero leading columns, forcing elimination to
    // pivot with the zone equations instead of assuming an invertible 2x2.
    inner.dfdy[0][2] = 1.0; inner.dfdy[1][3] = 1.0;
    outer.dfdy[0][0] = 1.0; outer.dfdy[1][1] = 1.0;
    std::vector<ZoneResidual> zones(points - 1);
    const double h = 1.0 / static_cast<double>(points);
    for (std::size_t i = 0; i < zones.size(); ++i) {
      auto& z = zones[i];
      for (std::size_t k = 0; k < NVAR; ++k) {
        const double unit = std::pow(10.0, 80.0 * std::sin(static_cast<double>(4 * i + k)));
        for (std::size_t v = 0; v < NVAR; ++v) {
          const double angle = static_cast<double>(7 * i + 3 * k + v);
          z.dfdy_lo[k][v] = unit * (-(k == v ? 1.0 : 0.0) + 0.1 * h * std::sin(angle));
          z.dfdy_hi[k][v] = unit * ((k == v ? 1.0 : 0.0) + 0.2 * h * std::cos(angle));
        }
      }
    }
    std::vector<std::array<double, NVAR>> expected(points);
    for (std::size_t i = 0; i < points; ++i)
      for (std::size_t v = 0; v < NVAR; ++v)
        expected[i][v] = std::sin(0.31 * static_cast<double>(4 * i + v)) + 0.2;
    for (std::size_t k = 0; k < 2; ++k)
      for (std::size_t v = 0; v < NVAR; ++v) {
        inner.f[k] -= inner.dfdy[k][v] * expected.front()[v];
        outer.f[k] -= outer.dfdy[k][v] * expected.back()[v];
      }
    for (std::size_t i = 0; i < zones.size(); ++i)
      for (std::size_t k = 0; k < NVAR; ++k)
        for (std::size_t v = 0; v < NVAR; ++v)
          zones[i].f[k] -= zones[i].dfdy_lo[k][v] * expected[i][v]
                           + zones[i].dfdy_hi[k][v] * expected[i + 1][v];
    const auto solved = solve_henyey(inner, zones, outer);
    double error = 0.0;
    for (std::size_t i = 0; i < points; ++i)
      for (std::size_t v = 0; v < NVAR; ++v)
        error = std::max(error, std::abs(solved.dy[i][v] - expected[i][v]));
    std::printf("       %zu mesh points\n", points);
    check(error < 1e-10, "known solution recovered despite row scales spanning 160 decades", error);
    check(solved.backward_error < 1e-12, "correction satisfies the original linear equations", solved.backward_error);
    if (points <= 31) {
      const std::size_t n = NVAR * points;
      std::vector<std::vector<double>> a(n, std::vector<double>(n + 1));
      for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t v = 0; v < NVAR; ++v) {
          a[k][v] = inner.dfdy[k][v]; a[n - 2 + k][n - NVAR + v] = outer.dfdy[k][v];
        }
        a[k][n] = -inner.f[k]; a[n - 2 + k][n] = -outer.f[k];
      }
      for (std::size_t i = 0; i < zones.size(); ++i)
        for (std::size_t k = 0; k < NVAR; ++k) {
          const auto row = 2 + NVAR * i + k;
          for (std::size_t v = 0; v < NVAR; ++v) {
            a[row][NVAR * i + v] = zones[i].dfdy_lo[k][v];
            a[row][NVAR * (i + 1) + v] = zones[i].dfdy_hi[k][v];
          }
          a[row][n] = -zones[i].f[k];
        }
      const auto dense = dense_solve(a);
      error = 0.0;
      for (std::size_t i = 0; i < points; ++i)
        for (std::size_t v = 0; v < NVAR; ++v)
          error = std::max(error, std::abs(solved.dy[i][v] - dense[NVAR * i + v]));
      check(error < 1e-11, "block result agrees with independent dense elimination", error);
    }
  }
  bool rejected = false;
  try { (void)solve_henyey({}, {}, {}); } catch (const std::runtime_error&) { rejected = true; }
  check(rejected, "singular systems are rejected");
  BoundaryBlock bad{}; bad.dfdy[0][0] = std::numeric_limits<double>::quiet_NaN();
  rejected = false;
  try { (void)solve_henyey(bad, {}, {}); } catch (const std::domain_error&) { rejected = true; }
  check(rejected, "non-finite coefficients are rejected");
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "ALL PASS", failures);
  return failures ? 1 : 0;
}
