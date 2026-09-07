#include "ember/henyey.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
using Row = std::array<double, 2 * NVAR + 1>;
using Active = std::array<Row, NVAR + 2>;
using Pivots = std::array<Row, NVAR>;
constexpr std::size_t rhs = 2 * NVAR;

void equilibrate(Row& row) {
  double scale = 0.0;
  for (std::size_t j = 0; j < rhs; ++j) scale = std::max(scale, std::abs(row[j]));
  for (double v : row)
    if (!std::isfinite(v)) throw std::domain_error("solve_henyey: non-finite matrix or residual");
  if (!(scale > 0.0)) throw std::runtime_error("solve_henyey: singular system (empty equation)");
  for (double& v : row) v /= scale;
}

void eliminate(Active& a, std::size_t rows) {
  for (std::size_t k = 0; k < NVAR; ++k) {
    std::size_t pivot = k;
    for (std::size_t i = k + 1; i < rows; ++i)
      if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    const double value = a[pivot][k];
    if (value == 0.0 || !std::isfinite(value))
      throw std::runtime_error("solve_henyey: singular or non-finite elimination pivot");
    std::swap(a[k], a[pivot]);
    for (std::size_t j = k + 1; j <= rhs; ++j) a[k][j] /= value;
    a[k][k] = 1.0;
    for (std::size_t i = k + 1; i < rows; ++i) {
      const double factor = a[i][k];
      a[i][k] = 0.0;
      for (std::size_t j = k + 1; j <= rhs; ++j) a[i][j] -= factor * a[k][j];
    }
  }
}

std::array<double, NVAR> back_substitute(const Pivots& a, const std::array<double, NVAR>& next) {
  std::array<double, NVAR> x{};
  for (std::size_t k = NVAR; k-- > 0;) {
    double value = a[k][rhs];
    for (std::size_t j = k + 1; j < NVAR; ++j) value -= a[k][j] * x[j];
    for (std::size_t j = 0; j < NVAR; ++j) value -= a[k][NVAR + j] * next[j];
    if (!std::isfinite(value)) throw std::runtime_error("solve_henyey: non-finite correction");
    x[k] = value;
  }
  return x;
}
} // namespace

HenyeyCorrection solve_henyey(const BoundaryBlock& inner, const std::vector<ZoneResidual>& zones,
                              const BoundaryBlock& outer) {
  std::array<Row, 2> carried{};
  for (std::size_t k = 0; k < 2; ++k) {
    std::copy(inner.dfdy[k].begin(), inner.dfdy[k].end(), carried[k].begin());
    carried[k][rhs] = -inner.f[k];
  }
  std::vector<Pivots> saved;
  saved.reserve(zones.size());
  for (const auto& zone : zones) {
    Active active{};
    active[0] = carried[0]; active[1] = carried[1];
    for (std::size_t k = 0; k < NVAR; ++k) {
      auto& row = active[k + 2];
      std::copy(zone.dfdy_lo[k].begin(), zone.dfdy_lo[k].end(), row.begin());
      std::copy(zone.dfdy_hi[k].begin(), zone.dfdy_hi[k].end(), row.begin() + NVAR);
      row[rhs] = -zone.f[k];
    }
    for (auto& row : active) equilibrate(row);
    eliminate(active, active.size());
    Pivots pivots{};
    std::copy_n(active.begin(), NVAR, pivots.begin());
    saved.push_back(pivots);
    for (std::size_t k = 0; k < 2; ++k) {
      carried[k] = {};
      std::copy_n(active[NVAR + k].begin() + NVAR, NVAR, carried[k].begin());
      carried[k][rhs] = active[NVAR + k][rhs];
    }
  }
  Active last{};
  last[0] = carried[0]; last[1] = carried[1];
  for (std::size_t k = 0; k < 2; ++k) {
    std::copy(outer.dfdy[k].begin(), outer.dfdy[k].end(), last[k + 2].begin());
    last[k + 2][rhs] = -outer.f[k];
  }
  for (std::size_t k = 0; k < NVAR; ++k) equilibrate(last[k]);
  eliminate(last, NVAR);
  Pivots final_pivots{};
  std::copy_n(last.begin(), NVAR, final_pivots.begin());
  HenyeyCorrection result{};
  result.dy.resize(zones.size() + 1);
  result.dy.back() = back_substitute(final_pivots, {});
  for (std::size_t i = zones.size(); i-- > 0;)
    result.dy[i] = back_substitute(saved[i], result.dy[i + 1]);

  // Check the original equations, not the eliminated ones. Long double
  // provides wider accumulation on platforms where it exceeds double.
  auto check = [&](double f, const auto& lo, const auto& hi, std::size_t i, std::size_t j) {
    long double residual = f, norm = std::abs(f);
    for (std::size_t v = 0; v < NVAR; ++v) {
      const long double a = static_cast<long double>(lo[v]) * result.dy[i][v];
      const long double b = static_cast<long double>(hi[v]) * result.dy[j][v];
      residual += a + b; norm += std::abs(a) + std::abs(b);
    }
    const double error = norm > 0.0L ? static_cast<double>(std::abs(residual) / norm) : 0.0;
    if (!std::isfinite(error)) throw std::runtime_error("solve_henyey: non-finite linear residual");
    result.backward_error = std::max(result.backward_error, error);
  };
  const std::array<double, NVAR> zero{};
  for (std::size_t k = 0; k < 2; ++k) {
    check(inner.f[k], inner.dfdy[k], zero, 0, 0);
    check(outer.f[k], outer.dfdy[k], zero, zones.size(), zones.size());
  }
  for (std::size_t i = 0; i < zones.size(); ++i)
    for (std::size_t k = 0; k < NVAR; ++k)
      check(zones[i].f[k], zones[i].dfdy_lo[k], zones[i].dfdy_hi[k], i, i + 1);
  if (result.backward_error > 1e-10)
    throw std::runtime_error("solve_henyey: correction fails the original linear equations");
  return result;
}

} // namespace ember
