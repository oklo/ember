#include "ember/henyey.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

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
void prescribed_boundary_values(const BoundaryBlock& boundary,std::array<double,NVAR>& x) {
  std::array<bool,NVAR> known{};
  for(std::size_t pass=0;pass<2;++pass)for(std::size_t row=0;row<2;++row) {
    std::size_t count=0,column=0;
    long double residual=boundary.f[row];
    for(std::size_t v=0;v<NVAR;++v)if(boundary.dfdy[row][v]!=0) {
      if(known[v])residual+=static_cast<long double>(boundary.dfdy[row][v])*x[v];
      else {++count;column=v;}
    }
    if(count==1) {
      x[column]=static_cast<double>(-residual/boundary.dfdy[row][column]);known[column]=true;
    }
  }
}
HenyeyCorrection eliminate_system(const BoundaryBlock& inner, const std::vector<ZoneResidual>& zones,
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
  // A prescribed boundary variable has an explicit scalar equation. Recover
  // its value from that equation before propagating inward. Otherwise tiny
  // elimination roundoff in an exactly zero correction has relative
  // backward error one, regardless of its absolute magnitude. The original
  // full-system check and its tolerance remain unchanged.
  prescribed_boundary_values(outer,result.dy.back());
  for (std::size_t i = zones.size(); i-- > 0;)
    result.dy[i] = back_substitute(saved[i], result.dy[i + 1]);
  prescribed_boundary_values(inner,result.dy.front());
  return result;
}
} // namespace

HenyeyCorrection solve_henyey(const BoundaryBlock& inner, const std::vector<ZoneResidual>& zones,
                              const BoundaryBlock& outer) {
  auto result = eliminate_system(inner, zones, outer);
  BoundaryBlock residual_inner = inner, residual_outer = outer;
  auto residual_zones = zones;
  const char* worst_block="";std::size_t worst_point=0,worst_equation=0;
  double worst_residual=0,worst_norm=0;
  for (unsigned refinement = 0; ; ++refinement) {
    result.backward_error = 0.0;
    // Check the original equations, not the eliminated ones. Long double
    // provides wider accumulation on platforms where it exceeds double.
    auto check = [&](double f, const auto& lo, const auto& hi, std::size_t i, std::size_t j,
                     const char* block,std::size_t equation) {
      long double residual = f, norm = std::abs(f);
      for (std::size_t v = 0; v < NVAR; ++v) {
        const long double a = static_cast<long double>(lo[v]) * result.dy[i][v];
        const long double b = static_cast<long double>(hi[v]) * result.dy[j][v];
        residual += a + b; norm += std::abs(a) + std::abs(b);
      }
      const double error = norm > 0.0L ? static_cast<double>(std::abs(residual) / norm) : 0.0;
      if (!std::isfinite(error)) throw std::runtime_error("solve_henyey: non-finite linear residual");
      if(error>result.backward_error) {
        result.backward_error=error;worst_block=block;worst_point=i;worst_equation=equation;
        worst_residual=static_cast<double>(residual);worst_norm=static_cast<double>(norm);
      }
      return static_cast<double>(residual);
    };
    const std::array<double, NVAR> zero{};
    for (std::size_t k = 0; k < 2; ++k) {
      residual_inner.f[k] = check(inner.f[k], inner.dfdy[k], zero, 0, 0,"inner",k);
      residual_outer.f[k] = check(outer.f[k], outer.dfdy[k], zero, zones.size(), zones.size(),"outer",k);
    }
    for (std::size_t i = 0; i < zones.size(); ++i)
      for (std::size_t k = 0; k < NVAR; ++k)
        residual_zones[i].f[k] = check(zones[i].f[k], zones[i].dfdy_lo[k], zones[i].dfdy_hi[k], i, i + 1,"zone",k);
    if (result.backward_error <= 1e-10) return result;
    if (refinement == 3) {
      std::ostringstream message;
      message << std::setprecision(4)<<"solve_henyey: correction fails the original linear equations (backward error="
              << result.backward_error << ", "<<worst_block<<" point "<<worst_point<<" equation "<<worst_equation
              <<", residual="<<worst_residual<<", norm="<<worst_norm<<')';
      throw std::runtime_error(message.str());
    }
    // Iterative refinement solves A*delta=-(f+A*dy) with the original
    // coefficients. Keep the acceptance threshold; improve the correction.
    const auto delta = eliminate_system(residual_inner, residual_zones, residual_outer);
    for (std::size_t i = 0; i < result.dy.size(); ++i)
      for (std::size_t v = 0; v < NVAR; ++v) result.dy[i][v] += delta.dy[i][v];
  }
}

} // namespace ember
