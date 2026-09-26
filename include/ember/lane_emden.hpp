#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ember::example {
struct LaneEmden {
  double xi{}, theta{}, mass{};  // mass = -xi^2 dtheta/dxi
};

inline std::vector<LaneEmden> lane_emden(const std::vector<double>& targets, double index = 3.0) {
  if (!std::isfinite(index) || index <= 0) throw std::invalid_argument("Lane-Emden index must be positive");
  // Independent IVP: dtheta/dxi=-mass/xi^2, dmass/dxi=xi^2 theta^index.
  // A regular series starts away from the coordinate singularity. RK4 uses
  // a fixed maximum step much smaller than the stellar relaxation mesh.
  double xi = 1e-6;
  std::array<double, 2> y{1.0 - xi * xi / 6.0, xi * xi * xi / 3.0};
  std::vector<LaneEmden> out;
  for (double target : targets) {
    if (target < xi || !std::isfinite(target)) throw std::invalid_argument("Lane-Emden targets must increase");
    auto derivative = [index](double x, const std::array<double, 2>& state) {
      if (state[0] <= 0) throw std::domain_error("Lane-Emden step crosses the polytrope surface");
      const double density = index == 3 ? state[0] * state[0] * state[0] : std::pow(state[0], index);
      return std::array{-state[1] / (x * x), x * x * density};
    };
    auto add = [](const auto& a, const auto& b, double scale) {
      return std::array{a[0] + scale * b[0], a[1] + scale * b[1]};
    };
    while (xi < target) {
      const double h = std::min({2e-4, 0.05 * xi, target - xi});
      const auto k1 = derivative(xi, y);
      const auto k2 = derivative(xi + h / 2.0, add(y, k1, h / 2.0));
      const auto k3 = derivative(xi + h / 2.0, add(y, k2, h / 2.0));
      const auto k4 = derivative(xi + h, add(y, k3, h));
      for (std::size_t k = 0; k < 2; ++k) y[k] += h * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6.0;
      xi += h;
    }
    if (!(y[0] > 0.0)) throw std::domain_error("Lane-Emden point is outside the positive polytrope");
    out.push_back({xi, y[0], y[1]});
  }
  return out;
}

} // namespace ember::example
