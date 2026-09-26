#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember::interp {

std::size_t locate(std::span<const double> x, double xq) {
  if (x.size() < 2) return 0;
  const auto it = std::upper_bound(x.begin(), x.end(), xq);
  const std::size_t i = static_cast<std::size_t>(it - x.begin());
  if (i == 0) return 0;
  if (i >= x.size()) return x.size() - 2;
  return i - 1;
}

namespace {
// Fritsch-Carlson limiter: the one-sided slope that keeps the cubic monotone.
struct Slope { double value, dp; };
Slope limited(Slope dl, Slope dc) {
  if (dl.value * dc.value <= 0.0) return {0.0, 0.0};
  const Slope s{0.5 * (dl.value + dc.value), 0.5 * (dl.dp + dc.dp)};
  const Slope small = std::abs(dl.value) <= std::abs(dc.value) ? dl : dc;
  return std::abs(s.value) > 3.0 * std::abs(small.value)
      ? Slope{3.0 * small.value, 3.0 * small.dp} : s;
}
// Density interpolation usually needs no independent parameter response.
// Compile that response out while retaining the same value/slope arithmetic.
template<bool parameter_response>
auto hermite_impl(std::span<const double> x, std::span<const double> y,
                  std::span<const double> dy_dp, double xq) {
  const std::size_t n = x.size();
  if (n == 0 || y.size() != n || (!dy_dp.empty() && dy_dp.size() != n))
    throw std::invalid_argument("hermite: incompatible or empty input arrays");
  auto dp = [&](std::size_t i) { return dy_dp.empty() ? 0.0 : dy_dp[i]; };
  if (n == 1) {
    if constexpr (parameter_response) return ParametricValue{y[0], 0.0, dp(0)};
    else return Value{y[0], 0.0};
  }
  const std::size_t k = locate(x, xq);
  const double h  = x[k + 1] - x[k];
  auto secant = [&](std::size_t j) {
    const double width = x[j + 1] - x[j];
    if constexpr (parameter_response)
      return Slope{(y[j + 1] - y[j]) / width, (dp(j + 1) - dp(j)) / width};
    else return Slope{(y[j + 1] - y[j]) / width, 0.0};
  };
  const auto dc = secant(k);
  const auto s1 = k > 0 ? limited(secant(k - 1), dc) : dc;
  const auto s2 = k + 2 < n ? limited(dc, secant(k + 1)) : dc;
  const double t  = (xq - x[k]) / h;
  const double t2 = t * t, t3 = t2 * t;
  const double h00 =  2 * t3 - 3 * t2 + 1, h10 =      t3 - 2 * t2 + t;
  const double h01 = -2 * t3 + 3 * t2,     h11 =      t3 -     t2;
  const double d00 =  6 * t2 - 6 * t,      d10 =  3 * t2 - 4 * t + 1;
  const double d01 = -6 * t2 + 6 * t,      d11 =  3 * t2 - 2 * t;
  const Value value{
           h00 * y[k] + h10 * h * s1.value + h01 * y[k + 1] + h11 * h * s2.value,
          (d00 * y[k] + d01 * y[k + 1]) / h + d10 * s1.value + d11 * s2.value};
  if constexpr (parameter_response)
    return ParametricValue{value.y, value.dydx,
          h00 * dp(k) + h10 * h * s1.dp + h01 * dp(k + 1) + h11 * h * s2.dp};
  else return value;
}
} // namespace

Value hermite(std::span<const double> x, std::span<const double> y, double xq) {
  return hermite_impl<false>(x, y, {}, xq);
}

ParametricValue hermite(std::span<const double> x, std::span<const double> y,
                       std::span<const double> dy_dp, double xq) {
  return hermite_impl<true>(x, y, dy_dp, xq);
}

} // namespace ember::interp
