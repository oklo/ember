#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>

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
double limited(double dl, double dc) {
  if (dl * dc <= 0.0) return 0.0;
  const double s = 0.5 * (dl + dc);
  const double cap = 3.0 * std::min(std::abs(dl), std::abs(dc));
  return (std::abs(s) > cap) ? std::copysign(cap, s) : s;
}
} // namespace

Value hermite(std::span<const double> x, std::span<const double> y, double xq) {
  const std::size_t n = x.size();
  if (n == 1) return {y[0], 0.0};
  const std::size_t k = locate(x, xq);
  const double h  = x[k + 1] - x[k];
  const double dc = (y[k + 1] - y[k]) / h;
  const double s1 = (k > 0)       ? limited((y[k] - y[k - 1]) / (x[k] - x[k - 1]), dc) : dc;
  const double s2 = (k + 2 < n)   ? limited(dc, (y[k + 2] - y[k + 1]) / (x[k + 2] - x[k + 1])) : dc;
  const double t  = (xq - x[k]) / h;
  const double t2 = t * t, t3 = t2 * t;
  const double h00 =  2 * t3 - 3 * t2 + 1, h10 =      t3 - 2 * t2 + t;
  const double h01 = -2 * t3 + 3 * t2,     h11 =      t3 -     t2;
  const double d00 =  6 * t2 - 6 * t,      d10 =  3 * t2 - 4 * t + 1;
  const double d01 = -6 * t2 + 6 * t,      d11 =  3 * t2 - 2 * t;
  return { h00 * y[k] + h10 * h * s1 + h01 * y[k + 1] + h11 * h * s2,
          (d00 * y[k] + d01 * y[k + 1]) / h + d10 * s1 + d11 * s2 };
}

} // namespace ember::interp
