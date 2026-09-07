#pragma once
#include <array>
#include <cmath>
#include <cstddef>

namespace ember::detail {

// First-order chain rule for the small, fixed number of zone variables.
// Physics values and their analytic partials are seeded at module boundaries.
// N=0 evaluates the same equations without derivative work for the numerical
// reference. No physics call is differentiated by finite differences here.
template<std::size_t N> struct Differential {
  double value{};
  std::array<double, N> d{};
  Differential() = default;
  Differential(double v) : value(v) {}
  static Differential variable(double v, std::size_t index) {
    Differential out(v);
    if constexpr (N > 0) out.d[index] = 1.0;
    return out;
  }
  friend Differential operator+(const Differential& a, const Differential& b) {
    Differential out(a.value + b.value);
    for (std::size_t i = 0; i < N; ++i) out.d[i] = a.d[i] + b.d[i];
    return out;
  }
  friend Differential operator-(const Differential& a, const Differential& b) {
    Differential out(a.value - b.value);
    for (std::size_t i = 0; i < N; ++i) out.d[i] = a.d[i] - b.d[i];
    return out;
  }
  friend Differential operator-(const Differential& a) { return Differential(0.0) - a; }
  friend Differential operator*(const Differential& a, const Differential& b) {
    Differential out(a.value * b.value);
    for (std::size_t i = 0; i < N; ++i) out.d[i] = a.d[i] * b.value + a.value * b.d[i];
    return out;
  }
  friend Differential operator/(const Differential& a, const Differential& b) {
    Differential out(a.value / b.value);
    for (std::size_t i = 0; i < N; ++i) out.d[i] = (a.d[i] - out.value * b.d[i]) / b.value;
    return out;
  }
};
template<std::size_t N> Differential<N> log(const Differential<N>& x) {
  Differential<N> out(std::log(x.value));
  for (std::size_t i = 0; i < N; ++i) out.d[i] = x.d[i] / x.value;
  return out;
}
template<std::size_t N> Differential<N> exp(const Differential<N>& x) {
  Differential<N> out(std::exp(x.value));
  for (std::size_t i = 0; i < N; ++i) out.d[i] = out.value * x.d[i];
  return out;
}

} // namespace ember::detail
