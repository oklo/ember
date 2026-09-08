#pragma once
#include "differential.hpp"

namespace ember::detail {
// Two independent coordinates and their symmetric Hessian. Used to
// differentiate the actual EOS interpolant, including active slope limiters.
struct Jet2 {
  double value{};
  std::array<double, 2> d{};
  std::array<std::array<double, 2>, 2> h{};
  Jet2() = default;
  Jet2(double v) : value(v) {}
  static Jet2 variable(double v, std::size_t i) { Jet2 a(v); a.d[i] = 1; return a; }
  Differential<2> first() const { Differential<2> a(value); a.d = d; return a; }
  Differential<2> partial(std::size_t i) const { Differential<2> a(d[i]); a.d = h[i]; return a; }
  friend Jet2 operator+(const Jet2& a, const Jet2& b) {
    Jet2 c(a.value + b.value);
    for (int i = 0; i < 2; ++i) {
      c.d[i] = a.d[i] + b.d[i];
      for (int j = 0; j < 2; ++j) c.h[i][j] = a.h[i][j] + b.h[i][j];
    }
    return c;
  }
  friend Jet2 operator-(const Jet2& a) {
    Jet2 c(-a.value);
    for (int i = 0; i < 2; ++i) {
      c.d[i] = -a.d[i];
      for (int j = 0; j < 2; ++j) c.h[i][j] = -a.h[i][j];
    }
    return c;
  }
  friend Jet2 operator-(const Jet2& a, const Jet2& b) { return a + (-b); }
  friend Jet2 operator*(const Jet2& a, const Jet2& b) {
    Jet2 c(a.value * b.value);
    for (int i = 0; i < 2; ++i) {
      c.d[i] = a.d[i] * b.value + a.value * b.d[i];
      for (int j = 0; j < 2; ++j)
        c.h[i][j] = a.h[i][j] * b.value + a.d[i] * b.d[j]
                  + a.d[j] * b.d[i] + a.value * b.h[i][j];
    }
    return c;
  }
};
inline Jet2 unary(const Jet2& a, double value, double first, double second) {
  Jet2 b(value);
  for (int i = 0; i < 2; ++i) {
    b.d[i] = first * a.d[i];
    for (int j = 0; j < 2; ++j) b.h[i][j] = first * a.h[i][j] + second * a.d[i] * a.d[j];
  }
  return b;
}
inline Jet2 operator/(const Jet2& a, const Jet2& b) {
  const double v = 1 / b.value;
  return a * unary(b, v, -v * v, 2 * v * v * v);
}
inline Jet2 exp(const Jet2& a) { const double v = std::exp(a.value); return unary(a, v, v, v); }
inline Jet2 log(const Jet2& a) { return unary(a, std::log(a.value), 1 / a.value, -1 / (a.value * a.value)); }
} // namespace ember::detail
