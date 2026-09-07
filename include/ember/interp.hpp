#pragma once
#include <span>
#include <vector>

namespace ember::interp {

// Monotone cubic Hermite (Fritsch-Carlson) on a possibly non-uniform grid,
// returning the value and its derivative.
//
// Monotone, not merely C1: an opacity table crosses the grain condensation
// edge with neighbouring cells differing by decades, and an unlimited cubic
// overshoots there.  A Newton solver handed an overshoot sees a spurious
// negative opacity gradient and diverges - which is not hypothetical; it is
// how the Fortran ancestor died the first time this table was installed.
struct Value { double y; double dydx; };

Value hermite(std::span<const double> x, std::span<const double> y, double xq);

// Bilinear-in-index location: index of the interval containing xq, clamped.
std::size_t locate(std::span<const double> x, double xq);

} // namespace ember::interp
