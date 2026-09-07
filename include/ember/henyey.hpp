#pragma once
#include "ember/structure.hpp"
#include <array>
#include <vector>

namespace ember {

struct BoundaryBlock {
  std::array<double, 2> f{};
  std::array<std::array<double, NVAR>, 2> dfdy{};
};

struct HenyeyCorrection {
  std::vector<std::array<double, NVAR>> dy;
  double backward_error{};  // max |J dy + f| / (|f| + sum |J_ij dy_j|)
};

// Solve J dy = -f for two inner rows, four rows per zone, and two outer rows.
// Inputs should use sensible variable units; individual equations are also
// equilibrated internally. Eliminates four variables at each step with row
// pivoting across the six active equations, passing two constraints outward.
// Storage and work are linear in the number of mesh points. Singular or
// non-finite systems throw; no diagonal regularization hides a failed solve.
HenyeyCorrection solve_henyey(const BoundaryBlock& inner,
                              const std::vector<ZoneResidual>& zones,
                              const BoundaryBlock& outer);

} // namespace ember
