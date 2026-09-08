#pragma once

namespace ember::detail {

// Shared backward gravitational heating convention. Scalar
// may be double or a chain-rule value; the previous state is held fixed.
template<class Scalar>
Scalar gravitational_heating(const Scalar& E, const Scalar& P, const Scalar& rho,
                             double previous_E, double previous_rho, double dt) {
  return -(E - previous_E + P * (1.0 / rho - 1.0 / previous_rho)) / dt;
}

} // namespace ember::detail
