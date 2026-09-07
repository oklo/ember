#pragma once

namespace ember::detail {

// Shared backward, left-endpoint gravitational heating convention. Scalar
// may be double or a chain-rule value; the previous state is held fixed.
template<class Scalar>
Scalar gravitational_heating(const Scalar& E, const Scalar& P, const Scalar& rho,
                             double previous_E, double previous_rho, double dt) {
  return -(E - previous_E - P / (rho * rho) * (rho - previous_rho)) / dt;
}

} // namespace ember::detail
