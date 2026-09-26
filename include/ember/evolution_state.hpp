#pragma once
#include "ember/model.hpp"
#include <array>

namespace ember {
struct EvolutionState {
  Model model;
  double next_dt{};
  std::size_t accepted{}, rejected{};
  // Accepted H1, He3 and metal face rates, in g/s. These recover the previous
  // material heat when the finite-convection coefficients are lagged.
  std::vector<std::array<double, 3>> metal_heat_rates{};
};
} // namespace ember
