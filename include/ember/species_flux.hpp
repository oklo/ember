#pragma once
#include <array>

namespace ember {
using SpeciesVector = std::array<double,2>;
using SpeciesMatrix = std::array<SpeciesVector,2>;
struct SpeciesFaceResponse {
  SpeciesVector rate{}; // outward g/s, H1 and He3
  SpeciesMatrix dleft{}, dright{}; // independent fractions replace He4
};
using CNSpeciesVector=std::array<double,5>; // H1, He3, C12, C13, N14 mass fractions
using CNSpeciesMatrix=std::array<CNSpeciesVector,5>;
struct CNSpeciesFaceResponse {
  CNSpeciesVector rate{}; // outward g/s; physical He4 flux is minus the sum
  CNSpeciesMatrix dleft{},dright{};
};
enum class CNMicroscopicApproximation {
  unselected, zero_catalyst_drift, helium_velocity
};
} // namespace ember
