#pragma once
#include "ember/composition.hpp"
#include <optional>

namespace ember {

struct OpacityState {
  double kappa{};      // Rosseland mean, radiative + conductive  [cm^2/g]
  double dlnk_dlnT{};
  double dlnk_dlnRho{};
};

class Opacity {
public:
  struct DensityRange { double min, max; };  // at a specified T and composition
  virtual ~Opacity() = default;
  virtual OpacityState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  // Optional domain information for constrained solves (e.g. an atmosphere's
  // top pressure). Absence means no declared bounds, not permission to ignore
  // errors from eval(). A table implementation should expose its actual range.
  virtual std::optional<DensityRange> density_range(double, const Composition&) const {
    return std::nullopt;
  }
};

} // namespace ember
