#pragma once
#include "ember/composition.hpp"

namespace ember {

struct OpacityState {
  double kappa{};      // Rosseland mean, radiative + conductive  [cm^2/g]
  double dlnk_dlnT{};
  double dlnk_dlnRho{};
};

class Opacity {
public:
  virtual ~Opacity() = default;
  virtual OpacityState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
};

} // namespace ember
