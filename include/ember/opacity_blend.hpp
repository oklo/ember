#pragma once
#include "ember/opacity.hpp"

namespace ember {
// Smooth log-opacity blend across a specified log10(T) interval. Both
// sources must cover every state inside it. Sources must outlive this object.
class BlendedOpacity final : public Opacity {
public:
  BlendedOpacity(const Opacity& low, const Opacity& high,
                 double logT_start = 4.0, double logT_end = 4.5);
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override { return "temperature-blended radiative opacity"; }
private:
  const Opacity& low_;
  const Opacity& high_;
  double start_, end_;
};
} // namespace ember
