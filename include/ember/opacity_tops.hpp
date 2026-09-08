#pragma once
#include "ember/opacity_table.hpp"
#include "ember/opacity_blend.hpp"

namespace ember {
// Original TOPS/ATOMIC cells at fixed X=.7, Z=.02. Two complete rectangles
// exclude every density that the server reported as clamped. Native log rho
// interpolation; no resampling to the OPAL log R grid.
class TopsOpacity final : public Opacity {
public:
  explicit TopsOpacity(const std::filesystem::path& directory)
      : low_(directory / "tops_gs98_x070_z020_low.dat", "TOPS ATOMIC low rectangle", TabulatedOpacity::DensityAxis::logRho),
        high_(directory / "tops_gs98_x070_z020_high.dat", "TOPS ATOMIC high rectangle", TabulatedOpacity::DensityAxis::logRho),
        blend_(low_, high_, 5.6, 5.7) {}
  OpacityState eval(double T, double rho, const Composition& c) const override { return blend_.eval(T, rho, c); }
  std::optional<DensityRange> density_range(double T, const Composition& c) const override { return blend_.density_range(T, c); }
  const char* name() const override { return "TOPS ATOMIC GS98 X=.7 Z=.02"; }
private:
  TabulatedOpacity low_, high_;
  BlendedOpacity blend_;
};
} // namespace ember
