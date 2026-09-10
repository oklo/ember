#pragma once
#include "ember/opacity.hpp"
#include <filesystem>
#include <vector>
#include <string>

namespace ember {

// Rectangular, fixed-Z Rosseland table in (X, log T, log density coordinate).
// Monotone Hermite in the density coordinate and log T; linear in X.
// A single X plane requires that exact composition. No extrapolation.
class TabulatedOpacity : public Opacity {
public:
  enum class DensityAxis { logR, logRho };
  TabulatedOpacity(const std::filesystem::path& file, std::string label,
                   DensityAxis axis = DensityAxis::logR);

  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override { return label_.c_str(); }

  struct Range { double logT_min, logT_max, logD_min, logD_max, X_min, X_max; DensityAxis density_axis; };
  Range range() const;
  double metallicity() const { return Z_; }
  // Geometric bounds only. eval/density_range also require the table's Z.
  bool covers(double T, double rho, double X) const;

private:
  std::string label_;
  std::vector<double> X_, logT_, logD_;
  std::vector<double> k_;   // [ix][it][ir] flattened, log10 kappa
  double Z_{};
  DensityAxis axis_;
};

} // namespace ember
