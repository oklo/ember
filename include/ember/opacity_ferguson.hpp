#pragma once
#include "ember/opacity.hpp"
#include <filesystem>
#include <vector>

namespace ember {

// Ferguson et al. (2005) low-temperature Rosseland means, tabulated in
// (X, log T, log R) with log R = log rho - 3 (log T - 6).
//
// The table edges are errors, not extrapolation.  Beyond them the caller must
// hand over to another source, because a table asked for a number it does not
// have is exactly how a cool giant lost its Hayashi limit in the Fortran line.
class FergusonOpacity final : public Opacity {
public:
  explicit FergusonOpacity(const std::filesystem::path& file);

  OpacityState eval(double T, double rho, const Composition&) const override;
  const char* name() const override { return "Ferguson+2005 (GS98)"; }

  struct Range { double logT_min, logT_max, logR_min, logR_max, X_min, X_max; };
  Range range() const;
  bool covers(double T, double rho, double X) const;

private:
  std::vector<double> X_, logT_, logR_;
  std::vector<double> k_;   // [ix][it][ir] flattened, log10 kappa
  double Z_{};
  double at(std::size_t ix, std::size_t it, std::size_t ir) const {
    return k_[(ix * logT_.size() + it) * logR_.size() + ir];
  }
};

} // namespace ember
