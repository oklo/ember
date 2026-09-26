#pragma once
#include "ember/opacity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {

// Continue the hydrogen dependence of a retained opacity using a second
// physical calculation. Below the anchor, return the retained source exactly.
// Above it, k = k_retained(anchor) * k_source(X) / k_source(anchor).
// This is an explicit physical approximation, not table extrapolation.
class HydrogenOpacityExtension final : public Opacity {
public:
  HydrogenOpacityExtension(const Opacity& retained, const Opacity& dependence,
                           double anchor, double maximum)
      : retained_(retained), dependence_(dependence), anchor_(anchor), maximum_(maximum) {
    if (!std::isfinite(anchor + maximum) || anchor < 0 || maximum > 1 || anchor >= maximum)
      throw std::invalid_argument("HydrogenOpacityExtension: invalid composition interval");
    if (retained.includes_conduction() || dependence.includes_conduction())
      throw std::invalid_argument("HydrogenOpacityExtension: radiative sources required");
  }

  OpacityState eval(double T, double rho, const Composition& c) const override {
    check(c);
    if (c.h1() <= anchor_) return retained_.eval(T, rho, c);
    const auto at_anchor = anchored(c);
    const auto a = retained_.eval(T, rho, at_anchor);
    const auto b = dependence_.eval(T, rho, c);
    const auto d = dependence_.eval(T, rho, at_anchor);
    OpacityState result = a;
    result.kappa = a.kappa * (b.kappa / d.kappa);
    result.dlnk_dlnT = a.dlnk_dlnT + b.dlnk_dlnT - d.dlnk_dlnT;
    result.dlnk_dlnRho = a.dlnk_dlnRho + b.dlnk_dlnRho - d.dlnk_dlnRho;
    result.dlnk_dX = b.dlnk_dX;
    result.dlnk_dZ = a.dlnk_dZ + b.dlnk_dZ - d.dlnk_dZ;
    result.dlnk_dY3 = 0.; // elemental source; isotope mapping belongs outside
    if (!std::isfinite(result.kappa) || !(result.kappa > 0))
      throw std::domain_error("HydrogenOpacityExtension: invalid source ratio");
    return result;
  }

  std::optional<DensityRange> density_range(double T, const Composition& c) const override {
    check(c);
    if (c.h1() <= anchor_) return retained_.density_range(T, c);
    const auto at_anchor = anchored(c);
    const auto a = retained_.density_range(T, at_anchor);
    const auto b = dependence_.density_range(T, c);
    const auto d = dependence_.density_range(T, at_anchor);
    if (!a || !b || !d)
      throw std::domain_error("HydrogenOpacityExtension: missing source density range");
    DensityRange range{std::max({a->min, b->min, d->min}),
                       std::min({a->max, b->max, d->max})};
    if (!(range.min < range.max))
      throw std::domain_error("HydrogenOpacityExtension: no common source density interval");
    return range;
  }

  const char* name() const override { return "opacity with source-based hydrogen dependence"; }

private:
  const Opacity& retained_;
  const Opacity& dependence_;
  double anchor_, maximum_;
  void check(const Composition& c) const {
    if (c.basis != AbundanceBasis::atomic_mass || c.X[1] != 0 ||
        !std::isfinite(c.h1()) || c.h1() < 0 || c.h1() > maximum_)
      throw std::domain_error("HydrogenOpacityExtension: elemental composition outside range");
  }
  Composition anchored(const Composition& c) const {
    auto out = c;
    out.X[2] += out.X[0] - anchor_;
    out.X[0] = anchor_;
    return out;
  }
};

} // namespace ember
