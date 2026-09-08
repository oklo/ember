#include "ember/opacity_blend.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
BlendedOpacity::BlendedOpacity(const Opacity& low, const Opacity& high, double start, double end)
    : low_(low), high_(high), start_(start), end_(end) {
  if (!std::isfinite(start) || !std::isfinite(end) || !(start < end))
    throw std::invalid_argument("BlendedOpacity: invalid transition temperatures");
}

OpacityState BlendedOpacity::eval(double T, double rho, const Composition& comp) const {
  const double lt = std::log10(T);
  if (!std::isfinite(lt)) throw std::domain_error("BlendedOpacity: invalid temperature");
  if (lt <= start_) return low_.eval(T, rho, comp);
  if (lt >= end_) return high_.eval(T, rho, comp);
  const auto a = low_.eval(T, rho, comp), b = high_.eval(T, rho, comp);
  const double u = (lt - start_) / (end_ - start_);
  const double w = u * u * (3.0 - 2.0 * u);
  const double dw = 6.0 * u * (1.0 - u) / ((end_ - start_) * std::log(10.0));
  const double la = std::log(a.kappa), lb = std::log(b.kappa);
  return {std::exp((1.0 - w) * la + w * lb),
          (1.0 - w) * a.dlnk_dlnT + w * b.dlnk_dlnT + dw * (lb - la),
          (1.0 - w) * a.dlnk_dlnRho + w * b.dlnk_dlnRho};
}

std::optional<Opacity::DensityRange> BlendedOpacity::density_range(double T, const Composition& comp) const {
  const double lt = std::log10(T);
  if (!std::isfinite(lt)) throw std::domain_error("BlendedOpacity: invalid temperature");
  if (lt <= start_) return low_.density_range(T, comp);
  if (lt >= end_) return high_.density_range(T, comp);
  const auto a = low_.density_range(T, comp), b = high_.density_range(T, comp);
  if (!a) return b;
  if (!b) return a;
  const DensityRange r{std::max(a->min, b->min), std::min(a->max, b->max)};
  if (!(r.min < r.max)) throw std::domain_error("BlendedOpacity: no shared density range");
  return r;
}
} // namespace ember
