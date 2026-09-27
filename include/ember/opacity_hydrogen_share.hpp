#pragma once
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_metal_lower_continuation.hpp"
#include "ember/interp.hpp"
#include <fstream>
#include <iomanip>

namespace ember {

namespace hydrogen_share {
inline void check(const Composition& c) {
  if (c.basis != AbundanceBasis::atomic_mass || c.X[1] != 0 ||
      std::abs(c.sum() - 1) > 1e-10 || !std::isfinite(c.Z()) || c.Z() >= 1)
    throw std::domain_error("hydrogen-share opacity: normalized elemental composition required");
  for (double x : c.X)
    if (!std::isfinite(x) || x < 0)
      throw std::domain_error("hydrogen-share opacity: invalid abundance");
  if (c.h1() / (1 - c.Z()) > 1 + 2e-14)
    throw std::domain_error("hydrogen-share opacity: hydrogen exceeds available mass");
}

// Change Z at fixed H/(H+He). Every comparison retains nonnegative helium.
// This changes the composition interpolation approximation, not source data.
inline Composition at(const Composition& c, double z) {
  auto out = c;
  const double old = c.Z(), u = std::clamp(c.h1() / (1 - old), 0., 1.);
  if (old == 0) {
    out = solar_scaled(u * (1 - z), z);
    out.basis = c.basis;
    out.metal_inventory = c.metal_inventory;
  } else {
    for (std::size_t j = 3; j < METAL_END; ++j) out.X[j] *= z / old;
    out.X[0] = u * (1 - z);
    out.X[2] = (1 - u) * (1 - z);
  }
  return out;
}
} // namespace hydrogen_share

// Candidate interpolation of actual fixed-Z tables at common H/(1-Z).
// The ordinary MixtureOpacity and its existing manifest interpretation remain
// unchanged. Selecting this class requires separate physical comparisons.
class HydrogenShareMixtureOpacity final : public Opacity {
public:
  explicit HydrogenShareMixtureOpacity(const std::filesystem::path& path) {
    std::ifstream in(path);
    std::string magic, axis, label;
    int version;
    std::size_t count;
    in >> magic >> version >> count >> axis >> label;
    if (!in || magic != "EMBER_OPACITY_MIXTURE" || version != 1 ||
        count < 2 || count > 100 || (axis != "logR" && axis != "logRho"))
      throw std::runtime_error("HydrogenShareMixtureOpacity: invalid manifest");
    for (std::size_t i = 0; i < count; ++i) {
      double z;
      std::string file;
      in >> z >> std::quoted(file);
      if (!in || !std::isfinite(z) || z < 0 || z >= 1 ||
          (i && z <= z_.back()) || file.empty())
        throw std::runtime_error("HydrogenShareMixtureOpacity: invalid metal axis");
      auto table = std::make_unique<TabulatedOpacity>(path.parent_path() / file,
          label.c_str(), axis == "logR" ? TabulatedOpacity::DensityAxis::logR
                                        : TabulatedOpacity::DensityAxis::logRho);
      if (std::abs(table->metallicity() - z) > 1e-12)
        throw std::runtime_error("HydrogenShareMixtureOpacity: source Z mismatch");
      z_.push_back(z);
      tables_.push_back(std::move(table));
    }
    if (in >> magic)
      throw std::runtime_error("HydrogenShareMixtureOpacity: trailing manifest data");
  }

  OpacityState eval(double T, double rho, const Composition& c) const override {
    const auto [i, f] = interval(c);
    const auto ca = hydrogen_share::at(c, z_[i]), cb = hydrogen_share::at(c, z_[i+1]);
    const auto a = tables_[i]->eval(T, rho, ca), b = tables_[i+1]->eval(T, rho, cb);
    const double la = std::log(a.kappa), lb = std::log(b.kappa), den = 1 - c.Z();
    return {std::exp((1-f)*la + f*lb),
            (1-f)*a.dlnk_dlnT + f*b.dlnk_dlnT,
            (1-f)*a.dlnk_dlnRho + f*b.dlnk_dlnRho,
            ((1-f)*a.dlnk_dX*(1-z_[i]) + f*b.dlnk_dX*(1-z_[i+1]))/den,
            (lb-la)/(z_[i+1]-z_[i]) +
              ((1-f)*a.dlnk_dX*ca.h1() + f*b.dlnk_dX*cb.h1())/den};
  }

  std::optional<DensityRange> density_range(double T, const Composition& c) const override {
    const auto [i, f] = interval(c);
    (void)f;
    const auto a = tables_[i]->density_range(T, hydrogen_share::at(c, z_[i]));
    const auto b = tables_[i+1]->density_range(T, hydrogen_share::at(c, z_[i+1]));
    if (!a || !b) throw std::logic_error("HydrogenShareMixtureOpacity: missing density bounds");
    DensityRange r{std::max(a->min,b->min), std::min(a->max,b->max)};
    if (!(r.min < r.max)) throw std::domain_error("HydrogenShareMixtureOpacity: no density overlap");
    return r;
  }
  const char* name() const override { return "source opacity at fixed hydrogen share"; }
private:
  std::vector<double> z_;
  std::vector<std::unique_ptr<TabulatedOpacity>> tables_;
  std::pair<std::size_t,double> interval(const Composition& c) const {
    hydrogen_share::check(c);
    if (c.Z() < z_.front()-2e-14 || c.Z() > z_.back()+2e-14)
      throw std::domain_error("HydrogenShareMixtureOpacity: Z outside source family");
    const double z = std::clamp(c.Z(),z_.front(),z_.back());
    const auto i = interp::locate(z_,z);
    return {i,(z-z_[i])/(z_[i+1]-z_[i])};
  }
};

// Bounded lower-Z continuation at fixed H/(1-Z). The source family must use
// physical composition coordinates too; this wrapper cannot repair missing
// temperature, density or source-composition coverage.
class LowerMetalHydrogenShareOpacity final : public Opacity {
public:
  using Method = LowerMetalOpacityContinuation::Method;
  LowerMetalHydrogenShareOpacity(const Opacity& source, double anchor,
      double minimum, Method method, double step = .01)
      : source_(source), anchor_(anchor), minimum_(minimum), step_(step), method_(method) {
    if (!std::isfinite(anchor+minimum+step) || minimum <= 0 || anchor <= minimum ||
        step <= 0 || anchor+step >= 1 || source.includes_conduction())
      throw std::invalid_argument("LowerMetalHydrogenShareOpacity: invalid radiative interval");
  }
  OpacityState eval(double T, double rho, const Composition& c) const override {
    check(c);
    if (c.Z() >= anchor_-2e-14) return source_.eval(T,rho,c);
    const auto ca = hydrogen_share::at(c,anchor_), cb = hydrogen_share::at(c,anchor_+step_);
    const auto a = source_.eval(T,rho,ca), b = source_.eval(T,rho,cb);
    const double f = (anchor_-c.Z())/step_, den = 1-c.Z();
    const double ax = a.dlnk_dX*(1-anchor_)/den, bx = b.dlnk_dX*(1-anchor_-step_)/den;
    const double az = a.dlnk_dX*ca.h1()/den, bz = b.dlnk_dX*cb.h1()/den;
    auto out = a;
    if (method_ == Method::linear_kappa) {
      const double wa = (1+f)*a.kappa, wb = -f*b.kappa, k = wa+wb;
      if (!(k > 0) || !std::isfinite(k))
        throw std::domain_error("LowerMetalHydrogenShareOpacity: nonpositive linear opacity");
      out.kappa = k;
      out.dlnk_dlnT = (wa*a.dlnk_dlnT + wb*b.dlnk_dlnT)/k;
      out.dlnk_dlnRho = (wa*a.dlnk_dlnRho + wb*b.dlnk_dlnRho)/k;
      out.dlnk_dX = (wa*ax + wb*bx)/k;
      out.dlnk_dZ = ((b.kappa-a.kappa)/step_ + wa*az + wb*bz)/k;
    } else {
      const double response = std::log(b.kappa/a.kappa);
      out.kappa = a.kappa*std::exp(-f*response);
      out.dlnk_dlnT = (1+f)*a.dlnk_dlnT - f*b.dlnk_dlnT;
      out.dlnk_dlnRho = (1+f)*a.dlnk_dlnRho - f*b.dlnk_dlnRho;
      out.dlnk_dX = (1+f)*ax - f*bx;
      out.dlnk_dZ = response/step_ + (1+f)*az - f*bz;
    }
    out.dlnk_dY3 = 0.;
    if (!(out.kappa > 0) || !std::isfinite(out.kappa) ||
        !std::isfinite(out.dlnk_dlnT) || !std::isfinite(out.dlnk_dlnRho) ||
        !std::isfinite(out.dlnk_dX) || !std::isfinite(out.dlnk_dZ))
      throw std::domain_error("LowerMetalHydrogenShareOpacity: invalid response");
    return out;
  }
  std::optional<DensityRange> density_range(double T, const Composition& c) const override {
    check(c);
    if (c.Z() >= anchor_-2e-14) return source_.density_range(T,c);
    const auto a = source_.density_range(T,hydrogen_share::at(c,anchor_));
    const auto b = source_.density_range(T,hydrogen_share::at(c,anchor_+step_));
    if (!a || !b) throw std::domain_error("LowerMetalHydrogenShareOpacity: missing density bounds");
    DensityRange r{std::max(a->min,b->min),std::min(a->max,b->max)};
    if (!(r.min < r.max)) throw std::domain_error("LowerMetalHydrogenShareOpacity: no density overlap");
    return r;
  }
  const char* name() const override { return "lower-metal opacity at fixed hydrogen share"; }
private:
  const Opacity& source_;
  double anchor_, minimum_, step_;
  Method method_;
  void check(const Composition& c) const {
    hydrogen_share::check(c);
    if (c.Z() < minimum_-2e-14)
      throw std::domain_error("LowerMetalHydrogenShareOpacity: Z outside continuation interval");
  }
};
} // namespace ember
