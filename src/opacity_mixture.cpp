#include "ember/opacity_mixture.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace ember {
MixtureOpacity::MixtureOpacity(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::string magic, axis;
  int version;
  std::size_t count;
  in >> magic >> version >> count >> axis >> name_;
  if (!in || magic != "EMBER_OPACITY_MIXTURE" || version != 1 || count < 2 || count > 100 ||
      (axis != "logR" && axis != "logRho"))
    throw std::runtime_error("MixtureOpacity: invalid manifest");
  for (std::size_t i = 0; i < count; ++i) {
    double z;
    std::string file;
    in >> z >> std::quoted(file);
    if (!in || !std::isfinite(z) || z <= 0 || z >= 1 || (i && z <= z_.back()) || file.empty())
      throw std::runtime_error("MixtureOpacity: invalid metallicity axis");
    auto table = std::make_unique<TabulatedOpacity>(path.parent_path() / file, name_.c_str(),
                                                    axis == "logR" ? TabulatedOpacity::DensityAxis::logR
                                                                   : TabulatedOpacity::DensityAxis::logRho);
    // The reader checks every X plane's common Z; query the recorded range.
    if (std::abs(table->metallicity() - z) > 1e-12)
      throw std::runtime_error("MixtureOpacity: source Z mismatch");
    z_.push_back(z);
    tables_.push_back(std::move(table));
  }
  if (in >> magic)
    throw std::runtime_error("MixtureOpacity: trailing manifest data");
}
std::pair<std::size_t, double> MixtureOpacity::interval(const Composition& c) const {
  double z = c.Z();
  if (std::abs(c.sum() - 1) > 1e-10 || c.X[1] != 0)
    throw std::domain_error("MixtureOpacity: normalized elemental abundances required");
  for (double v : c.X)
    if (!std::isfinite(v) || v < 0)
      throw std::domain_error("MixtureOpacity: invalid abundance");
  constexpr double roundoff = 2e-14;
  if (c.basis != AbundanceBasis::atomic_mass || !std::isfinite(z) || z < z_.front() - roundoff ||
      z > z_.back() + roundoff)
    throw std::domain_error("MixtureOpacity: atomic metallicity outside source family");
  z = std::clamp(z, z_.front(), z_.back());
  const auto i = interp::locate(z_, z);
  return {i, (z - z_[i]) / (z_[i + 1] - z_[i])};
}
Composition MixtureOpacity::source(const Composition& c, std::size_t i) const {
  auto out = c;
  const double z = c.Z();
  for (std::size_t j = 3; j < NSPEC; ++j)
    out.X[j] *= z_[i] / z;
  out.X[2] += z - z_[i];
  return out;
}
OpacityState MixtureOpacity::eval(double T, double rho, const Composition& c) const {
  const auto [i, f] = interval(c);
  const auto a = tables_[i]->eval(T, rho, source(c, i)), b = tables_[i + 1]->eval(T, rho, source(c, i + 1));
  const double la = std::log(a.kappa), lb = std::log(b.kappa);
  return {std::exp((1 - f) * la + f * lb), (1 - f) * a.dlnk_dlnT + f * b.dlnk_dlnT,
          (1 - f) * a.dlnk_dlnRho + f * b.dlnk_dlnRho, (1 - f) * a.dlnk_dX + f * b.dlnk_dX,
          (lb - la) / (z_[i + 1] - z_[i])};
}
std::optional<Opacity::DensityRange> MixtureOpacity::density_range(double T, const Composition& c) const {
  const auto [i, f] = interval(c);
  (void)f;
  const auto a = tables_[i]->density_range(T, source(c, i)),
             b = tables_[i + 1]->density_range(T, source(c, i + 1));
  if (!a || !b)
    throw std::logic_error("MixtureOpacity: missing table bounds");
  DensityRange r{std::max(a->min, b->min), std::min(a->max, b->max)};
  if (r.min >= r.max)
    throw std::domain_error("MixtureOpacity: no source density overlap");
  return r;
}
ElementalOpacity::Mapping ElementalOpacity::map(const Composition& c) {
  if (std::abs(c.sum() - 1) > 1e-10)
    throw std::domain_error("ElementalOpacity: normalized abundances required");
  for (double x : c.X)
    if (!std::isfinite(x) || x < 0)
      throw std::domain_error("ElementalOpacity: invalid abundance");
  Mapping m{};
  m.c.basis = AbundanceBasis::atomic_mass;
  for (std::size_t j = 0; j < NSPEC; ++j)
    m.c.X[j] = nuclides[j].A * c.X[j] / c.abundance_weight(j);
  if(c.metal_inventory==MetalInventory::gs98) {
    const double scale=c.basis==AbundanceBasis::baryon_mass?gs98_atomic_mass_scale():1.;
    // Convert the aggregate GS98 baryonic metal mass to source grams.
    // The source retains its fixed atomic-weight GS98 pattern; individual
    // representative-isotope number ratios are therefore approximate.
    for(std::size_t j=3;j<NSPEC;++j)m.c.X[j]=c.X[j]*scale;
  }
  m.c.X[2] += nuclides[2].A * c.X[1] / c.abundance_weight(1);
  m.c.X[1] = 0;
  m.scale = m.c.sum();
  for (double& x : m.c.X)
    x /= m.scale;
  m.ds = {nuclides[0].A / c.abundance_weight(0) - nuclides[2].A / c.abundance_weight(2),
          nuclides[2].A / c.abundance_weight(1) - nuclides[2].A / c.abundance_weight(2)};
  for (std::size_t j = 0; j < 2; ++j) {
    m.dx[j] = ((j == 0 ? nuclides[0].A / c.abundance_weight(0) : 0) - m.c.X[0] * m.ds[j]) / m.scale;
    m.dz[j] = -m.c.Z() * m.ds[j] / m.scale;
  }
  return m;
}
OpacityState ElementalOpacity::eval(double T, double rho, const Composition& c) const {
  const auto m = map(c);
  auto out = table_.eval(T, rho * m.scale, m.c);
  if (!std::isfinite(out.dlnk_dX) || !std::isfinite(out.dlnk_dZ))
    throw std::domain_error("ElementalOpacity: source X and Z derivatives required");
  const auto derivative = [&](std::size_t j) {
    return (1 + out.dlnk_dlnRho) * m.ds[j] / m.scale + out.dlnk_dX * m.dx[j] + out.dlnk_dZ * m.dz[j];
  };
  const double dx = derivative(0), d3 = derivative(1);
  out.kappa *= m.scale;
  out.dlnk_dX = dx;
  out.dlnk_dY3 = d3;
  // Unavailable: arbitrary metal-pattern replacement in the input basis.
  out.dlnk_dZ = std::numeric_limits<double>::quiet_NaN();
  return out;
}
std::optional<Opacity::DensityRange> ElementalOpacity::density_range(double T, const Composition& c) const {
  const auto m = map(c);
  auto r = table_.density_range(T, m.c);
  if (r) {
    r->min /= m.scale;
    r->max /= m.scale;
  }
  return r;
}
} // namespace ember
