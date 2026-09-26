#include "ember/atmosphere_metal_response.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
bool positive(double v) { return std::isfinite(v) && v > 0; }
constexpr double ln10 = 2.302585092994045684;
}

MetalResponseAtmosphere::MetalResponseAtmosphere(
    const Eos &eos, const Atmosphere &reference, const std::filesystem::path &path,
    Approximation, Options options)
    : eos_(eos), reference_(reference), options_(options) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("MetalResponseAtmosphere: cannot open " + path.string());
  read(in);
}
MetalResponseAtmosphere::MetalResponseAtmosphere(
    const Eos &eos, const Atmosphere &reference, std::istream &in,
    Approximation, Options options)
    : eos_(eos), reference_(reference), options_(options) { read(in); }

void MetalResponseAtmosphere::read(std::istream &in) {
  const auto label = [&](const char *expected) {
    std::string actual;
    if (!(in >> actual) || actual != expected)
      throw std::runtime_error(std::string("MetalResponseAtmosphere: expected ") + expected);
  };
  label("EMBER_METAL_ATMOSPHERE_RESPONSE");
  int version{};
  if (!(in >> version) || version != 1) throw std::runtime_error("metal atmosphere version");
  label("source"); in >> std::quoted(source_);
  label("approximation"); in >> std::quoted(approximation_);
  label("basis"); label("baryon_mass");
  label("reference_Z"); in >> reference_z_;
  label("source_Z"); in >> source_z_;
  label("maximum_helium3"); in >> maximum_helium3_;
  label("tau"); in >> tau_;
  const double excess = options_.maximum_reference_metal_excess;
  if (!in || source_.empty() || approximation_.empty() || !positive(tau_) ||
      !std::isfinite(source_z_) || source_z_ < 0 || !positive(reference_z_) ||
      source_z_ >= reference_z_ || reference_z_ >= 1 ||
      !std::isfinite(maximum_helium3_) || maximum_helium3_ < 0 || maximum_helium3_ >= 1 ||
      !std::isfinite(excess) || excess < 0 || excess > .01*(reference_z_-source_z_))
    throw std::runtime_error("MetalResponseAtmosphere: invalid source or continuation domain");
  std::size_t count = 1;
  const std::array labels{"hydrogen", "log_teff", "log_g"};
  for (std::size_t k = 0; k < 3; ++k) {
    label(labels[k]);
    std::size_t n{};
    if (!(in >> n) || n < 2 || n > 10000 || count > 1000000/n)
      throw std::runtime_error("MetalResponseAtmosphere: invalid axis size");
    count *= n;
    axes_[k].resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      auto &v = axes_[k][i]; in >> v;
      if (!in || !std::isfinite(v) || (i && v <= axes_[k][i-1]) ||
          (k == 0 ? v < 0 || v >= 1 : !positive(std::pow(10., v))))
        throw std::runtime_error("MetalResponseAtmosphere: invalid axis");
    }
  }
  if (axes_[0].back()+maximum_helium3_+reference_z_+excess >= 1)
    throw std::runtime_error("MetalResponseAtmosphere: reference mixture has no helium-4");
  label("data");
  valid_.resize(count); delta_logT_.resize(count); delta_logPg_.resize(count);
  std::size_t present_count{};
  for (std::size_t i = 0; i < count; ++i) {
    int present{}; in >> present;
    if (!in || (present != 0 && present != 1)) throw std::runtime_error("invalid response mask");
    valid_[i] = present;
    if (!present) continue;
    in >> delta_logT_[i] >> delta_logPg_[i];
    if (!in || !std::isfinite(delta_logT_[i]) || !std::isfinite(delta_logPg_[i]))
      throw std::runtime_error("invalid atmosphere metal response");
    ++present_count;
  }
  std::string extra;
  if (!present_count || in >> extra) throw std::runtime_error("invalid response data extent");
}

std::optional<std::array<std::size_t, 3>>
MetalResponseAtmosphere::stencil(const std::array<double, 3> &q) const {
  std::array<std::size_t, 3> preferred{};
  for (std::size_t k = 0; k < 3; ++k) {
    if (!std::isfinite(q[k]) || q[k] < axes_[k].front() || q[k] > axes_[k].back())
      return std::nullopt;
    preferred[k] = interp::locate(axes_[k], q[k]);
  }
  for (unsigned alternative = 0; alternative < 8; ++alternative) {
    auto base = preferred;
    bool eligible = true;
    for (unsigned k = 0; k < 3; ++k) if (alternative & (1U << k)) {
      if (!base[k] || q[k] != axes_[k][base[k]]) { eligible = false; break; }
      --base[k];
    }
    if (!eligible) continue;
    bool complete = true;
    for (unsigned corner = 0; corner < 8; ++corner) {
      std::size_t index{};
      for (unsigned k = 0; k < 3; ++k)
        index = index*axes_[k].size()+base[k]+bool(corner & (1U << k));
      if (!valid_[index]) { complete = false; break; }
    }
    if (complete) return base;
  }
  return std::nullopt;
}

std::array<double, 4> MetalResponseAtmosphere::interpolate(
    const std::vector<double> &values, const std::array<double, 3> &q) const {
  const auto selected = stencil(q);
  if (!selected) throw std::domain_error("MetalResponseAtmosphere: incomplete source coverage");
  std::array<double, 3> u{}, width{};
  const auto &base = *selected;
  for (std::size_t k = 0; k < 3; ++k) {
    width[k] = axes_[k][base[k]+1]-axes_[k][base[k]];
    u[k] = (q[k]-axes_[k][base[k]])/width[k];
  }
  std::array<double, 4> result{};
  for (unsigned corner = 0; corner < 8; ++corner) {
    std::size_t index{};
    std::array<double, 3> weight{};
    for (unsigned k = 0; k < 3; ++k) {
      const bool high = corner & (1U << k);
      index = index*axes_[k].size()+base[k]+high;
      weight[k] = high ? u[k] : 1-u[k];
    }
    result[0] += values[index]*weight[0]*weight[1]*weight[2];
    for (unsigned k = 0; k < 3; ++k) {
      double w = (corner & (1U << k) ? 1. : -1.)/width[k];
      for (unsigned j = 0; j < 3; ++j) if (j != k) w *= weight[j];
      result[k+1] += values[index]*w;
    }
  }
  result[2] /= ln10; result[3] /= ln10;
  return result;
}

AtmosphereState MetalResponseAtmosphere::eval(double teff, double g, const Composition &c) const {
  if (!positive(teff) || !positive(g) || c.basis != AbundanceBasis::baryon_mass ||
      c.metal_inventory != MetalInventory::gs98 || std::abs(c.sum()-1) > 1e-10)
    throw std::domain_error("MetalResponseAtmosphere: invalid physical query");
  for (double v : c.X) if (!std::isfinite(v) || v < 0 || v > 1)
    throw std::domain_error("MetalResponseAtmosphere: invalid composition");
  const double z = c.Z();
  // Z is a sum of five stored mass fractions. Permit only its summation
  // roundoff at a source endpoint; keep the actual value in the response.
  const double roundoff = 8*std::numeric_limits<double>::epsilon()*reference_z_;
  if (z < source_z_-roundoff || z > reference_z_+options_.maximum_reference_metal_excess+roundoff ||
      c.X[1] > maximum_helium3_)
    throw std::domain_error("MetalResponseAtmosphere: composition outside response domain");
  const std::array q{c.X[0], std::log10(teff), std::log10(g)};
  const auto t = interpolate(delta_logT_, q), pg = interpolate(delta_logPg_, q);
  auto source_c = solar_scaled(c.X[0], reference_z_);
  source_c.basis = AbundanceBasis::baryon_mass;
  source_c.metal_inventory = MetalInventory::gs98;
  source_c.X[1] = c.X[1]; source_c.X[2] -= c.X[1];
  const auto base = reference_.eval(teff, g, source_c);
  if (base.tau != tau_ || !positive(base.Pgas))
    throw std::domain_error("MetalResponseAtmosphere: inconsistent reference boundary");
  const double f = (z-reference_z_)/(source_z_-reference_z_);
  AtmosphereState result{};
  result.T = base.T*std::exp(f*t[0]);
  result.Pgas = base.Pgas*std::exp(f*pg[0]);
  const double base_pr = constants::a_rad*std::pow(base.T, 4)/3;
  const double pr = constants::a_rad*std::pow(result.T, 4)/3;
  result.P = result.Pgas+pr; result.tau = tau_;
  if (!positive(result.T) || !positive(result.P) || !positive(result.Pgas) || result.P == pr)
    throw std::domain_error("MetalResponseAtmosphere: unrepresentable pressure or temperature");
  result.dlnT_dlnTeff = base.dlnT_dlnTeff+f*t[2];
  result.dlnT_dlng = base.dlnT_dlng+f*t[3];
  const double base_pgt = (base.P*base.dlnP_dlnTeff-4*base_pr*base.dlnT_dlnTeff)/base.Pgas;
  const double base_pgg = (base.P*base.dlnP_dlng-4*base_pr*base.dlnT_dlng)/base.Pgas;
  result.dlnP_dlnTeff = (result.Pgas*(base_pgt+f*pg[2])+4*pr*result.dlnT_dlnTeff)/result.P;
  result.dlnP_dlng = (result.Pgas*(base_pgg+f*pg[3])+4*pr*result.dlnT_dlng)/result.P;
  result.rho = eos_.rho_from_PT(result.T, result.P, c, base.rho);
  return result;
}
} // namespace ember
