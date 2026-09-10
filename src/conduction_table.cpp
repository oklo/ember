#include "ember/conduction_table.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>

namespace ember {
namespace {
void check_state(double T, double rho, const Composition& c) {
  if (!std::isfinite(T) || !std::isfinite(rho) || T <= 0 || rho <= 0 || std::abs(c.sum() - 1) > 1e-10)
    throw std::domain_error("Conduction: invalid state/composition");
  for (double v : c.X)
    if (!std::isfinite(v) || v < 0)
      throw std::domain_error("Conduction: invalid abundance");
}
} // namespace

TabulatedConduction::TabulatedConduction(const std::filesystem::path& file) {
  std::ifstream in(file);
  std::string magic;
  int version;
  std::size_t nz, nt, nr;
  in >> magic >> version >> nz >> nt >> nr >> label_;
  if (!in || magic != "EMBER_CONDUCTIVITY" || version != 1 || nz < 2 || nz > 100 || nt < 4 || nt > 1000 ||
      nr < 4 || nr > 1000)
    throw std::runtime_error("TabulatedConduction: invalid header");
  logT_.resize(nt);
  logRho_.resize(nr);
  logZ_.resize(nz);
  mass_.resize(nz);
  logK_.resize(nz * nt * nr);
  for (auto& v : logT_)
    in >> v;
  for (auto& v : logRho_)
    in >> v;
  for (std::size_t z = 0; z < nz; ++z) {
    double charge;
    in >> charge >> mass_[z];
    if (!in || !std::isfinite(charge) || charge <= 0 || !std::isfinite(mass_[z]) || mass_[z] < charge)
      throw std::runtime_error("TabulatedConduction: invalid source ion");
    logZ_[z] = std::log10(charge);
    for (std::size_t r = 0; r < nr; ++r)
      for (std::size_t t = 0; t < nt; ++t)
        in >> logK_[(z * nr + r) * nt + t];
  }
  if (!in)
    throw std::runtime_error("TabulatedConduction: truncated data");
  for (const auto* axis : {&logT_, &logRho_, &logZ_})
    if (!std::all_of(axis->begin(), axis->end(), [](double x) { return std::isfinite(x); }) ||
        std::adjacent_find(axis->begin(), axis->end(), std::greater_equal<double>{}) != axis->end())
      throw std::runtime_error("TabulatedConduction: invalid axis");
  for (double v : logK_)
    if (!std::isfinite(v) || std::abs(v) > 100)
      throw std::runtime_error("TabulatedConduction: invalid conductivity");
  if (in >> magic)
    throw std::runtime_error("TabulatedConduction: trailing data");
}
TabulatedConduction::Conductivity TabulatedConduction::pure(std::size_t z, double lt, double le) const {
  // Electron number density is rho*Ye/m_u. The tables use integer A.
  const double lr = le + std::log10(mass_[z]) - logZ_[z];
  if (lt < logT_.front() || lt > logT_.back() || lr < logRho_.front() || lr > logRho_.back())
    throw std::domain_error("TabulatedConduction: T/electron density outside source support");
  const std::size_t it = interp::locate(logT_, lt), ir = interp::locate(logRho_, lr);
  const std::size_t t0 = it ? it - 1 : 0, t1 = std::min(it + 2, logT_.size() - 1);
  const std::size_t r0 = ir ? ir - 1 : 0, r1 = std::min(ir + 2, logRho_.size() - 1);
  const auto tt = std::span(logT_).subspan(t0, t1 - t0 + 1), rr = std::span(logRho_).subspan(r0, r1 - r0 + 1);
  std::array<double, 4> vals{}, dR{}, row{};
  for (std::size_t t = t0; t <= t1; ++t) {
    for (std::size_t r = r0; r <= r1; ++r)
      row[r - r0] = logK_[(z * logRho_.size() + r) * logT_.size() + t];
    const auto a = interp::hermite(rr, std::span(row).first(rr.size()), lr);
    vals[t - t0] = a.y;
    dR[t - t0] = a.dydx;
  }
  const auto a = interp::hermite(tt, std::span(vals).first(tt.size()), std::span(dR).first(tt.size()), lt);
  return {a.y, a.dydx, a.dydp};
}
TabulatedConduction::Conductivity TabulatedConduction::ion(
    double charge, double lt, double le, std::span<std::optional<Conductivity>> pure_cache) const {
  const double z = std::log10(charge);
  if (z < logZ_.front() || z > logZ_.back())
    throw std::domain_error("TabulatedConduction: unsupported ion charge");
  const auto source = [&](std::size_t i) {
    auto& value = pure_cache[i];
    if (!value) value = pure(i, lt, le);
    return *value;
  };
  const auto exact = std::lower_bound(logZ_.begin(), logZ_.end(), z);
  if (exact != logZ_.end() && *exact == z)
    return source(static_cast<std::size_t>(exact - logZ_.begin()));
  const auto i = interp::locate(logZ_, z);
  const double f = (z - logZ_[i]) / (logZ_[i + 1] - logZ_[i]);
  const auto a = source(i), b = source(i + 1);
  return {(1 - f) * a.logK + f * b.logK, (1 - f) * a.dT + f * b.dT, (1 - f) * a.drho + f * b.drho};
}
OpacityState TabulatedConduction::eval(double T, double rho, const Composition& c) const {
  check_state(T, rho, c);
  const double Ye = c.mu_elec_inv(), lt = std::log10(T), le = std::log10(rho * Ye);
  // Several elements interpolate between the same source ions, and both
  // helium isotopes use the same charge. All share exactly this T and electron
  // density. Reuse each source lookup within this call, preserving the order
  // of mixture sums and derivatives. No mutable state survives the call.
  std::array<std::optional<Conductivity>, 100> pure_cache{};
  double resistance = 0, dT = 0, drho = 0;
  std::array<double, 2> dc{};
  const std::array<double, 2> de{1 / c.abundance_weight(0) - 2 / c.abundance_weight(2),
                                 2 / c.abundance_weight(1) - 2 / c.abundance_weight(2)};
  auto add_ion=[&](double charge,double ej,std::array<double,2> dej) {
    const double f=ej/Ye;
    const auto k = ion(charge, lt, le, std::span(pure_cache).first(logZ_.size()));
    const double ik = std::pow(10., -k.logK), r = f * ik;
    resistance += r;
    dT -= r * k.dT;
    drho -= r * k.drho;
    for (std::size_t v = 0; v < 2; ++v) {
      dc[v] += ik * ((dej[v] - f * de[v]) / Ye - f * k.drho * de[v] / Ye);
    }
  };
  for (std::size_t j = 0; j < (c.metal_inventory==MetalInventory::gs98?3:NSPEC); ++j) {
    if(c.X[j]==0 && j>=3)continue;
    const double charge=nuclides[j].Z;
    std::array<double,2> dej{};
    for(std::size_t v=0;v<2;++v)dej[v]=(j==v?charge/c.abundance_weight(j):j==2?-charge/c.abundance_weight(j):0);
    add_ion(charge,c.X[j]*charge/c.abundance_weight(j),dej);
  }
  if(c.metal_inventory==MetalInventory::gs98 && c.Z()>0)
    for(const auto& e:gs98_metals)add_ion(e.charge,c.Z()*e.fraction*e.charge/e.mass_number/
        (c.basis==AbundanceBasis::baryon_mass?1.:gs98_atomic_mass_scale()),{});
  OpacityState out{};
  out.kappa = 16 * constants::sigma_SB * T * T * T / (3 * rho) * resistance;
  out.dlnk_dlnT = 3 + dT / resistance;
  out.dlnk_dlnRho = -1 + drho / resistance;
  out.dlnk_dX = dc[0] / resistance;
  out.dlnk_dY3 = dc[1] / resistance;
  return out;
}

std::optional<Opacity::DensityRange> TabulatedConduction::density_range(double T,
                                                                        const Composition& c) const {
  check_state(T, 1, c);
  const double lt = std::log10(T), Ye = c.mu_elec_inv();
  if (lt < logT_.front() || lt > logT_.back())
    throw std::domain_error("TabulatedConduction: temperature outside source support");
  Opacity::DensityRange range{0, std::numeric_limits<double>::infinity()};
  const auto include = [&](std::size_t i) {
    const double scale = std::pow(10., logZ_[i]) / (mass_[i] * Ye);
    range.min = std::max(range.min, std::pow(10., logRho_.front()) * scale);
    range.max = std::min(range.max, std::pow(10., logRho_.back()) * scale);
  };
  auto include_charge=[&](double charge) {
    const double z = std::log10(charge);
    if (z < logZ_.front() || z > logZ_.back())
      throw std::domain_error("TabulatedConduction: unsupported ion charge");
    const auto exact = std::lower_bound(logZ_.begin(), logZ_.end(), z);
    if (exact != logZ_.end() && *exact == z) {
      include(static_cast<std::size_t>(exact - logZ_.begin()));
    } else {
      const auto i = interp::locate(logZ_, z);
      include(i);
      include(i + 1);
    }
  };
  for(std::size_t j=0;j<(c.metal_inventory==MetalInventory::gs98?3:NSPEC);++j)
    if(c.X[j]!=0 || j<3)include_charge(nuclides[j].Z);
  if(c.metal_inventory==MetalInventory::gs98 && c.Z()>0)
    for(const auto& e:gs98_metals)include_charge(e.charge);
  if (range.min >= range.max)
    throw std::domain_error("TabulatedConduction: no common source density support");
  return range;
}
HotConduction::HotConduction(const Conduction& c, double start, double end)
    : source_(c), start_(std::log(start)), end_(std::log(end)) {
  if (!std::isfinite(start_) || !std::isfinite(end_) || start_ >= end_)
    throw std::invalid_argument("HotConduction: invalid join");
}
OpacityState HotConduction::eval(double T, double rho, const Composition& c) const {
  check_state(T, rho, c);
  const double t = std::log(T);
  if (t <= start_)
    return {std::numeric_limits<double>::infinity(), 0, 0, 0, 0, 0};
  auto out = source_.eval(T, rho, c);
  if (t >= end_)
    return out;
  const double u = (t - start_) / (end_ - start_), w = u * u * (3 - 2 * u),
               dw = 6 * u * (1 - u) / (end_ - start_);
  out.kappa /= w;
  out.dlnk_dlnT -= dw / w;
  return out;
}
std::optional<Opacity::DensityRange> HotConduction::density_range(double T, const Composition& c) const {
  check_state(T, 1, c);
  return std::log(T) <= start_ ? std::nullopt : source_.density_range(T, c);
}
} // namespace ember
