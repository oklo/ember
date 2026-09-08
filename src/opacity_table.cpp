#include "ember/opacity_table.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace ember {
namespace {
bool inside(double value, double lo, double hi) {
  // Decimal table axes round on the (logT,logR)->(T,rho)->log round trip.
  // Admit only floating-point roundoff, then snap to the boundary in eval.
  const double tolerance = 16 * std::numeric_limits<double>::epsilon()
                         * std::max({1.0, std::abs(lo), std::abs(hi)});
  return std::isfinite(value) && value >= lo - tolerance && value <= hi + tolerance;
}
}

TabulatedOpacity::TabulatedOpacity(const std::filesystem::path& file, std::string label, DensityAxis axis)
    : label_(std::move(label)), axis_(axis) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("TabulatedOpacity: cannot open " + file.string());
  std::size_t nx = 0, nt = 0, nr = 0;
  in >> nx >> nt >> nr;
  if (!in || nx < 1 || nt < 4 || nr < 4 || nx > 1000 || nt > 10000 || nr > 10000
      || nx * nt * nr > 10000000)
    throw std::runtime_error("TabulatedOpacity: bad header in " + file.string());
  std::string rest; std::getline(in, rest);

  logD_.resize(nr); for (auto& v : logD_) in >> v;
  logT_.resize(nt); for (auto& v : logT_) in >> v;
  X_.resize(nx);
  k_.resize(nx * nt * nr);
  for (std::size_t ix = 0; ix < nx; ++ix) {
    double z{}; in >> X_[ix] >> z;
    if (ix == 0) Z_ = z;
    if (!in || !std::isfinite(z) || z < 0.0 || z > 1.0 || z != Z_)
      throw std::runtime_error("TabulatedOpacity: inconsistent or invalid metallicity");
    for (std::size_t it = 0; it < nt; ++it)
      for (std::size_t ir = 0; ir < nr; ++ir)
        in >> k_[(ix * nt + it) * nr + ir];
  }
  if (!in) throw std::runtime_error("TabulatedOpacity: truncated table " + file.string());
  auto valid_axis = [](const auto& axis) {
    return std::all_of(axis.begin(), axis.end(), [](double v) { return std::isfinite(v); })
        && std::adjacent_find(axis.begin(), axis.end(), std::greater_equal<double>{}) == axis.end();
  };
  if (!valid_axis(X_) || !valid_axis(logT_) || !valid_axis(logD_))
    throw std::runtime_error("TabulatedOpacity: axes not monotone in " + file.string());
  if (X_.front() < 0.0 || X_.back() + Z_ > 1.0 + 1e-12
      || !std::all_of(k_.begin(), k_.end(), [](double v) {
           return std::isfinite(v) && v != 9.999 && std::isfinite(std::pow(10.0, v))
               && std::pow(10.0, v) > 0.0;
         }))
    throw std::runtime_error("TabulatedOpacity: invalid composition, opacity, or missing cell");
  std::string extra;
  if (in >> extra) throw std::runtime_error("TabulatedOpacity: unexpected trailing data");
}

TabulatedOpacity::Range TabulatedOpacity::range() const {
  return {logT_.front(), logT_.back(), logD_.front(), logD_.back(), X_.front(), X_.back(), axis_};
}

bool TabulatedOpacity::covers(double T, double rho, double X) const {
  const double lt = std::log10(T);
  const double lr = std::log10(rho) - (axis_ == DensityAxis::logR ? 3.0 * (lt - 6.0) : 0.0);
  return inside(lt, logT_.front(), logT_.back())
      && inside(lr, logD_.front(), logD_.back())
      && X  >= X_.front()    && X  <= X_.back();
}

std::optional<Opacity::DensityRange> TabulatedOpacity::density_range(double T, const Composition& comp) const {
  const double raw_lt = std::log10(T), X = comp.h1();
  if (!inside(raw_lt, logT_.front(), logT_.back())
      || !std::isfinite(X) || X < X_.front() || X > X_.back()
      || !std::isfinite(comp.Z()) || std::abs(comp.Z() - Z_) > 1e-10)
    throw std::domain_error("TabulatedOpacity: temperature or composition outside table");
  const double lt = std::clamp(raw_lt, logT_.front(), logT_.back());
  const double shift = axis_ == DensityAxis::logR ? 3.0 * (lt - 6.0) : 0.0;
  return DensityRange{std::pow(10.0, logD_.front() + shift), std::pow(10.0, logD_.back() + shift)};
}

OpacityState TabulatedOpacity::eval(double T, double rho, const Composition& comp) const {
  double lt = std::log10(T);
  double lr = std::log10(rho) - (axis_ == DensityAxis::logR ? 3.0 * (lt - 6.0) : 0.0);
  const double X  = comp.h1();
  if (!std::isfinite(comp.Z()) || std::abs(comp.Z() - Z_) > 1e-10)
    throw std::domain_error(label_ + ": composition metallicity does not match fixed-Z table");
  if (!covers(T, rho, X))
    throw std::domain_error(label_ + ": (logT=" + std::to_string(lt) +
                            (axis_ == DensityAxis::logR ? ", logR=" : ", logrho=") + std::to_string(lr) + ", X=" + std::to_string(X) +
                            ") outside table");
  lt = std::clamp(lt, logT_.front(), logT_.back());
  lr = std::clamp(lr, logD_.front(), logD_.back());

  // Interpolate in the density coordinate along four bracketing isotherms, then across them in
  // log T, then linearly in X.  X is frozen through a Newton solve, so only T
  // and rho need the smoother treatment; the chain rule below therefore only
  // has to carry the first two.
  const std::size_t ix = X_.size() == 1 ? 0 : interp::locate(X_, X);
  const std::size_t it = interp::locate(logT_, lt);
  const std::size_t t0 = std::min(it > 0 ? it - 1 : 0, logT_.size() - 4);

  double v[2]{}, dvdT[2]{}, dvdR[2]{};
  for (std::size_t m = 0; m < 2; ++m) {
    const std::size_t xi = std::min(ix + m, X_.size() - 1);
    double col[4], dcol[4], tt[4];
    for (std::size_t j = 0; j < 4; ++j) {
      const std::size_t ti = t0 + j;
      tt[j] = logT_[ti];
      // one isotherm, interpolated in the density coordinate
      const auto row = std::span<const double>(&k_[(xi * logT_.size() + ti) * logD_.size()],
                                               logD_.size());
      const auto r = interp::hermite(logD_, row, lr);
      col[j] = r.y; dcol[j] = r.dydx;
    }
    const auto a = interp::hermite(std::span<const double>(tt, 4),
                                   std::span<const double>(col, 4),
                                   std::span<const double>(dcol, 4), lt);
    v[m] = a.y; dvdT[m] = a.dydx;
    dvdR[m] = a.dydp;
  }
  double fx = 0.0;
  if (X_.size() > 1 && ix + 1 < X_.size())
    fx = std::clamp((X - X_[ix]) / (X_[ix + 1] - X_[ix]), 0.0, 1.0);
  const double logk   = (1 - fx) * v[0]    + fx * v[1];
  const double dl_dlt = (1 - fx) * dvdT[0] + fx * dvdT[1];  // at fixed density coordinate
  const double dl_dlr = (1 - fx) * dvdR[0] + fx * dvdR[1];

  // log R = log rho - 3 log T + 18, so at fixed rho a change in log T moves
  // log R by -3.  Forgetting that term is the classic error with these tables.
  OpacityState s{};
  s.kappa       = std::pow(10.0, logk);
  s.dlnk_dlnRho = dl_dlr;
  s.dlnk_dlnT   = dl_dlt - (axis_ == DensityAxis::logR ? 3.0 * dl_dlr : 0.0);
  return s;
}

} // namespace ember
