#include "ember/opacity_ferguson.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ember {

FergusonOpacity::FergusonOpacity(const std::filesystem::path& file) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("FergusonOpacity: cannot open " + file.string());
  std::size_t nx = 0, nt = 0, nr = 0;
  in >> nx >> nt >> nr;
  if (!in || nx < 2 || nt < 4 || nr < 4)
    throw std::runtime_error("FergusonOpacity: bad header in " + file.string());
  std::string rest; std::getline(in, rest);

  logR_.resize(nr); for (auto& v : logR_) in >> v;
  logT_.resize(nt); for (auto& v : logT_) in >> v;
  X_.resize(nx);
  k_.resize(nx * nt * nr);
  for (std::size_t ix = 0; ix < nx; ++ix) {
    double z; in >> X_[ix] >> z;
    if (ix == 0) Z_ = z;
    for (std::size_t it = 0; it < nt; ++it)
      for (std::size_t ir = 0; ir < nr; ++ir)
        in >> k_[(ix * nt + it) * nr + ir];
  }
  if (!in) throw std::runtime_error("FergusonOpacity: truncated table " + file.string());
  if (!std::is_sorted(X_.begin(), X_.end()) ||
      !std::is_sorted(logT_.begin(), logT_.end()) ||
      !std::is_sorted(logR_.begin(), logR_.end()))
    throw std::runtime_error("FergusonOpacity: axes not monotone in " + file.string());
}

FergusonOpacity::Range FergusonOpacity::range() const {
  return {logT_.front(), logT_.back(), logR_.front(), logR_.back(), X_.front(), X_.back()};
}

bool FergusonOpacity::covers(double T, double rho, double X) const {
  const double lt = std::log10(T);
  const double lr = std::log10(rho) - 3.0 * (lt - 6.0);
  return lt >= logT_.front() && lt <= logT_.back()
      && lr >= logR_.front() && lr <= logR_.back()
      && X  >= X_.front()    && X  <= X_.back();
}

std::optional<Opacity::DensityRange> FergusonOpacity::density_range(double T, const Composition& comp) const {
  const double lt = std::log10(T), X = comp.h1();
  if (!std::isfinite(lt) || lt < logT_.front() || lt > logT_.back()
      || !std::isfinite(X) || X < X_.front() || X > X_.back())
    throw std::domain_error("FergusonOpacity: temperature or composition outside table");
  return DensityRange{std::pow(10.0, logR_.front() + 3.0 * (lt - 6.0)),
                      std::pow(10.0, logR_.back() + 3.0 * (lt - 6.0))};
}

OpacityState FergusonOpacity::eval(double T, double rho, const Composition& comp) const {
  const double lt = std::log10(T);
  const double lr = std::log10(rho) - 3.0 * (lt - 6.0);
  const double X  = comp.h1();
  if (!covers(T, rho, X))
    throw std::domain_error("FergusonOpacity: (T=" + std::to_string(T) +
                            ", rho=" + std::to_string(rho) + ", X=" + std::to_string(X) +
                            ") outside table");

  // Interpolate in log R along four bracketing isotherms, then across them in
  // log T, then linearly in X.  X is frozen through a Newton solve, so only T
  // and rho need the smoother treatment; the chain rule below therefore only
  // has to carry the first two.
  const std::size_t ix = interp::locate(X_, X);
  const std::size_t it = interp::locate(logT_, lt);
  const std::size_t t0 = std::min(it > 0 ? it - 1 : 0, logT_.size() - 4);

  double v[2]{}, dvdT[2]{}, dvdR[2]{};
  for (std::size_t m = 0; m < 2; ++m) {
    const std::size_t xi = std::min(ix + m, X_.size() - 1);
    double col[4], dcol[4], tt[4];
    for (std::size_t j = 0; j < 4; ++j) {
      const std::size_t ti = t0 + j;
      tt[j] = logT_[ti];
      // one isotherm, interpolated in log R
      const auto row = std::span<const double>(&k_[(xi * logT_.size() + ti) * logR_.size()],
                                               logR_.size());
      const auto r = interp::hermite(logR_, row, lr);
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
  const double dl_dlt = (1 - fx) * dvdT[0] + fx * dvdT[1];  // at fixed log R
  const double dl_dlr = (1 - fx) * dvdR[0] + fx * dvdR[1];

  // log R = log rho - 3 log T + 18, so at fixed rho a change in log T moves
  // log R by -3.  Forgetting that term is the classic error with these tables.
  OpacityState s{};
  s.kappa       = std::pow(10.0, logk);
  s.dlnk_dlnRho = dl_dlr;
  s.dlnk_dlnT   = dl_dlt - 3.0 * dl_dlr;
  return s;
}

} // namespace ember
