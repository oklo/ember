#include "ember/atmosphere_grid.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>

namespace ember {
namespace {
bool positive(double v) { return std::isfinite(v) && v > 0; }
constexpr double ln10 = 2.302585092994045684;
} // namespace

CompositionAtmosphereGrid::CompositionAtmosphereGrid(
    const Eos &eos, const std::filesystem::path &path, Mixture mixture)
    : eos_(eos) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("CompositionAtmosphereGrid: cannot open " +
                             path.string());
  read(input, mixture);
}
CompositionAtmosphereGrid::CompositionAtmosphereGrid(const Eos &eos,
                                                     std::istream &input,
                                                     Mixture mixture)
    : eos_(eos) {
  read(input, mixture);
}

void CompositionAtmosphereGrid::read(std::istream &in, Mixture mixture) {
  const auto label = [&](const char *expected) {
    std::string s;
    if (!(in >> s) || s != expected)
      throw std::runtime_error(
          std::string("CompositionAtmosphereGrid: expected ") + expected);
  };
  label("EMBER_COMPOSITION_ATMOSPHERE");
  int version{};
  if (!(in >> version) || version != 1)
    throw std::runtime_error("CompositionAtmosphereGrid: version");
  label("source");
  in >> std::quoted(source_);
  label("approximation");
  in >> std::quoted(approximation_);
  if (!in || source_.empty() || approximation_.empty())
    throw std::runtime_error("CompositionAtmosphereGrid: missing provenance");
  if (approximation_ != "none" && mixture != Mixture::allow_documented_proxy)
    throw std::invalid_argument("CompositionAtmosphereGrid: explicit "
                                "documented proxy selection required");
  label("basis");
  label("baryon_mass");
  label("tau");
  in >> tau_;
  if (!in || !positive(tau_))
    throw std::runtime_error("CompositionAtmosphereGrid: invalid depth");
  label("metals");
  for (auto &m : metals_) {
    in >> m;
    if (!in || !std::isfinite(m) || m < 0 || m >= 1)
      throw std::runtime_error("CompositionAtmosphereGrid: invalid metals");
  }
  const double z = std::accumulate(metals_.begin(), metals_.end(), 0.);
  if (z >= 1)
    throw std::runtime_error(
        "CompositionAtmosphereGrid: invalid total metallicity");
  const std::array labels{"hydrogen", "helium3", "log_teff", "log_g"};
  std::size_t cells = 1;
  for (std::size_t k = 0; k < axes_.size(); ++k) {
    label(labels[k]);
    std::size_t n{};
    if (!(in >> n) || n < 2 || n > 10000 || cells > 1000000 / n)
      throw std::runtime_error("CompositionAtmosphereGrid: invalid grid size");
    cells *= n;
    auto &axis = axes_[k];
    axis.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      in >> axis[i];
      if (!in || !std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]) ||
          (k < 2 ? axis[i] < 0 || axis[i] > 1
                 : !positive(std::pow(10., axis[i]))))
        throw std::runtime_error("CompositionAtmosphereGrid: invalid axis");
    }
  }
  if (axes_[0].back() + axes_[1].back() + z > 1 + 1e-12)
    throw std::runtime_error(
        "CompositionAtmosphereGrid: grid includes negative He4");
  label("data");
  logT_.resize(cells);
  logPg_.resize(cells);
  for (std::size_t i = 0; i < cells; ++i) {
    in >> logT_[i] >> logPg_[i];
    if (!in || !positive(std::pow(10., logT_[i])) ||
        !positive(std::pow(10., logPg_[i])))
      throw std::runtime_error(
          "CompositionAtmosphereGrid: invalid or missing source state");
  }
  std::string extra;
  if (in >> extra)
    throw std::runtime_error("CompositionAtmosphereGrid: trailing data");
}

CompositionAtmosphereGrid::Support CompositionAtmosphereGrid::support() const {
  return {{axes_[0].front(), axes_[0].back()},
          {axes_[1].front(), axes_[1].back()},
          {std::pow(10., axes_[2].front()), std::pow(10., axes_[2].back())},
          {std::pow(10., axes_[3].front()), std::pow(10., axes_[3].back())}};
}

bool CompositionAtmosphereGrid::covers(double Teff, double g,
                                       const Composition &c) const {
  if (!positive(Teff) || !positive(g) || c.basis != AbundanceBasis::baryon_mass)
    return false;
  for (double x : c.X)
    if (!std::isfinite(x) || x < 0 || x > 1)
      return false;
  if (std::abs(c.sum() - 1) > 1e-10)
    return false;
  for (std::size_t i = 0; i < metals_.size(); ++i)
    if (std::abs(c.X[i + 3] - metals_[i]) > 1e-12)
      return false;
  const std::array q{c.X[0], c.X[1], std::log10(Teff), std::log10(g)};
  for (std::size_t i = 0; i < q.size(); ++i)
    if (q[i] < axes_[i].front() || q[i] > axes_[i].back())
      return false;
  return true;
}
std::array<double, 4>
CompositionAtmosphereGrid::coordinates(double Teff, double g,
                                       const Composition &c) const {
  if (!covers(Teff, g, c))
    throw std::domain_error("CompositionAtmosphereGrid: Teff, gravity or "
                            "composition outside source grid");
  return {c.X[0], c.X[1], std::log10(Teff), std::log10(g)};
}
std::array<double, 5>
CompositionAtmosphereGrid::interpolate(const std::vector<double> &f,
                                       const std::array<double, 4> &q) const {
  std::array<std::size_t, 4> base{};
  std::array<double, 4> u{}, width{};
  for (std::size_t k = 0; k < q.size(); ++k) {
    base[k] = interp::locate(axes_[k], q[k]);
    width[k] = axes_[k][base[k] + 1] - axes_[k][base[k]];
    u[k] = (q[k] - axes_[k][base[k]]) / width[k];
  }
  std::array<double, 5> result{};
  for (unsigned corner = 0; corner < 16; ++corner) {
    std::size_t index = 0;
    std::array<double, 4> w{};
    for (unsigned k = 0; k < 4; ++k) {
      const bool high = corner & (1U << k);
      index = index * axes_[k].size() + base[k] + high;
      w[k] = high ? u[k] : 1 - u[k];
    }
    result[0] += f[index] * w[0] * w[1] * w[2] * w[3];
    for (unsigned k = 0; k < 4; ++k) {
      double derivative = (corner & (1U << k) ? 1. : -1.) / width[k];
      for (unsigned j = 0; j < 4; ++j)
        if (j != k)
          derivative *= w[j];
      result[k + 1] += f[index] * derivative;
    }
  }
  result[0] = std::pow(10., result[0]);
  result[1] *= ln10;
  result[2] *= ln10;
  return result;
}
AtmosphereState CompositionAtmosphereGrid::eval(double Teff, double g,
                                                const Composition &c) const {
  const auto q = coordinates(Teff, g, c);
  const auto t = interpolate(logT_, q), pg = interpolate(logPg_, q);
  const double pr = constants::a_rad * std::pow(t[0], 4) / 3;
  AtmosphereState s{};
  s.T = t[0];
  s.Pgas = pg[0];
  s.P = pg[0] + pr;
  s.tau = tau_;
  if (!positive(s.P) || s.P == pr)
    throw std::domain_error(
        "CompositionAtmosphereGrid: unrepresentable total pressure");
  s.rho = eos_.rho_from_PT(s.T, s.P, c);
  s.dlnT_dlnTeff = t[3];
  s.dlnT_dlng = t[4];
  s.dlnP_dlnTeff = (s.Pgas * pg[3] + 4 * pr * t[3]) / s.P;
  s.dlnP_dlng = (s.Pgas * pg[4] + 4 * pr * t[4]) / s.P;
  return s;
}
CompositionAtmosphereGrid::CompositionResponse
CompositionAtmosphereGrid::composition_response(double Teff, double g,
                                                const Composition &c) const {
  const auto q = coordinates(Teff, g, c);
  const auto t = interpolate(logT_, q), pg = interpolate(logPg_, q);
  const double pr = constants::a_rad * std::pow(t[0], 4) / 3;
  const double p = pg[0] + pr;
  if (!positive(p) || p == pr)
    throw std::domain_error("CompositionAtmosphereGrid: pressure overflow");
  return {t[1], t[2], (pg[0] * pg[1] + 4 * pr * t[1]) / p,
          (pg[0] * pg[2] + 4 * pr * t[2]) / p};
}
} // namespace ember
