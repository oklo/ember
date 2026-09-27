#include "ember/atmosphere_grid.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <limits>
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
  if (!(in >> version) || (version < 1 || version > 4))
    throw std::runtime_error("CompositionAtmosphereGrid: version");
  helium_fraction_coordinates_ = version == 3;
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
  if (helium_fraction_coordinates_ || version == 4) {
    label("metal_tolerance");
    if (!(in >> metal_tolerance_) || !std::isfinite(metal_tolerance_) ||
        metal_tolerance_ < 0 || metal_tolerance_ > 1e-12)
      throw std::runtime_error("CompositionAtmosphereGrid: invalid metal tolerance");
  }
  const std::array labels{"hydrogen",
                         helium_fraction_coordinates_ ? "helium3_fraction" : "helium3",
                         "log_teff", "log_g"};
  std::size_t cells = 1;
  for (std::size_t k = 0; k < axes_.size(); ++k) {
    label(labels[k]);
    std::size_t n{};
    if (!(in >> n) || n < ((version == 4 && k == 1) ? 1u : 2u) || n > 10000 || cells > 1000000 / n)
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
  if (helium_fraction_coordinates_ ? axes_[0].back() + z >= 1
                                   : axes_[0].back() + axes_[1].back() + z > 1 + 1e-12)
    throw std::runtime_error(
        "CompositionAtmosphereGrid: invalid helium support");
  label("data");
  logT_.resize(cells);
  logPg_.resize(cells);
  valid_.assign(cells, true);
  std::size_t accepted = 0;
  for (std::size_t i = 0; i < cells; ++i) {
    if (version >= 2) {
      int present{};
      if (!(in >> present) || (present != 0 && present != 1))
        throw std::runtime_error("CompositionAtmosphereGrid: invalid source mask");
      if (!present) {
        valid_[i] = false;
        has_missing_states_ = true;
        continue;
      }
    }
    in >> logT_[i] >> logPg_[i];
    if (!in || !positive(std::pow(10., logT_[i])) ||
        !positive(std::pow(10., logPg_[i])))
      throw std::runtime_error(
          "CompositionAtmosphereGrid: invalid or missing source state");
    ++accepted;
  }
  if (!accepted)
    throw std::runtime_error("CompositionAtmosphereGrid: no source states");
  if(version == 4) {
    label("cells");
    std::size_t count=1;
    for(const auto& axis:axes_)count*=std::max(std::size_t{1},axis.size()-1);
    cell_valid_.resize(count);bool any=false;
    for(std::size_t i=0;i<count;++i) {
      int present{};
      if(!(in>>present) || (present!=0 && present!=1))
        throw std::runtime_error("CompositionAtmosphereGrid: invalid cell mask");
      cell_valid_[i]=present;any|=present;has_missing_states_|=!present;
      if(!present)continue;
      std::array<std::size_t,4> base{};auto value=i;
      for(int k=3;k>=0;--k) {const auto n=std::max(std::size_t{1},axes_[k].size()-1);base[k]=value%n;value/=n;}
      for(unsigned corner=0;corner<16;++corner) {
        std::size_t node=0;
        for(unsigned k=0;k<4;++k)node=node*axes_[k].size()+base[k]+(axes_[k].size()>1 && (corner&(1u<<k)));
        if(!valid_[node])throw std::runtime_error("CompositionAtmosphereGrid: supported cell has a missing vertex");
      }
    }
    if(!any)throw std::runtime_error("CompositionAtmosphereGrid: no supported cells");
  }
  std::string extra;
  if (in >> extra)
    throw std::runtime_error("CompositionAtmosphereGrid: trailing data");
}

CompositionAtmosphereGrid::Support CompositionAtmosphereGrid::support() const {
  const auto helium3 = helium_fraction_coordinates_
      ? std::array{axes_[1].front() * (1 - axes_[0].back() - reference_metallicity()),
                   axes_[1].back() * (1 - axes_[0].front() - reference_metallicity())}
      : std::array{axes_[1].front(), axes_[1].back()};
  return {{axes_[0].front(), axes_[0].back()},
          helium3,
          {std::pow(10., axes_[2].front()), std::pow(10., axes_[2].back())},
          {std::pow(10., axes_[3].front()), std::pow(10., axes_[3].back())}};
}

std::size_t CompositionAtmosphereGrid::check_temperature_extension(
    const CompositionAtmosphereGrid &old) const {
  if (source_ != old.source_ || approximation_ != old.approximation_ ||
      tau_ != old.tau_ || metals_ != old.metals_ ||
      helium_fraction_coordinates_ != old.helium_fraction_coordinates_ ||
      metal_tolerance_ != old.metal_tolerance_ || cell_valid_.empty()!=old.cell_valid_.empty())
    throw std::runtime_error("atmosphere extension: physics or matching depth changed");
  for (std::size_t k : {0u, 1u, 3u})
    if (axes_[k] != old.axes_[k])
      throw std::runtime_error("atmosphere extension: composition or gravity axis changed");
  if (axes_[2].size() <= old.axes_[2].size() ||
      !std::equal(old.axes_[2].begin(), old.axes_[2].end(), axes_[2].begin()))
    throw std::runtime_error("atmosphere extension: original temperature axis changed");
  std::size_t added = 0;
  for (std::size_t h = 0; h < axes_[0].size(); ++h)
    for (std::size_t he = 0; he < axes_[1].size(); ++he)
      for (std::size_t t = 0; t < axes_[2].size(); ++t)
        for (std::size_t g = 0; g < axes_[3].size(); ++g) {
          const auto i = ((h * axes_[1].size() + he) * axes_[2].size() + t) * axes_[3].size() + g;
          if (t >= old.axes_[2].size()) {
            added += valid_[i];
            continue;
          }
          const auto j = ((h * axes_[1].size() + he) * old.axes_[2].size() + t) * axes_[3].size() + g;
          if (valid_[i] != old.valid_[j] ||
              (valid_[i] && (logT_[i] != old.logT_[j] || logPg_[i] != old.logPg_[j])))
            throw std::runtime_error("atmosphere extension: original source values or mask changed");
        }
  for(std::size_t i=0;i<old.cell_valid_.size();++i) {
    std::array<std::size_t,4> base{};auto value=i;
    for(int k=3;k>=0;--k) {const auto n=std::max(std::size_t{1},old.axes_[k].size()-1);base[k]=value%n;value/=n;}
    std::size_t j=0;
    for(std::size_t k=0;k<4;++k)j=j*std::max(std::size_t{1},axes_[k].size()-1)+base[k];
    if(cell_valid_[j]!=old.cell_valid_[i])throw std::runtime_error("atmosphere extension: original cell support changed");
  }
  if (!added)
    throw std::runtime_error("atmosphere extension: no source states added");
  return added;
}

bool CompositionAtmosphereGrid::covers(double Teff, double g,
                                       const Composition &c) const {
  if (!positive(Teff) || !positive(g) || c.basis != AbundanceBasis::baryon_mass)
    return false;
  if(c[Species::H2]!=0)return false; // Requires an explicitly selected isotope mapping.
  for (double x : c.X)
    if (!std::isfinite(x) || x < 0 || x > 1)
      return false;
  if (std::abs(c.sum() - 1) > 1e-10)
    return false;
  for (std::size_t i = 0; i < metals_.size(); ++i)
    if (std::abs(c.X[i + 3] - metals_[i]) > metal_tolerance_)
      return false;
  const double helium = c.X[1] + c.X[2];
  if (helium_fraction_coordinates_ && !(helium > 0))
    return false;
  const std::array q{c.X[0], helium_fraction_coordinates_ ? c.X[1] / helium : c.X[1],
                     std::log10(Teff), std::log10(g)};
  for (std::size_t i = 0; i < q.size(); ++i)
    if (q[i] < axes_[i].front() || q[i] > axes_[i].back())
      return false;
  return !has_missing_states_ || stencil(q).has_value();
}

std::optional<std::array<std::size_t, 4>>
CompositionAtmosphereGrid::stencil(const std::array<double, 4> &q) const {
  std::array<std::size_t, 4> preferred{};
  for (std::size_t k = 0; k < q.size(); ++k)
    preferred[k] = axes_[k].size()==1 ? 0 : interp::locate(axes_[k], q[k]);
  if (!has_missing_states_ && cell_valid_.empty())
    return preferred;
  // At an exact knot either incident cell provides a one-sided derivative.
  // Prefer the ordinary upper-side cell. If absent, use a complete lower-
  // side cell, preserving closed edges of the existing supported domain.
  // No tolerance, extrapolation, or omitted derivative corner is allowed.
  for (unsigned alternative = 0; alternative < 16; ++alternative) {
    auto base = preferred;
    bool eligible = true;
    for (unsigned k = 0; k < 4; ++k) {
      if (!(alternative & (1U << k)))
        continue;
      if (!base[k] || q[k] != axes_[k][base[k]]) {
        eligible = false;
        break;
      }
      --base[k];
    }
    if (!eligible)
      continue;
    if(!cell_valid_.empty()) {
      std::size_t cell=0;
      for(unsigned k=0;k<4;++k)cell=cell*std::max(std::size_t{1},axes_[k].size()-1)+base[k];
      if(!cell_valid_[cell])continue;
    }
    bool complete = true;
    for (unsigned corner = 0; corner < 16; ++corner) {
      std::size_t index = 0;
      for (unsigned k = 0; k < 4; ++k)
        index = index * axes_[k].size() + base[k] + (axes_[k].size()>1 && (corner & (1U << k)));
      if (!valid_[index]) {
        complete = false;
        break;
      }
    }
    if (complete)
      return base;
  }
  return std::nullopt;
}
std::array<double, 4>
CompositionAtmosphereGrid::coordinates(double Teff, double g,
                                       const Composition &c) const {
  if (!covers(Teff, g, c))
    throw std::domain_error("CompositionAtmosphereGrid: Teff, gravity or "
                            "composition outside source grid");
  return {c.X[0], helium_fraction_coordinates_ ? c.X[1] / (c.X[1] + c.X[2]) : c.X[1],
          std::log10(Teff), std::log10(g)};
}
std::array<double, 5>
CompositionAtmosphereGrid::interpolate(const std::vector<double> &f,
                                       const std::array<double, 4> &q) const {
  const auto selected = stencil(q);
  if (!selected)
    throw std::domain_error("CompositionAtmosphereGrid: incomplete source stencil");
  const auto &base = *selected;
  std::array<double, 4> u{}, width{};
  for (std::size_t k = 0; k < q.size(); ++k) {
    if(axes_[k].size()==1) {width[k]=1;u[k]=0;continue;}
    width[k] = axes_[k][base[k] + 1] - axes_[k][base[k]];
    u[k] = (q[k] - axes_[k][base[k]]) / width[k];
  }
  std::array<double, 5> result{};
  for (unsigned corner = 0; corner < 16; ++corner) {
    std::size_t index = 0;
    std::array<double, 4> w{};
    for (unsigned k = 0; k < 4; ++k) {
      const bool high = corner & (1U << k);
      index = index * axes_[k].size() + base[k] + (axes_[k].size()>1 && high);
      w[k] = high ? u[k] : 1 - u[k];
    }
    result[0] += f[index] * w[0] * w[1] * w[2] * w[3];
    for (unsigned k = 0; k < 4; ++k) {
      if(axes_[k].size()==1)continue;
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
  auto t = interpolate(logT_, q), pg = interpolate(logPg_, q);
  const double pr = constants::a_rad * std::pow(t[0], 4) / 3;
  const double p = pg[0] + pr;
  if (!positive(p) || p == pr)
    throw std::domain_error("CompositionAtmosphereGrid: pressure overflow");
  if (helium_fraction_coordinates_) {
    // H replaces He4 at fixed He3: df3/dXH=f3/Y. He3 replaces He4 at
    // fixed H: df3/dX3=1/Y. Return derivatives in physical mass fractions.
    const double helium = c.X[1] + c.X[2];
    t[1] += t[2] * q[1] / helium;
    pg[1] += pg[2] * q[1] / helium;
    t[2] /= helium;
    pg[2] /= helium;
  }
  if(axes_[1].size()==1)t[2]=pg[2]=std::numeric_limits<double>::quiet_NaN();
  return {t[1], t[2], (pg[0] * pg[1] + 4 * pr * t[1]) / p,
          (pg[0] * pg[2] + 4 * pr * t[2]) / p};
}
} // namespace ember
