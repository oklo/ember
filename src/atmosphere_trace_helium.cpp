#include "ember/atmosphere_trace_helium.hpp"
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

TraceHeliumAtmosphereGrid::TraceHeliumAtmosphereGrid(
    const PressureDensity &eos, const std::filesystem::path &path, Approximation approximation, double maximum_helium3)
    : eos_(eos) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("TraceHeliumAtmosphereGrid: cannot open " +
                             path.string());
  read(input, approximation, maximum_helium3);
}
TraceHeliumAtmosphereGrid::TraceHeliumAtmosphereGrid(const PressureDensity &eos,
                                                     std::istream &input,
                                                     Approximation approximation, double maximum_helium3)
    : eos_(eos) {
  read(input, approximation, maximum_helium3);
}

void TraceHeliumAtmosphereGrid::read(std::istream &in, Approximation approximation, double maximum_helium3) {
  const auto label = [&](const char *expected) {
    std::string s;
    if (!(in >> s) || s != expected)
      throw std::runtime_error(
          std::string("TraceHeliumAtmosphereGrid: expected ") + expected);
  };
  label("EMBER_TRACE_HELIUM_ATMOSPHERE");
  int version{};
  if (!(in >> version) || version < 1 || version > 3)
    throw std::runtime_error("TraceHeliumAtmosphereGrid: version");
  hydrogen_share_ = version >= 2;
  irregular_ = version == 3;
  label("source");
  in >> std::quoted(source_);
  label("approximation");
  in >> std::quoted(approximation_);
  if (!in || source_.empty() || approximation_.empty())
    throw std::runtime_error("TraceHeliumAtmosphereGrid: missing provenance");
  if (approximation != Approximation::neglect_atmospheric_helium3 ||
      !std::isfinite(maximum_helium3) || maximum_helium3 <= 0 || maximum_helium3 >= 1)
    throw std::invalid_argument("TraceHeliumAtmosphereGrid: explicit isotope approximation and bound required");
  label("basis");
  label("baryon_mass");
  label("tau");
  in >> tau_;
  if (!in || !positive(tau_))
    throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid depth");
  label("maximum_helium3"); in >> maximum_helium3_;
  if (!in || !std::isfinite(maximum_helium3_) || maximum_helium3_ <= 0 ||
      maximum_helium3_ > maximum_helium3)
    throw std::invalid_argument("TraceHeliumAtmosphereGrid: isotope bound exceeds caller selection");
  label("metal_pattern");
  for (auto &m : metal_pattern_) {
    in >> m;
    if (!in || !std::isfinite(m) || m < 0 || m >= 1)
      throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid metals");
  }
  const double z = std::accumulate(metal_pattern_.begin(), metal_pattern_.end(), 0.);
  if (std::abs(z-1) > 1e-12)
    throw std::runtime_error(
        "TraceHeliumAtmosphereGrid: invalid total metallicity");
  if (irregular_) {
    label("irregular_scale");
    for (auto& v : offset_) {
      if (!(in >> v) || !std::isfinite(v))
        throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid irregular offset");
    }
    for (auto& v : width_) { in >> v; if (!in || !positive(v)) throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid irregular scale"); }
    std::size_t n{}; label("vertices");
    if (!(in >> n) || n < 5 || n > 100000) throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid vertex count");
    vertex_.resize(n); vertex_value_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      std::array<double,4> raw{}; for (auto& v : raw) in >> v;
      in >> vertex_value_[i][0] >> vertex_value_[i][1];
      if (!in || !(raw[0] > 0 && raw[0] <= 1) || !(raw[1] >= 0 && raw[1] < 1) || !std::isfinite(raw[2]) || !std::isfinite(raw[3])
          || !positive(std::pow(10., raw[2])) || !positive(std::pow(10., raw[3]))
          || !positive(std::pow(10., vertex_value_[i][0])) || !positive(std::pow(10., vertex_value_[i][1])))
        throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid irregular vertex");
      for (int k = 0; k < 4; ++k) {
        vertex_[i][k] = (raw[k] - offset_[k]) / width_[k];
        if (!std::isfinite(vertex_[i][k]))
          throw std::runtime_error("TraceHeliumAtmosphereGrid: unrepresentable scaled vertex");
      }
    }
    std::size_t m{}; label("simplices");
    if (!(in >> m) || m < 1 || m > 10000000) throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid simplex count");
    simplex_.resize(m); inverse_.resize(m); gradient_.resize(m);
    for (std::size_t j = 0; j < m; ++j) {
      auto& s = simplex_[j];
      for (auto& i : s) if (!(in >> i) || i >= n) throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid simplex index");
      // A = [v1-v0 ... v4-v0] as columns; invert by Gauss-Jordan with partial pivoting.
      std::array<double,32> a{};
      for (int r = 0; r < 4; ++r) { for (int c = 0; c < 4; ++c) a[r*8+c] = vertex_[s[c+1]][r] - vertex_[s[0]][r]; a[r*8+4+r] = 1; }
      double scale = 0; for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) scale = std::max(scale, std::abs(a[r*8+c]));
      for (int c = 0; c < 4; ++c) {
        int piv = c; for (int r = c+1; r < 4; ++r) if (std::abs(a[r*8+c]) > std::abs(a[piv*8+c])) piv = r;
        if (!(std::abs(a[piv*8+c]) > 1e-12*scale)) throw std::runtime_error("TraceHeliumAtmosphereGrid: degenerate irregular simplex");
        if (piv != c) for (int k = 0; k < 8; ++k) std::swap(a[c*8+k], a[piv*8+k]);
        const double d = a[c*8+c]; for (int k = 0; k < 8; ++k) a[c*8+k] /= d;
        for (int r = 0; r < 4; ++r) if (r != c) { const double f = a[r*8+c]; for (int k = 0; k < 8; ++k) a[r*8+k] -= f*a[c*8+k]; }
      }
      for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) inverse_[j][r*4+c] = a[r*8+4+c];
      // Gradient in scaled space: g_k = sum_r (f_{r+1}-f_0) * Ainv[r][k]; raw derivative divides by the width.
      for (int q = 0; q < 2; ++q) for (int k = 0; k < 4; ++k) {
        double g = 0; for (int r = 0; r < 4; ++r) g += (vertex_value_[s[r+1]][q] - vertex_value_[s[0]][q]) * inverse_[j][r*4+k];
        gradient_[j][q*4+k] = g / width_[k];
      }
    }
    std::string extra;
    if (in >> extra) throw std::runtime_error("TraceHeliumAtmosphereGrid: trailing data");
    return;
  }
  const std::array labels{hydrogen_share_ ? "hydrogen_share" : "hydrogen", "metallicity", "log_teff", "log_g"};
  std::size_t cells = 1;
  for (std::size_t k = 0; k < axes_.size(); ++k) {
    label(labels[k]);
    std::size_t n{};
    if (!(in >> n) || n < 2 || n > 10000 || cells > 1000000 / n)
      throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid grid size");
    cells *= n;
    auto &axis = axes_[k];
    axis.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      in >> axis[i];
      if (!in || !std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]) ||
          (k < 2 ? axis[i] < 0 || axis[i] > 1
                 : !positive(std::pow(10., axis[i]))))
        throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid axis");
    }
  }
  if (hydrogen_share_ ? axes_[1].back() >= 1 : axes_[0].back() + axes_[1].back() > 1 + 1e-12)
    throw std::runtime_error(
        "TraceHeliumAtmosphereGrid: grid includes negative He4");
  label("data");
  logT_.resize(cells);
  logPg_.resize(cells);
  valid_.assign(cells, true);
  std::size_t accepted = 0;
  for (std::size_t i = 0; i < cells; ++i) {
    {
      int present{};
      if (!(in >> present) || (present != 0 && present != 1))
        throw std::runtime_error("TraceHeliumAtmosphereGrid: invalid source mask");
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
          "TraceHeliumAtmosphereGrid: invalid or missing source state");
    ++accepted;
  }
  if (!accepted)
    throw std::runtime_error("TraceHeliumAtmosphereGrid: no source states");
  std::string extra;
  if (in >> extra)
    throw std::runtime_error("TraceHeliumAtmosphereGrid: trailing data");
}

std::array<double,5> TraceHeliumAtmosphereGrid::interpolate_value(std::size_t which,const std::array<double,4>& q) const {
  // which: 0 = temperature, 1 = gas pressure (explicit index, not identified by storage).
  if (which > 1) throw std::invalid_argument("TraceHeliumAtmosphereGrid: unknown interpolated quantity");
  if (irregular_) return interpolate_irregular(which, q);
  return interpolate(which == 0 ? logT_ : logPg_, q);
}
std::optional<std::size_t> TraceHeliumAtmosphereGrid::locate(const std::array<double,4>& raw) const {
  std::array<double,4> x{}; for (int k = 0; k < 4; ++k) x[k] = (raw[k] - offset_[k]) / width_[k];
  // Closed simplices; the lowest-index containing simplex is used, so shared faces are deterministic.
  constexpr double tol = 1e-12;
  for (std::size_t j = 0; j < simplex_.size(); ++j) {
    const auto& v0 = vertex_[simplex_[j][0]];
    double sum = 0; bool inside = true;
    for (int r = 0; r < 4 && inside; ++r) {
      double l = 0; for (int k = 0; k < 4; ++k) l += inverse_[j][r*4+k] * (x[k] - v0[k]);
      inside = l >= -tol; sum += l;
    }
    if (inside && 1 - sum >= -tol) return j;
  }
  return std::nullopt;
}
std::array<double,5> TraceHeliumAtmosphereGrid::interpolate_irregular(std::size_t q,const std::array<double,4>& raw) const {
  const auto j = locate(raw);
  if (!j) throw std::domain_error("TraceHeliumAtmosphereGrid: outside the irregular source hull");
  const auto& s = simplex_[*j];
  double value = vertex_value_[s[0]][q];
  for (int k = 0; k < 4; ++k) value += gradient_[*j][q*4+k] * (raw[k] - (vertex_[s[0]][k]*width_[k] + offset_[k]));
  // Same convention as the rectangular reader: value, then d ln f/dh, d ln f/dZ, d ln f/d ln Teff, d ln f/d ln g.
  return {std::pow(10., value), gradient_[*j][q*4+0]*ln10, gradient_[*j][q*4+1]*ln10, gradient_[*j][q*4+2], gradient_[*j][q*4+3]};
}

TraceHeliumAtmosphereGrid::Support TraceHeliumAtmosphereGrid::support() const {
  if (irregular_) {
    Support out{{1,0},{1,0},{1e300,0},{1e300,0}};
    for (const auto& v : vertex_) {
      std::array<double,4> r{}; for (int k = 0; k < 4; ++k) r[k] = v[k]*width_[k] + offset_[k];
      out.hydrogen = {std::min(out.hydrogen[0], r[0]*(1-r[1])), std::max(out.hydrogen[1], r[0]*(1-r[1]))};
      out.metals = {std::min(out.metals[0], r[1]), std::max(out.metals[1], r[1])};
      out.teff = {std::min(out.teff[0], std::pow(10., r[2])), std::max(out.teff[1], std::pow(10., r[2]))};
      out.gravity = {std::min(out.gravity[0], std::pow(10., r[3])), std::max(out.gravity[1], std::pow(10., r[3]))};
    }
    return out;
  }
  return {{axes_[0].front()*(hydrogen_share_ ? 1-axes_[1].back() : 1),
           axes_[0].back()*(hydrogen_share_ ? 1-axes_[1].front() : 1)},
          {axes_[1].front(), axes_[1].back()},
          {std::pow(10., axes_[2].front()), std::pow(10., axes_[2].back())},
          {std::pow(10., axes_[3].front()), std::pow(10., axes_[3].back())}};
}

bool TraceHeliumAtmosphereGrid::covers(double Teff, double g, const Composition& c) const {
  if (!positive(Teff) || !positive(g) || c.basis != AbundanceBasis::baryon_mass ||
      c.metal_inventory != MetalInventory::gs98 || c.X[1] > maximum_helium3_)
    return false;
  for (double x : c.X) if (!std::isfinite(x) || x < 0 || x > 1) return false;
  if (std::abs(c.sum()-1) > 1e-10) return false;
  const double z=c.Z();
  for (std::size_t j=0;j<metal_pattern_.size();++j)
    if (std::abs(c.X[j+3]-z*metal_pattern_[j]) > 1e-12*std::max(z,1e-30)) return false;
  if (irregular_) return locate({c.X[0]/(1-z), z, std::log10(Teff), std::log10(g)}).has_value();
  // Only the floating-point summation error in Z may be rounded to an edge.
  const double roundoff=8*std::numeric_limits<double>::epsilon()*axes_[1].back();
  if (z < axes_[1].front()-roundoff || z > axes_[1].back()+roundoff) return false;
  const double zs=std::clamp(z,axes_[1].front(),axes_[1].back());
  const std::array q{c.X[0]/(hydrogen_share_ ? 1-zs : 1),zs,std::log10(Teff),std::log10(g)};
  for (std::size_t j=0;j<4;++j)
    if (q[j] < axes_[j].front() || q[j] > axes_[j].back()) return false;
  return !has_missing_states_ || stencil(q).has_value();
}

std::optional<std::array<std::size_t, 4>>
TraceHeliumAtmosphereGrid::stencil(const std::array<double, 4> &q) const {
  std::array<std::size_t, 4> preferred{};
  for (std::size_t k = 0; k < q.size(); ++k)
    preferred[k] = interp::locate(axes_[k], q[k]);
  if (!has_missing_states_)
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
    bool complete = true;
    for (unsigned corner = 0; corner < 16; ++corner) {
      std::size_t index = 0;
      for (unsigned k = 0; k < 4; ++k)
        index = index * axes_[k].size() + base[k] + bool(corner & (1U << k));
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
TraceHeliumAtmosphereGrid::coordinates(double Teff, double g,
                                       const Composition &c) const {
  if (!covers(Teff, g, c))
    throw std::domain_error("TraceHeliumAtmosphereGrid: Teff, gravity or "
                            "composition outside source grid");
  if (irregular_) { const double z=c.Z(); return {c.X[0]/(1-z), z, std::log10(Teff), std::log10(g)}; }
  const double zs=std::clamp(c.Z(),axes_[1].front(),axes_[1].back());
  return {c.X[0]/(hydrogen_share_ ? 1-zs : 1), zs, std::log10(Teff), std::log10(g)};
}
std::array<double, 5>
TraceHeliumAtmosphereGrid::interpolate(const std::vector<double> &f,
                                       const std::array<double, 4> &q) const {
  const auto selected = stencil(q);
  if (!selected)
    throw std::domain_error("TraceHeliumAtmosphereGrid: incomplete source stencil");
  const auto &base = *selected;
  std::array<double, 4> u{}, width{};
  for (std::size_t k = 0; k < q.size(); ++k) {
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
AtmosphereState TraceHeliumAtmosphereGrid::eval(double Teff, double g,
                                                const Composition &c) const {
  const auto q = coordinates(Teff, g, c);
  const auto t = interpolate_value(0, q), pg = interpolate_value(1, q);
  const double pr = constants::a_rad * std::pow(t[0], 4) / 3;
  AtmosphereState s{};
  s.T = t[0];
  s.Pgas = pg[0];
  s.P = pg[0] + pr;
  s.tau = tau_;
  if (!positive(s.P) || s.P == pr)
    throw std::domain_error(
        "TraceHeliumAtmosphereGrid: unrepresentable total pressure");
  s.rho = eos_.rho_from_PT(s.T, s.P, c);
  s.dlnT_dlnTeff = t[3];
  s.dlnT_dlng = t[4];
  s.dlnP_dlnTeff = (s.Pgas * pg[3] + 4 * pr * t[3]) / s.P;
  s.dlnP_dlng = (s.Pgas * pg[4] + 4 * pr * t[4]) / s.P;
  return s;
}
TraceHeliumAtmosphereGrid::CompositionResponse
TraceHeliumAtmosphereGrid::composition_response(double Teff, double g,
                                                const Composition &c) const {
  const auto q = coordinates(Teff, g, c);
  const auto t = interpolate_value(0, q), pg = interpolate_value(1, q);
  const double pr = constants::a_rad * std::pow(t[0], 4) / 3;
  const double p = pg[0] + pr;
  if (!positive(p) || p == pr)
    throw std::domain_error("TraceHeliumAtmosphereGrid: pressure overflow");
  const double dx=hydrogen_share_ ? 1/(1-q[1]) : 1;
  const double dz=hydrogen_share_ ? q[0]/(1-q[1]) : 0;
  const double tx=t[1]*dx, tz=t[2]+t[1]*dz;
  const double px=pg[1]*dx, pz=pg[2]+pg[1]*dz;
  return {tx, tz, (pg[0]*px+4*pr*tx)/p, (pg[0]*pz+4*pr*tz)/p};
}
} // namespace ember
