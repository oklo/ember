#include "ember/eos_cms19.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include "jet2.hpp"
#include <algorithm>
#include <fstream>
#include <limits>
#include <vector>

namespace ember {
namespace {
using detail::Jet2;
using detail::exp;
using detail::log;
constexpr double ln10 = 2.302585092994045684;
constexpr double unavailable = std::numeric_limits<double>::quiet_NaN();

// Interior monotone cubic Hermite. The four-node stencil always brackets
// the query with one extra node on each side; no endpoint extrapolation.
Jet2 cubic(const std::array<double, 4>& x, const std::array<Jet2, 4>& y, const Jet2& q) {
  std::array<double, 3> h{}; std::array<Jet2, 3> d{};
  for (int i = 0; i < 3; ++i) { h[i] = x[i + 1] - x[i]; d[i] = (y[i + 1] - y[i]) / h[i]; }
  auto slope = [&](int i) -> Jet2 {
    if (d[i].value * d[i + 1].value <= 0) return 0;
    const double w1 = 2 * h[i + 1] + h[i], w2 = h[i + 1] + 2 * h[i];
    return (w1 + w2) / (w1 / d[i] + w2 / d[i + 1]);
  };
  const auto u = (q - x[1]) / h[1], u2 = u * u, u3 = u2 * u;
  return (2 * u3 - 3 * u2 + 1) * y[1] + (u3 - 2 * u2 + u) * h[1] * slope(0)
       + (-2 * u3 + 3 * u2) * y[2] + (u3 - u2) * h[1] * slope(1);
}

struct Table {
  std::vector<double> t, p;
  struct Cell { double r, s; bool valid; };
  std::vector<Cell> cells;
  std::vector<std::pair<double, double>> bounds;
  std::string label;
  Table(const std::filesystem::path& file, bool helium) : label(helium ? "He" : "H") {
    std::ifstream input(file);
    std::size_t nt{}, np{}; input >> nt >> np;
    if (!input || nt != 121 || np != 441) throw std::runtime_error("CMS19: invalid table header");
    std::string title; std::getline(input, title);
    if (title.find(std::string("CMS19 ") + (helium ? "HE:" : "H:")) == std::string::npos)
      throw std::runtime_error("CMS19: source species label disagrees with constructor argument");
    t.resize(nt); p.resize(np);
    for (double& v : t) input >> v;
    for (double& v : p) input >> v;
    for (std::size_t i = 0; i < nt; ++i)
      if (!std::isfinite(t[i]) || std::abs(t[i] - (2 + .05 * static_cast<double>(i))) > 1e-10)
        throw std::runtime_error("CMS19: unexpected temperature axis");
    for (std::size_t j = 0; j < np; ++j)
      if (!std::isfinite(p[j]) || std::abs(p[j] - (-9 + .05 * static_cast<double>(j))) > 1e-10)
        throw std::runtime_error("CMS19: unexpected pressure axis");
    cells.resize(nt * np);
    for (std::size_t i = 0; i < nt; ++i) for (std::size_t j = 0; j < np; ++j) {
      auto& c = cells[i * np + j]; input >> c.r >> c.s;
      if (!input || !std::isfinite(c.r) || !std::isfinite(c.s))
        throw std::runtime_error("CMS19: truncated or nonfinite source cells");
      // Conservative stellar-fluid support, evaluated on every stencil node.
      // Source rectangular corners are NOT a physical validity statement.
      const double T = std::pow(10., t[i]), rho = std::pow(10., c.r);
      const double A = helium ? 4.0026 : 1.0078, Z = helium ? 2. : 1.;
      const double quantum = 7.71e3 * std::sqrt(rho) / T * Z * std::pow(A, -5. / 3.);
      const double Tm_ocp = std::pow(10., 3.1 + c.r / 3.) * Z * Z / std::cbrt(A);
      // The H2 melting curve peaks below the temperature subset used here.
      // For neutral He, apply the Simon curve; the ionized domain uses OCP.
      // The paper places closure of helium's electronic gap around 10 g/cm3;
      // extending the neutral Simon curve into that plasma would spuriously
      // exclude million-K stellar material. Use OCP in the ionized domain.
      const double Tm = helium && T < 1e6 && rho < 10
          ? 61 * std::pow(std::pow(10., p[j] + 1), .639) : Tm_ocp;
      c.valid = c.r >= -8 && c.r <= 4 && std::abs(c.s) < 30
             && quantum * quantum / 24 < .7 && T > Tm;
      c.r *= ln10; c.s = (c.s + 10) * ln10; // g/cm3 and erg/g/K
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("CMS19: trailing table tokens");
    for (double& v : t) v *= ln10;
    for (double& v : p) v = (v + 10) * ln10; // GPa -> dyn/cm2
    bounds.resize(nt - 1, {unavailable, unavailable});
    for (std::size_t i = 1; i + 2 < nt; ++i) {
      std::size_t first = np, last = 0;
      bool ended = false;
      for (std::size_t j = 1; j + 2 < np; ++j) {
        bool ok = true;
        for (std::size_t a = i - 1; a <= i + 2; ++a)
          for (std::size_t b = j - 1; b <= j + 2; ++b) ok &= cells[a * np + b].valid;
        // Retain a single contiguous physical pressure interval. A later
        // disconnected source island is not a second allowed branch.
        if (!ok) { if (first != np) ended = true; continue; }
        if (ended) continue;
        first = std::min(first, j); last = j;
      }
      if (first != np) bounds[i] = {p[first] + 1e-11, p[last + 1] - 1e-11};
    }
  }
  std::pair<double, double> range(double logT) const {
    if (!std::isfinite(logT) || logT < 3.2 * ln10 || logT > 7.3 * ln10)
      throw std::domain_error("CMS19: temperature outside supported stellar subset (logT=3.2..7.3)");
    const auto result = bounds[interp::locate(t, logT)];
    if (!std::isfinite(result.first)) throw std::domain_error("CMS19: no supported fluid pressure interval");
    return result;
  }
  std::pair<Jet2, Jet2> eval(const Jet2& logT, const Jet2& logP) const {
    const auto [lo, hi] = range(logT.value);
    if (!(logP.value >= lo && logP.value <= hi))
      throw std::domain_error("CMS19 " + label + ": pressure outside supported density/fluid stencil");
    const auto i = interp::locate(t, logT.value), j = interp::locate(p, logP.value);
    std::array<Jet2, 4> rr{}, ss{};
    std::array<double, 4> tt{}, pp{};
    for (std::size_t b = 0; b < 4; ++b) pp[b] = p[j - 1 + b];
    for (std::size_t a = 0; a < 4; ++a) {
      tt[a] = t[i - 1 + a];
      std::array<Jet2, 4> r{}, s{};
      for (std::size_t b = 0; b < 4; ++b) {
        const auto& c = cells[(i - 1 + a) * p.size() + j - 1 + b];
        if (!c.valid) throw std::domain_error("CMS19: unsupported physical stencil");
        r[b] = c.r; s[b] = c.s;
      }
      rr[a] = cubic(pp, r, logP); ss[a] = cubic(pp, s, logP);
    }
    return {cubic(tt, rr, logT), exp(cubic(tt, ss, logT))};
  }
};
} // namespace

struct Cms19Eos::Impl {
  Table h, he;
  Metals metals;
  Impl(const std::filesystem::path& H, const std::filesystem::path& He, Metals m) : h(H, false), he(He, true), metals(m) {}
  double hydrogen(const Composition& c) const {
    for (double x : c.X) if (!std::isfinite(x) || x < 0) throw std::domain_error("CMS19: invalid composition");
    if (std::abs(c.sum() - 1) > 1e-10 || c[Species::He3] != 0)
      throw std::domain_error("CMS19: normalized composition with He3=0 required");
    if (c.Z() > 1e-12 && metals != Metals::helium_proxy)
      throw std::domain_error("CMS19: metals require explicit helium-proxy approximation");
    return c.h1();
  }
  std::pair<double, double> range(double t, double X) const {
    auto a = X == 0 ? he.range(t) : h.range(t);
    if (X > 0 && X < 1) {
      const auto b = he.range(t); a = {std::max(a.first, b.first), std::min(a.second, b.second)};
    }
    if (!(a.first < a.second)) throw std::domain_error("CMS19: H/He support does not overlap");
    return a;
  }
  struct Material { Jet2 r, p, S; };
  Material material(double lt, double lp, double X) const {
    const auto t = Jet2::variable(lt, 0), p = Jet2::variable(lp, 1);
    Jet2 volume{}, entropy{};
    if (X > 0) { const auto [r, S] = h.eval(t, p); volume = volume + X * exp(-r); entropy = entropy + X * S; }
    if (X < 1) { const auto [r, S] = he.eval(t, p); volume = volume + (1 - X) * exp(-r); entropy = entropy + (1 - X) * S; }
    const double nH = X / 1.0078, nHe = (1 - X) / 4.0026, n = nH + nHe;
    double mixing = 0;
    for (double amount : {nH, nHe}) if (amount > 0) mixing -= constants::R_gas * amount * std::log(amount / n);
    // CMS19's ideal mixing entropy neglects free-electron mixing. This is
    // the published approximation, not a chemical species abundance model.
    const auto Pr = (constants::a_rad / 3) * exp(4 * t);
    return {-log(volume), log(exp(p) + Pr), entropy + mixing + 4 * Pr * volume / exp(t)};
  }
  double pressure(double T, double rho, double X) const {
    if (!std::isfinite(T) || !std::isfinite(rho) || !(T > 0) || !(rho > 0))
      throw std::domain_error("CMS19: positive finite T and rho required");
    const double t = std::log(T), r = std::log(rho);
    auto [lo, hi] = range(t, X);
    if (r < material(t, lo, X).r.value || r > material(t, hi, X).r.value)
      throw std::domain_error("CMS19: density outside supported fluid table at logT="
                              + std::to_string(t / ln10) + ", logrho=" + std::to_string(r / ln10));
    double p = std::clamp(r + t + std::log(constants::R_gas * (2 * X + .75 * (1 - X))), lo, hi);
    for (int iteration = 0; iteration < 80; ++iteration) {
      const auto a = material(t, p, X).r;
      if (!(a.d[1] > 0)) throw std::domain_error("CMS19: nonpositive compressibility");
      const double f = a.value - r;
      if (std::abs(f) < 2e-13) return p;
      if (f > 0) hi = p; else lo = p;
      const double next = p - f / a.d[1];
      p = next > lo && next < hi ? next : .5 * (lo + hi);
    }
    throw std::runtime_error("CMS19: pressure inversion did not converge");
  }
  EosResponse response(double T, double rho, double X) const {
    const auto m = material(std::log(T), pressure(T, rho, X), X);
    const auto rt = m.r.partial(0), rp = m.r.partial(1), pt = m.p.partial(0), pp = m.p.partial(1);
    const auto st = m.S.partial(0), sp = m.S.partial(1);
    const auto cp = st - sp * pt / pp;
    const auto cv = st - sp * rt / rp;
    const auto delta = -rt + rp * pt / pp;
    const auto grad = -sp / (pp * cp);
    EosResponse out{}; auto& e = out.state;
    e.P = std::exp(m.p.value); e.S = m.S.value; e.E = unavailable;
    e.mu = e.free_e = unavailable;
    e.chiRho = pp.value / rp.value; e.chiT = pt.value - pp.value * rt.value / rp.value;
    e.cp = cp.value; e.cv = cv.value; e.delta = delta.value; e.grad_ad = grad.value;
    e.Gamma1 = e.chiRho * e.cp / e.cv;
    for (double v : {e.P, e.S, e.cp, e.cv, e.delta, e.grad_ad, e.chiRho, e.chiT, e.Gamma1})
      if (!std::isfinite(v)) throw std::domain_error("CMS19: nonfinite thermodynamic response");
    if (!(e.P > 0 && e.cp > 0 && e.cv > 0 && e.delta > 0 && e.grad_ad > 0 && e.chiRho > 0))
      throw std::domain_error("CMS19: unstable or nonphysical interpolated entropy response");
    auto transform = [&](const auto& d, double& dT, double& dr) {
      dT = d.d[0] - d.d[1] * rt.value / rp.value; dr = d.d[1] / rp.value;
    };
    out.dE_dlnRho = unavailable;
    transform(cp, out.dcp_dlnT, out.dcp_dlnRho);
    transform(delta, out.ddelta_dlnT, out.ddelta_dlnRho);
    transform(grad, out.dgrad_ad_dlnT, out.dgrad_ad_dlnRho);
    return out;
  }
};

Cms19Eos::Cms19Eos(const std::filesystem::path& h, const std::filesystem::path& he, Metals m)
    : impl_(std::make_unique<Impl>(h, he, m)) {}
Cms19Eos::~Cms19Eos() = default;
EosState Cms19Eos::eval(double T, double rho, const Composition& c) const {
  return impl_->response(T, rho, impl_->hydrogen(c)).state;
}
EosResponse Cms19Eos::eval_with_derivatives(double T, double rho, const Composition& c) const {
  return impl_->response(T, rho, impl_->hydrogen(c));
}
double Cms19Eos::rho_from_PT(double T, double P, const Composition& c, double guess) const {
  if (!(T > 0 && P > 0) || !std::isfinite(T) || !std::isfinite(P) || !std::isfinite(guess) || guess < 0)
    throw std::domain_error("CMS19: invalid pressure inversion input");
  const double Pg = P - constants::a_rad * std::pow(T, 4) / 3;
  if (!(Pg > 32 * std::numeric_limits<double>::epsilon() * P))
    throw std::domain_error("CMS19: total pressure does not resolve positive gas pressure");
  return std::exp(impl_->material(std::log(T), std::log(Pg), impl_->hydrogen(c)).r.value);
}
std::optional<Eos::DensityRange> Cms19Eos::density_range(double T, const Composition& c) const {
  const double t = std::log(T), X = impl_->hydrogen(c);
  const auto [lo, hi] = impl_->range(t, X);
  return DensityRange{std::exp(impl_->material(t, lo, X).r.value), std::exp(impl_->material(t, hi, X).r.value)};
}
} // namespace ember
