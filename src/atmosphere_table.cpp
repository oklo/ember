#include "ember/atmosphere_table.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace ember {
namespace {
bool positive(double x) { return std::isfinite(x) && x > 0.0; }
}

TabulatedAtmosphere::TabulatedAtmosphere(const Eos& eos, const std::filesystem::path& path, Mixture mixture)
    : eos_(eos), mixture_(mixture) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("TabulatedAtmosphere: cannot open " + path.string());
  read(in);
}
TabulatedAtmosphere::TabulatedAtmosphere(const Eos& eos, std::istream& in, Mixture mixture)
    : eos_(eos), mixture_(mixture) { read(in); }

void TabulatedAtmosphere::read(std::istream& in) {
  auto label = [&](const char* expected) {
    std::string word;
    if (!(in >> word) || word != expected)
      throw std::runtime_error(std::string("TabulatedAtmosphere: expected ") + expected);
  };
  label("EMBER_ATMOSPHERE");
  int version = 0; in >> version;
  if (!in || (version != 1 && version != 2)) throw std::runtime_error("TabulatedAtmosphere: unsupported version");
  label("source"); in >> std::quoted(source_);
  if (!in || source_.empty()) throw std::runtime_error("TabulatedAtmosphere: missing provenance");
  label("tau"); in >> tau_;
  if (!in || !positive(tau_)) throw std::runtime_error("TabulatedAtmosphere: invalid optical depth");
  if (version == 2) {
    label("composition_proxy"); in >> std::quoted(composition_proxy_);
    if (!in || composition_proxy_.empty())
      throw std::runtime_error("TabulatedAtmosphere: missing composition approximation");
    if (mixture_ != Mixture::allow_documented_proxy)
      throw std::invalid_argument("TabulatedAtmosphere: explicit solar-mixture proxy selection required");
  }
  label("composition");
  for (double& x : composition_.X) {
    in >> x;
    if (!in || !std::isfinite(x) || x < 0.0 || x > 1.0)
      throw std::runtime_error("TabulatedAtmosphere: invalid composition");
  }
  if (std::abs(composition_.sum() - 1.0) > 1e-10)
    throw std::runtime_error("TabulatedAtmosphere: composition must sum to one");
  auto axis = [&](const char* key, std::vector<double>& values) {
    label(key);
    std::size_t n = 0; in >> n;
    if (!in || n < 2 || n > 10000)
      throw std::runtime_error("TabulatedAtmosphere: invalid axis size");
    values.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      in >> values[i];
      if (!in || !std::isfinite(values[i]) || !positive(std::pow(10.0, values[i]))
          || (i > 0 && !(values[i] > values[i - 1])))
        throw std::runtime_error("TabulatedAtmosphere: axes must be finite and strictly increasing");
    }
  };
  axis("log_teff", logTeff_); axis("log_g", logg_);
  const std::size_t n = logTeff_.size() * logg_.size();
  if (n > 1000000) throw std::runtime_error("TabulatedAtmosphere: grid is too large");
  logT_.resize(n); logPg_.resize(n);
  label("data");
  for (std::size_t i = 0; i < n; ++i) {
    in >> logT_[i] >> logPg_[i];
    if (!in || !positive(std::pow(10.0, logT_[i])) || !positive(std::pow(10.0, logPg_[i])))
      throw std::runtime_error("TabulatedAtmosphere: invalid or truncated data");
  }
  std::string extra;
  if (in >> extra) throw std::runtime_error("TabulatedAtmosphere: unexpected trailing data");
}

bool TabulatedAtmosphere::covers(double Teff, double gravity, const Composition& comp) const {
  if (!positive(Teff) || !positive(gravity)) return false;
  for (std::size_t i = 0; i < NSPEC; ++i)
    if (!std::isfinite(comp.X[i]) || comp.X[i] < 0.0
        || std::abs(comp.X[i] - composition_.X[i]) > 1e-10) return false;
  const double t = std::log10(Teff), g = std::log10(gravity);
  return t >= logTeff_.front() && t <= logTeff_.back() && g >= logg_.front() && g <= logg_.back();
}

AtmosphereState TabulatedAtmosphere::eval(double Teff, double gravity, const Composition& comp) const {
  if (!covers(Teff, gravity, comp))
    throw std::domain_error("TabulatedAtmosphere: Teff, gravity, or composition outside table");
  const double t = std::log10(Teff), g = std::log10(gravity);
  const std::size_t it = interp::locate(logTeff_, t), ig = interp::locate(logg_, g);
  const double dt = logTeff_[it + 1] - logTeff_[it], dg = logg_[ig + 1] - logg_[ig];
  const double u = (t - logTeff_[it]) / dt, v = (g - logg_[ig]) / dg;
  auto interpolate = [&](const std::vector<double>& field) {
    const std::size_t j = it * logg_.size() + ig;
    const double f00 = field[j], f01 = field[j + 1];
    const double f10 = field[j + logg_.size()], f11 = field[j + logg_.size() + 1];
    const double f = (1.0 - u) * ((1.0 - v) * f00 + v * f01)
                   + u * ((1.0 - v) * f10 + v * f11);
    const double dT = ((1.0 - v) * (f10 - f00) + v * (f11 - f01)) / dt;
    const double dG = ((1.0 - u) * (f01 - f00) + u * (f11 - f10)) / dg;
    return std::array{std::pow(10.0, f), dT, dG};
  };
  const auto temp = interpolate(logT_), gas = interpolate(logPg_);
  const double Pr = constants::a_rad * std::pow(temp[0], 4) / 3.0;
  AtmosphereState s{};
  s.T = temp[0]; s.Pgas = gas[0]; s.P = s.Pgas + Pr; s.tau = tau_;
  if (!positive(s.P) || s.P == Pr)
    throw std::domain_error("TabulatedAtmosphere: total pressure is not representable");
  s.rho = eos_.rho_from_PT(s.T, s.P, comp);
  s.dlnT_dlnTeff = temp[1]; s.dlnT_dlng = temp[2];
  s.dlnP_dlnTeff = (s.Pgas / s.P) * gas[1] + (4.0 * Pr / s.P) * temp[1];
  s.dlnP_dlng = (s.Pgas / s.P) * gas[2] + (4.0 * Pr / s.P) * temp[2];
  return s;
}

} // namespace ember
