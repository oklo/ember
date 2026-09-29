#pragma once
#include "ember/atmosphere_grid.hpp"
#include "ember/atmosphere_metal_chain.hpp"
#include "ember/constants.hpp"
#include <fstream>
#include <iomanip>
#include <memory>
#include <cmath>

namespace ember {

// Use neighboring physical hydrogen coordinates from two source grids in one
// continuous interpolation. Outside the interval, return the corresponding
// grid exactly. No source temperature, gravity or composition is extrapolated.
class HydrogenIntervalAtmosphere final : public Atmosphere {
public:
  HydrogenIntervalAtmosphere(const PressureDensity& eos, const Atmosphere& lower,
      const std::filesystem::path& path, unsigned depth=0) : eos_(eos), lower_(&lower) {
    if (depth > 8)
      throw std::runtime_error("HydrogenIntervalAtmosphere: excessive interval nesting");
    std::ifstream in(path);
    std::string label, reference, chain;
    int version{};
    double z{};
    const auto expect = [&](const char* wanted) {
      if (!(in >> label) || label != wanted)
        throw std::runtime_error("HydrogenIntervalAtmosphere: invalid manifest field");
    };
    expect("EMBER_HYDROGEN_ATMOSPHERE_INTERVAL");
    if (!(in >> version) || (version != 1 && version != 2))
      throw std::runtime_error("HydrogenIntervalAtmosphere: invalid version");
    if (version == 2) {
      std::string lower_interval, lower_chain;
      double lower_z{};
      expect("lower_interval"); in >> std::quoted(lower_interval);
      expect("lower_metal_chain"); in >> std::quoted(lower_chain);
      expect("lower_reference_Z"); in >> lower_z;
      if (!in || lower_interval.empty() || lower_chain.empty() ||
          !std::isfinite(lower_z) || lower_z<=0 || lower_z>=1)
        throw std::runtime_error("HydrogenIntervalAtmosphere: invalid retained lower boundary");
      retained_lower_ = std::make_unique<HydrogenIntervalAtmosphere>(
          eos_,lower,path.parent_path()/lower_interval,depth+1);
      continued_lower_ = std::make_unique<MetalAtmosphereChain>(
          eos_,*retained_lower_,path.parent_path()/lower_chain,lower_z);
      lower_ = continued_lower_.get();
    }
    expect("reference"); in >> std::quoted(reference);
    expect("chain"); in >> std::quoted(chain);
    expect("reference_Z"); in >> z;
    expect("hydrogen"); in >> low_ >> high_;
    if (!in || !std::isfinite(z+low_+high_) || z <= 0 || z >= 1 ||
        low_ < 0 || low_ >= high_ || high_ >= 1 || reference.empty() || chain.empty())
      throw std::runtime_error("HydrogenIntervalAtmosphere: invalid source interval");
    if (retained_lower_ && low_ < retained_lower_->high_)
      throw std::runtime_error("HydrogenIntervalAtmosphere: overlapping hydrogen intervals");
    if (in >> label)
      throw std::runtime_error("HydrogenIntervalAtmosphere: trailing manifest data");
    reference_ = std::make_unique<CompositionAtmosphereGrid>(eos_,
        path.parent_path()/reference, CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    upper_ = std::make_unique<MetalAtmosphereChain>(eos_,*reference_,path.parent_path()/chain,z);
  }

  AtmosphereState eval(double teff, double gravity, const Composition& c) const override {
    if (c.basis != AbundanceBasis::baryon_mass || !std::isfinite(c.h1()))
      throw std::domain_error("HydrogenIntervalAtmosphere: physical composition required");
    if (c.h1() <= low_) return lower_->eval(teff,gravity,c);
    if (c.h1() >= high_) return upper_->eval(teff,gravity,c);
    auto a_c = c, b_c = c;
    a_c.X[2] += c.h1()-low_; a_c.X[0] = low_;
    b_c.X[2] += c.h1()-high_; b_c.X[0] = high_;
    const auto a = lower_->eval(teff,gravity,a_c), b = upper_->eval(teff,gravity,b_c);
    if (a.tau != b.tau || !(a.Pgas > 0) || !(b.Pgas > 0))
      throw std::domain_error("HydrogenIntervalAtmosphere: inconsistent source boundaries");
    const double f = (c.h1()-low_)/(high_-low_);
    const auto blend = [&](double x,double y) { return (1-f)*x+f*y; };
    AtmosphereState out{};
    out.T = std::exp(blend(std::log(a.T),std::log(b.T)));
    out.Pgas = std::exp(blend(std::log(a.Pgas),std::log(b.Pgas)));
    out.tau = a.tau;
    const double ar = constants::a_rad*std::pow(a.T,4)/3;
    const double br = constants::a_rad*std::pow(b.T,4)/3;
    const double pr = constants::a_rad*std::pow(out.T,4)/3;
    out.P = out.Pgas+pr;
    out.dlnT_dlnTeff = blend(a.dlnT_dlnTeff,b.dlnT_dlnTeff);
    out.dlnT_dlng = blend(a.dlnT_dlng,b.dlnT_dlng);
    const auto pressure_derivative = [&](double adp,double adt,double bdp,double bdt,double dt) {
      const double dpg = blend((a.P*adp-4*ar*adt)/a.Pgas,(b.P*bdp-4*br*bdt)/b.Pgas);
      return (out.Pgas*dpg+4*pr*dt)/out.P;
    };
    out.dlnP_dlnTeff = pressure_derivative(a.dlnP_dlnTeff,a.dlnT_dlnTeff,
        b.dlnP_dlnTeff,b.dlnT_dlnTeff,out.dlnT_dlnTeff);
    out.dlnP_dlng = pressure_derivative(a.dlnP_dlng,a.dlnT_dlng,
        b.dlnP_dlng,b.dlnT_dlng,out.dlnT_dlng);
    out.rho = eos_.rho_from_PT(out.T,out.P,c,a.rho);
    return out;
  }
  const char* name() const override { return "continuous hydrogen-source atmosphere interpolation"; }
private:
  const PressureDensity& eos_;
  const Atmosphere* lower_;
  std::unique_ptr<HydrogenIntervalAtmosphere> retained_lower_;
  std::unique_ptr<MetalAtmosphereChain> continued_lower_;
  std::unique_ptr<CompositionAtmosphereGrid> reference_;
  std::unique_ptr<MetalAtmosphereChain> upper_;
  double low_{}, high_{};
};
} // namespace ember
