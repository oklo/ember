#pragma once
#include "ember/atmosphere_gravity_interval.hpp"
#include "ember/atmosphere_helium_blend.hpp"
#include "ember/atmosphere_hydrogen_dominated.hpp"
#include "ember/atmosphere_metal_interval.hpp"
#include "ember/atmosphere_trace_helium.hpp"
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>

namespace ember {

// Composition and gravity intervals among sources sharing one optical depth.
// The full-composition source remains exact outside the declared overlaps.
class HydrogenEnvelopeAtmosphere final : public Atmosphere {
public:
  struct Specification {
    std::filesystem::path trace, low_gravity, middle_gravity, high_gravity;
    double maximum_helium3{}, maximum_helium{};
    std::array<double,2> gravity_low{}, gravity_high{}, trace_helium_join{},
        composition_helium_join{}, metal_join{}, pure_metal_join{};
  };

  static Specification read(const std::filesystem::path& path) {
    std::ifstream in(path); std::string label, filename; int version{};
    const auto expect=[&](const char* name) {
      if (!(in>>label) || label!=name)
        throw std::runtime_error("hydrogen envelope atmosphere: invalid manifest field");
    };
    const auto file=[&](const char* name) {
      expect(name);
      if (!(in>>std::quoted(filename)) || filename.empty())
        throw std::runtime_error("hydrogen envelope atmosphere: missing source file");
      return path.parent_path()/filename;
    };
    expect("EMBER_HYDROGEN_ENVELOPE_ATMOSPHERE");
    if (!(in>>version) || version!=1)
      throw std::runtime_error("hydrogen envelope atmosphere: invalid manifest version");
    Specification s;
    s.trace=file("trace_helium");
    s.low_gravity=file("low_gravity");
    s.middle_gravity=file("middle_gravity");
    s.high_gravity=file("high_gravity");
    expect("maximum_helium3"); in>>s.maximum_helium3;
    expect("maximum_helium"); in>>s.maximum_helium;
    const auto interval=[&](const char* name,std::array<double,2>& q) {
      expect(name); in>>q[0]>>q[1];
    };
    interval("gravity_low",s.gravity_low);
    interval("gravity_high",s.gravity_high);
    interval("trace_helium_join",s.trace_helium_join);
    interval("composition_helium_join",s.composition_helium_join);
    interval("metal_join",s.metal_join);
    interval("pure_metal_join",s.pure_metal_join);
    if (!in || !std::isfinite(s.maximum_helium+s.maximum_helium3))
      throw std::runtime_error("hydrogen envelope atmosphere: incomplete or invalid manifest");
    if (in>>label)
      throw std::runtime_error("hydrogen envelope atmosphere: trailing manifest data");
    return s;
  }

  HydrogenEnvelopeAtmosphere(const Eos& eos,const Atmosphere& full,
      const std::filesystem::path& path):HydrogenEnvelopeAtmosphere(eos,full,read(path)) {}

  HydrogenEnvelopeAtmosphere(const Eos& eos,const Atmosphere& full,const Specification& s)
      : trace_(eos,s.trace,TraceHeliumAtmosphereGrid::Approximation::neglect_atmospheric_helium3,s.maximum_helium3),
        low_(eos,s.low_gravity,HydrogenDominatedAtmosphereGrid::Approximation::neglect_trace_atmospheric_helium,s.maximum_helium),
        middle_(eos,s.middle_gravity,HydrogenDominatedAtmosphereGrid::Approximation::neglect_trace_atmospheric_helium,s.maximum_helium),
        high_(eos,s.high_gravity,HydrogenDominatedAtmosphereGrid::Approximation::neglect_trace_atmospheric_helium,s.maximum_helium),
        low_g_(eos,low_,middle_,s.gravity_low[0],s.gravity_low[1]),
        pure_(eos,low_g_,high_,s.gravity_high[0],s.gravity_high[1]),
        trace_join_(eos,pure_,trace_,s.trace_helium_join[0],s.trace_helium_join[1]),
        metal_join_(eos,trace_join_,full,s.metal_join[0],s.metal_join[1]),
        composition_join_(eos,metal_join_,full,s.composition_helium_join[0],s.composition_helium_join[1]),
        selected_(eos,pure_,composition_join_,s.pure_metal_join[0],s.pure_metal_join[1]) {
    if (s.gravity_high[0]<s.gravity_low[1] ||
        s.pure_metal_join[1]!=s.metal_join[0] ||
        low_.tau_match()!=middle_.tau_match() || low_.tau_match()!=high_.tau_match() ||
        low_.tau_match()!=trace_.tau_match())
      throw std::invalid_argument("hydrogen envelope atmosphere: incompatible source intervals");
  }

  AtmosphereState eval(double t,double g,const Composition& c) const override {
    return selected_.eval(t,g,c);
  }
  const char* name() const override {return "composition-dependent hydrogen-envelope atmosphere";}
private:
  TraceHeliumAtmosphereGrid trace_;
  HydrogenDominatedAtmosphereGrid low_, middle_, high_;
  GravityIntervalAtmosphere low_g_, pure_;
  HeliumBlendAtmosphere trace_join_;
  MetalIntervalAtmosphere metal_join_;
  HeliumBlendAtmosphere composition_join_;
  MetalIntervalAtmosphere selected_;
};
} // namespace ember
