#pragma once
#include "ember/atmosphere.hpp"
#include <array>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>
namespace ember {
// Actual zero-He3 gas sources on (XH, Z, log10 Teff, log10 g). The explicit
// approximation neglects atmospheric isotope dependence only; density is
// inverted with the caller's actual stellar composition. No extrapolation
// or synthetic isotope source rows are used.
class TraceHeliumAtmosphereGrid final : public Atmosphere {
public:
  enum class Approximation { neglect_atmospheric_helium3 };
  TraceHeliumAtmosphereGrid(const PressureDensity&, const std::filesystem::path&, Approximation,
                           double maximum_helium3);
  TraceHeliumAtmosphereGrid(const PressureDensity&, std::istream&, Approximation,
                           double maximum_helium3);
  AtmosphereState eval(double, double, const Composition&) const override;
  bool covers(double, double, const Composition&) const;
  const char* name() const override { return source_.c_str(); }
  const std::string& approximation() const { return approximation_; }
  double tau_match() const { return tau_; }
  struct Support { std::array<double,2> hydrogen, metals, teff, gravity; };
  Support support() const;
  bool has_missing_states() const { return has_missing_states_; }
  struct CompositionResponse { double dlnT_dXH, dlnT_dZ, dlnP_dXH, dlnP_dZ; };
  CompositionResponse composition_response(double, double, const Composition&) const;
private:
  void read(std::istream&, Approximation, double);
  std::optional<std::array<std::size_t,4>> stencil(const std::array<double,4>&) const;
  std::array<double,5> interpolate(const std::vector<double>&,const std::array<double,4>&) const;
  std::array<double,4> coordinates(double,double,const Composition&) const;
  const PressureDensity& eos_;
  std::string source_, approximation_;
  double tau_{}, maximum_helium3_{};
  std::array<double,NMETALS> metal_pattern_{};
  std::array<std::vector<double>,4> axes_;
  std::vector<double> logT_, logPg_;
  std::vector<bool> valid_;
  bool has_missing_states_{};
};
} // namespace ember
