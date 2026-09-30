#pragma once
#include "ember/atmosphere.hpp"
#include <array>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>
namespace ember {
// Actual zero-He3 gas sources on (XH, Z, log10 Teff, log10 g), or with
// hydrogen share XH/(1-Z) in versions 2 and 3. Version 3 stores irregular
// source points and their simplices. The explicit
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
  // Bounding box only; hydrogen bounds are physical mass fractions.
  // covers() checks the actual joint domain, including irregular hulls.
  struct Support { std::array<double,2> hydrogen, metals, teff, gravity; };
  Support support() const;
  bool has_missing_states() const { return irregular_ || has_missing_states_; }
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
  bool hydrogen_share_{}; // Version 2 stores X/(1-Z), keeping every metal anchor physical.
  // Version 3 (irregular): accepted source vertices in (X/(1-Z), Z, log10 Teff, log10 g) with offline Delaunay
  // simplices. Values are piecewise linear in the affinely scaled coordinates; each simplex has a constant gradient.
  // Queries outside every stored simplex are rejected. Old rectangular formats are unchanged.
  bool irregular_{};
  std::array<double,4> offset_{}, width_{};
  std::vector<std::array<double,4>> vertex_;          // scaled coordinates
  std::vector<std::array<double,2>> vertex_value_;    // log10 T, log10 Pgas
  std::vector<std::array<std::size_t,5>> simplex_;
  std::vector<std::array<double,16>> inverse_;        // inverse of [v1-v0 ... v4-v0] (row-major)
  std::vector<std::array<double,8>> gradient_;        // d log10(T,Pg) / d raw coordinate, per simplex
  std::optional<std::size_t> locate(const std::array<double,4>& raw) const;
  std::array<double,5> interpolate_irregular(std::size_t value,const std::array<double,4>& raw) const;
  std::array<double,5> interpolate_value(std::size_t which,const std::array<double,4>& q) const;
};
} // namespace ember
