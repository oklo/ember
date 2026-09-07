#pragma once
#include "ember/atmosphere.hpp"
#include <filesystem>
#include <istream>
#include <string>
#include <vector>

namespace ember {

// A rectangular (log10 Teff, log10 g) grid at one composition and optical
// depth. Stores log10 T and log10 Pgas, adding LTE radiation pressure before
// EOS inversion. Bilinear interpolation has exact cellwise derivatives and
// cannot overshoot the cell's corner values. Table edges and composition
// mismatches throw. There is no implicit grey filling or fallback.
// See docs/ATMOSPHERE.md for the versioned, whitespace-delimited format.
class TabulatedAtmosphere final : public Atmosphere {
public:
  TabulatedAtmosphere(const Eos&, const std::filesystem::path&);
  TabulatedAtmosphere(const Eos&, std::istream&);
  AtmosphereState eval(double Teff, double gravity, const Composition&) const override;
  bool covers(double Teff, double gravity, const Composition&) const;
  const char* name() const override { return source_.c_str(); }
  const Composition& composition() const { return composition_; }
  double tau_match() const { return tau_; }
private:
  void read(std::istream&);
  const Eos& eos_;
  std::string source_;
  Composition composition_{};
  double tau_{};
  std::vector<double> logTeff_, logg_, logT_, logPg_;
};

} // namespace ember
