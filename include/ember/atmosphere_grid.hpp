#pragma once
#include "ember/atmosphere.hpp"
#include <array>
#include <filesystem>
#include <istream>
#include <string>
#include <vector>

namespace ember {

// Non-grey source states on (XH, XHe3, log10 Teff, log10 g), at fixed
// baryonic metal abundances and Rosseland matching depth. Interpolate log T
// and log Pgas, then invert the caller's EOS at the actual composition.
// A source grid's declared mixture/EOS/isotope approximations require an
// explicit opt-in. No extrapolation, abundance clipping or grey fallback.
class CompositionAtmosphereGrid final : public Atmosphere {
public:
  enum class Mixture { exact, allow_documented_proxy };
  CompositionAtmosphereGrid(const Eos &, const std::filesystem::path &,
                            Mixture = Mixture::exact);
  CompositionAtmosphereGrid(const Eos &, std::istream &,
                            Mixture = Mixture::exact);
  AtmosphereState eval(double Teff, double gravity,
                       const Composition &) const override;
  bool covers(double Teff, double gravity, const Composition &) const;
  const char *name() const override { return source_.c_str(); }
  const std::string &approximation() const { return approximation_; }
  double tau_match() const { return tau_; }
  struct Support {
    std::array<double, 2> hydrogen, helium3, teff, gravity;
  };
  Support support() const;
  // Derivatives with H1 or He3 replacing He4; same interpolant as eval().
  struct CompositionResponse {
    double dlnT_dXH, dlnT_dX3, dlnP_dXH, dlnP_dX3;
  };
  CompositionResponse composition_response(double Teff, double gravity,
                                           const Composition &) const;

private:
  void read(std::istream &, Mixture);
  // Value (linear units), then dln(value)/d(XH, X3, lnTeff, lng).
  std::array<double, 5> interpolate(const std::vector<double> &,
                                    const std::array<double, 4> &) const;
  std::array<double, 4> coordinates(double, double, const Composition &) const;
  const Eos &eos_;
  std::string source_, approximation_;
  double tau_{};
  std::array<double, NSPEC - 3> metals_{};
  std::array<std::vector<double>, 4> axes_;
  std::vector<double> logT_, logPg_;
};

} // namespace ember
