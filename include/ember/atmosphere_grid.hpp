#pragma once
#include "ember/atmosphere.hpp"
#include <array>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>

namespace ember {

// Non-grey source states on (XH, XHe3, log10 Teff, log10 g), at fixed
// baryonic metal abundances and Rosseland matching depth. Interpolate log T
// and log Pgas, then invert the caller's EOS at the actual composition.
// Version 3 uses XHe3/(XHe3+XHe4) as its second table coordinate so a
// rectangular grid can include both hydrogen-rich and helium-rich mixtures.
// Version 4 permits a fixed He3 axis and explicit supported-cell masks. A fixed
// axis supplies no derivative or physical coverage away from its one value.
// A source grid's declared mixture/EOS/isotope approximations require an
// explicit opt-in. No extrapolation, abundance clipping or grey fallback.
class CompositionAtmosphereGrid final : public Atmosphere {
public:
  enum class Mixture { exact, allow_documented_proxy };
  CompositionAtmosphereGrid(const PressureDensity &, const std::filesystem::path &,
                            Mixture = Mixture::exact);
  CompositionAtmosphereGrid(const PressureDensity &, std::istream &,
                            Mixture = Mixture::exact);
  AtmosphereState eval(double Teff, double gravity,
                       const Composition &) const override;
  bool covers(double Teff, double gravity, const Composition &) const;
  const char *name() const override { return source_.c_str(); }
  const std::string &approximation() const { return approximation_; }
  double tau_match() const { return tau_; }
  const std::array<double,NMETALS>& reference_metals() const { return metals_; }
  double reference_metallicity() const {
    double z=0; for (double x : metals_) z+=x; return z;
  }
  struct Support {
    std::array<double, 2> hydrogen, helium3, teff, gravity;
  };
  // Outer bounds only: a version-2 table can contain missing source states.
  // covers() also checks the complete interpolation/derivative stencil.
  Support support() const;
  bool has_missing_states() const { return has_missing_states_; }
  // Require identical physics labels, composition/gravity axes and every old
  // source value/mask. Only higher-temperature rows may be appended. Returns
  // the number of added source states; throws for any other table change.
  std::size_t check_temperature_extension(const CompositionAtmosphereGrid &) const;
  // Derivatives with H1 or He3 replacing He4; same interpolant as eval().
  // A fixed He3 axis returns NaN for the unmeasured He3 derivatives.
  struct CompositionResponse {
    double dlnT_dXH, dlnT_dX3, dlnP_dXH, dlnP_dX3;
  };
  CompositionResponse composition_response(double Teff, double gravity,
                                           const Composition &) const;

private:
  void read(std::istream &, Mixture);
  std::optional<std::array<std::size_t, 4>>
  stencil(const std::array<double, 4> &) const;
  // Value (linear units), then dln(value)/d(XH, X3, lnTeff, lng).
  std::array<double, 5> interpolate(const std::vector<double> &,
                                    const std::array<double, 4> &) const;
  std::array<double, 4> coordinates(double, double, const Composition &) const;
  const PressureDensity &eos_;
  std::string source_, approximation_;
  double tau_{};
  std::array<double, NMETALS> metals_{};
  std::array<std::vector<double>, 4> axes_;
  std::vector<double> logT_, logPg_;
  std::vector<bool> valid_, cell_valid_;
  bool has_missing_states_{};
  bool helium_fraction_coordinates_{};
  double metal_tolerance_{1e-12};
};

} // namespace ember
