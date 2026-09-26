#pragma once
#include "ember/composition.hpp"
#include <optional>
#include <limits>

namespace ember {

struct OpacityState {
  double kappa{};      // Rosseland mean, radiative + conductive  [cm^2/g]
  double dlnk_dlnT{};
  double dlnk_dlnRho{};
  double dlnk_dX{std::numeric_limits<double>::quiet_NaN()}; // fixed Z, He4 replaced by H1
  double dlnk_dZ{std::numeric_limits<double>::quiet_NaN()}; // fixed H1, scaled metals replace He4
  double dlnk_dY3{std::numeric_limits<double>::quiet_NaN()}; // fixed H1/metals, He3 replaces He4
};

class Opacity {
public:
  struct DensityRange { double min, max; };  // at a specified T and composition
  virtual ~Opacity() = default;
  virtual OpacityState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  // Wrappers must propagate this property. Microscopic heat transport supplies
  // its own conductivity and cannot be added to an already combined opacity.
  virtual bool includes_conduction() const {return false;}
  // Optional domain information for constrained solves (e.g. an atmosphere's
  // top pressure). Absence means no declared bounds, not permission to ignore
  // errors from eval(). A table implementation should expose its actual range.
  virtual std::optional<DensityRange> density_range(double, const Composition&) const {
    return std::nullopt;
  }
};

} // namespace ember
