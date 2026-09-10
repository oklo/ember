#pragma once
#include "ember/atmosphere_table.hpp"

namespace ember {
struct ConvectiveAtmosphereOptions {
  double tau_match{100}, tau_start{.01}, tau_top{.001};
  double alpha{1.9}, henyey_y{1.0 / 3.0};
  double tolerance{2e-8};
  double sensitivity_tolerance{2e-8}; // Jacobian integration, distinct from boundary values
  std::size_t max_steps{10000};
};

// Plane-parallel grey radiative/convective column. Henyey element cooling
// includes finite optical thickness; y=1/3 gives the interior MLT limit.
// Radiation uses the grey diffusion closure with an Eddington top. This is
// an explicit atmosphere approximation, not a non-grey model calculation.
class ConvectiveAtmosphere final : public Atmosphere {
public:
  ConvectiveAtmosphere(const Eos&, const Opacity&, ConvectiveAtmosphereOptions = {});
  AtmosphereState eval(double Teff, double gravity, const Composition&) const override;
  const char* name() const override { return "grey radiative/convective column; Henyey element cooling"; }

private:
  const Eos& eos_;
  const Opacity& opacity_;
  ConvectiveAtmosphereOptions options_;
};

// Differential composition correction: preserve the COND T and Pgas at
// a declared reference composition; multiply by column(c)/column(reference).
// This assumes the non-grey correction transfers between compositions.
// No new non-grey composition data are implied; compare with the raw column.
class CompositionCorrectedAtmosphere final : public Atmosphere {
public:
  CompositionCorrectedAtmosphere(const Eos&, const TabulatedAtmosphere&, const Atmosphere& column,
                                 Composition reference);
  AtmosphereState eval(double Teff, double gravity, const Composition&) const override;
  const char* name() const override { return "COND-anchored grey convective composition correction"; }

private:
  const Eos& eos_;
  const TabulatedAtmosphere& table_;
  const Atmosphere& column_;
  Composition reference_;
};
} // namespace ember
