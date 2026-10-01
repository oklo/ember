#pragma once
#include <array>

namespace ember {

enum class FingeringRegime { stable, overturning, fingering };

struct FingeringFlux {
  FingeringRegime regime{FingeringRegime::stable};
  double growth_rate{};            // dimensionless lambda
  double wavenumber_squared{};     // dimensionless l^2
  double thermal_nusselt_excess{}; // Nu_T - 1; retain small excess explicitly
  double chemical_nusselt_excess{};// Nu_mu - 1, excludes microscopic diffusion
  bool within_calibrated_ratio{};  // tau <= Pr; fit is poorer outside this range
};

// Brown, Garaud & Stellmach (2013), equations 19, 20, 32 and 33, C=7.
// https://arxiv.org/abs/1212.1688
// Pr = viscosity / thermal_diffusivity; tau = chemical / thermal diffusivity.
// R0 is the stabilizing thermal buoyancy / destabilizing composition buoyancy.
// CALLER must first establish a stable thermal gradient and an inverse
// composition gradient. Positive R0 alone does not establish those signs.
// Requires finite Pr>0, 0<tau<1, R0>0. R0<=1 is returned as overturning;
// R0>=1/tau is stable. Neither case is assigned a fingering diffusivity.
//
// Additional composition diffusivity is kappa_mu*(Nu_mu-1). Adding Nu_mu
// instead would count microscopic diffusion twice. Thermal transport is
// returned too; it must be assessed using the superadiabatic gradient.
// This local, nonrotating, nonmagnetic homogeneous-fingering model does not
// prescribe convective boundary mixing, layers, or semiconvection.
FingeringFlux brown_fingering_flux(double prandtl, double diffusivity_ratio,
                                 double density_ratio);

struct FingeringResponse {
  FingeringFlux value;
  // Rows: lambda, l^2, Nu_T-1, Nu_mu-1. Columns: ln Pr, ln tau, ln R0.
  // These are absolute derivatives with respect to logarithmic parameters.
  std::array<std::array<double,3>,4> partials{};
  // Only the open fingering regime has these derivatives. At a stability
  // boundary or in another regime, zero storage is not a derivative claim.
  bool derivatives_defined{};
};

// Analytic differentiation of the dispersion relation and its stationary
// fastest-mode condition. Intended for a coupled composition/heat response;
// this component does not by itself select a stellar transport prescription.
FingeringResponse brown_fingering_response(double prandtl, double diffusivity_ratio,
                                         double density_ratio);

struct TwoCompositionFingering {
  FingeringRegime regime{FingeringRegime::stable};
  double growth_rate{},wavenumber_squared{},thermal_nusselt_excess{};
  // Extra diffusivity / thermal diffusivity, in the two chemical diffusion
  // eigenmodes. Does not include microscopic diffusion.
  std::array<double,2> mixing_over_thermal{};
};

// Two independently diffusing composition fields. Positive driving is
// destabilizing; sum(driving) is inverse R0. Diffusivity ratios are ordered
// from slow to fast, with 0<tau[0]<=tau[1]<1. A destabilizing slow field and
// a stabilizing fast field are supported even if their net buoyancy is stable.
// For the opposite ordering, both spectral branches are compared; a dominant
// oscillatory mode is explicitly refused. Buoyancy signs alone cannot classify it.
// The linear dispersion relation retains both fields exactly. Applying the
// Brown C=7 saturation hypothesis to it is an UNCALIBRATED extension, not a
// published multicomponent transport fit. Equal diffusivities recover Brown.
TwoCompositionFingering two_composition_fingering(double prandtl,
    std::array<double,2> diffusivity_ratios,std::array<double,2> driving);

} // namespace ember
