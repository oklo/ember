#pragma once
#include "ember/eos.hpp"

namespace ember {

// Local, optically thick Bohm-Vitense MLT (a=1/8, b=1/2, c=24), using the
// Schwarzschild criterion. Composition mixing and Ledoux terms are separate.
struct ConvectionState {
  bool unstable{};
  double grad{};                     // background dlnT/dlnP
  double grad_element{};             // rising element's dlnT/dlnP
  // Keep the small differences explicitly: subtracting the full gradients
  // loses them in nearly adiabatic or nearly radiative layers.
  double superadiabaticity{};         // grad - grad_ad (zero in stable layers)
  double element_contrast{};          // grad - grad_element
  double convective_excess{};         // grad_rad - grad; proportional to F_conv

  // Independent partials at fixed other arguments, for Jacobian assembly.
  double dgrad_dgrad_rad{1.0};
  double dgrad_dgrad_ad{};
  double dgrad_dlnU{};
};

// U is the dimensionless radiative-loss parameter, finite and strictly
// positive. Small U means efficient convection. Negative grad_rad is allowed
// (inward luminosity); grad_ad must be finite and non-negative.
ConvectionState mixing_length_gradient(double grad_rad, double grad_ad, double U);

// Finite-face composition buoyancy B. Both compositions are compared at the
// SAME midpoint T/P; thermal and density gradients do not masquerade as a
// composition gradient. The ideal gas limit is B=dln(mu)/dln(P).
struct CompositionBuoyancy {
  double B{},dB_dlnT{},dB_dlnP{},dB_ddelta{},dB_dpressure_contrast{};
};
CompositionBuoyancy composition_buoyancy(const Eos&,double T,double P,double delta,
    double pressure_contrast,const Composition& lo,const Composition& hi,double rho_guess=0);

// Cox/Giuli Ledoux MLT: replace the buoyancy-neutral gradient by grad_ad+B.
// The derivative with respect to B equals dgrad_dgrad_ad. Stable layers keep
// the diffusive gradient; semi-convection/thermohaline transport are separate.
ConvectionState ledoux_mixing_length_gradient(double grad_rad,double grad_ad,double B,double U);

// CGS inputs. H_P = P/(rho*g), l = alpha*H_P, and
// U = 3*a_rad*c*T^3/(cp*rho^2*kappa*l^2) * sqrt(8*H_P/(g*delta)).
// Requires positive finite inputs; evaluate away from the central singularity
// g=0. This diffusion treatment needs an optically thick interior; an
// optically thin atmosphere requires a different radiative-loss prescription.
double mixing_length_U(double T, double rho, double kappa, double gravity,
                       const EosState& eos, double alpha);

} // namespace ember
