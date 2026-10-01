#pragma once
#include <array>

namespace ember {
struct LiquidViscosity {
  double electron{}, ion_classical{}; // dynamic viscosity, g/(cm s)
  double electron_ion_frequency{}, electron_electron_frequency{}; // 1/s
  double coulomb_log{}, gamma{}, ion_quantum_parameter{}, temperature_over_fermi{};
};
struct LiquidViscosityResponse {
  LiquidViscosity value;
  // Rows: electron and classical-ion viscosity. Columns: ln T, ln rho,
  // ln ion charge, ln mass number. Absolute derivatives, not logarithmic.
  std::array<std::array<double,4>,2> partials{};
};

// Fully ionized, nonmagnetic one-component liquid with degenerate electrons.
// Electron-ion viscosity: Chugunov & Yakovlev (2005), effective-potential fit.
// Electron-electron term: their nonrelativistic expression; limited here to
// p_F/(m_e c)<=0.5. Ion term: classical OCP fit of Daligault et al. (2014).
// Requires T/T_F<=0.05 and Gamma<=175. The ion term is explicitly classical;
// a returned ion_quantum_parameter near unity requires a sensitivity check.
// No ionization, mixture rule, magnetic field or crystal rheology is inferred.
LiquidViscosityResponse ocp_liquid_viscosity(double temperature,double density,
                                          double charge,double mass_number);
} // namespace ember
