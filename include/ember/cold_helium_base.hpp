#pragma once
#include "ember/eos_helmholtz.hpp"
#include <functional>

namespace ember {
// Dense liquid H/He/GS98 mixture: relativistic degenerate electrons through T^2,
// classical Coulomb ions, electron polarization and exchange-correlation,
// and the Baiko-Chugunov quantum-ion correction at a common electron density.
// Crystallization is optional and assumes no element fractionation.
// See docs/DENSE_EOS.md for domains and approximations.
struct ColdHeliumOptions {
  // Optional metal-depleted, fully ionized H/He transition below 300 kK.
  // The existing helium-core join and phase response are unchanged.
  bool dense_transition=false;
  // Default liquid-only host limit. Metals can move a mixture phase boundary
  // above the pure-helium melting temperature; this is not a phase criterion.
  double max_mixture_hydrogen=.20,max_mixture_metals=.16,max_helium_gamma=130;
  double minimum_density=1e3;
  double maximum_electron_temperature_ratio=.05;
  double maximum_hydrogen_quantum=4;
  // Dilute hydrogen (mass fraction <= trace_hydrogen) may reach a larger T_p,H/T: its per-ion
  // quantum model spread changes the native diffusion force by <=3.8e-3 (99th percentile)
  // and the local heat capacity by ~X (docs/DENSE_EOS.md).
  double trace_hydrogen=1e-3,maximum_trace_hydrogen_quantum=8;
  // Further dilute-H range assessed with the optional mixture phase model.
  double ultratrace_hydrogen=1e-6,maximum_ultratrace_hydrogen_quantum=12;
  double join_cold=5e5,join_hot=8e5;
  // The narrow X=0.01--0.05 overlap introduces a non-convex composition
  // potential between the differently aligned sources. This broader overlap
  // retains positive H/He3 curvature in the assessed liquid transition.
  double hydrogen_join_full=.005,hydrogen_join_zero=.20;
  // Optional same-composition crystallization (ion_phase.hpp, ion_mixture_phase_difference_jets):
  // a soft minimum of the liquid and a solid solution of the per-species bcc branches, width
  // phase_width k per ion. It replaces max_helium_gamma by the supercooled-liquid limit below and adds
  // nothing where Gamma_He < minimum_solid_gamma. Metal fractionation and precipitation are not included.
  // Classical-liquid continuation coupling (detail/ion_ocp_components.hpp); 0 keeps FITION9 at all couplings.
  // Applies to every species' liquid in both the base and the phase difference.
  double liquid_continuation_gamma=0;
  bool mixture_phase=false;
  double phase_width=.005,minimum_solid_gamma=90,max_phase_helium_gamma=250;
};
// Material F/T and its thermal, density and composition derivatives. Radiation and
// ideal mixing are supplied once by the surrounding EOS.
using ColdHeliumTable=std::function<std::array<HelmholtzJet,10>(double,double,const Composition&,std::size_t)>;
std::array<HelmholtzJet,10> cold_helium_material_jets(double T,double rho,const Composition&,
    std::size_t channels=10,const ColdHeliumOptions& = {});
// G=a(rho,X)+b(rho,X)/T matches the two potentials at both temperature anchors.
// Its density/composition dependence carries model differences below the join;
// it is not merely an energy-zero convention. The cold heat capacity is retained.
std::array<HelmholtzJet,10> cold_helium_alignment_jets(const ColdHeliumTable&,double T,double rho,
    const Composition&,std::size_t channels=10,const ColdHeliumOptions& = {});
std::array<HelmholtzJet,10> cold_helium_joined_jets(const ColdHeliumTable&,double T,double rho,
    const Composition&,std::size_t channels=10,const ColdHeliumOptions& = {});
void cold_helium_validate(double T,double rho,const Composition&,const ColdHeliumOptions& = {});
// Product of C2 temperature and hydrogen weights; zero returns the original EOS.
double cold_helium_join_weight(double T,const Composition&,const ColdHeliumOptions& = {});
// H/He transition: temperature weight 200--300 kK, density weight 300--600 g/cm3.
// Valid at 100 kK or warmer, rho <= 6000, Z <= 1e-8, and T_p,H/T <= 2.5.
// Both source anchors and any fractionally weighted source must remain supported.
double dense_hhe_transition_weight(double T,double rho);
void dense_hhe_transition_validate(double T,double rho,const Composition&,const ColdHeliumOptions&);
std::array<HelmholtzJet,10> dense_hhe_transition_jets(const ColdHeliumTable&,double T,double rho,
    const Composition&,std::size_t channels,const ColdHeliumOptions&);
// Solid weight including the existing temperature/composition joins, followed
// by derivatives in ln T, ln rho, XH, X3 and Z. Zero when the phase is off.
std::array<double,6> cold_helium_solid_response(double T,double rho,const Composition&,
    const ColdHeliumOptions& = {});
// Electron components for independent checks: 0 ideal, 1 exchange-correlation.
HelmholtzJet cold_helium_component_jet(int,double T,double rho,const Composition&,const ColdHeliumOptions& = {});
} // namespace ember
