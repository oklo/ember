#pragma once
#include "ember/eos_helmholtz.hpp"
#include <functional>

namespace ember {
// Dense liquid H/He/GS98 mixture: relativistic degenerate electrons through T^2,
// classical Coulomb ions, electron polarization and exchange-correlation,
// and the Baiko-Chugunov quantum-ion correction at a common electron density.
// Crystallization is not included. See docs/DENSE_EOS.md for domains and approximations.
struct ColdHeliumOptions {
  // Liquid host limit: helium coupling Gamma_He = 2^(5/3) Gamma_e stays below both assessed
  // melting references (OCP ~175; electron-screened He 141-147 at 1e4-1e5 g/cm3).
  // Trace metals are more strongly coupled but are not the host phase.
  double max_mixture_hydrogen=.05,max_mixture_metals=.16,max_helium_gamma=130;
  double minimum_density=1e3;
  double maximum_electron_temperature_ratio=.05;
  double maximum_hydrogen_quantum=4;
  // Dilute hydrogen (mass fraction <= trace_hydrogen) may reach a larger T_p,H/T: its per-ion
  // quantum model spread changes the native diffusion force by <=2.4e-3 (99th percentile)
  // and the local heat capacity by ~X (docs/DENSE_EOS.md).
  double trace_hydrogen=1e-3,maximum_trace_hydrogen_quantum=8;
  double join_cold=5e5,join_hot=8e5;
  double hydrogen_join_full=.01,hydrogen_join_zero=.05;
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
// Electron components for independent checks: 0 ideal, 1 exchange-correlation.
HelmholtzJet cold_helium_component_jet(int,double T,double rho,const Composition&,const ColdHeliumOptions& = {});
} // namespace ember
