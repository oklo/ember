#pragma once
#include "ember/eos_helmholtz.hpp"
#include "ember/composition.hpp"

namespace ember {
// Fully ionized one-component-plasma ion free energies for a crystal and a
// liquid, from the Ioffe fits in the pinned eos22.f (not selected by any EOS).
//   liquid: FITION9 (classical OCP) + ideal ions + LIQUBC (Baiko & Chugunov
//           2022 quantum liquid, with density-dependent coefficients);
//   solid:  FHARM12 (quantum harmonic bcc lattice, HLfit12 with zero point)
//           + ANHBC (Baiko & Chugunov 2022 quantum anharmonic crystal).
// Optionally each branch adds its separately fitted electron-ion screening
// term (FSCRliq8 / FSCRsol8). Electron exchange-correlation is identical in
// both phases and is NOT included here (the base electron model must supply it once). Only the free energy is ported; every
// response differentiates it through total degree three (Taylor3).
enum class IonScreening { ocp, fitted_screening };

// F/T per unit mass of a pure ion species (mass number A, charge Z) at
// temperature exp(logT) and electron density exp(logne) [cm^-3].
// Entries [i][j] are d^i/dlnT^i d^j/dln(ne)^j, i+j<=3.
HelmholtzJet ion_ocp_liquid_jet(double logT,double logne,double A,double Z,IonScreening);
HelmholtzJet ion_ocp_solid_jet(double logT,double logne,double A,double Z,IonScreening);

// Phase correction relative to the liquid for He4, using a soft minimum of
// the per-ion potentials f=F/(N_i k T):
//   f = -c ln(exp(-f_l/c) + exp(-f_s/c)),   returned as f - f_l per mass.
// The soft minimum preserves concavity of F/T in 1/T, so it cannot create a
// negative heat capacity; its only artefact is a bounded -c ln2 offset at the
// melting point. Other species (H, He3, metals) carry no phase term and are
// linearly mixed at a common electron density, like ion_quantum_liquid_jets.
// This is a reference for nearly pure He only: it throws unless the non-He4
// mass fraction is <= max_non_helium. It is not a mixture phase diagram and
// is not added to any EOS by default.
struct IonPhaseOptions {
  IonScreening screening=IonScreening::fitted_screening;
  double width=.005;            // c, per ion in units of k
  double max_non_helium=1e-3;   // composition gate (mass fraction)
  // Below this coupling the forced-solid formula re-crosses spuriously (Gamma~40-90). The cut sits near the
  // largest solid-liquid separation (Gamma~80-100); near the cut the solid weight must be < exp(-20).
  double minimum_solid_gamma=90;
};
std::array<HelmholtzJet,10> ion_phase_difference_jets(
    double T,double rho,const Composition&,std::size_t channels=10,const IonPhaseOptions& = {});

// Same-composition (no-fractionation) phase term for a He-dominated H/He3/He4/GS98 mixture.
// Liquid and solid solutions are each linear mixtures of the per-species branches above at the
// common electron density, with identical ideal mixing. The soft minimum acts on the whole
// mixture with a width c per ion (c times the ion number per mass), so
//   F/T = F_L/T + W g(D/W),  D = sum_s X_s R/A_s (f_s,s - f_l,s),  W = c R sum_s X_s/A_s,
// g(y) = -ln(1 + exp(-y)); returned as F/T - F_L/T per mass with all composition channels.
// No phase term below the helium coupling cut; the solid weight must be negligible near it.
// Hydrogen above trace_hydrogen is outside the assessed lattice treatment and throws if the
// solid weight is not negligible there. Phase separation of metals is not represented.
struct MixturePhaseOptions {
  IonScreening screening=IonScreening::fitted_screening;
  double width=.005;
  double minimum_solid_gamma=90;   // helium coupling Gamma_He below which the mixture is liquid
  double trace_hydrogen=1e-3;
  double liquid_continuation_gamma=0;   // as ColdHeliumOptions: must match the liquid base
};
std::array<HelmholtzJet,10> ion_mixture_phase_difference_jets(
    double T,double rho,const Composition&,std::size_t channels=10,const MixturePhaseOptions& = {});
// Diagnostic: solid-solution weight 1/(1+exp(D/W)) and D/W at a state (value channel only).
std::array<double,2> ion_mixture_phase_weight(double T,double rho,const Composition&,const MixturePhaseOptions& = {});
} // namespace ember
