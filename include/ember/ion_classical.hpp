#pragma once
#include "ember/eos_helmholtz.hpp"
namespace ember {
// Classical ion-ion liquid free energy, Potekhin & Chabrier (2000), using
// the DeWitt–Slattery coefficients also used by the nuclear screening fit.
// Linear mixing at common electron density includes H/He isotope and GS98
// composition derivatives. Above individual Gamma=200 the free energy is
// continued linearly in T; this sets that species' heat-capacity term to zero.
// No ideal ions, electrons, polarization, nonlinear mixture correction,
// quantum term or phase selection is included. Never add this to an EOS
// already containing ion-ion interactions without removing its existing term.
// Channels follow ion_quantum_liquid_jets: F/T, three gradients and six Hessians.
std::array<HelmholtzJet,10> ion_classical_liquid_jets(
    double T,double rho,const Composition&,std::size_t channels=10);
}
