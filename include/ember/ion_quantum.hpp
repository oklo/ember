#pragma once
#include "ember/eos_helmholtz.hpp"

namespace ember {
// Quantum addition to an otherwise classical liquid-ion free energy.
// Baiko & Chugunov (2022), equations 33–34, linearly mixed at a common
// electron density. The leading Wigner–Kirkwood mixture limit is exact.
// Beyond that limit linear mixing is an approximation, not a mixture fit.
// Returned channels: F/T, H/He3/Z gradients, symmetric Hessian, with all
// thermal derivatives through total degree three. Composition replaces He4.
// No electron quantum, radiation, screening, or classical ion term is added.
std::array<HelmholtzJet,10> ion_quantum_liquid_jets(
    double T,double rho,const Composition&,std::size_t channels=10);
namespace detail {
// The same correction for one pure species (unit mass fraction), F/T per mass, at ln T and ln n_e,
// without domain guards (callers apply their own). Used by the cold He liquid mixture.
HelmholtzJet bc22_liquid_quantum_per_mass(double logT,double logne,double A,double Z);
}
} // namespace ember
