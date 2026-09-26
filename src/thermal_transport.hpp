#pragma once
#include "ember/constants.hpp"
#include "ember/structure.hpp"
#include "differential.hpp"
#include <stdexcept>

namespace ember::detail {
inline void check_thermal_transport(const Physics& p) {
  if(p.microscopic && (!p.opacity || p.opacity->includes_conduction()))
    throw std::invalid_argument("microscopic transport requires radiative-only opacity; conduction comes from the same microscopic evaluator");
}
template<std::size_t N> Differential<N> logarithmic_temperature(
    const Differential<N>& lnlo,const Differential<N>& lnhi) {
  const auto d=lnhi-lnlo;
  if(std::abs(d.value)<1e-3) {
    const auto z=d*d;
    return exp(.5*(lnlo+lnhi))*(1+z*(1./24+z*(1./1920+z/322560.)));
  }
  return (exp(lnhi)-exp(lnlo))/d;
}
template<std::size_t N> struct ThermalTransport {
  Differential<N> opacity,gradient;
};
// Volume-face models use the same thermal difference with or without a
// species-transport provider. Only older nodal models keep the old midpoint
// expression. Adding a vanishing material flux must not change the stencil.
// Tlog converts delta ln T into delta T exactly. When hydrostatic balance
// holds, the radiative solution reproduces the mass-coordinate face flux.
// MLT still uses the physical local thermal diffusivity at the midpoint.
template<std::size_t N> ThermalTransport<N> thermal_transport(
    const Differential<N>& T,const Differential<N>& rho,const Differential<N>& P,double mass,
    const Differential<N>& opacity,const Differential<N>& luminosity,bool conservative_face,
    const Differential<N>& carried,const Differential<N>& conductivity,
    const Differential<N>& lnlo,const Differential<N>& lnhi) {
  using constants::a_rad;using constants::c;using constants::G;
  if(!conservative_face)return {opacity,3*opacity*luminosity*P/(16*M_PI*a_rad*c*G*mass*T*T*T*T)};
  const auto thermal=4*a_rad*c*T*T*T/(3*rho*opacity)+conductivity;
  const auto k=4*a_rad*c*T*T*T/(3*rho*thermal);
  const auto grad=(luminosity-carried)*P/(4*M_PI*G*mass*rho*thermal*logarithmic_temperature(lnlo,lnhi));
  return {k,grad};
}
} // namespace ember::detail
