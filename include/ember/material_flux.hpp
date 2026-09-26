#pragma once
#include "ember/material_transport.hpp"

namespace ember {
struct MaterialHeatSplit {
  std::array<double,2> transported_enthalpy{}; // erg/g, includes EOS and kinetic terms
  double conductivity{};                     // erg/(cm s K), at zero species flux
};
struct MaterialFaceFlux {
  std::array<double,2> species_rate{}; // g/s, positive outward
  double material_luminosity{}, carried_luminosity{}, conductive_luminosity{};
  double entropy_production{}; // erg/(s K), integrated over this face interval
};

// Transform a reduced-heat mobility to U=(XH,X3,E/e0), W=(Phi_H,Phi_3,-e0/T)/s0.
// Input reduced heat is Q'=Q-h_EOS.j, scaled by local_energy_scale.
// The returned local matrix includes s0 but no spatial/area factors.
MaterialMatrix full_material_mobility(const MaterialMatrix& reduced,
    const std::array<double,2>& eos_exchange_enthalpy,double local_energy_scale,
    double energy_scale,double entropy_scale);

MaterialHeatSplit split_material_heat(const MaterialMatrix& full,double temperature,
    double energy_scale,double entropy_scale);

// Exact discrete face flux on the baryonic mass coordinate. The same area and
// density enter both species and heat transport. Phi=F/T composition gradients
// are supplied at the endpoints. No extra electron or gravitational force is
// added: these neutral-matter exchange directions already impose zero baryonic
// mass flux and zero electric current. The caller supplies matched full mobility.
// Density, area and local mobility are face values. For this discrete face,
// thermal conductivity is evaluated at sqrt(Tlo*Thi), so conductive heat equals
// -kappa0*area^2*rho*(Thi-Tlo)/dm exactly. This is a flux, not a heating source.
MaterialFaceFlux material_face_flux(const MaterialMatrix& full,
    const std::array<double,2>& phi_lo,const std::array<double,2>& phi_hi,
    double temperature_lo,double temperature_hi,double density,double area,
    double mass_interval,double energy_scale,double entropy_scale);
} // namespace ember
