#pragma once
#include "ember/eos_helmholtz.hpp"

namespace ember {
struct MetalCompositionPotentialResponse {
  double phi{};
  // Physical mass-fraction directions H1, He3, GS98 metals replacing He4.
  std::array<double,3> gradient{},dgradient_dlnT{},dgradient_dlnRho{};
  std::array<std::array<double,3>,3> hessian{};
};
struct MetalCompositionHeatResponse {
  double material_delta{};
  std::array<double,3> exchange_enthalpy{},radiation_enthalpy{};
  // Coordinates ln T, ln rho, XH, X3, Z.
  std::array<double,5> delta_partials{};
  std::array<std::array<double,5>,3> enthalpy_partials{},radiation_enthalpy_partials{};
};
struct IsobaricCompositionResponse {
  // H1, He3 and GS98 metals replace He4, at fixed temperature and pressure.
  std::array<double,3> dlnRho{};
  std::array<std::array<double,3>,3> potential_hessian{}; // Gibbs free energy / T
};

// Material potential shared by structure, composition forces and heat transport.
// Coordinates replace He4 with H1, He3 or a fixed GS98 metal mixture. Providers
// return F/T before analytic ionic mixing and radiation; both are restored here
// exactly once. This convention permits table and phase-equilibrium providers
// to use the same transport equations.
class MaterialEos : public Eos {
 public:
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  MetalCompositionPotentialResponse composition_potential(double,double,const Composition&,
      std::array<bool,3> active={true,true,true},bool hessian=true) const;
  IsobaricCompositionResponse isobaric_composition_response(double,double,const Composition&,
      std::array<bool,3> active={true,true,true}) const;
  MetalCompositionHeatResponse composition_heat(double,double,const Composition&,
      std::array<bool,3> active={true,true,true},bool derivatives=true,
      bool composition_derivatives=true) const;
  // Channels: value, H, He3, Z, HH, HHe3, HZ, He3He3, He3Z, ZZ.
  // Thermal/density derivatives through total order three. Valid channel
  // counts are 1, 4 and 10. An unrequested composition entry may be NaN.
  virtual std::array<HelmholtzJet,10> material_jets(double,double,
      const Composition&,std::size_t channels) const=0;
  // Reuse of nearby responses must still enforce the provider's support.
  virtual void validate_composition_domain(double,double,const Composition&) const=0;
};
} // namespace ember
