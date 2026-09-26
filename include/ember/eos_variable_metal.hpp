#pragma once
#include "ember/eos_helmholtz.hpp"
#include <memory>

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

// Free-energy interpolation in Z, u=XH/(1-Z), and v=X3/(1-Z-XH).
// The source grid remains physical as the helium reservoir shrinks. Source
// ionic mixing is removed before interpolation, and exact isotope mixing
// is restored afterwards. Radiation is added once by helmholtz_response.
// Composition forces and transported enthalpies differentiate this same
// potential. This class is independently selectable; it is not production
// selection of any newly generated source family.
// Manifest v1 uses a global degree <=3 polynomial in Z. Version 2 keeps
// the first four Z planes identical and appends C2 quintic intervals.
// Endpoint derivatives come from the preceding four source planes; adding
// a higher-Z plane never changes the already covered composition interval.
class VariableMetalHelmholtzEos final : public Eos {
 public:
  explicit VariableMetalHelmholtzEos(const std::filesystem::path&,
      HelmholtzTableEos::Mixture=HelmholtzTableEos::Mixture::exact);
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  MetalCompositionPotentialResponse composition_potential(double,double,const Composition&,
      std::array<bool,3> active={true,true,true}) const;
  MetalCompositionHeatResponse composition_heat(double,double,const Composition&,
      std::array<bool,3> active={true,true,true},bool derivatives=true) const;
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  const char* name() const override {return "FreeEOS variable GS98 metals, H and helium-isotope potential";}
 private:
  using WeightedTable=HelmholtzTableEos::WeightedCompositionTable;
  struct Weights {std::array<WeightedTable,80> tables{};std::size_t count{};};
  Weights weights(const Composition&,std::size_t channels) const;
  std::array<HelmholtzJet,10> jets(double,double,const Composition&,std::size_t) const;
  std::size_t index(std::size_t z,std::size_t u,std::size_t v) const {
    return (z*u_.size()+u)*v_.size()+v;
  }
  std::vector<double> z_,u_,v_;
  bool extend_metals_{};
  // Five source weights, each polynomial in the interval coordinate.
  using MetalExtension=std::array<std::array<double,6>,5>;
  std::vector<MetalExtension> metal_extensions_;
  std::vector<std::unique_ptr<HelmholtzTableEos>> tables_,slopes_;
  std::array<double,5> metal_pattern_{};
};
} // namespace ember
