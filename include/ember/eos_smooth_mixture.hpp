#pragma once
#include "ember/eos_mixture.hpp"

namespace ember {
// Composition derivatives of MATERIAL Phi=F/T at fixed T, rho and metals.
// Directions replace He4 with H1 or He3. Radiation has zero composition
// response. Units: erg/g/K per unit mass fraction (squared for Hessian).
struct CompositionPotentialResponse {
  double phi{};
  std::array<double,2> gradient{}, dgradient_dlnT{}, dgradient_dlnRho{};
  std::array<std::array<double,2>,2> hessian{};
};
struct CompositionHeatResponse {
  double material_delta{};
  std::array<double,2> exchange_enthalpy{}; // erg/g, independent species replace He4
  // Coordinates ln T, ln rho, XH, X3, with the metal inventory fixed.
  std::array<double,4> delta_partials{};
  std::array<std::array<double,4>,2> enthalpy_partials{};
  // Additional composition enthalpy of radiation in a pressure-balanced LTE
  // fluid parcel, separate from the microscopic material convention.
  std::array<double,2> radiation_enthalpy{};
  std::array<std::array<double,4>,2> radiation_enthalpy_partials{};
};

// Explicit, independently selectable composition interpolation. Requires
// three or four He3 planes; uses a quadratic or cubic polynomial there,
// respectively, and a C2 not-a-knot H spline.
// H slopes are formed once from contiguous valid source nodes. Analytic
// isotope/ion mixing is removed before interpolation and restored exactly.
// This class does not change the default EOS used by the evolution executable.
class SmoothMetalHelmholtzEos final : public Eos {
public:
  explicit SmoothMetalHelmholtzEos(const std::filesystem::path&,
      HelmholtzTableEos::Mixture = HelmholtzTableEos::Mixture::exact);
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  CompositionPotentialResponse composition_potential(double,double,const Composition&) const;
  // Finite chemical part psi_k = Phi_k - (R_gas/A_k)*ln(X_k).
  // Returned gradient/Hessian are psi and dpsi/dX; phi remains full F/T.
  // Defined at XH=0 or X3=0, with positive reference He4. The omitted ideal
  // logarithm must be treated analytically in a finite species flux.
  CompositionPotentialResponse regular_composition_potential(double,double,const Composition&) const;
  // Select independent directions [H1 replacing He4, He3 replacing He4].
  // Only selected species and the He4 reference must be positive. Derivatives
  // in unselected directions are NaN, not zero; phi always remains defined.
  // Selecting neither returns only phi. No absent species is given a floor.
  CompositionPotentialResponse active_composition_potential(
      double,double,const Composition&,std::array<bool,2> active) const;
  // Material partial enthalpy differences and analytic derivatives. Radiation
  // is excluded from delta. Unselected species and their composition partials
  // are NaN, as in active_composition_potential; no abundance floor is used.
  // Separate radiation fields give the extra isobaric parcel coefficient.
  CompositionHeatResponse composition_heat(double,double,const Composition&,
      std::array<bool,2> active,bool derivatives=true) const;
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  std::optional<DensityRange> density_range_near(double,const Composition&,double) const override;
  const char* name() const override {
    return "FreeEOS GS98 smooth H/He3 potential; trace K omission and isotope approximation";
  }
private:
  using WeightedTable=HelmholtzTableEos::WeightedTable;
  struct WeightedTables {
    std::array<WeightedTable,16> tables{};
    std::size_t count{};
    std::span<const WeightedTable> span() const {return {tables.data(),count};}
  };
  WeightedTables weights(const Composition&,unsigned dx,unsigned dy) const;
  HelmholtzJet residual_jet(double,double,const Composition&,unsigned,unsigned) const;
  MetalHelmholtzEos residual_;
  std::vector<std::unique_ptr<HelmholtzTableEos>> slopes_;
};
} // namespace ember
