#pragma once
#include "ember/eos_helmholtz.hpp"
#include "ember/eos_additive_volume.hpp"
#include "ember/cold_helium_base.hpp"
#include <optional>
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
  static constexpr const char* cold_model_identifier=
      "additive_volume.thermal_potential_join.He_proxy_metals.v2";
  // Optional local low-Z potential; C2 blend back to cubic by z_[2].
  enum class LowMetalInterpolation { cubic, quadratic };
  explicit VariableMetalHelmholtzEos(const std::filesystem::path&,
      HelmholtzTableEos::Mixture=HelmholtzTableEos::Mixture::exact,
      LowMetalInterpolation=LowMetalInterpolation::cubic,
      const std::filesystem::path& cold_potential={},bool quantum_ions=false,
      std::optional<ColdHeliumOptions> cold_helium={});
  // Optional cold dense-He liquid-mixture join (cold_helium_base.hpp), applied to the material potential
  // itself so that eval, composition forces, transported enthalpies and domain checks share it. It requires
  // quantum_ions (the base carries the same Baiko-Chugunov term, so it is present exactly once on each side),
  // is exclusive with cold_potential, and never replaces a failed source query below its window.
  static constexpr const char* cold_helium_identifier=
      "cold_helium.liquid_mixture.full_anchor.bc22_linear.X005_200.v2";
  // Opt-in same-composition phase variant (ColdHeliumOptions::mixture_phase).
  static constexpr const char* cold_helium_phase_identifier=
      "cold_helium.mixture_softmin.same_composition.common_liquid.width0.005.full_anchor.bc22_linear.X005_200.v3";
  static constexpr const char* dense_transition_identifier=
      "dense_hhe.liquid.bc22.T200_300.rho300_600.full_anchor.v2";
  // Convert a text family and its planes to one relocatable binary input.
  // Stored doubles, masks and logarithmic coordinates remain bit-identical.
  static void pack_binary(const std::filesystem::path& source,const std::filesystem::path& destination);
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  // Unrequested derivative entries remain NaN. Thermal-only callers do not
  // need second derivatives with respect to composition.
  MetalCompositionPotentialResponse composition_potential(double,double,const Composition&,
      std::array<bool,3> active={true,true,true},bool hessian=true) const;
  MetalCompositionHeatResponse composition_heat(double,double,const Composition&,
      std::array<bool,3> active={true,true,true},bool derivatives=true,
      bool composition_derivatives=true) const;
  // Source support for all composition derivative channels, without
  // evaluating the free-energy polynomial (used by optional response reuse).
  void validate_composition_domain(double T,double rho,const Composition&) const;
  // Material F/T jets (radiation excluded) of the same potential used by eval(); exposed so that an
  // external component can be joined to this table at the potential level.
  std::array<HelmholtzJet,10> material_jets(double T,double rho,const Composition& c,std::size_t channels) const {
    return jets(T,rho,c,channels);
  }
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  std::optional<DensityRange> density_range_near(double,const Composition&,double) const override;
  const char* name() const override {return cold_helium_?(cold_helium_->mixture_phase?
      "FreeEOS variable GS98 metals with cold dense-He mixture join and same-composition crystallization":
      "FreeEOS variable GS98 metals with cold dense-He liquid-mixture join"):cold_?
      "FreeEOS with cold additive-volume H/He potential; cold metal helium proxy":
      "FreeEOS variable GS98 metals, H and helium-isotope potential";}
 private:
  using WeightedTable=HelmholtzTableEos::WeightedCompositionTable;
  struct Weights {std::array<WeightedTable,80> tables{};std::size_t count{};};
  Weights weights(const Composition&,std::size_t channels) const;
  std::optional<DensityRange> density_range_impl(double,const Composition&,std::optional<double>) const;
  std::array<HelmholtzJet,10> jets(double,double,const Composition&,std::size_t) const;
  std::array<HelmholtzJet,10> source_jets(double,double,const Composition&,std::size_t) const;
  std::array<HelmholtzJet,10> raw_source_jets(double,double,const Composition&,std::size_t) const;
  void validate_source_domain(double,double,const Composition&) const;
  void validate_raw_source_domain(double,double,const Composition&) const;
  std::size_t index(std::size_t z,std::size_t u,std::size_t v) const {
    return (z*u_.size()+u)*v_.size()+v;
  }
  std::vector<double> z_,u_,v_;
  bool extend_metals_{};
  LowMetalInterpolation low_metal_interpolation_{};
  // Optional physical model, never a fallback on a failed source query.
  // In the cold component only, metals use a helium electronic/caloric proxy;
  // actual species inventories and analytic metal mixing entropy are retained.
  std::unique_ptr<AdditiveVolumePotential> cold_;
  bool quantum_ions_{};
  std::optional<ColdHeliumOptions> cold_helium_;
  // Five source weights, each polynomial in the interval coordinate.
  using MetalExtension=std::array<std::array<double,6>,5>;
  std::vector<MetalExtension> metal_extensions_;
  std::vector<std::unique_ptr<HelmholtzTableEos>> tables_,slopes_;
  std::array<double,5> metal_pattern_{};
  // Cache ownership must survive address reuse without retaining the tables.
  std::shared_ptr<const char> jet_cache_identity_{std::make_shared<const char>(0)};
};
} // namespace ember
