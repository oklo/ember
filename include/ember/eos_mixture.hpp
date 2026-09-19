#pragma once
#include "ember/eos_helmholtz.hpp"
#include <memory>

namespace ember {
// Fixed GS98 Z=.02, two-dimensional baryonic H/He3 material potentials.
// Source element counts already include He3 and metals; radiation and ideal
// helium isotope entropy are added once, after potential interpolation.
class MetalHelmholtzEos final : public Eos {
public:
  explicit MetalHelmholtzEos(const std::filesystem::path&,
      HelmholtzTableEos::Mixture = HelmholtzTableEos::Mixture::exact);
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  // Require identical original potentials and masks; only hotter rows may be added.
  std::size_t check_temperature_extension(const MetalHelmholtzEos& original) const;
  // Require identical original potentials and masks; only denser columns may be added.
  std::size_t check_density_extension(const MetalHelmholtzEos& original) const;
  const char* name() const override {return "FreeEOS GS98 baryonic H/He3 potential; trace K omission and isotope approximation";}
private:
  std::size_t check_extension(const MetalHelmholtzEos& original, bool density) const;
  struct Coordinates {std::size_t x,y;double u,v;};
  Coordinates coordinates(const Composition&) const;
  const HelmholtzTableEos& table(std::size_t x,std::size_t y) const {return *tables_[x*y_.size()+y];}
  std::vector<double> x_,y_;
  std::vector<std::unique_ptr<HelmholtzTableEos>> tables_;
  Composition metals_;
};
} // namespace ember
