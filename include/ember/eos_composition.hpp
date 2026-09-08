#pragma once
#include "ember/eos_helmholtz.hpp"
#include <memory>

namespace ember {

// Interpolate MATERIAL potentials in source hydrogen fraction. Helium isotope
// number densities are preserved by a density/composition transformation;
// ideal isotope mixing/translation entropy is added analytically. Metals
// remain the explicitly selected He4 proxy. Supports either abundance basis.
class CompositionHelmholtzEos final : public Eos {
public:
  explicit CompositionHelmholtzEos(const std::filesystem::path& manifest,
      HelmholtzTableEos::Mixture mixture = HelmholtzTableEos::Mixture::exact);
  EosState eval(double T,double rho,const Composition& c) const override {
    return eval_with_derivatives(T,rho,c).state;
  }
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  const char* name() const override { return "FreeEOS composition potential; explicit He4 metal proxy and isotope approximation"; }
private:
  struct Coordinates {
    double scale{},x{},fraction{},isotope_phi{};
    std::array<double,2> dscale{},dx{};
    std::size_t interval{};
  };
  Coordinates coordinates(const Composition&) const;
  std::vector<double> x_;
  std::vector<std::unique_ptr<HelmholtzTableEos>> tables_;
  Composition metals_;
};
} // namespace ember
