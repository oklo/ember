#pragma once
#include "ember/opacity_hydrogen.hpp"
#include "ember/opacity_hydrogen_dense.hpp"
#include "ember/opacity_hydrogen_continuation.hpp"
#include "ember/opacity_hydrogen_share.hpp"
#include "ember/opacity_metal_extension.hpp"
#include "ember/opacity_cold_dense.hpp"

namespace ember {

// Own the radiative tables and their composition/temperature joins. Conduction
// is deliberately supplied separately by the selected heat transport law.
class RadiativeOpacity final : public Opacity {
public:
  struct Tables {
    std::filesystem::path low,warm,bridge,hot;
    std::filesystem::path cold_dense{};
    double cold_dense_scale{1};
  };
  struct Extension {
    std::filesystem::path hydrogen_response;
    double minimum_Z{},maximum_Z{.16},maximum_X{1};
    double dense_hydrogen_maximum_logR{1.8};
    double dense_hydrogen_maximum_logT{6.1};
  };
  explicit RadiativeOpacity(const Tables& tables) { build(tables,nullptr); }
  RadiativeOpacity(const Tables& tables,const Extension& extension) { build(tables,&extension); }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    return mapped_->eval(T,rho,c);
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    return mapped_->density_range(T,c);
  }
  const char* name() const override {return "joined radiative tables with elemental isotope mapping";}
private:
  std::vector<std::unique_ptr<Opacity>> sources_;
  const Opacity* mapped_{};
  template<class T,class... Args> const Opacity& add(Args&&... args) {
    auto source=std::make_unique<T>(std::forward<Args>(args)...);
    const auto& result=*source;sources_.push_back(std::move(source));return result;
  }
  void build(const Tables& tables,const Extension* extension) {
    const auto& warm=add<MixtureOpacity>(tables.warm);
    const auto& bridge=add<MixtureOpacity>(tables.bridge);
    const auto& hot=add<MixtureOpacity>(tables.hot);
    if(!extension) {
      const auto& low=add<MixtureOpacity>(tables.low);
      const auto& mid=add<BlendedOpacity>(warm,bridge,5.05,5.10);
      const auto& upper=add<BlendedOpacity>(mid,hot,5.6,5.7);
      const auto& raw=add<BlendedOpacity>(low,upper,4.4,4.47);
      map(raw,tables);return;
    }
    const auto& e=*extension;
    if(e.hydrogen_response.empty() || !std::isfinite(e.minimum_Z+e.maximum_Z+e.maximum_X)
        || e.minimum_Z<0 || e.minimum_Z>=.01 || e.maximum_Z<=.03 || e.maximum_Z>.16
        || e.maximum_X<.95 || e.maximum_X>1)
      throw std::invalid_argument("RadiativeOpacity: invalid bounded composition extension");
    const auto& low=add<HydrogenShareMixtureOpacity>(tables.low,true);
    const auto& dependence=add<HydrogenShareMixtureOpacity>(e.hydrogen_response);
    const auto& warm_x=add<HydrogenOpacityContinuation>(warm,.75,.05,e.maximum_X);
    const auto& bridge_x=add<HydrogenOpacityContinuation>(bridge,.75,.05,e.maximum_X);
    const auto& hot_ratio=add<HydrogenOpacityExtension>(hot,dependence,.75,e.maximum_X);
    const auto& hot_slope=add<HydrogenOpacityContinuation>(hot,.75,.05,e.maximum_X);
    const auto& hot_x=add<DenseHydrogenOpacity>(hot_ratio,hot_slope,e.dense_hydrogen_maximum_logR,e.dense_hydrogen_maximum_logT);
    const auto& mid=add<BlendedOpacity>(warm_x,bridge_x,5.05,5.10);
    const auto& upper=add<BlendedOpacity>(mid,hot_x,5.6,5.7);
    const auto& lower=add<LowerMetalHydrogenShareOpacity>(upper,.01,e.minimum_Z,
        LowerMetalOpacityContinuation::Method::linear_kappa);
    const auto& raw=add<BlendedOpacity>(low,lower,4.4,4.47);
    // The linear-kappa choice uses the .02 and .03 source planes. It does not
    // consult a second high-Z family, so do not load that unused dataset.
    const auto& extended=add<MetalOpacityExtension>(raw,raw,.03,e.maximum_Z,
        MetalOpacityExtension::Method::linear_kappa);
    map(extended,tables);
  }
  void map(const Opacity& raw,const Tables& tables) {
    const auto& source=tables.cold_dense.empty()?raw:
        add<ColdDenseOpacity>(raw,tables.cold_dense,tables.cold_dense_scale);
    mapped_=&add<ElementalOpacity>(source);
  }
};
} // namespace ember
