#pragma once
#include "ember/atmosphere_table.hpp"
#include "ember/opacity.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
// Deliberate, bounded initial-evolution approximations. These wrappers do
// not claim isotope-resolved opacity or composition-dependent COND data.
// Their names and limitations must be recorded with an evolutionary run.
class NominalAbundanceOpacity final : public Opacity {
public:
  explicit NominalAbundanceOpacity(const Opacity& opacity) : opacity_(opacity) {}
  OpacityState eval(double T,double rho,const Composition& c) const override {
    return opacity_.eval(T,rho,nominal(c));
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    return opacity_.density_range(T,nominal(c));
  }
  const char* name() const override { return "nominal X/Z opacity proxy; He3 treated as He4; Y3<=.005"; }
private:
  static Composition nominal(const Composition& c) {
    if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10 || c.X[1]>.005)
      throw std::domain_error("NominalAbundanceOpacity: requires normalized baryon abundances and Y3<=.005");
    for(double v:c.X) if(!std::isfinite(v) || v<0) throw std::domain_error("NominalAbundanceOpacity: invalid abundance");
    auto proxy=c;proxy.X[2]+=proxy.X[1];proxy.X[1]=0;proxy.basis=AbundanceBasis::atomic_mass;
    return proxy;
  }
  const Opacity& opacity_;
};

class FrozenCompositionAtmosphere final : public Atmosphere {
public:
  FrozenCompositionAtmosphere(const Eos& eos,const TabulatedAtmosphere& table)
      : eos_(eos),table_(table) {}
  AtmosphereState eval(double Teff,double gravity,const Composition& c) const override {
    if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10 || c.X[1]>.005
        || std::abs(c.X[0]-table_.composition().X[0])>.005)
      throw std::domain_error("FrozenCompositionAtmosphere: initial-evolution proxy limited to |X-Xref|<=.005 and Y3<=.005");
    for(std::size_t i=0;i<NSPEC;++i)
      if(!std::isfinite(c.X[i]) || c.X[i]<0 || (i>=3 && std::abs(c.X[i]-table_.composition().X[i])>1e-10))
        throw std::domain_error("FrozenCompositionAtmosphere: changed metals or invalid abundances");
    auto state=table_.eval(Teff,gravity,table_.composition());
    state.rho=eos_.rho_from_PT(state.T,state.P,c);
    return state;
  }
  const char* name() const override { return "COND tau100, frozen GN93 T/P boundary for bounded initial evolution"; }
private:
  const Eos& eos_;
  const TabulatedAtmosphere& table_;
};
} // namespace ember
