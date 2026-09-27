#pragma once
#include "ember/atmosphere_grid.hpp"
#include "ember/nuclear_cn.hpp"
#include <cmath>
#include <limits>

namespace ember {
// Explicit boundary approximation for minute changes of a fixed metal pattern.
// T and P use the source metal abundance, with He4 absorbing the difference in
// the lookup alone. H and He3 keep their actual values. The returned density is
// inverted with the actual composition. No stellar or nuclear inventory changes.
// The caller must justify the selected bound for its physical regime; the hard
// ceiling below is a scope limit, not a certified error bound on an atmosphere.
// This cannot represent settling or a changed elemental pattern.
class FixedMetalAtmosphere final : public Atmosphere {
public:
  FixedMetalAtmosphere(const Eos& eos,const CompositionAtmosphereGrid& source,
                      double maximum_delta_Z)
      :eos_(eos),source_(source),maximum_delta_Z_(maximum_delta_Z),
       reference_Z_(source.reference_metallicity()) {
    if(!std::isfinite(maximum_delta_Z_) || maximum_delta_Z_<=0 || reference_Z_<=0 ||
        maximum_delta_Z_>1e-6*reference_Z_)
      throw std::invalid_argument("fixed-metal atmosphere: invalid or excessive approximation bound");
  }
  AtmosphereState eval(double Teff,double gravity,const Composition& c) const override {
    if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98 ||
        !std::isfinite(c.sum()) || std::abs(c.sum()-1)>1e-10 || c[Species::H2]!=0)
      throw std::domain_error("fixed-metal atmosphere: invalid composition or unmapped deuterium");
    for(double x:c.X)if(!std::isfinite(x) || x<0 || x>1)
      throw std::domain_error("fixed-metal atmosphere: invalid abundance");
    if(c.cn_molality) {
      if(c.cn_mass_convention!=CNMassConvention::explicit_metal_mass)
        throw std::domain_error("fixed-metal atmosphere: active CN requires explicit metal mass");
      (void)cn_physical_ledger(c,*c.cn_molality);
    }
    const double z=c.Z(),delta=z-reference_Z_;
    const double roundoff=16*std::numeric_limits<double>::epsilon()*reference_Z_;
    if(std::abs(delta)>maximum_delta_Z_+roundoff)
      throw std::domain_error("fixed-metal atmosphere: metal change exceeds selected approximation bound");
    if(delta==0)return source_.eval(Teff,gravity,c);
    auto lookup=c;
    // Retain the actual relative pattern. The strict source reader checks
    // it against its recorded mixture with its original component tolerance.
    for(std::size_t j=METAL_BEGIN;j<METAL_END;++j)lookup.X[j]*=reference_Z_/z;
    lookup.X[2]+=delta;
    if(lookup.X[2]<0)throw std::domain_error("fixed-metal atmosphere: negative reference helium");
    // Synthetic material lookup only; its CN ledger would no longer match Z.
    lookup.cn_molality.reset();lookup.cn_mass_convention=CNMassConvention::fixed_metal_proxy;
    auto result=source_.eval(Teff,gravity,lookup);
    result.rho=eos_.rho_from_PT(result.T,result.P,c,result.rho);
    return result;
  }
  const char* name() const override {return "bounded fixed-metal T/P boundary; actual-composition density";}
private:
  const Eos& eos_;
  const CompositionAtmosphereGrid& source_;
  double maximum_delta_Z_,reference_Z_;
};
} // namespace ember
