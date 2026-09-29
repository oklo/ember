#pragma once
#include "ember/atmosphere_grid.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
// Use He4 source atmospheres at the same H/He number ratio and g/kappa.
// Neglect isotope-dependent cross sections, line broadening and molecular
// partition shifts. Trace metals are replaced only inside the atmosphere
// source query; the caller's composition supplies the matching density.
class HeliumIsotopeAtmosphere final : public Atmosphere {
public:
  HeliumIsotopeAtmosphere(const PressureDensity& eos,const CompositionAtmosphereGrid& source,
                         double maximum_helium3,double maximum_metals)
      :eos_(eos),source_(source),maximum_helium3_(maximum_helium3),maximum_metals_(maximum_metals) {
    if(!std::isfinite(maximum_helium3) || maximum_helium3<=0 || maximum_helium3>=1
        || !std::isfinite(maximum_metals) || maximum_metals<0 || maximum_metals>1e-10
        || source.reference_metallicity()>maximum_metals
        || source.support().helium3!=std::array<double,2>{0.,0.})
      throw std::invalid_argument("helium isotope atmosphere: invalid source or approximation bounds");
  }
  AtmosphereState eval(double teff,double gravity,const Composition& c) const override {
    if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98
        || !std::isfinite(c.sum()) || std::abs(c.sum()-1)>1e-10 || c[Species::H2]!=0
        || c[Species::He3]>maximum_helium3_ || c.Z()>maximum_metals_)
      throw std::domain_error("helium isotope atmosphere: outside declared composition approximation");
    for(double x:c.X)if(!std::isfinite(x) || x<0 || x>1)
      throw std::domain_error("helium isotope atmosphere: invalid composition");
    const double scale=1+c[Species::He3]/3;
    Composition mapped{};mapped.basis=AbundanceBasis::baryon_mass;mapped.metal_inventory=MetalInventory::gs98;
    mapped[Species::H1]=c[Species::H1]/scale;
    const auto& metals=source_.reference_metals();
    for(std::size_t k=0;k<NMETALS;++k)mapped.X[METAL_BEGIN+k]=metals[k];
    mapped[Species::He4]=1-mapped[Species::H1]-source_.reference_metallicity();
    auto s=source_.eval(teff,gravity/scale,mapped);
    // At fixed composition the logarithmic T/g derivatives are unchanged.
    s.rho=eos_.rho_from_PT(s.T,s.P,c,s.rho/scale);
    return s;
  }
  const char* name()const override{return "mixed-helium atmosphere with isotope mapping";}
private:
  const PressureDensity& eos_;
  const CompositionAtmosphereGrid& source_;
  double maximum_helium3_,maximum_metals_;
};
} // namespace ember
