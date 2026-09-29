#pragma once
#include "ember/atmosphere.hpp"
#include "ember/nuclear_cn.hpp"
#include <cmath>

namespace ember {
// Explicit trace-isotope approximation for an existing atmosphere family.
// Its T/P boundary uses the total hydrogen isotope mass as ordinary H;
// density is then inverted with the actual deuterium-aware interior EOS.
// This keeps the source Z and helium inventory fixed, but omits the O(XD)
// atmosphere change from isotope number counts, HD and line/band shifts.
// It is separate from the exact fuel inventory used by the nuclear solver.
class TraceDeuteriumAtmosphere final : public Atmosphere {
 public:
  TraceDeuteriumAtmosphere(const PressureDensity& eos,const Atmosphere& source)
      :eos_(eos),source_(source){}
  AtmosphereState eval(double Teff,double gravity,const Composition& c) const override {
    const double D=c[Species::H2];
    if(D==0)return source_.eval(Teff,gravity,c);
    if(c.basis!=AbundanceBasis::baryon_mass || D<0 || D>1e-4 || !std::isfinite(c.sum())
        || std::abs(c.sum()-1)>1e-10)
      throw std::domain_error("trace deuterium atmosphere: invalid or unsupported composition");
    if(c.cn_molality) {
      if(c.cn_mass_convention!=CNMassConvention::explicit_metal_mass)
        throw std::domain_error("trace deuterium atmosphere: active CN requires explicit physical metal mass");
      (void)cn_physical_ledger(c,*c.cn_molality);
    }
    for(double x:c.X)if(!std::isfinite(x) || x<0)
      throw std::domain_error("trace deuterium atmosphere: invalid abundance");
    auto proxy=c;proxy.X[0]+=D;proxy[Species::H2]=0;
    auto result=source_.eval(Teff,gravity,proxy);
    result.rho=eos_.rho_from_PT(result.T,result.P,c,result.rho);return result;
  }
  const char* name() const override {return "trace-D atmosphere: total H mass in T/P lookup, actual isotope EOS density";}
 private:
  const PressureDensity& eos_;const Atmosphere& source_;
};
} // namespace ember
