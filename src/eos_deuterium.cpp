#include "ember/eos_deuterium.hpp"
#include "ember/constants.hpp"
#include "ember/nuclear_cn.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
DeuteriumApproxEos::Mapping DeuteriumApproxEos::map(const Composition& c) {
  Mapping m{c};
  const double D=c[Species::H2];
  if(D==0)return m;
  if(c.basis!=AbundanceBasis::baryon_mass || !std::isfinite(c.sum()) || std::abs(c.sum()-1)>1e-10)
    throw std::domain_error("deuterium EOS: normalized baryon fractions required");
  for(double x:c.X)if(!std::isfinite(x) || x<0)
    throw std::domain_error("deuterium EOS: nonnegative finite abundances required");
  if(D>1e-4)throw std::domain_error("deuterium EOS: trace approximation requires XD <= 1e-4");
  if(c.cn_molality) {
    if(c.cn_mass_convention!=CNMassConvention::explicit_metal_mass)
      throw std::domain_error("deuterium EOS: active CN requires explicit physical metal mass");
    (void)cn_physical_ledger(c,*c.cn_molality);
  }
  m.f=1-D/2;
  m.c[Species::H1]+=D/2;m.c[Species::H2]=0;
  for(double& x:m.c.X)x/=m.f;
  // The source gram has mass f times the physical gram. Carry its actual
  // catalyst numbers with the same normalization as every other species.
  if(m.c.cn_molality)for(double& y:*m.c.cn_molality)y/=m.f;
  const auto xlog=[](double x){return x>0?x*std::log(x):0.;};
  const double p=c[Species::H1],d=D/2;
  m.entropy=constants::R_gas*(xlog(p+d)-xlog(p)-xlog(d)+1.5*d*std::log(2.));
  return m;
}
EosState DeuteriumApproxEos::convert(EosState e,const Mapping& m) {
  e.E*=m.f;e.S=e.S*m.f+m.entropy;e.cv*=m.f;e.cp*=m.f;
  // Electron/nucleus ratio is unchanged when both number densities match.
  e.mu/=m.f;return e;
}
EosState DeuteriumApproxEos::eval(double T,double rho,const Composition& c) const {
  const auto m=map(c);return convert(source_.eval(T,m.f*rho,m.c),m);
}
EosResponse DeuteriumApproxEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  const auto m=map(c);auto r=source_.eval_with_derivatives(T,m.f*rho,m.c);
  r.state=convert(r.state,m);r.dE_dlnRho*=m.f;r.dcp_dlnT*=m.f;r.dcp_dlnRho*=m.f;return r;
}
EosCompositionResponse DeuteriumApproxEos::composition_response(double T,double rho,const Composition& c) const {
  const auto m=map(c);auto r=source_.composition_response(T,m.f*rho,m.c);
  // H1 or He3 replaces He4 at fixed D: dX_source=dX/f. Specific energy
  // contributes the inverse factor f, so its composition derivative cancels.
  for(double& d:r.dP)d/=m.f;return r;
}
std::optional<Eos::DensityRange> DeuteriumApproxEos::density_range(double T,const Composition& c) const {
  const auto m=map(c);auto r=source_.density_range(T,m.c);
  if(r){r->min/=m.f;r->max/=m.f;}return r;
}
double DeuteriumApproxEos::rho_from_PT(double T,double P,const Composition& c,double guess) const {
  const auto m=map(c);return source_.rho_from_PT(T,P,m.c,guess*m.f)/m.f;
}
} // namespace ember
