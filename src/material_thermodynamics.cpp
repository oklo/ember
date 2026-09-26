#include "ember/material_transport.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
MaterialPoint native_material_point(const SmoothMetalHelmholtzEos& eos,double T,double rho,
    const Composition& c,double energy_scale,double entropy_scale,bool derivatives) {
  if(!(T>0) || !std::isfinite(T) || !(rho>0) || !std::isfinite(rho)
      || !(energy_scale>0) || !std::isfinite(energy_scale)
      || !(entropy_scale>0) || !std::isfinite(entropy_scale))
    throw std::domain_error("native material thermodynamics: invalid state or normalization");
  const auto e=eos.eval(T,rho,c);
  MaterialPoint p;p.conserved={c.X[0],c.X[1],e.E/energy_scale};p.entropy=e.S/entropy_scale;
  if(!derivatives)return p;
  if(!(e.cv>0) || !std::isfinite(e.cv))throw std::domain_error("native material thermodynamics: nonpositive heat capacity");
  const auto ec=eos.composition_response(T,rho,c).dE;
  const auto phi=eos.composition_potential(T,rho,c);
  double error=0,scale=entropy_scale*T;
  for(std::size_t i=0;i<2;++i) {
    error=std::max(error,std::abs(ec[i]+T*phi.dgradient_dlnT[i]));scale=std::max(scale,std::abs(ec[i]));
  }
  if(!std::isfinite(error) || error>1e-9*scale)
    throw std::domain_error("native material thermodynamics: chemical energy identity failed");
  MaterialMatrix h{};
  for(std::size_t i=0;i<2;++i) {
    p.potential[i]=phi.gradient[i]/entropy_scale;
    for(std::size_t j=0;j<2;++j)h[i][j]=(phi.hessian[i][j]+ec[i]*ec[j]/(e.cv*T*T))/entropy_scale;
    h[i][2]=h[2][i]=-energy_scale*ec[i]/(e.cv*T*T*entropy_scale);
    p.primitive_from_conserved[i][i]=1;p.primitive_from_conserved[2][i]=-ec[i]/(e.cv*T);
  }
  h[2][2]=energy_scale*energy_scale/(e.cv*T*T*entropy_scale);
  p.capacity=material_positive_inverse(h,3);p.potential[2]=-energy_scale/(T*entropy_scale);
  p.primitive_from_conserved[2][2]=energy_scale/(e.cv*T);
  return p;
}
} // namespace ember
