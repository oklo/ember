#include "ember/material_flux.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
void positive(double x) {
  if(!(x>0) || !std::isfinite(x))throw std::domain_error("material flux: positive finite scale required");
}
void validate(const MaterialMatrix& a) {
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
    if(!std::isfinite(a[i][j]) || a[i][j]!=a[j][i])
      throw std::domain_error("material flux: finite reciprocal mobility required");
  material_positive_inverse(a,3);
}
}

MaterialMatrix full_material_mobility(const MaterialMatrix& reduced,
    const std::array<double,2>& enthalpy,double local_energy_scale,
    double energy_scale,double entropy_scale) {
  positive(local_energy_scale);positive(energy_scale);positive(entropy_scale);validate(reduced);
  for(double h:enthalpy)if(!std::isfinite(h))throw std::domain_error("material flux: non-finite exchange enthalpy");
  MaterialMatrix g{},full{};g[0][0]=g[1][1]=1;
  g[2]={enthalpy[0]/energy_scale,enthalpy[1]/energy_scale,local_energy_scale/energy_scale};
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<=i;++j) {
    double sum=0;
    for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b)sum+=g[i][a]*reduced[a][b]*g[j][b];
    full[i][j]=full[j][i]=entropy_scale*sum;
  }
  validate(full);return full;
}

MaterialHeatSplit split_material_heat(const MaterialMatrix& full,double T,double e0,double s0) {
  positive(T);positive(e0);positive(s0);validate(full);
  const auto inverse=material_positive_inverse(full,2);
  std::array<double,2> ratio{};MaterialHeatSplit result;
  for(std::size_t i=0;i<2;++i) {
    for(std::size_t j=0;j<2;++j)ratio[i]+=inverse[i][j]*full[j][2];
    result.transported_enthalpy[i]=e0*ratio[i];
  }
  // Cholesky residual after scaling the diagonal: no large dimensional
  // energy terms are subtracted. A nonpositive residual is rejected.
  const double sA=std::sqrt(full[0][0]),sB=std::sqrt(full[1][1]),sE=std::sqrt(full[2][2]);
  const double l10=full[1][0]/sA/sB,l20=full[2][0]/sA/sE;
  const double l11=std::sqrt(1-l10*l10);
  const double l21=(full[2][1]/sB/sE-l10*l20)/l11;
  const double residual=1-l20*l20-l21*l21;positive(residual);
  result.conductivity=(e0/T)*(e0/T)/s0*full[2][2]*residual;
  positive(result.conductivity);
  for(double h:result.transported_enthalpy)if(!std::isfinite(h))
    throw std::domain_error("material flux: non-finite transported enthalpy");
  return result;
}

MaterialFaceFlux material_face_flux(const MaterialMatrix& full,
    const std::array<double,2>& phi_lo,const std::array<double,2>& phi_hi,
    double Ta,double Tb,double rho,double area,double dm,double e0,double s0) {
  positive(Ta);positive(Tb);positive(rho);positive(area);positive(dm);positive(e0);positive(s0);
  const auto split=split_material_heat(full,std::sqrt(Ta)*std::sqrt(Tb),e0,s0);
  MaterialVector difference{},rate{};
  for(std::size_t i=0;i<2;++i) {
    if(!std::isfinite(phi_lo[i]) || !std::isfinite(phi_hi[i]))
      throw std::domain_error("material flux: non-finite chemical potential");
    difference[i]=(phi_hi[i]-phi_lo[i])/s0;
  }
  difference[2]=(e0/Ta)/s0*((Tb-Ta)/Tb);
  const double geometry=area*(area/dm)*rho;positive(geometry);
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)rate[i]-=geometry*full[i][j]*difference[j];
  MaterialFaceFlux result;result.species_rate={rate[0],rate[1]};result.material_luminosity=e0*rate[2];
  for(std::size_t i=0;i<2;++i)result.carried_luminosity+=split.transported_enthalpy[i]*rate[i];
  result.conductive_luminosity=-split.conductivity*geometry*(Tb-Ta);
  for(std::size_t i=0;i<3;++i)result.entropy_production-=s0*rate[i]*difference[i];
  for(double v:{rate[0],rate[1],result.material_luminosity,result.carried_luminosity,
                result.conductive_luminosity,result.entropy_production})
    if(!std::isfinite(v))throw std::domain_error("material flux: non-finite face flux");
  return result;
}
} // namespace ember
