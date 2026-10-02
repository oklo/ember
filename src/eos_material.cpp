#include "ember/eos_material.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
constexpr std::size_t second_channel(std::size_t i,std::size_t j) {
  if(i>j)std::swap(i,j);
  return 4+i*3-i*(i-1)/2+j-i;
}
double xlogx(double x){return x>0?x*std::log(x):0.;}
bool potassium(const Gs98Metal& m){return m.charge==19 && m.mass_number==39;}
double metal_mixing(double Z) {
  double sum=0;for(const auto& m:gs98_metals)if(!potassium(m))sum+=xlogx(Z*m.fraction/m.mass_number);
  return sum;
}
double full_mixing(const Composition& c) {
  const double n3=c.X[1]/3;
  return constants::R_gas*(xlogx(c.X[0])+xlogx(n3)+xlogx(c.X[2]/4)+metal_mixing(c.Z())
      -1.5*n3*std::log(nuclides[1].A/nuclides[2].A));
}
bool active_any(std::array<bool,3> active){return active[0] || active[1] || active[2];}
void check_active(const Composition& c,std::array<bool,3> active) {
  if((active[0] && c.X[0]<=0) || (active[1] && c.X[1]<=0) || (active[2] && c.Z()<=0)
      || (active_any(active) && c.X[2]<=0))
    throw std::domain_error("variable EOS: positive active species and reference He4 required");
}
} // namespace
EosResponse MaterialEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  auto f=material_jets(T,rho,c,1)[0];f[0][0]+=full_mixing(c);return helmholtz_response(T,rho,f);
}
EosCompositionResponse MaterialEos::composition_response(double T,double rho,const Composition& c) const {
  const auto j=material_jets(T,rho,c,4);EosCompositionResponse out;
  for(std::size_t k=0;k<2;++k){out.dP[k]=rho*T*j[1+k][0][1];out.dE[k]=-T*j[1+k][1][0];}
  return out;
}
MetalCompositionPotentialResponse MaterialEos::composition_potential(double T,double rho,
    const Composition& c,std::array<bool,3> active,bool hessian) const {
  check_active(c,active);const auto j=material_jets(T,rho,c,active_any(active)?(hessian?10:4):1);
  MetalCompositionPotentialResponse out;out.phi=j[0][0][0]+full_mixing(c);
  const double nan=std::numeric_limits<double>::quiet_NaN();
  out.gradient.fill(nan);out.dgradient_dlnT.fill(nan);out.dgradient_dlnRho.fill(nan);for(auto& h:out.hessian)h.fill(nan);
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    out.gradient[k]=j[1+k][0][0];out.dgradient_dlnT[k]=j[1+k][1][0];out.dgradient_dlnRho[k]=j[1+k][0][1];
    if(hessian)for(std::size_t l=0;l<3;++l)if(active[l])out.hessian[k][l]=j[second_channel(k,l)][0][0];
  }
  auto add=[&](double n,std::array<double,3> b) {
    if(n<=0)return;
    for(std::size_t k=0;k<3;++k)if(active[k]) {
      out.gradient[k]+=constants::R_gas*b[k]*(std::log(n)+1);
      if(hessian)for(std::size_t l=0;l<3;++l)if(active[l])out.hessian[k][l]+=constants::R_gas*b[k]*b[l]/n;
    }
  };
  add(c.X[0],{1,0,0});add(c.X[1]/3,{0,1./3,0});add(c.X[2]/4,{-.25,-.25,-.25});
  for(const auto& m:gs98_metals)if(!potassium(m))add(c.Z()*m.fraction/m.mass_number,{0,0,m.fraction/m.mass_number});
  if(active[1])out.gradient[1]-=.5*constants::R_gas*std::log(nuclides[1].A/nuclides[2].A);
  return out;
}
IsobaricCompositionResponse MaterialEos::isobaric_composition_response(
    double T,double rho,const Composition& c,std::array<bool,3> active) const {
  const auto p=composition_potential(T,rho,c,active,true);
  const auto e=eval(T,rho,c);
  // dP/dX = rho*T*d(phi)/dlnrho/dX; eliminate the density change required
  // to keep P fixed. Radiation contributes no composition or density term.
  const double scale=rho*T/(e.P*e.chiRho);
  if(!(scale>0) || !std::isfinite(scale))
    throw std::domain_error("isobaric composition: nonpositive compressibility");
  IsobaricCompositionResponse result;
  result.dlnRho.fill(std::numeric_limits<double>::quiet_NaN());
  result.potential_hessian=p.hessian;
  for(std::size_t j=0;j<3;++j)if(active[j]) {
    result.dlnRho[j]=-scale*p.dgradient_dlnRho[j];
    for(std::size_t k=0;k<3;++k)if(active[k])
      result.potential_hessian[k][j]-=scale*p.dgradient_dlnRho[k]*p.dgradient_dlnRho[j];
  }
  return result;
}
MetalCompositionHeatResponse MaterialEos::composition_heat(double T,double rho,
    const Composition& c,std::array<bool,3> active,bool derivatives,bool composition_derivatives) const {
  check_active(c,active);const auto j=material_jets(T,rho,c,active_any(active)?(derivatives && composition_derivatives?10:4):1);const auto& f=j[0];
  (void)helmholtz_response(T,rho,f);
  const double denominator=f[0][1]+f[0][2],delta=(f[0][1]+f[1][1])/denominator;
  const double radiation=4*constants::a_rad*T*T*T/(3*rho*denominator);
  if(!(denominator>0) || !std::isfinite(delta+radiation))throw std::domain_error("variable EOS: invalid material heat response");
  MetalCompositionHeatResponse out;out.material_delta=delta;
  const double nan=std::numeric_limits<double>::quiet_NaN();out.exchange_enthalpy.fill(nan);out.radiation_enthalpy.fill(nan);out.delta_partials.fill(nan);
  for(auto& h:out.enthalpy_partials)h.fill(nan);for(auto& h:out.radiation_enthalpy_partials)h.fill(nan);
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    const auto& g=j[1+k];out.exchange_enthalpy[k]=T*(delta*g[0][1]-g[1][0]);out.radiation_enthalpy[k]=T*radiation*g[0][1];
  }
  if(!derivatives)return out;
  auto& d=out.delta_partials;
  std::array<double,5> dr{radiation*(3-(f[1][1]+f[1][2])/denominator),radiation*(-1-(f[0][2]+f[0][3])/denominator),nan,nan,nan};
  d[0]=(f[1][1]+f[2][1]-delta*(f[1][1]+f[1][2]))/denominator;
  d[1]=(f[0][2]+f[1][2]-delta*(f[0][2]+f[0][3]))/denominator;
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    const auto& g=j[1+k];
    if(composition_derivatives) {
      d[2+k]=(g[0][1]+g[1][1]-delta*(g[0][1]+g[0][2]))/denominator;
      dr[2+k]=-radiation*(g[0][1]+g[0][2])/denominator;
    }
    auto& h=out.enthalpy_partials[k];h[0]=out.exchange_enthalpy[k]+T*(d[0]*g[0][1]+delta*g[1][1]-g[2][0]);
    h[1]=T*(d[1]*g[0][1]+delta*g[0][2]-g[1][1]);
    auto& hr=out.radiation_enthalpy_partials[k];hr[0]=out.radiation_enthalpy[k]+T*(dr[0]*g[0][1]+radiation*g[1][1]);
    hr[1]=T*(dr[1]*g[0][1]+radiation*g[0][2]);
  }
  if(composition_derivatives)for(std::size_t k=0;k<3;++k)if(active[k])for(std::size_t l=k;l<3;++l)if(active[l]) {
    const auto& second=j[second_channel(k,l)];const double fixed=T*(delta*second[0][1]-second[1][0]);
    out.enthalpy_partials[k][2+l]=fixed+T*d[2+l]*j[1+k][0][1];out.enthalpy_partials[l][2+k]=fixed+T*d[2+k]*j[1+l][0][1];
    const double radfixed=T*radiation*second[0][1];
    out.radiation_enthalpy_partials[k][2+l]=radfixed+T*dr[2+l]*j[1+k][0][1];out.radiation_enthalpy_partials[l][2+k]=radfixed+T*dr[2+k]*j[1+l][0][1];
  }
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    if(!std::isfinite(out.exchange_enthalpy[k]+out.radiation_enthalpy[k]))throw std::domain_error("variable EOS: nonfinite enthalpy");
    for(std::size_t l=0;l<5;++l)if(l<2 || (composition_derivatives && active[l-2]))
      if(!std::isfinite(out.enthalpy_partials[k][l]+out.radiation_enthalpy_partials[k][l]+d[l]))
        throw std::domain_error("variable EOS: nonfinite enthalpy derivative");
  }
  return out;
}
} // namespace ember
