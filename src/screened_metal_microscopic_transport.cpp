#include "ember/metal_microscopic_transport.hpp"
#include "ember/constants.hpp"
#include "ember/eos_component.hpp"
#include "ember/detail/differential.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ember {
namespace {
void require(bool ok,const char* why) {
  if(!ok)throw std::domain_error(std::string("screened microscopic transport: ")+why);
}
template<std::size_t N,bool abundances=true> MetalMicroscopicFaceResponse evaluate(
    const VariableMetalHelmholtzEos& eos,const ScreenedCollisionTransport& collision,
    bool ions,double minimum_T,std::array<bool,3> active,double mlo,double mhi,const Point& lo,const Composition& a,
    const Point& hi,const Composition& b,const MetalSpeciesVector* total_rate=nullptr,
    bool radiation_with_redistribution=false) {
  using D=detail::Differential<N>;using detail::exp;using detail::log;
  constexpr std::size_t stride=abundances?7:NVAR;
  require(std::isfinite(mlo) && std::isfinite(mhi) && mlo>=0 && mhi>mlo,"invalid mass interval");
  for(const auto* c:{&a,&b}) {
    require(c->basis==AbundanceBasis::baryon_mass && c->metal_inventory==MetalInventory::gs98,
        "baryonic GS98 mixture required");
    require(std::abs(c->sum()-1)<1e-12,"unnormalized composition");
    for(double x:c->X)require(std::isfinite(x) && x>=0,"invalid abundance");
    require(c->X[2]>0,"positive reference He4 required");
    require(!c->cn_molality || c->cn_mass_convention==CNMassConvention::explicit_metal_mass,
        "CN inventory requires physical helium and total metal mass");
    for(std::size_t k=0;k<3;++k)
      require(active[k]?(k==2?c->Z():c->X[k])>0:(k==2?c->Z():c->X[k])==0,"selected species must be positive; removed species must stay absent");
  }
  if(total_rate)for(std::size_t k=0;k<3;++k)
    require(std::isfinite((*total_rate)[k]) && (active[k] || (*total_rate)[k]==0),
        "total rate must be finite and removed species must have zero rate");
  std::array<D,2> radius,rho,T,lt;
  std::array<std::array<D,3>,2> fraction,phi;
  for(std::size_t e=0;e<2;++e) {
    const auto& p=e?hi:lo;const auto& c=e?b:a;const auto o=stride*e;
    for(double v:{p.lnr,p.lnrho,p.lnT,p.L})require(std::isfinite(v),"nonfinite stellar point");
    radius[e]=exp(D::variable(p.lnr,o));rho[e]=exp(D::variable(p.lnrho,o+1));
    lt[e]=D::variable(p.lnT,o+2);T[e]=exp(lt[e]);
    require(T[e].value>=minimum_T,"temperature below the declared hot domain");
    const auto chemical=eos.composition_potential(T[e].value,rho[e].value,c,active);
    for(std::size_t k=0;k<3;++k) {
      if(!active[k])continue;
      if constexpr(abundances)fraction[e][k]=D::variable(k==2?c.Z():c.X[k],o+4+k);
      else fraction[e][k]=k==2?c.Z():c.X[k];
      phi[e][k]=chemical.gradient[k];
      if constexpr(N>0) {
        phi[e][k].d[o+1]=chemical.dgradient_dlnRho[k];phi[e][k].d[o+2]=chemical.dgradient_dlnT[k];
        if constexpr(abundances)for(std::size_t j=0;j<3;++j)
          if(active[j])phi[e][k].d[o+4+j]=chemical.hessian[k][j];
      }
    }
  }
  const D rb=.5*(radius[0]+radius[1]),rhob=.5*(rho[0]+rho[1]);
  // This temperature makes K*(Thi-Tlo) identical to the inverse-T force
  // in the full material mobility. Radius and density match structure.cpp.
  const D ltb=.5*(lt[0]+lt[1]),Tb=exp(ltb);
  const D xb=.5*(fraction[0][0]+fraction[1][0]),yb=.5*(fraction[0][1]+fraction[1][1]);
  const D zb=.5*(fraction[0][2]+fraction[1][2]);
  Composition c=a;for(std::size_t k=0;k<NSPEC;++k)c.X[k]=.5*(a.X[k]+b.X[k]);
  const double ne=c.mu_elec_inv()*constants::NA*rhob.value;
  const ElectronGas electrons;
  const auto electron=N>0?electrons.eval_with_derivatives(Tb.value,rhob.value,c):electrons.eval(Tb.value,rhob.value,c);
  D stiffness(electron.dP_dlnRho/ne);
  require(std::isfinite(stiffness.value) && stiffness.value>0,"invalid electron stiffness");
  const std::array<D,6> coords{ltb,log(rhob),xb,yb,zb,D(0)};
  if constexpr(N>0)for(std::size_t j=0;j<N;++j) {
    const double dne=rhob.d[j]/rhob.value+(.5*xb.d[j]+yb.d[j]/6+(gs98_ion_moment(1)-.5)*zb.d[j])/c.mu_elec_inv();
    stiffness.d[j]=electron.d2P_dlnTdlnRho/ne*ltb.d[j]
      +(electron.d2P_dlnRho2-electron.dP_dlnRho)/ne*dne;
  }
  D length;
  if constexpr(N>0) {
    const auto screening=ScreenedCollisionTransport::screening_derivatives(Tb.value,rhob.value,xb.value,yb.value,c.Z(),stiffness.value,ions);
    length.value=screening.value;
    for(std::size_t j=0;j<N;++j) {
      for(std::size_t k=0;k<5;++k)length.d[j]+=screening.partials[k]*coords[k].d[j];
      length.d[j]+=screening.partials[5]*stiffness.d[j]/stiffness.value;
    }
  }else length=ScreenedCollisionTransport::screening_length(Tb.value,rhob.value,xb.value,yb.value,c.Z(),stiffness.value,ions);
  BulkMetalCollisionDerivatives kinetic;
  if constexpr(N>0)kinetic=collision.bulk_metal_derivatives(Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value);
  else kinetic.value=collision.bulk_metal_eval(Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value);
  auto coefficient=[&](double value,auto partial) {
    D result(value);
    if constexpr(N>0)for(std::size_t j=0;j<N;++j) {
      for(std::size_t k=0;k<5;++k)result.d[j]+=partial(kinetic.partials[k])*coords[k].d[j];
      result.d[j]+=partial(kinetic.partials[5])*length.d[j]/length.value;
    }
    return result;
  };
  const auto heat=eos.composition_heat(Tb.value,rhob.value,c,active,N>0);
  std::array<D,3> enthalpy,eos_enthalpy,kinetic_enthalpy,radiation_enthalpy,force,rate;
  for(std::size_t k=0;k<3;++k) {
    if(!active[k])continue;
    D h(heat.exchange_enthalpy[k]);
    if constexpr(N>0)for(std::size_t j=0;j<N;++j)for(std::size_t v=0;v<5;++v)
      if(v<2 || active[v-2])h.d[j]+=heat.enthalpy_partials[k][v]*coords[v].d[j];
    eos_enthalpy[k]=h;
    if(total_rate && radiation_with_redistribution) {
      D hr(heat.radiation_enthalpy[k]);
      if constexpr(N>0)for(std::size_t j=0;j<N;++j)for(std::size_t v=0;v<5;++v)
        if(v<2 || active[v-2])hr.d[j]+=heat.radiation_enthalpy_partials[k][v]*coords[v].d[j];
      radiation_enthalpy[k]=hr;
    }
    kinetic_enthalpy[k]=coefficient(kinetic.value.transport_enthalpy[k],[&](const auto& p){return p.transport_enthalpy[k];});
    enthalpy[k]=eos_enthalpy[k]+kinetic_enthalpy[k];
    force[k]=phi[1][k]-phi[0][k]+enthalpy[k]*(T[1]-T[0])/(T[0]*T[1]);
  }
  const D area=4*M_PI*rb*rb,geometry=area*area*rhob/(mhi-mlo);
  D carried;
  for(std::size_t k=0;k<3;++k) {
    for(std::size_t j=0;j<3;++j)
      rate[k]-=geometry*coefficient(kinetic.value.mobility[k][j],[&](const auto& p){return p.mobility[k][j];})*force[j];
    if(total_rate)carried+=eos_enthalpy[k]*(*total_rate)[k]+kinetic_enthalpy[k]*rate[k];
    else carried+=enthalpy[k]*rate[k];
    if(total_rate && radiation_with_redistribution)
      carried+=radiation_enthalpy[k]*((*total_rate)[k]-rate[k]);
  }
  const D conductivity=coefficient(kinetic.value.conductivity,[](const auto& p){return p.conductivity;});
  MetalMicroscopicFaceResponse result;result.carried_luminosity=carried.value;result.conductivity=conductivity.value;
  for(std::size_t k=0;k<3;++k) {
    result.species.rate[k]=rate[k].value;
    if(total_rate) {
      const double h=eos_enthalpy[k].value+(radiation_with_redistribution?radiation_enthalpy[k].value:0);
      if(k<2)result.total_rate_enthalpy[k]=h;else result.total_metal_rate_enthalpy=h;
    }
    if constexpr(N>0 && abundances)for(std::size_t j=0;j<3;++j) {
      result.species.dleft[k][j]=rate[k].d[4+j];result.species.dright[k][j]=rate[k].d[stride+4+j];
    }
  }
  if constexpr(N>0)for(std::size_t j=0;j<NVAR;++j) {
    result.dcarried_lo[j]=carried.d[j];result.dcarried_hi[j]=carried.d[stride+j];
    result.dconductivity_lo[j]=conductivity.d[j];result.dconductivity_hi[j]=conductivity.d[stride+j];
  }
  return result;
}
} // namespace
ScreenedMetalMicroscopicTransport::ScreenedMetalMicroscopicTransport(const VariableMetalHelmholtzEos& eos,
    ScreenedCollisionTransport collision,bool ions,double minimum_T,std::array<bool,3> active,
    bool radiation_with_redistribution):eos_(eos),collisions_(std::move(collision)),
    include_ions_(ions),minimum_temperature_(minimum_T),active_(active),
    radiation_with_redistribution_(radiation_with_redistribution) {
  require(std::isfinite(minimum_T) && minimum_T>0,"positive temperature bound required");
}
MetalMicroscopicFaceResponse ScreenedMetalMicroscopicTransport::metal_eval(std::size_t,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  if(derivatives)return evaluate<14>(eos_,collisions_,include_ions_,minimum_temperature_,active_,mlo,mhi,lo,a,hi,b);
  return evaluate<0>(eos_,collisions_,include_ions_,minimum_temperature_,active_,mlo,mhi,lo,a,hi,b);
}
MicroscopicHeatResponse ScreenedMetalMicroscopicTransport::heat(std::size_t,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  auto present=active_;
  for(std::size_t k=0;k<3;++k)if((k==2?a.Z():a.X[k])==0 && (k==2?b.Z():b.X[k])==0)present[k]=false;
  // A species missing at only one endpoint is still rejected. Removing it
  // there would change the physical face rather than restrict thermal
  // derivatives to a fixed-composition subspace.
  if(derivatives)return evaluate<2*NVAR,false>(eos_,collisions_,include_ions_,minimum_temperature_,present,mlo,mhi,lo,a,hi,b);
  return evaluate<0,false>(eos_,collisions_,include_ions_,minimum_temperature_,present,mlo,mhi,lo,a,hi,b);
}
MicroscopicHeatResponse ScreenedMetalMicroscopicTransport::heat_with_total_metal_rate(std::size_t,
    double mlo,double mhi,const Point& lo,const Composition& a,const Point& hi,const Composition& b,
    const MetalSpeciesVector& total_rate,bool derivatives) const {
  auto present=active_;
  for(std::size_t k=0;k<3;++k)if((k==2?a.Z():a.X[k])==0 && (k==2?b.Z():b.X[k])==0)present[k]=false;
  if(derivatives)return evaluate<2*NVAR,false>(eos_,collisions_,include_ions_,minimum_temperature_,
      present,mlo,mhi,lo,a,hi,b,&total_rate,radiation_with_redistribution_);
  return evaluate<0,false>(eos_,collisions_,include_ions_,minimum_temperature_,
      present,mlo,mhi,lo,a,hi,b,&total_rate,radiation_with_redistribution_);
}
} // namespace ember
