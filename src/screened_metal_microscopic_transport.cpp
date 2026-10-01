#include "ember/metal_microscopic_transport.hpp"
#include "ember/constants.hpp"
#include "ember/eos_component.hpp"
#include "ember/detail/differential.hpp"
#include <atomic>
#include <cmath>
#include <iomanip>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ember {
namespace {
void require(bool ok,const char* why) {
  if(!ok)throw std::domain_error(std::string("screened microscopic transport: ")+why);
}
} // namespace
struct ScreenedMetalMicroscopicTransport::EosCache {
  static constexpr std::size_t capacity=8192;
  struct Coordinates {std::array<double,5> x{};std::array<bool,3> active{};bool valid{};};
  struct PotentialSlot {std::mutex mutex;Coordinates at;MetalCompositionPotentialResponse exact;};
  struct HeatSlot {std::mutex mutex;Coordinates at;bool radiation{};MetalCompositionHeatResponse exact;};
  double radius;bool verify;
  std::unique_ptr<PotentialSlot[]> points=std::make_unique<PotentialSlot[]>(capacity);
  std::unique_ptr<HeatSlot[]> faces=std::make_unique<HeatSlot[]>(capacity);
  std::atomic<std::size_t> hits{0},exact{0},verified{0};
  std::mutex worst_mutex;double worst_potential{},worst_enthalpy{};
  EosCache(double r,bool v):radius(r),verify(v) {}
  static std::array<double,5> coordinates(double T,double rho,const Composition& c) {
    return {std::log(T),std::log(rho),c.X[0],c.X[1],c.Z()};
  }
  bool inside(const Coordinates& a,const std::array<double,5>& x,std::array<bool,3> active,std::array<double,5>& d) const {
    if(!a.valid || a.active!=active)return false;
    const double he=1-a.x[2]-a.x[3]-a.x[4],now_he=1-x[2]-x[3]-x[4];
    if(!(he>0 && now_he>0) || std::abs(now_he-he)>radius*he)return false;
    // Absolute in ln T and ln rho; relative in each mass fraction, because the
    // composition potential varies like ln X for a trace species.
    for(std::size_t k=0;k<5;++k) {
      d[k]=x[k]-a.x[k];
      if(k>=2 && !active[k-2] && d[k]!=0)return false;
      if(std::abs(d[k])>(k<2?radius:radius*std::abs(a.x[k])))return false;
    }
    return true;
  }
  MetalCompositionPotentialResponse potential(const VariableMetalHelmholtzEos& eos,std::size_t point,
      double T,double rho,const Composition& c,std::array<bool,3> active) {
    if(point>=capacity)return eos.composition_potential(T,rho,c,active,true);
    const auto x=coordinates(T,rho,c);std::array<double,5> d{};
    auto& slot=points[point];std::unique_lock lock(slot.mutex);
    if(inside(slot.at,x,active,d)) {
      eos.validate_composition_domain(T,rho,c);
      auto out=slot.exact;lock.unlock();
      for(std::size_t k=0;k<3;++k)if(active[k]) {
        out.gradient[k]+=out.dgradient_dlnT[k]*d[0]+out.dgradient_dlnRho[k]*d[1];
        for(std::size_t j=0;j<3;++j)if(active[j])out.gradient[k]+=out.hessian[k][j]*d[2+j];
      }
      hits.fetch_add(1,std::memory_order_relaxed);
      if(verify) {
        const auto e=eos.composition_potential(T,rho,c,active,true);double scale=0,err=0;
        for(std::size_t k=0;k<3;++k)if(active[k]){scale=std::max(scale,std::abs(e.gradient[k]));err=std::max(err,std::abs(out.gradient[k]-e.gradient[k]));}
        verified.fetch_add(1,std::memory_order_relaxed);
        std::lock_guard w(worst_mutex);worst_potential=std::max(worst_potential,scale>0?err/scale:0);
      }
      return out;
    }
    slot.exact=eos.composition_potential(T,rho,c,active,true);slot.at={x,active,true};
    exact.fetch_add(1,std::memory_order_relaxed);
    return slot.exact;
  }
  MetalCompositionHeatResponse heat(const VariableMetalHelmholtzEos& eos,std::size_t face,
      double T,double rho,const Composition& c,std::array<bool,3> active,bool radiation) {
    if(face>=capacity)return eos.composition_heat(T,rho,c,active,true,true);
    const auto x=coordinates(T,rho,c);std::array<double,5> d{};
    auto& slot=faces[face];std::unique_lock lock(slot.mutex);
    if(inside(slot.at,x,active,d)) {
      eos.validate_composition_domain(T,rho,c);
      auto out=slot.exact;lock.unlock();
      for(std::size_t v=0;v<5;++v)if(d[v]!=0) {
        out.material_delta+=out.delta_partials[v]*d[v];
        for(std::size_t k=0;k<3;++k)if(active[k]) {
          out.exchange_enthalpy[k]+=out.enthalpy_partials[k][v]*d[v];
          out.radiation_enthalpy[k]+=out.radiation_enthalpy_partials[k][v]*d[v];
        }
      }
      hits.fetch_add(1,std::memory_order_relaxed);
      if(verify) {
        const auto e=eos.composition_heat(T,rho,c,active,true,true);double scale=0,err=0;
        for(std::size_t k=0;k<3;++k)if(active[k]) {
          scale=std::max({scale,std::abs(e.exchange_enthalpy[k]),radiation?std::abs(e.radiation_enthalpy[k]):0.});
          err=std::max({err,std::abs(out.exchange_enthalpy[k]-e.exchange_enthalpy[k]),
              radiation?std::abs(out.radiation_enthalpy[k]-e.radiation_enthalpy[k]):0.});
        }
        verified.fetch_add(1,std::memory_order_relaxed);
        std::lock_guard w(worst_mutex);worst_enthalpy=std::max(worst_enthalpy,scale>0?err/scale:0);
      }
      return out;
    }
    slot.exact=eos.composition_heat(T,rho,c,active,true,true);slot.at={x,active,true};slot.radiation=radiation;
    exact.fetch_add(1,std::memory_order_relaxed);
    return slot.exact;
  }
};
void ScreenedMetalMicroscopicTransport::use_eos_taylor(double radius,bool verify) {
  require(std::isfinite(radius) && radius>=0 && radius<=.01,"EOS reuse radius must lie in [0,0.01]");
  eos_cache_=radius>0?std::make_shared<EosCache>(radius,verify):nullptr;
}
ScreenedMetalMicroscopicTransport::EosReuse ScreenedMetalMicroscopicTransport::eos_reuse_statistics() const {
  if(!eos_cache_)return {};
  std::lock_guard w(eos_cache_->worst_mutex);
  return {eos_cache_->hits.load(),eos_cache_->exact.load(),eos_cache_->verified.load(),
          eos_cache_->worst_potential,eos_cache_->worst_enthalpy};
}
void ScreenedMetalMicroscopicTransport::use_phase_mobility(const ColdHeliumOptions& phase,double remaining) {
  require(phase.mixture_phase,"solid mobility requires the same-composition phase EOS");
  require(std::isfinite(remaining)&&remaining>=0&&remaining<=1,"solid mobility fraction must lie in [0,1]");
  phase_options_=phase;solid_mobility_=remaining;
}

namespace {
template<std::size_t N,bool abundances=true> MetalMicroscopicFaceResponse evaluate(
    const VariableMetalHelmholtzEos& eos,const ScreenedCollisionTransport& collision,
    const CollisionTaylorCache* cache,ScreenedMetalMicroscopicTransport::EosCache* eos_cache,
    std::size_t face,bool ions,double minimum_T,std::array<bool,3> active,double mlo,double mhi,const Point& lo,const Composition& a,
    const Point& hi,const Composition& b,const MetalSpeciesVector* total_rate=nullptr,
    bool radiation_with_redistribution=false,const ColdHeliumOptions* phase=nullptr,double solid_mobility=1) {
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
    if(T[e].value<minimum_T) {
      std::ostringstream message;
      message << std::setprecision(4) << "screened microscopic transport: temperature below the declared hot domain"
        << " (face=" << face << ", endpoint=" << e << ", T=" << T[e].value
        << ", minimum=" << minimum_T << ", rho=" << rho[e].value
        << ", X=" << c.X[0] << ", Y3=" << c.X[1] << ", Z=" << c.Z() << ')';
      throw std::domain_error(message.str());
    }
    const auto chemical=eos_cache?eos_cache->potential(eos,face+e,T[e].value,rho[e].value,c,active)
        :eos.composition_potential(T[e].value,rho[e].value,c,active,N>0 && abundances);
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
  const auto c=mean_composition(a,b);
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
  if constexpr(N>0)kinetic=cache?cache->derivatives(collision,face,Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value)
      :collision.bulk_metal_derivatives(Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value);
  else kinetic.value=cache?cache->value(collision,face,Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value)
      :collision.bulk_metal_eval(Tb.value,rhob.value,xb.value,yb.value,c.Z(),length.value);
  auto coefficient=[&](double value,auto partial) {
    D result(value);
    if constexpr(N>0)for(std::size_t j=0;j<N;++j) {
      for(std::size_t k=0;k<5;++k)result.d[j]+=partial(kinetic.partials[k])*coords[k].d[j];
      result.d[j]+=partial(kinetic.partials[5])*length.d[j]/length.value;
    }
    return result;
  };
  const auto heat=eos_cache?eos_cache->heat(eos,face,Tb.value,rhob.value,c,active,radiation_with_redistribution)
      :eos.composition_heat(Tb.value,rhob.value,c,active,N>0,abundances);
  std::array<D,3> enthalpy,eos_enthalpy,kinetic_enthalpy,radiation_enthalpy,force,rate;
  for(std::size_t k=0;k<3;++k) {
    if(!active[k])continue;
    D h(heat.exchange_enthalpy[k]);
    if constexpr(N>0)for(std::size_t j=0;j<N;++j)for(std::size_t v=0;v<(abundances?5:2);++v)
      if(v<2 || active[v-2])h.d[j]+=heat.enthalpy_partials[k][v]*coords[v].d[j];
    eos_enthalpy[k]=h;
    if(total_rate && radiation_with_redistribution) {
      D hr(heat.radiation_enthalpy[k]);
      if constexpr(N>0)for(std::size_t j=0;j<N;++j)for(std::size_t v=0;v<(abundances?5:2);++v)
        if(v<2 || active[v-2])hr.d[j]+=heat.radiation_enthalpy_partials[k][v]*coords[v].d[j];
      radiation_enthalpy[k]=hr;
    }
    kinetic_enthalpy[k]=coefficient(kinetic.value.transport_enthalpy[k],[&](const auto& p){return p.transport_enthalpy[k];});
    enthalpy[k]=eos_enthalpy[k]+kinetic_enthalpy[k];
    force[k]=phi[1][k]-phi[0][k]+enthalpy[k]*(T[1]-T[0])/(T[0]*T[1]);
  }
  const D area=4*M_PI*rb*rb,geometry=area*area*rhob/(mhi-mlo);
  D mobility_factor(1);
  if(phase && solid_mobility<1) {
    const auto w=cold_helium_solid_response(Tb.value,rhob.value,c,*phase);
    mobility_factor.value=1-(1-solid_mobility)*w[0];
    if constexpr(N>0)for(std::size_t j=0;j<N;++j)for(std::size_t v=0;v<5;++v)
      mobility_factor.d[j]-=(1-solid_mobility)*w[1+v]*coords[v].d[j];
  }
  D carried;
  for(std::size_t k=0;k<3;++k) {
    for(std::size_t j=0;j<3;++j)
      rate[k]-=geometry*coefficient(kinetic.value.mobility[k][j],[&](const auto& p){return p.mobility[k][j];})*force[j];
    rate[k]=rate[k]*mobility_factor;
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
MetalMicroscopicFaceResponse ScreenedMetalMicroscopicTransport::metal_eval(std::size_t face,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  if(derivatives)return evaluate<14>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,active_,mlo,mhi,lo,a,hi,b,nullptr,false,phase_options_?&*phase_options_:nullptr,solid_mobility_);
  return evaluate<0>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,active_,mlo,mhi,lo,a,hi,b,nullptr,false,phase_options_?&*phase_options_:nullptr,solid_mobility_);
}
MicroscopicHeatResponse ScreenedMetalMicroscopicTransport::heat(std::size_t face,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  auto present=active_;
  for(std::size_t k=0;k<3;++k)if((k==2?a.Z():a.X[k])==0 && (k==2?b.Z():b.X[k])==0)present[k]=false;
  // A species missing at only one endpoint is still rejected. Removing it
  // there would change the physical face rather than restrict thermal
  // derivatives to a fixed-composition subspace.
  if(derivatives)return evaluate<2*NVAR,false>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,present,mlo,mhi,lo,a,hi,b,nullptr,false,phase_options_?&*phase_options_:nullptr,solid_mobility_);
  return evaluate<0,false>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,present,mlo,mhi,lo,a,hi,b,nullptr,false,phase_options_?&*phase_options_:nullptr,solid_mobility_);
}
MicroscopicHeatResponse ScreenedMetalMicroscopicTransport::heat_with_total_metal_rate(std::size_t face,
    double mlo,double mhi,const Point& lo,const Composition& a,const Point& hi,const Composition& b,
    const MetalSpeciesVector& total_rate,bool derivatives) const {
  auto present=active_;
  for(std::size_t k=0;k<3;++k)if((k==2?a.Z():a.X[k])==0 && (k==2?b.Z():b.X[k])==0)present[k]=false;
  if(derivatives)return evaluate<2*NVAR,false>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,
      present,mlo,mhi,lo,a,hi,b,&total_rate,radiation_with_redistribution_,phase_options_?&*phase_options_:nullptr,solid_mobility_);
  return evaluate<0,false>(eos_,collisions_,taylor_.get(),eos_cache_.get(),face,include_ions_,minimum_temperature_,
      present,mlo,mhi,lo,a,hi,b,&total_rate,radiation_with_redistribution_,phase_options_?&*phase_options_:nullptr,solid_mobility_);
}
} // namespace ember
