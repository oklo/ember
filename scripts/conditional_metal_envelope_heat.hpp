#pragma once
// Whole-star integration experiment, not an accepted cool transport law.
// Species exchange requires the native hot domain at every region boundary.
// EOS enthalpy follows the total species flux throughout the star. Below the
// stated temperature join the reduced microscopic heat is set to zero and
// the retained material conductivity is used. No error bound for the omitted
// cool kinetic term is implied. A new cold radiative region rejects the step.
#include "ember/metal_microscopic_transport.hpp"
#include "ember/opacity.hpp"
#include "ember/conduction.hpp"
#include "ember/constants.hpp"
#include "../src/differential.hpp"
#include <span>
#include <atomic>
#include <stdexcept>

class ConditionalMetalEnvelopeHeat final:public ember::MetalMicroscopicTransport {
 public:
  ConditionalMetalEnvelopeHeat(const ember::VariableMetalHelmholtzEos& eos,
      const ember::ScreenedCollisionTransport& collisions,const ember::Conduction& conduction,
      bool ions,double lower=2e6,double upper=3e6,bool radiation_with_redistribution=false):eos_(eos),conduction_(conduction),
      hot_(eos,collisions,ions,lower,{true,true,true},radiation_with_redistribution),lower_(lower),upper_(upper),
      radiation_with_redistribution_(radiation_with_redistribution) {
    if(!(upper>lower && lower>0))throw std::invalid_argument("conditional heat: invalid temperature join");
  }
  mutable std::atomic<std::size_t> maximum_species_face{};
  // Explicit fully convective control: material conduction and composition
  // enthalpy remain, while reduced microscopic heat is omitted. Any region
  // boundary needing microscopic species exchange rejects rather than
  // extrapolating the collision table or switching physics automatically.
  bool fully_convective_only{};
  // Used only to reproduce the final step's diagnostics after evolution.
  std::span<const ember::MetalSpeciesVector> prescribed;
  bool requires_positive_species_guess()const override{return true;}
  const char* name()const override{return fully_convective_only?
    "fully convective material conduction and composition enthalpy":
    "conditional cool-heat omission with total species enthalpy";}
  ember::MetalMicroscopicFaceResponse metal_eval(std::size_t face,double mlo,double mhi,const ember::Point& lo,
      const ember::Composition& a,const ember::Point& hi,const ember::Composition& b,bool derivatives)const override {
    auto recorded=maximum_species_face.load(std::memory_order_relaxed);
    while(recorded<face && !maximum_species_face.compare_exchange_weak(recorded,face,std::memory_order_relaxed)) {}
    if(fully_convective_only)throw std::domain_error("fully convective transport: radiative species boundary requires microscopic transport");
    return hot_.metal_eval(face,mlo,mhi,lo,a,hi,b,derivatives);
  }
  ember::MicroscopicHeatResponse heat(std::size_t face,double mlo,double mhi,const ember::Point& lo,
      const ember::Composition& a,const ember::Point& hi,const ember::Composition& b,bool derivatives)const override {
    if(!prescribed.empty()) {
      if(face>=prescribed.size())throw std::out_of_range("conditional heat: missing prescribed metal rate");
      return heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,prescribed[face],derivatives);
    }
    if(derivatives)return evaluate<8>(face,mlo,mhi,lo,a,hi,b,nullptr);
    return evaluate<0>(face,mlo,mhi,lo,a,hi,b,nullptr);
  }
  ember::MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t face,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,const ember::Point& hi,const ember::Composition& b,
      const ember::MetalSpeciesVector& total,bool derivatives)const override {
    if(derivatives)return evaluate<8>(face,mlo,mhi,lo,a,hi,b,&total);
    return evaluate<0>(face,mlo,mhi,lo,a,hi,b,&total);
  }
 private:
  const ember::VariableMetalHelmholtzEos& eos_;
  const ember::Conduction& conduction_;
  ember::ScreenedMetalMicroscopicTransport hot_;
  double lower_,upper_;
  bool radiation_with_redistribution_;
  template<std::size_t N> ember::MicroscopicHeatResponse evaluate(std::size_t face,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,const ember::Point& hi,const ember::Composition& b,
      const ember::MetalSpeciesVector* total)const {
    using D=ember::detail::Differential<N>;using ember::detail::exp;using ember::detail::log;
    const D ltlo=D::variable(lo.lnT,2),lthi=D::variable(hi.lnT,6);
    const auto gate=[&](const D& lt)->D {
      if(lt.value<=std::log(lower_))return D(0);
      if(lt.value>=std::log(upper_))return D(1);
      const D z=(lt-std::log(lower_))/std::log(upper_/lower_);
      return z*z*z*(10+z*(-15+6*z));
    };
    const D weight=fully_convective_only?D(0):gate(ltlo)*gate(lthi);
    ember::MicroscopicHeatResponse hot;
    if(weight.value>0) {
      hot=total?hot_.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,*total,N>0):
                hot_.heat(face,mlo,mhi,lo,a,hi,b,N>0);
      if(weight.value==1)return hot;
    }
    const D lt=.5*(ltlo+lthi),T=exp(lt);
    const D rho=.5*(exp(D::variable(lo.lnrho,1))+exp(D::variable(hi.lnrho,5))),lrho=log(rho);
    auto c=a;for(std::size_t k=0;k<ember::NSPEC;++k)c.X[k]=.5*(a.X[k]+b.X[k]);
    const auto material=conduction_.eval(T.value,rho.value,c);
    if(!(material.kappa>0))throw std::domain_error("conditional heat: invalid material opacity");
    D cold_K(4*ember::constants::a_rad*ember::constants::c*T.value*T.value*T.value/(3*rho.value*material.kappa));
    if constexpr(N>0)for(std::size_t v=0;v<N;++v)
      cold_K.d[v]=cold_K.value*((3-material.dlnk_dlnT)*lt.d[v]-(1+material.dlnk_dlnRho)*lrho.d[v]);
    D cold_Q;ember::MetalSpeciesVector h{};
    if(total) {
      const std::array<bool,3> active{c.X[0]>0,c.X[1]>0,c.Z()>0};
      for(std::size_t k=0;k<3;++k)if(!active[k] && (*total)[k]!=0)
        throw std::domain_error("conditional heat: nonzero rate of absent species");
      const auto e=eos_.composition_heat(T.value,rho.value,c,active,N>0);h=e.exchange_enthalpy;
      for(std::size_t k=0;k<3;++k) {
        if(!active[k]){h[k]=0;continue;}
        if(radiation_with_redistribution_)h[k]+=e.radiation_enthalpy[k];
        D enthalpy(h[k]);
        if constexpr(N>0)for(std::size_t v=0;v<N;++v)
          enthalpy.d[v]=e.enthalpy_partials[k][0]*lt.d[v]+e.enthalpy_partials[k][1]*lrho.d[v];
        if constexpr(N>0)if(radiation_with_redistribution_)for(std::size_t v=0;v<N;++v)
          enthalpy.d[v]+=e.radiation_enthalpy_partials[k][0]*lt.d[v]+e.radiation_enthalpy_partials[k][1]*lrho.d[v];
        cold_Q+=enthalpy*(*total)[k];
      }
    }
    auto thermal=[&](double value,const auto& left,const auto& right) {
      D d(value);
      if constexpr(N>0)for(std::size_t v=0;v<ember::NVAR;++v){d.d[v]=left[v];d.d[4+v]=right[v];}
      return d;
    };
    const D Q=weight*thermal(hot.carried_luminosity,hot.dcarried_lo,hot.dcarried_hi)+(1-weight)*cold_Q;
    const D K=weight*thermal(hot.conductivity,hot.dconductivity_lo,hot.dconductivity_hi)+(1-weight)*cold_K;
    ember::MicroscopicHeatResponse out;out.carried_luminosity=Q.value;out.conductivity=K.value;
    for(std::size_t k=0;k<2;++k)out.total_rate_enthalpy[k]=weight.value*hot.total_rate_enthalpy[k]+(1-weight.value)*h[k];
    out.total_metal_rate_enthalpy=weight.value*hot.total_metal_rate_enthalpy+(1-weight.value)*h[2];
    if constexpr(N>0)for(std::size_t v=0;v<ember::NVAR;++v) {
      out.dcarried_lo[v]=Q.d[v];out.dcarried_hi[v]=Q.d[v+4];
      out.dconductivity_lo[v]=K.d[v];out.dconductivity_hi[v]=K.d[v+4];
    }
    return out;
  }
};
