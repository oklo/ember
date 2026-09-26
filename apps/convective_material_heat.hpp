#pragma once
#include "ember/metal_microscopic_transport.hpp"
#include "ember/conduction.hpp"
#include "ember/constants.hpp"
#include "../src/differential.hpp"
#include <span>

namespace ember::driver {
// Leading, well-mixed limit for a closed, wholly convective star. Total
// composition enthalpy is retained. The caller must check mixing times and
// assess the omitted microscopic drift/kinetic heat. This is NOT a collision
// law for cool radiative layers: any species boundary explicitly rejects.
class ConvectiveMaterialHeat final:public MetalMicroscopicTransport {
 public:
  ConvectiveMaterialHeat(const VariableMetalHelmholtzEos& eos,const Conduction& conduction)
      :eos_(eos),conduction_(conduction) {}
  // Accepted rates for diagnostic queries. Evolution supplies its own rates
  // explicitly and converges them with the structure.
  std::span<const MetalSpeciesVector> diagnostic_rates;
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool)const override {
    throw std::domain_error("convective material heat: radiative species boundary needs microscopic transport");
  }
  MicroscopicHeatResponse heat(std::size_t face,double a,double b,const Point& p,
      const Composition& x,const Point& q,const Composition& y,bool derivatives)const override {
    if(!diagnostic_rates.empty() && face>=diagnostic_rates.size())
      throw std::out_of_range("convective material heat: missing diagnostic face rate");
    return heat_with_total_metal_rate(face,a,b,p,x,q,y,
        diagnostic_rates.empty()?MetalSpeciesVector{}:diagnostic_rates[face],derivatives);
  }
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,
      const Point& p,const Composition& x,const Point& q,const Composition& y,
      const MetalSpeciesVector& rate,bool derivatives)const override {
    using D=detail::Differential<8>;using detail::exp;using detail::log;
    if(x[Species::H2]!=0 || y[Species::H2]!=0)
      throw std::domain_error("convective material heat: initial D needs its fourth heat rate");
    const D lt=.5*(D::variable(p.lnT,2)+D::variable(q.lnT,6)),T=exp(lt);
    const D rho=.5*(exp(D::variable(p.lnrho,1))+exp(D::variable(q.lnrho,5))),lr=log(rho);
    auto c=x;for(std::size_t k=0;k<NSPEC;++k)c.X[k]=.5*(x.X[k]+y.X[k]);
    if(x.cn_molality.has_value()!=y.cn_molality.has_value())
      throw std::domain_error("convective material heat: inconsistent CN inventory");
    if(c.cn_molality)for(std::size_t k=0;k<3;++k)
      (*c.cn_molality)[k]=.5*((*x.cn_molality)[k]+(*y.cn_molality)[k]);
    const auto k=conduction_.eval(T.value,rho.value,c);
    D K(4*constants::a_rad*constants::c*std::pow(T.value,3)/(3*rho.value*k.kappa));
    if(derivatives)for(std::size_t v=0;v<8;++v)
      K.d[v]=K.value*((3-k.dlnk_dlnT)*lt.d[v]-(1+k.dlnk_dlnRho)*lr.d[v]);
    const std::array<bool,3> active{c.X[0]>0,c.X[1]>0,c.Z()>0};
    const auto h=eos_.composition_heat(T.value,rho.value,c,active,derivatives);
    D Q;MicroscopicHeatResponse out;
    for(std::size_t i=0;i<3;++i) {
      if(!active[i] && rate[i]!=0)throw std::domain_error("convective material heat: absent transported species");
      const double value=active[i]?h.exchange_enthalpy[i]+h.radiation_enthalpy[i]:0;
      D H(value);
      if(derivatives && active[i])for(std::size_t v=0;v<8;++v)
        H.d[v]=(h.enthalpy_partials[i][0]+h.radiation_enthalpy_partials[i][0])*lt.d[v]
              +(h.enthalpy_partials[i][1]+h.radiation_enthalpy_partials[i][1])*lr.d[v];
      Q+=H*rate[i];
      if(i<2)out.total_rate_enthalpy[i]=value;else out.total_metal_rate_enthalpy=value;
    }
    out.conductivity=K.value;out.carried_luminosity=Q.value;
    if(derivatives)for(std::size_t i=0;i<4;++i) {
      out.dcarried_lo[i]=Q.d[i];out.dcarried_hi[i]=Q.d[4+i];
      out.dconductivity_lo[i]=K.d[i];out.dconductivity_hi[i]=K.d[4+i];
    }
    return out;
  }
  const char* name()const override{return "whole-star convective limit with material conduction and composition enthalpy";}
 private:
  const VariableMetalHelmholtzEos& eos_;const Conduction& conduction_;
};
} // namespace ember::driver
