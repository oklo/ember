#pragma once
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace ember {
// Join two actual source families by interpolating log matching T/Pgas
// between fixed metal abundances. Each anchor must lie within its source
// domain. This is interpolation in composition, not time-based smoothing.
class MetalIntervalAtmosphere final : public Atmosphere {
public:
 MetalIntervalAtmosphere(const Eos& eos,const Atmosphere& lower,const Atmosphere& upper,
                         double low,double high):eos_(eos),lower_(lower),upper_(upper),low_(low),high_(high){
  if(!std::isfinite(low+high)||low<0||low>=high||high>=1)
   throw std::invalid_argument("MetalIntervalAtmosphere: invalid interval");
 }
 AtmosphereState eval(double teff,double gravity,const Composition& c)const override{
  if(c.basis!=AbundanceBasis::baryon_mass||c.metal_inventory!=MetalInventory::gs98||
     !std::isfinite(c.sum())||std::abs(c.sum()-1)>1e-10)
   throw std::domain_error("MetalIntervalAtmosphere: normalized physical composition required");
  for(double x:c.X)if(!std::isfinite(x)||x<0||x>1)
   throw std::domain_error("MetalIntervalAtmosphere: invalid composition");
  const double z=c.Z();
  const double eps=8*std::numeric_limits<double>::epsilon();
  // Preserve endpoint-only support under summation roundoff, without
  // changing the physical composition passed to either source family.
  if(z<=low_+eps*low_)return lower_.eval(teff,gravity,c);
  if(z>=high_-eps*high_)return upper_.eval(teff,gravity,c);
  const auto anchor=[&](double value){
   auto q=solar_scaled(c.h1(),value);q.X[1]=c.X[1];q.X[2]-=c.X[1];
   q.basis=AbundanceBasis::baryon_mass;q.metal_inventory=MetalInventory::gs98;
   if(q.X[2]<0)throw std::domain_error("MetalIntervalAtmosphere: source anchor has negative helium");
   return q;
  };
  const auto a=lower_.eval(teff,gravity,anchor(low_)),b=upper_.eval(teff,gravity,anchor(high_));
  if(a.tau!=b.tau||!(a.Pgas>0)||!(b.Pgas>0))
   throw std::domain_error("MetalIntervalAtmosphere: inconsistent matching states");
  const double f=(z-low_)/(high_-low_);
  const auto mix=[&](double x,double y){return (1-f)*x+f*y;};
  AtmosphereState out{};out.T=std::exp(mix(std::log(a.T),std::log(b.T)));
  out.Pgas=std::exp(mix(std::log(a.Pgas),std::log(b.Pgas)));out.tau=a.tau;
  const double ar=constants::a_rad*std::pow(a.T,4)/3,br=constants::a_rad*std::pow(b.T,4)/3;
  const double pr=constants::a_rad*std::pow(out.T,4)/3;out.P=out.Pgas+pr;
  out.dlnT_dlnTeff=mix(a.dlnT_dlnTeff,b.dlnT_dlnTeff);out.dlnT_dlng=mix(a.dlnT_dlng,b.dlnT_dlng);
  const auto dp=[&](double adp,double adt,double bdp,double bdt,double dt){
   const double dpg=mix((a.P*adp-4*ar*adt)/a.Pgas,(b.P*bdp-4*br*bdt)/b.Pgas);
   return (out.Pgas*dpg+4*pr*dt)/out.P;
  };
  out.dlnP_dlnTeff=dp(a.dlnP_dlnTeff,a.dlnT_dlnTeff,b.dlnP_dlnTeff,b.dlnT_dlnTeff,out.dlnT_dlnTeff);
  out.dlnP_dlng=dp(a.dlnP_dlng,a.dlnT_dlng,b.dlnP_dlng,b.dlnT_dlng,out.dlnT_dlng);
  out.rho=eos_.rho_from_PT(out.T,out.P,c,a.rho);return out;
 }
 const char* name()const override{return "interpolated metal-source atmosphere interval";}
private:
 const Eos& eos_;const Atmosphere& lower_;const Atmosphere& upper_;double low_,high_;
};
} // namespace ember
