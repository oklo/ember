#pragma once
#include "ember/opacity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {

// Selected bounded approximation, with independent half/double-opacity
// stellar controls in docs/results/dense_hydrogen_temperature_controls_v2.json.
// Blend the established hydrogen ratio with the retained source's measured
// log-opacity hydrogen slope before the ratio source reaches log R = 1.5.
class DenseHydrogenOpacity final: public Opacity {
 const Opacity &ratio_, &slope_;
 static constexpr double start=1.25,end=1.45,maximum=1.8;
 static double logR(double T,double rho){return std::log10(rho)-3*(std::log10(T)-6);}
 static void domain(double T,double r){
   if(!(T>=std::pow(10.,5.6) && T<=std::pow(10.,6.1) && r<=maximum))
     throw std::domain_error("DenseHydrogenOpacity: bounded approximation requires logT=5.6..6.1, logR<=1.8");
 }
 public:
 DenseHydrogenOpacity(const Opacity& ratio,const Opacity& slope):ratio_(ratio),slope_(slope){}
 OpacityState eval(double T,double rho,const Composition& c)const override {
   const double r=logR(T,rho);
   if(c.h1()<=.75 || r<=start)return ratio_.eval(T,rho,c);
   domain(T,r);
   auto b=slope_.eval(T,rho,c);
   if(r>=end)return b;
   const auto a=ratio_.eval(T,rho,c);const double u=(r-start)/(end-start);
   const double w=u*u*(3-2*u),dw=6*u*(1-u)/((end-start)*std::log(10.));
   const double contrast=std::log(b.kappa/a.kappa);
   return {std::exp(std::log(a.kappa)+w*contrast),
     (1-w)*a.dlnk_dlnT+w*b.dlnk_dlnT-3*dw*contrast,
     (1-w)*a.dlnk_dlnRho+w*b.dlnk_dlnRho+dw*contrast,
     (1-w)*a.dlnk_dX+w*b.dlnk_dX,
     (1-w)*a.dlnk_dZ+w*b.dlnk_dZ,
     (1-w)*a.dlnk_dY3+w*b.dlnk_dY3};
 }
 std::optional<DensityRange> density_range(double T,const Composition& c)const override {
   if(c.h1()<=.75)return ratio_.density_range(T,c);
   if(T<std::pow(10.,5.6) || T>std::pow(10.,6.1)) {
     auto a=ratio_.density_range(T,c);
     if(a)a->max=std::min(a->max,std::pow(10.,start)*std::pow(T/1e6,3));
     if(a && !(a->min<a->max))throw std::domain_error("DenseHydrogenOpacity: no bounded density interval");
     return a;
   }
   const auto a=ratio_.density_range(T,c),b=slope_.density_range(T,c);
   if(!a || !b)throw std::domain_error("DenseHydrogenOpacity: both source domains required");
   const double lo=std::pow(10.,start)*std::pow(T/1e6,3),hi=std::pow(10.,end)*std::pow(T/1e6,3);
   if(a->max<hi || b->min>lo || b->max<hi)
     throw std::domain_error("DenseHydrogenOpacity: missing source overlap");
   return DensityRange{a->min,std::min(b->max,std::pow(10.,maximum)*std::pow(T/1e6,3))};
 }
 const char* name()const override{return "bounded dense hydrogen source-slope approximation";}
};
} // namespace ember
