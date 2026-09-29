#pragma once
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>

namespace ember {
// Join two validated atmosphere families in a declared common region.
// The main-sequence family takes over as gravity rises or hydrogen decreases.
// The composition condition prevents a later expansion from returning the star
// to the initial, nearly solar-composition contraction table. This is a fixed
// function of physical coordinates, never an age-triggered boundary switch.
class AtmosphereOverlap final:public Atmosphere {
public:
  struct Options {
    double log_g_low{},log_g_high{},hydrogen_low{},hydrogen_high{};
  };
  AtmosphereOverlap(const PressureDensity& eos,const Atmosphere& contraction,
                    const Atmosphere& main_sequence,Options options)
      :eos_(eos),contraction_(contraction),main_(main_sequence),o_(options) {
    if(!std::isfinite(o_.log_g_low+o_.log_g_high+o_.hydrogen_low+o_.hydrogen_high) ||
        o_.log_g_high<=o_.log_g_low || o_.hydrogen_low<0 ||
        o_.hydrogen_high>1 || o_.hydrogen_high<=o_.hydrogen_low)
      throw std::invalid_argument("atmosphere overlap: invalid transition region");
  }
  AtmosphereState eval(double teff,double g,const Composition& c)const override {
    if(!(g>0) || !std::isfinite(g+c.h1()) || c.basis!=AbundanceBasis::baryon_mass)
      throw std::domain_error("atmosphere overlap: invalid physical coordinates");
    const double raw=(std::log10(g)-o_.log_g_low)/(o_.log_g_high-o_.log_g_low);
    const double u=std::clamp(raw,0.,1.);
    const auto smooth=[](double x){return x*x*x*(10+x*(-15+6*x));};
    const double h=smooth(std::clamp((c.h1()-o_.hydrogen_low)/(o_.hydrogen_high-o_.hydrogen_low),0.,1.));
    const double w=1-h*(1-smooth(u));
    if(w==0)return contraction_.eval(teff,g,c);
    if(w==1)return main_.eval(teff,g,c);
    const double dw=h*30*u*u*(1-u)*(1-u)/((o_.log_g_high-o_.log_g_low)*std::log(10.));
    const auto a=contraction_.eval(teff,g,c),b=main_.eval(teff,g,c);
    if(a.tau!=b.tau || !(a.T>0 && b.T>0 && a.Pgas>0 && b.Pgas>0))
      throw std::domain_error("atmosphere overlap: incompatible source boundaries");
    const auto mix=[&](double x,double y){return (1-w)*x+w*y;};
    AtmosphereState s;
    const double lt=std::log(b.T/a.T),lp=std::log(b.Pgas/a.Pgas);
    s.T=a.T*std::exp(w*lt);s.Pgas=a.Pgas*std::exp(w*lp);s.tau=a.tau;
    const double ar=constants::a_rad*std::pow(a.T,4)/3,br=constants::a_rad*std::pow(b.T,4)/3;
    const double pr=constants::a_rad*std::pow(s.T,4)/3;s.P=s.Pgas+pr;
    s.dlnT_dlnTeff=mix(a.dlnT_dlnTeff,b.dlnT_dlnTeff);
    s.dlnT_dlng=mix(a.dlnT_dlng,b.dlnT_dlng)+dw*lt;
    const auto pressure=[&](double ap,double at,double bp,double bt,double st,double extra) {
      const double dpg=mix((a.P*ap-4*ar*at)/a.Pgas,(b.P*bp-4*br*bt)/b.Pgas)+extra;
      return (s.Pgas*dpg+4*pr*st)/s.P;
    };
    s.dlnP_dlnTeff=pressure(a.dlnP_dlnTeff,a.dlnT_dlnTeff,b.dlnP_dlnTeff,b.dlnT_dlnTeff,s.dlnT_dlnTeff,0);
    s.dlnP_dlng=pressure(a.dlnP_dlng,a.dlnT_dlng,b.dlnP_dlng,b.dlnT_dlng,s.dlnT_dlng,dw*lp);
    s.rho=eos_.rho_from_PT(s.T,s.P,c,a.rho);
    return s;
  }
  const char* name()const override{return "smooth overlap of contraction and main-sequence atmospheres";}
private:
  const PressureDensity& eos_;const Atmosphere& contraction_;const Atmosphere& main_;Options o_;
};
} // namespace ember
