#pragma once
#include "ember/metal_microscopic_transport.hpp"
#include <cmath>
#include <span>
#include <utility>

namespace ember {
// Species always use the microscopic law, including its domain checks.
// Heat approaches the supplied material conduction/enthalpy law in the cool
// envelope. There the omitted kinetic heat must be assessed separately, and
// convection must mix any boundary outside the microscopic species domain.
class EnvelopeTransport final : public MetalMicroscopicTransport {
public:
  EnvelopeTransport(const MetalMicroscopicTransport& material,
      const MetalMicroscopicTransport& microscopic,double lower_T,double upper_T)
      :material_(material),microscopic_(microscopic),lower_(std::log(lower_T)),
       width_(std::log(upper_T/lower_T)) {
    if(!std::isfinite(lower_T) || !std::isfinite(upper_T) || !(upper_T>lower_T && lower_T>0))
      throw std::invalid_argument("envelope transport: invalid temperature overlap");
  }
  // Used for assessments at an accepted model; evolution passes rates directly.
  std::span<const MetalSpeciesVector> diagnostic_rates;
  bool requires_positive_species_guess() const override {
    return microscopic_.requires_positive_species_guess();
  }
  MetalMicroscopicFaceResponse metal_eval(std::size_t i,double a,double b,
      const Point& p,const Composition& x,const Point& q,const Composition& y,bool d) const override {
    return microscopic_.metal_eval(i,a,b,p,x,q,y,d);
  }
  MicroscopicHeatResponse heat(std::size_t i,double a,double b,
      const Point& p,const Composition& x,const Point& q,const Composition& y,bool d) const override {
    if(!diagnostic_rates.empty()) {
      if(i>=diagnostic_rates.size())throw std::out_of_range("envelope transport: missing diagnostic rate");
      return heat_with_total_metal_rate(i,a,b,p,x,q,y,diagnostic_rates[i],d);
    }
    return join(p,q,d,[&](const auto& source){return source.heat(i,a,b,p,x,q,y,d);});
  }
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t i,double a,double b,
      const Point& p,const Composition& x,const Point& q,const Composition& y,
      const MetalSpeciesVector& rate,bool d) const override {
    return join(p,q,d,[&](const auto& source){return source.heat_with_total_metal_rate(i,a,b,p,x,q,y,rate,d);});
  }
  const char* name() const override {return "hot microscopic transport with cool material conduction and mixing heat";}
private:
  const MetalMicroscopicTransport& material_;
  const MetalMicroscopicTransport& microscopic_;
  double lower_,width_;
  std::pair<double,double> weight(double logT) const {
    if(!std::isfinite(logT))throw std::domain_error("envelope transport: invalid temperature");
    const double u=(logT-lower_)/width_;
    if(u<=0)return {0,0};
    if(u>=1)return {1,0};
    return {u*u*u*(10+u*(-15+6*u)),30*u*u*(1-u)*(1-u)/width_};
  }
  template<class Evaluate> MicroscopicHeatResponse join(const Point& p,const Point& q,
      bool derivatives,const Evaluate& evaluate) const {
    const auto [a,da]=weight(p.lnT);const auto [b,db]=weight(q.lnT);
    const double w=a*b;
    if(w==0)return evaluate(material_);
    if(w==1)return evaluate(microscopic_);
    const auto cold=evaluate(material_),hot=evaluate(microscopic_);
    const auto mix=[&](double c,double h){return (1-w)*c+w*h;};
    MicroscopicHeatResponse r;
    r.carried_luminosity=mix(cold.carried_luminosity,hot.carried_luminosity);
    r.conductivity=mix(cold.conductivity,hot.conductivity);
    for(std::size_t k=0;k<2;++k)r.total_rate_enthalpy[k]=mix(cold.total_rate_enthalpy[k],hot.total_rate_enthalpy[k]);
    r.total_metal_rate_enthalpy=mix(cold.total_metal_rate_enthalpy,hot.total_metal_rate_enthalpy);
    if(derivatives)for(std::size_t k=0;k<NVAR;++k) {
      r.dcarried_lo[k]=mix(cold.dcarried_lo[k],hot.dcarried_lo[k]);
      r.dcarried_hi[k]=mix(cold.dcarried_hi[k],hot.dcarried_hi[k]);
      r.dconductivity_lo[k]=mix(cold.dconductivity_lo[k],hot.dconductivity_lo[k]);
      r.dconductivity_hi[k]=mix(cold.dconductivity_hi[k],hot.dconductivity_hi[k]);
      if(k==2) {
        r.dcarried_lo[k]+=da*b*(hot.carried_luminosity-cold.carried_luminosity);
        r.dcarried_hi[k]+=a*db*(hot.carried_luminosity-cold.carried_luminosity);
        r.dconductivity_lo[k]+=da*b*(hot.conductivity-cold.conductivity);
        r.dconductivity_hi[k]+=a*db*(hot.conductivity-cold.conductivity);
      }
    }
    return r;
  }
};
} // namespace ember
