#include "ember/convective_envelope_transport.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
ConvectiveEnvelopeTransport::ConvectiveEnvelopeTransport(const MetalMicroscopicTransport& hot,
    const MetalMicroscopicTransport& convective,double low,double high)
    :hot_(hot),convective_(convective),log_cool_(std::log(low)),log_hot_(std::log(high)) {
  if(!(std::isfinite(low+high) && low>0 && high>low))
    throw std::invalid_argument("convective envelope transport: invalid temperature interval");
}
MetalMicroscopicFaceResponse ConvectiveEnvelopeTransport::metal_eval(std::size_t face,
    double mlo,double mhi,const Point& lo,const Composition& a,const Point& hi,
    const Composition& b,bool derivatives)const {
  if(!(lo.lnT>=log_hot_ && hi.lnT>=log_hot_))
    throw std::domain_error("convective envelope transport: radiative species boundary outside the full hot domain");
  return hot_.metal_eval(face,mlo,mhi,lo,a,hi,b,derivatives);
}
MicroscopicHeatResponse ConvectiveEnvelopeTransport::heat(std::size_t face,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives)const {
  if(!diagnostic_rates.empty() && face>=diagnostic_rates.size())
    throw std::out_of_range("convective envelope transport: missing diagnostic rate");
  return heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,
      diagnostic_rates.empty()?MetalSpeciesVector{}:diagnostic_rates[face],derivatives);
}
MicroscopicHeatResponse ConvectiveEnvelopeTransport::heat_with_total_metal_rate(std::size_t face,
    double mlo,double mhi,const Point& lo,const Composition& a,const Point& hi,const Composition& b,
    const MetalSpeciesVector& rate,bool derivatives)const {
  const auto gate=[&](double lnT)->std::array<double,2> {
    if(!std::isfinite(lnT))throw std::domain_error("convective envelope transport: nonfinite temperature");
    if(lnT<=log_cool_)return {0,0};if(lnT>=log_hot_)return {1,0};
    const double x=(lnT-log_cool_)/(log_hot_-log_cool_);
    return {x*x*x*(10+x*(-15+6*x)),30*x*x*(1-x)*(1-x)/(log_hot_-log_cool_)};
  };
  const auto left=gate(lo.lnT),right=gate(hi.lnT);const double w=left[0]*right[0];
  if(w==0)return convective_.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,rate,derivatives);
  if(w==1)return hot_.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,rate,derivatives);
  const auto cool=convective_.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,rate,derivatives);
  const auto hot=hot_.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,rate,derivatives);
  MicroscopicHeatResponse out;
  out.conductivity=(1-w)*cool.conductivity+w*hot.conductivity;
  out.carried_luminosity=(1-w)*cool.carried_luminosity+w*hot.carried_luminosity;
  for(std::size_t k=0;k<2;++k)
    out.total_rate_enthalpy[k]=(1-w)*cool.total_rate_enthalpy[k]+w*hot.total_rate_enthalpy[k];
  out.total_metal_rate_enthalpy=(1-w)*cool.total_metal_rate_enthalpy+w*hot.total_metal_rate_enthalpy;
  if(derivatives)for(std::size_t k=0;k<NVAR;++k) {
    const double dl=k==2?left[1]*right[0]:0,dh=k==2?left[0]*right[1]:0;
    out.dconductivity_lo[k]=(1-w)*cool.dconductivity_lo[k]+w*hot.dconductivity_lo[k]+dl*(hot.conductivity-cool.conductivity);
    out.dconductivity_hi[k]=(1-w)*cool.dconductivity_hi[k]+w*hot.dconductivity_hi[k]+dh*(hot.conductivity-cool.conductivity);
    out.dcarried_lo[k]=(1-w)*cool.dcarried_lo[k]+w*hot.dcarried_lo[k]+dl*(hot.carried_luminosity-cool.carried_luminosity);
    out.dcarried_hi[k]=(1-w)*cool.dcarried_hi[k]+w*hot.dcarried_hi[k]+dh*(hot.carried_luminosity-cool.carried_luminosity);
  }
  return out;
}
} // namespace ember
