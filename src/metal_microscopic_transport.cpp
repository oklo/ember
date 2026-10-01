#include "ember/metal_microscopic_transport.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
template<class A> void finite(const A& a) {
  for(double x:a)if(!std::isfinite(x))
    throw std::domain_error("metal microscopic transport: nonfinite response");
}
void interval(double lo,double hi) {
  if(!(hi>lo && lo>0) || !std::isfinite(hi))
    throw std::invalid_argument("metal microscopic transport: invalid face mass interval");
}
void heat_check(const MicroscopicHeatResponse& r,bool derivatives) {
  finite(std::array{r.carried_luminosity,r.conductivity});
  if(r.conductivity<0)throw std::domain_error("metal microscopic transport: negative conductivity");
  if(derivatives) {
    finite(r.dcarried_lo);finite(r.dcarried_hi);
    finite(r.dconductivity_lo);finite(r.dconductivity_hi);
  }
}
}
MicroscopicFaceResponse MetalMicroscopicTransport::eval(std::size_t,double,double,
    const Point&,const Composition&,const Point&,const Composition&,bool) const {
  throw std::logic_error("metal microscopic transport requires the three-mass interface");
}
MicroscopicHeatResponse MetalMicroscopicTransport::heat(std::size_t face,double mlo,double mhi,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  return metal_eval(face,mlo,mhi,lo,a,hi,b,derivatives);
}
MetalMicroscopicFaceResponse metal_microscopic_face(const MetalMicroscopicTransport& transport,
    std::size_t face,double mlo,double mhi,const Point& lo,const Composition& a,
    const Point& hi,const Composition& b,bool derivatives) {
  interval(mlo,mhi);
  const auto r=transport.metal_eval(face,mlo,mhi,lo,a,hi,b,derivatives);
  heat_check(r,derivatives);finite(r.species.rate);
  if(derivatives) {
    for(const auto& row:r.species.dleft)finite(row);
    for(const auto& row:r.species.dright)finite(row);
  }
  return r;
}
MicroscopicHeatResponse microscopic_heat_with_total_metal_rate(const MetalMicroscopicTransport& transport,
    std::size_t face,double mlo,double mhi,const Point& lo,const Composition& a,
    const Point& hi,const Composition& b,const MetalSpeciesVector& total,bool derivatives) {
  interval(mlo,mhi);finite(total);
  const auto r=transport.heat_with_total_metal_rate(face,mlo,mhi,lo,a,hi,b,total,derivatives);
  heat_check(r,derivatives);finite(r.total_rate_enthalpy);
  finite(std::array{r.total_metal_rate_enthalpy});return r;
}
MetalSpeciesFaceResponse metal_species_face(const MetalMicroscopicTransport& transport,
    std::size_t face,double mlo,double mhi,const Point& lo,const Composition& a,
    const Point& hi,const Composition& b,bool derivatives) {
  interval(mlo,mhi);const auto r=transport.species(face,mlo,mhi,lo,a,hi,b,derivatives);
  finite(r.rate);
  if(derivatives){for(const auto& row:r.dleft)finite(row);for(const auto& row:r.dright)finite(row);}
  return r;
}
} // namespace ember
