#include "ember/microscopic_transport.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
void interval(double lo,double hi) {
  if(!(hi>lo) || !(lo>0) || !std::isfinite(hi))
    throw std::invalid_argument("microscopic transport: invalid face mass interval");
}
template<class Array> void finite(const Array& array) {
  for(double x:array)if(!std::isfinite(x))throw std::domain_error("microscopic transport: non-finite face response");
}
void check_heat(const MicroscopicHeatResponse& result,bool derivatives) {
  finite(std::array{result.carried_luminosity,result.conductivity});
  if(result.conductivity<0)throw std::domain_error("microscopic transport: negative conductivity");
  if(derivatives) {
    finite(result.dcarried_lo);finite(result.dcarried_hi);
    finite(result.dconductivity_lo);finite(result.dconductivity_hi);
  }
}
}
MicroscopicHeatResponse MicroscopicTransport::heat_with_total_species_rate(std::size_t,
    double,double,const Point&,const Composition&,const Point&,const Composition&,
    const SpeciesVector&,bool) const {
  throw std::logic_error("microscopic transport: total-species heat is not implemented by this provider");
}
MicroscopicHeatResponse microscopic_heat_with_total_species_rate(const MicroscopicTransport& transport,
    std::size_t face,double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,const SpeciesVector& total_rate,bool derivatives) {
  interval(mass_lo,mass_hi);finite(total_rate);
  if(!transport.uses_total_species_heat())
    throw std::invalid_argument("microscopic transport: total-species heat capability required");
  const auto result=transport.heat_with_total_species_rate(face,mass_lo,mass_hi,lo,comp_lo,
      hi,comp_hi,total_rate,derivatives);
  check_heat(result,derivatives);finite(result.total_rate_enthalpy);return result;
}
MicroscopicFaceResponse microscopic_face(const MicroscopicTransport& transport,std::size_t face,
    double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,bool derivatives) {
  interval(mass_lo,mass_hi);
  const auto result=transport.eval(face,mass_lo,mass_hi,lo,comp_lo,hi,comp_hi,derivatives);
  check_heat(result,derivatives);finite(result.species.rate);
  if(derivatives) {
    for(const auto& row:result.species.dleft)finite(row);
    for(const auto& row:result.species.dright)finite(row);
  }
  return result;
}
MicroscopicHeatResponse microscopic_heat(const MicroscopicTransport& transport,std::size_t face,
    double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,bool derivatives) {
  interval(mass_lo,mass_hi);
  const auto result=transport.heat(face,mass_lo,mass_hi,lo,comp_lo,hi,comp_hi,derivatives);
  check_heat(result,derivatives);return result;
}
} // namespace ember
