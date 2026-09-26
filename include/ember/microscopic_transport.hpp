#pragma once
#include "ember/model.hpp"
#include "ember/species_flux.hpp"

namespace ember {
struct MicroscopicHeatResponse {
  double carried_luminosity{}; // h dot species_rate, erg/s, positive outward
  double conductivity{}; // zero-species-flux conductivity, erg/(cm s K)
  // Derivatives with respect to [ln r, ln rho, ln T, L] at fixed composition.
  std::array<double,NVAR> dcarried_lo{},dcarried_hi{},dconductivity_lo{},dconductivity_hi{};
  // Only defined by heat_with_total_species_rate: exact linear derivative
  // of carried luminosity with respect to the prescribed H1/He3 rates.
  SpeciesVector total_rate_enthalpy{};
  // Additional derivative for providers that transport total metal mass.
  double total_metal_rate_enthalpy{};
};
struct MicroscopicFaceResponse : MicroscopicHeatResponse {
  SpeciesFaceResponse species;
};

// One consistent physical evaluator supplies species transport and the heat
// split. Temperature, density, composition and geometry refer to the trial
// state. The caller must not substitute saved coefficients for an evolution
// closure. Implementations declare/check their physical domain explicitly.
// Closed exterior species fluxes and a radiative outer boundary are assumed;
// no material heat may leave through that boundary.
class MicroscopicTransport {
public:
  virtual ~MicroscopicTransport()=default;
  virtual MicroscopicFaceResponse eval(std::size_t face,double mass_lo,double mass_hi,
      const Point& lo,const Composition& comp_lo,const Point& hi,const Composition& comp_hi,
      bool derivatives) const=0;
  // Thermal variations hold composition fixed. A provider can define this
  // response on a subspace with absent species even when derivatives that
  // introduce those species are undefined. No species Jacobian is returned.
  virtual MicroscopicHeatResponse heat(std::size_t face,double mass_lo,double mass_hi,
      const Point& lo,const Composition& comp_lo,const Point& hi,const Composition& comp_hi,
      bool derivatives) const {
    return eval(face,mass_lo,mass_hi,lo,comp_lo,hi,comp_hi,derivatives);
  }
  virtual const char* name() const=0;
  // Providers with this capability account for EOS enthalpy transported by
  // the TOTAL species rate (microscopic drift plus other redistribution).
  // Kinetic transported heat still follows microscopic drift only. Thermal
  // derivatives hold the supplied total rate fixed; its nonlocal dependence
  // must be converged by the evolution driver outside the structure solve.
  // The supplied rate enters carried heat linearly through total_rate_enthalpy;
  // the conductivity and kinetic drift contribution are independent of it.
  virtual bool uses_total_species_heat() const {return false;}
  virtual MicroscopicHeatResponse heat_with_total_species_rate(std::size_t face,
      double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
      const Point& hi,const Composition& comp_hi,const SpeciesVector& total_rate,
      bool derivatives) const;
  // An implicit step may need an interior initial iterate even when previous
  // cells are empty. This changes the guess only, never the conserved RHS.
  virtual bool requires_positive_species_guess() const {return false;}
};

MicroscopicFaceResponse microscopic_face(const MicroscopicTransport&,std::size_t face,
    double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,bool derivatives);
MicroscopicHeatResponse microscopic_heat(const MicroscopicTransport&,std::size_t face,
    double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,bool derivatives);
MicroscopicHeatResponse microscopic_heat_with_total_species_rate(const MicroscopicTransport&,
    std::size_t face,double mass_lo,double mass_hi,const Point& lo,const Composition& comp_lo,
    const Point& hi,const Composition& comp_hi,const SpeciesVector& total_rate,bool derivatives);
} // namespace ember
