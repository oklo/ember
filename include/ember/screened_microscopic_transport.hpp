#pragma once
#include "ember/collision_transport.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include "ember/microscopic_transport.hpp"

namespace ember {
// Conditional fully stripped transport, using the native material potential.
// The caller supplies a temperature lower bound and a screening choice. The
// bound restricts evaluation; it does not establish that the ions are stripped.
// Selected independent species must be positive at both evaluated endpoints.
// An implicit step may start from exactly zero OLD fractions using an interior
// guess; fluxes are evaluated only at the new-time candidate. Globally absent
// species can be removed explicitly, and must then stay exactly absent.
class ScreenedMicroscopicTransport final : public MicroscopicTransport {
 public:
  ScreenedMicroscopicTransport(const SmoothMetalHelmholtzEos& eos,
      ScreenedCollisionTransport collisions,bool include_ion_screening,
      double minimum_temperature,std::array<bool,2> active_species={true,true},
      bool radiation_with_redistribution=false);
  // The optional parcel term is h_rad*(J_total-J_micro). The material force
  // and microscopic energy convention remain unchanged. It assumes LTE and
  // pressure-balanced fluid parcels for the additional redistribution.
  MicroscopicFaceResponse eval(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  MicroscopicHeatResponse heat(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  bool uses_total_species_heat() const override {return true;}
  MicroscopicHeatResponse heat_with_total_species_rate(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,const SpeciesVector&,bool) const override;
  const char* name() const override {return "conditional fully stripped H/He species and heat transport";}
  bool requires_positive_species_guess() const override {return true;}
 private:
  const SmoothMetalHelmholtzEos& eos_;
  ScreenedCollisionTransport collisions_;
  bool include_ions_;
  double minimum_temperature_;
  std::array<bool,2> active_;
  bool radiation_with_redistribution_;
};
} // namespace ember
