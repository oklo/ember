#pragma once
#include "ember/collision_transport.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/metal_cn_transport.hpp"
#include "ember/microscopic_transport.hpp"

namespace ember {
struct MetalMicroscopicFaceResponse : MicroscopicHeatResponse {
  MetalSpeciesFaceResponse species;
};
// An explicit three-mass capability. Calling the two-mass interface rejects
// rather than silently dropping the metal flux or its transported enthalpy.
class MetalMicroscopicTransport : public MicroscopicTransport {
 public:
  virtual MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,
      const Point&,const Composition&,const Point&,const Composition&,bool) const=0;
  virtual MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,
      const Point&,const Composition&,const Point&,const Composition&,
      const MetalSpeciesVector&,bool) const=0;
  MicroscopicFaceResponse eval(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const final;
  MicroscopicHeatResponse heat(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  bool uses_total_species_heat() const final {return true;}
};

MetalMicroscopicFaceResponse metal_microscopic_face(const MetalMicroscopicTransport&,
    std::size_t,double,double,const Point&,const Composition&,const Point&,const Composition&,bool);
MicroscopicHeatResponse microscopic_heat_with_total_metal_rate(const MetalMicroscopicTransport&,
    std::size_t,double,double,const Point&,const Composition&,const Point&,const Composition&,
    const MetalSpeciesVector&,bool);

// Conditional fully stripped transport. All GS98 metals share one mass
// velocity, while their collision charges and heat variables remain distinct.
// The temperature bound limits evaluation; it does not prove full ionization.
class ScreenedMetalMicroscopicTransport final : public MetalMicroscopicTransport {
 public:
  ScreenedMetalMicroscopicTransport(const VariableMetalHelmholtzEos&,
      ScreenedCollisionTransport,bool include_ion_screening,double minimum_temperature,
      std::array<bool,3> active_species={true,true,true},bool radiation_with_redistribution=false);
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  MicroscopicHeatResponse heat(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,const MetalSpeciesVector&,bool) const override;
  const char* name() const override {return "conditional fully stripped H/He/metal species and heat transport";}
  bool requires_positive_species_guess() const override {return true;}
 private:
  const VariableMetalHelmholtzEos& eos_;
  ScreenedCollisionTransport collisions_;
  bool include_ions_;
  double minimum_temperature_;
  std::array<bool,3> active_;
  bool radiation_with_redistribution_;
};
} // namespace ember
