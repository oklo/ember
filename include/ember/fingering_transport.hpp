#pragma once
#include "ember/metal_microscopic_transport.hpp"
#include "ember/opacity.hpp"

namespace ember {
// Brown saturation extended to two composition fields, with the same sensible heat and
// composition mixing in every face evaluation. Native microscopic species
// fluxes remain separate; add_mixing_flux transports each actual CN isotope.
// The coefficient model currently supports fully ionized, degenerate liquid
// H/He with negligible metals. It is not an atmosphere or solid mixing law.
class BrownFingeringTransport final:public MetalMicroscopicTransport {
 public:
  BrownFingeringTransport(const VariableMetalHelmholtzEos&,const Opacity& radiation,
      const MetalMicroscopicTransport& base,const ScreenedCollisionTransport&);
  struct Face {
    double heat_luminosity{},mass_conductance{},diffusivity{};
    MetalSpeciesMatrix mixing_conductance{};
    double density_ratio{},prandtl{},diffusivity_ratio{},thermal_nusselt_excess{};
  };
  Face face(std::size_t,double,double,const Point&,const Composition&,
            const Point&,const Composition&) const;
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  MetalSpeciesFaceResponse species(std::size_t i,double ml,double mh,const Point& lo,
      const Composition& a,const Point& hi,const Composition& b,bool derivatives) const override {
    return base_.species(i,ml,mh,lo,a,hi,b,derivatives);
  }
  MicroscopicHeatResponse heat(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,const MetalSpeciesVector&,bool) const override;
  void add_mixing_flux(MetalCNFaceResponse&,std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const override;
  bool requires_positive_species_guess() const override {return base_.requires_positive_species_guess();}
  const char* name() const override {return "microscopic transport with Brown fingering heat and mixing";}
 private:
  void add_heat(MicroscopicHeatResponse&,std::size_t,double,double,const Point&,
      const Composition&,const Point&,const Composition&,bool) const;
  const VariableMetalHelmholtzEos& eos_;
  const Opacity& radiation_;
  const MetalMicroscopicTransport& base_;
  const ScreenedCollisionTransport& collisions_;
};
} // namespace ember
