#pragma once
#include "ember/metal_microscopic_transport.hpp"
#include <span>

namespace ember {
// Conditional hot-region diffusion with a separately assessed convective
// envelope heat approximation. Species transport is NEVER temperature-tapered:
// a required exchange face must lie above hot_temperature at both nodes.
// The cooler heat provider may be used only inside a mixed convective region;
// the evolution solver's region-boundary flux queries enforce that condition.
// The transition of conductivity/residual heat must be tested on stellar
// structures. It does not establish a partially ionized collision law.
class ConvectiveEnvelopeTransport final:public MetalMicroscopicTransport {
 public:
  ConvectiveEnvelopeTransport(const MetalMicroscopicTransport& hot,
      const MetalMicroscopicTransport& convective,double cool_temperature,
      double hot_temperature);
  std::span<const MetalSpeciesVector> diagnostic_rates;
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,
      const Point&,const Composition&,const Point&,const Composition&,bool)const override;
  MicroscopicHeatResponse heat(std::size_t,double,double,
      const Point&,const Composition&,const Point&,const Composition&,bool)const override;
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,
      const Point&,const Composition&,const Point&,const Composition&,
      const MetalSpeciesVector&,bool)const override;
  bool requires_positive_species_guess()const override{return hot_.requires_positive_species_guess();}
  const char* name()const override{return "hot microscopic exchange with assessed convective envelope heat";}
 private:
  const MetalMicroscopicTransport &hot_,&convective_;
  double log_cool_,log_hot_;
};
} // namespace ember
