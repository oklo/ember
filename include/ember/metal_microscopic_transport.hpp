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
  // Optional per-face first-order reuse of collision responses (see
  // CollisionTaylorCache). Null restores exact evaluation at every call.
  void use_collision_taylor(std::shared_ptr<const CollisionTaylorCache> cache) {taylor_=std::move(cache);}
  // Optional first-order reuse of the EOS composition potential (per mesh
  // point) and exchange/radiation enthalpies (per face) within radius r of an
  // exact anchor in ln T and ln rho (absolute) and XH, X3, Z (relative to
  // the anchor's fraction), using the EOS's own first
  // derivatives. Returned derivative arrays are the anchor's. r=0 disables.
  void use_eos_taylor(double radius,bool verify=false);
  // Optional phase scenario: multiply the entire species mobility matrix by
  // (1-w)+remaining*w, where w is the selected EOS solid weight. The same
  // factor applies to microscopic carried heat; zero-flux conduction remains.
  // remaining=0 is an immobile-solid limit, not a calibrated crystal law.
  void use_phase_mobility(const ColdHeliumOptions&,double remaining);
  struct EosReuse {std::size_t hits{},exact{},verified{};double worst_potential{},worst_enthalpy{};};
  EosReuse eos_reuse_statistics() const;
 private:
  const VariableMetalHelmholtzEos& eos_;
  ScreenedCollisionTransport collisions_;
  bool include_ions_;
  double minimum_temperature_;
  std::array<bool,3> active_;
  bool radiation_with_redistribution_;
  std::shared_ptr<const CollisionTaylorCache> taylor_;
 public:
  struct EosCache;
 private:
  std::shared_ptr<EosCache> eos_cache_;
  std::optional<ColdHeliumOptions> phase_options_;
  double solid_mobility_=1;
};
} // namespace ember
