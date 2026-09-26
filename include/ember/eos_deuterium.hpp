#pragma once
#include "ember/eos.hpp"

namespace ember {

// Trace-D material approximation, selected explicitly. Replace D nuclei by H
// nuclei at the same number density in an existing baryonic EOS. The nuclear
// network continues to use the actual proton and D abundances. This preserves
// ion/electron counts, but omits HD/D2 chemistry and isotope energy shifts.
// The translational isotope entropy correction uses unit nuclear partitions.
// No composition or material-domain extrapolation is performed.
// An explicit CN inventory is rescaled with the source mass; any fixed metal
// pattern approximation of the underlying EOS remains in effect.
class DeuteriumApproxEos final : public Eos {
 public:
  explicit DeuteriumApproxEos(const Eos& source) : source_(source) {}
  EosState eval(double,double,const Composition&) const override;
  EosResponse eval_with_derivatives(double,double,const Composition&) const override;
  EosCompositionResponse composition_response(double,double,const Composition&) const override;
  std::optional<DensityRange> density_range(double,const Composition&) const override;
  double rho_from_PT(double,double,const Composition&,double=0.) const override;
  bool has_internal_energy() const override {return source_.has_internal_energy();}
  const char* name() const override {return "trace deuterium: number-mapped EOS, isotope chemistry approximated";}
 private:
  struct Mapping {Composition c;double f{1.},entropy{};};
  static Mapping map(const Composition&);
  static EosState convert(EosState,const Mapping&);
  const Eos& source_;
};
} // namespace ember
