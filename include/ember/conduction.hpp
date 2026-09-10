#pragma once
#include "ember/composition.hpp"
#include "ember/opacity.hpp"
#include <memory>

namespace ember {

// Electron conduction, expressed as an opacity so it combines with the
// radiative one by the usual reciprocal sum, 1/k = 1/k_rad + 1/k_cond.
//
// In a white dwarf interior conduction is not a correction: degenerate
// electrons carry essentially all of the heat, and the star's cooling rate is
// set by where the conductive interior meets the radiative envelope.  The
// modern source is Cassisi et al. (2007), superseding Hubbard & Lampe (1969);
// Blouin et al. (2020) revise it further in the regime a cold massive remnant
// occupies.
class Conduction {
public:
  virtual ~Conduction() = default;
  virtual OpacityState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
  virtual std::optional<Opacity::DensityRange> density_range(double, const Composition&) const {
    return std::nullopt;
  }
};

// Radiative and conductive opacities combined.  Owns neither; both are shared
// so that a run can report which sources produced a given track.
class CombinedOpacity final : public Opacity {
public:
  CombinedOpacity(std::shared_ptr<Opacity> rad, std::shared_ptr<Conduction> cond)
      : rad_(std::move(rad)), cond_(std::move(cond)) {
    if(!rad_)throw std::invalid_argument("CombinedOpacity: radiative opacity is required");
  }
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T,const Composition&) const override;
  const char* name() const override { return "radiative + conductive"; }
private:
  std::shared_ptr<Opacity> rad_;
  std::shared_ptr<Conduction> cond_;
};

} // namespace ember
