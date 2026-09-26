#pragma once
#include "ember/atmosphere.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
// Two atmosphere responses share the same matching state at their Z join.
// The lower response must use the upper response as its reference so that
// the existing isotope response and all matching-state derivatives persist.
class PiecewiseMetalAtmosphere final : public Atmosphere {
public:
  PiecewiseMetalAtmosphere(const Atmosphere& upper,const Atmosphere& lower,double join)
      : upper_(upper),lower_(lower),join_(join) {
    if(!std::isfinite(join) || join<=0 || join>=1)
      throw std::invalid_argument("PiecewiseMetalAtmosphere: invalid composition join");
  }
  AtmosphereState eval(double T,double g,const Composition& c) const override {
    // Z sums five stored metal fractions. Treat only summation roundoff
    // at the shared source endpoint as the upper interval. This avoids
    // requiring an unused lower response stencil for that exact source.
    // Keep the physical composition unchanged for all material queries.
    const double roundoff=8*std::numeric_limits<double>::epsilon()*join_;
    return (c.Z()>=join_-roundoff?upper_:lower_).eval(T,g,c);
  }
  const char* name() const override{return "joined measured atmospheric metal responses";}
private:
  const Atmosphere& upper_;const Atmosphere& lower_;double join_;
};
} // namespace ember
