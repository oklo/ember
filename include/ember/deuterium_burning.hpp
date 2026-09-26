#pragma once
#include "ember/evolution.hpp"
#include "ember/deuterium.hpp"

namespace ember {
// Coupled backward-Euler H1, H2 and He3 evolution, with homogeneous mixing
// within each physically supplied convective region and local burning elsewhere.
std::vector<Composition> burn_deuterium_and_mix(const Model& thermal,const Model& previous,
    const PPDeuterium&,const MixingRegions&,double dt,double tolerance);
// The same reaction system with finite composition diffusion between regions.
// Diffusivities have units cm^2/s and apply equally to every H/He isotope.
// Metals must be spatially uniform. This is macroscopic mixing, not an
// isotope-specific microscopic settling prescription.
std::vector<Composition> burn_deuterium_and_transport(const Model& thermal,const Model& previous,
    const PPDeuterium&,const MixingRegions&,const std::vector<double>& diffusivity,
    double dt,double tolerance);
}
