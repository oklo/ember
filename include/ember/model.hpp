#pragma once
#include "ember/composition.hpp"
#include <cmath>
#include <vector>

namespace ember {

// Logarithmic radius, density and temperature, plus signed luminosity.
// Luminosity remains linear because contraction and expansion can reverse
// its sign within the star.
enum class Var : std::size_t { lnr = 0, lnrho, lnT, L, COUNT };
inline constexpr std::size_t NVAR = static_cast<std::size_t>(Var::COUNT);

// Temperature, density, radius and composition are stored at mass nodes.
// Total luminosity can instead be stored at each nodal volume's outer face,
// where the conservative composition update also carries material heat.
enum class LuminosityGrid { mass_nodes, volume_faces };

struct Point {
  double lnr{}, lnrho{}, lnT{}, L{};
  double& operator[](Var v) {
    switch (v) { case Var::lnr: return lnr; case Var::lnrho: return lnrho;
                 case Var::lnT: return lnT; default: return L; }
  }
  double operator[](Var v) const {
    switch (v) { case Var::lnr: return lnr; case Var::lnrho: return lnrho;
                 case Var::lnT: return lnT; default: return L; }
  }
};

// A stellar model on a Lagrangian mass mesh. Each accepted step retains
// the same mesh; explicit remapping must conserve its material inventories.
struct Model {
  double M{};                       // total mass, g
  double age{};                     // s
  std::vector<double> m;            // enclosed mass, g; ends at M - envelope_mass
  std::vector<Point>  y;            // state at point i
  std::vector<Composition> comp;    // composition at point i
  std::vector<double> Lsurf_hist;   // diagnostics
  LuminosityGrid luminosity_grid{LuminosityGrid::mass_nodes};

  // A homogeneous outer reservoir, represented thermally at the last node.
  // Its mass participates in composition, mixing and energy conservation.
  // A deep atmosphere supplies its pressure/temperature profile and true radius.
  // The base-state thermal approximation must be bounded for each application.
  double envelope_mass{};
  bool valid_outer_mass() const {
    return !m.empty() && std::isfinite(M) && M>0 && std::isfinite(envelope_mass)
        && envelope_mass>=0 && envelope_mass<M && m.back()==M-envelope_mass
        && (envelope_mass==0 || luminosity_grid==LuminosityGrid::volume_faces);
  }

  std::size_t size() const { return y.size(); }
  double r(std::size_t i)   const { return std::exp(y[i].lnr); }
  double rho(std::size_t i) const { return std::exp(y[i].lnrho); }
  double T(std::size_t i)   const { return std::exp(y[i].lnT); }
};

} // namespace ember
