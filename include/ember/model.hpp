#pragma once
#include "ember/composition.hpp"
#include <cmath>
#include <vector>

namespace ember {

// The four variables carried at each mesh point.
//
// Logarithms of r, rho and T; luminosity linear, because it passes through
// zero at the centre and changes sign wherever the star is contracting.
// Density spans eighteen decades between a giant's photosphere and a white
// dwarf's core, and no linear variable can be conditioned across that; the
// Fortran ancestor used rho^(1/3) as a compromise and paid for it at both
// ends.
enum class Var : std::size_t { lnr = 0, lnrho, lnT, L, COUNT };
inline constexpr std::size_t NVAR = static_cast<std::size_t>(Var::COUNT);

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

// A stellar model on a Lagrangian mass mesh.
//
// The mesh is stored as the mass interior to each point, and is allowed to
// change between steps: refining where a burning shell steepens and coarsening
// where an isothermal core does not need the points is not an optimisation but
// a requirement, and doing it badly ate the core resolution of every giant in
// the Fortran line.
struct Model {
  double M{};                       // total mass, g
  double age{};                     // s
  std::vector<double> m;            // mass interior to point i, g  (m[0] = 0)
  std::vector<Point>  y;            // state at point i
  std::vector<Composition> comp;    // composition at point i
  std::vector<double> Lsurf_hist;   // diagnostics

  std::size_t size() const { return y.size(); }
  double r(std::size_t i)   const { return std::exp(y[i].lnr); }
  double rho(std::size_t i) const { return std::exp(y[i].lnrho); }
  double T(std::size_t i)   const { return std::exp(y[i].lnT); }
};

} // namespace ember
