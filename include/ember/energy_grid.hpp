#pragma once
#include "ember/model.hpp"
#include <stdexcept>

namespace ember {

inline bool face_luminosities(const Model& m) {
  switch(m.luminosity_grid) {
    case LuminosityGrid::mass_nodes: return false;
    case LuminosityGrid::volume_faces: return true;
  }
  throw std::invalid_argument("unknown luminosity grid");
}

// The same nodal control volumes used for conservative abundance changes.
inline double nodal_volume_mass(const Model& m,std::size_t i) {
  if(m.size()<2 || m.m.size()!=m.size() || i>=m.size())
    throw std::invalid_argument("nodal volume: invalid mesh or index");
  const double inner=i ? .5*(m.m[i]-m.m[i-1]) : m.m[0];
  return inner+(i+1<m.size() ? .5*(m.m[i+1]-m.m[i]) : 0.);
}

inline double luminosity_mass(const Model& m,std::size_t i) {
  if(i>=m.size() || m.m.size()!=m.size())
    throw std::invalid_argument("luminosity position: invalid mesh or index");
  return face_luminosities(m) && i+1<m.size() ? .5*(m.m[i]+m.m[i+1]) : m.m[i];
}

// Energy row i relates luminosities i and i+1. Its mass differs from the
// interval used by hydrostatic/temperature rows when luminosities are faces.
inline double energy_interval_mass(const Model& m,std::size_t i) {
  if(i+1>=m.size() || m.m.size()!=m.size())
    throw std::invalid_argument("energy interval: invalid mesh or index");
  return face_luminosities(m) ? nodal_volume_mass(m,i+1) : m.m[i+1]-m.m[i];
}

inline double thermal_face_luminosity(const Model& m,std::size_t i) {
  if(i+1>=m.size())throw std::invalid_argument("thermal face: invalid index");
  return face_luminosities(m) ? m.y[i].L : .5*(m.y[i].L+m.y[i+1].L);
}
} // namespace ember
