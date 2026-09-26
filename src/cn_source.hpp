#pragma once
#include "ember/nuclear_cn.hpp"
#include "ember/species_flux.hpp"
#include <stdexcept>

namespace ember::detail {
using CNVector = std::array<double,4>;
using CNMatrix = std::array<CNVector,4>;
struct CNReducedSource { CNVector source{}; CNMatrix jacobian{}; };
struct CNFullSource { CNSpeciesVector source{}; CNSpeciesMatrix jacobian{}; };
inline CNFullSource cn_full_source(double T,double rho,const Composition& c,
    const PPChains& pp,const CNNetwork& cn) {
  if(!c.cn_molality)throw std::invalid_argument("CN source: missing catalysts");
  const auto p=pp.composition_response(T,rho,c);
  const auto r=cn.response(T,rho,c,*c.cn_molality);const auto& n=r.physical;
  constexpr std::array<std::size_t,5> species{0,1,3,4,5};
  CNFullSource out;
  for(std::size_t row=0;row<5;++row) {
    const auto s=species[row];out.source[row]=p.state.dXdt[s]+n.state.dXdt[s];
    for(std::size_t col=0;col<2;++col)
      out.jacobian[row][col]=p.d_dXdt_dX[s][col]-p.d_dXdt_dX[s][2]
        +n.d_dXdt_dX[s][col]-n.d_dXdt_dX[s][2];
    for(std::size_t col=0;col<3;++col)
      out.jacobian[row][col+2]=r.d_dXdt_dY[s][col]/mass_numbers[col+3];
  }
  return out;
}
// Independent variables H1, He3, Y12, Y13. He4 closes the lookup mass;
// Y14 closes catalyst number. Used by local burning and spatial transport.
inline CNReducedSource cn_reduced_source(double T,double rho,const Composition& c,
    const CNAbundances& y,const PPChains& pp,const CNNetwork& cn) {
  const auto p=pp.composition_response(T,rho,c);
  const auto r=cn.response(T,rho,c,y);const auto& n=r.physical;
  CNReducedSource out;
  out.source={p.state.dXdt[0]+n.state.dXdt[0],p.state.dXdt[1],
              n.state.dXdt[3]/12,n.state.dXdt[4]/13};
  for(std::size_t col=0;col<2;++col) {
    const CNVector derivative{
      p.d_dXdt_dX[0][col]-p.d_dXdt_dX[0][2]+n.d_dXdt_dX[0][col]-n.d_dXdt_dX[0][2],
      p.d_dXdt_dX[1][col]-p.d_dXdt_dX[1][2],
      (n.d_dXdt_dX[3][col]-n.d_dXdt_dX[3][2])/12,
      (n.d_dXdt_dX[4][col]-n.d_dXdt_dX[4][2])/13};
    for(std::size_t row=0;row<4;++row)out.jacobian[row][col]=derivative[row];
    const CNVector catalyst_derivative{
      r.d_dXdt_dY[0][col]-r.d_dXdt_dY[0][2],0,
      (r.d_dXdt_dY[3][col]-r.d_dXdt_dY[3][2])/12,
      (r.d_dXdt_dY[4][col]-r.d_dXdt_dY[4][2])/13};
    for(std::size_t row=0;row<4;++row)out.jacobian[row][col+2]=catalyst_derivative[row];
  }
  return out;
}
} // namespace ember::detail
