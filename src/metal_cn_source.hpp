#pragma once
#include "ember/metal_cn_transport.hpp"

namespace ember::detail {
struct MetalCNSource {MetalCNVector source{};MetalCNMatrix jacobian{};};
inline MetalCNSource metal_cn_source(double T,double rho,const Composition& c,
    const Nuclear& light,const CNNetwork& cn) {
  if(!c.cn_molality || c.cn_mass_convention!=CNMassConvention::explicit_metal_mass)
    throw std::invalid_argument("metal CN source: explicit physical inventory required");
  const auto p=light.composition_response(T,rho,c);
  const auto r=cn.response(T,rho,c,*c.cn_molality);const auto& n=r.physical;
  constexpr std::array<std::size_t,METAL_CN_SIZE> species{0,1,3,4,5,0,
      static_cast<std::size_t>(Species::H2)};
  MetalCNSource out;const double Z=c.Z();
  for(std::size_t row=0;row<METAL_CN_SIZE;++row) {
    if(row==5)continue; // inert metals
    const auto s=species[row];out.source[row]=p.state.dXdt[s]+n.state.dXdt[s];
    for(const auto col:{0UL,1UL,METAL_CN_D}) {
      const auto j=species[col];
      out.jacobian[row][col]=p.d_dXdt_dX[s][j]-p.d_dXdt_dX[s][2]
          +n.d_dXdt_dX[s][j]-n.d_dXdt_dX[s][2];
    }
    double metal_partial=-p.d_dXdt_dX[s][2]-n.d_dXdt_dX[s][2];
    if(Z>0)for(std::size_t k=3;k<METAL_END;++k)
      metal_partial+=c.X[k]/Z*(p.d_dXdt_dX[s][k]+n.d_dXdt_dX[s][k]);
    for(std::size_t col=2;col<6;++col) {
      out.jacobian[row][col]=metal_partial;
      if(col<5)out.jacobian[row][col]+=r.d_dXdt_dY[s][col-2]/mass_numbers[col+1];
    }
  }
  // The inert metal mass has no nuclear source.
  return out;
}
} // namespace ember::detail
