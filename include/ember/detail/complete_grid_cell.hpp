#pragma once
#include "ember/interp.hpp"
#include <algorithm>
#include <array>
#include <optional>
#include <vector>

namespace ember::detail {
// Choose a complete tensor cell. At an exact knot a complete cell on either
// side supplies a one-sided derivative; points inside a missing cell fail.
// The caller validates axis bounds and mask sizes before using this helper.
template<std::size_t N>
std::optional<std::array<std::size_t,N>> complete_grid_cell(
    const std::array<std::vector<double>,N>& axes,
    const std::array<double,N>& q,const std::vector<bool>& vertices,
    const std::vector<bool>& cells={}) {
  std::array<std::size_t,N> preferred{};
  for(std::size_t k=0;k<N;++k)
    preferred[k]=axes[k].size()==1?0:interp::locate(axes[k],q[k]);
  if(vertices.empty() && cells.empty())return preferred;
  for(unsigned alternative=0;alternative<(1U<<N);++alternative) {
    auto base=preferred;bool eligible=true;
    for(std::size_t k=0;k<N;++k) {
      if(!(alternative&(1U<<k)))continue;
      if(!base[k] || q[k]!=axes[k][base[k]]) {eligible=false;break;}
      --base[k];
    }
    if(!eligible)continue;
    if(!cells.empty()) {
      std::size_t index=0;
      for(std::size_t k=0;k<N;++k)
        index=index*std::max(std::size_t{1},axes[k].size()-1)+base[k];
      if(!cells[index])continue;
    }
    bool complete=true;
    for(unsigned corner=0;corner<(1U<<N) && complete;++corner) {
      std::size_t index=0;
      for(std::size_t k=0;k<N;++k)
        index=index*axes[k].size()+base[k]+(axes[k].size()>1 && (corner&(1U<<k)));
      complete=vertices.empty() || vertices[index];
    }
    if(complete)return base;
  }
  return std::nullopt;
}
} // namespace ember::detail
