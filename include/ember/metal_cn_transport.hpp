#pragma once
#include "ember/cn_transport.hpp"

namespace ember {
using MetalSpeciesVector=std::array<double,3>; // H1, He3, total metal mass
using MetalSpeciesMatrix=std::array<MetalSpeciesVector,3>;
struct MetalSpeciesFaceResponse {
  MetalSpeciesVector rate{};
  MetalSpeciesMatrix dleft{},dright{};
};
// Physical independent masses. He4 is one minus their sum. Separating inert
// metal mass from CN makes its conservation independent of nuclear captures.
inline constexpr std::size_t METAL_CN_SIZE=7;
inline constexpr std::size_t METAL_CN_D=6;
using MetalCNVector=std::array<double,METAL_CN_SIZE>; // H1, He3, C12, C13, N14, inert metals, initial D
using MetalCNMatrix=std::array<MetalCNVector,METAL_CN_SIZE>;
struct MetalCNFaceResponse {MetalCNVector rate{};MetalCNMatrix dleft{},dright{};};
using MetalCNFlux=std::function<MetalCNFaceResponse(std::size_t,const Composition&,const Composition&,bool)>;
struct MetalCNBoundaryFlux {std::size_t face{};MetalCNVector rate{};};
struct MetalCNTransportResult {
  std::vector<Composition> composition;
  std::vector<MetalCNBoundaryFlux> boundary_fluxes;
  MetalCNVector integrated_balance{};
  double residual{},abundance_correction{};
  std::size_t iterations{};
  std::vector<double> residual_history;
};
MetalCNVector metal_cn_abundances(const Composition&);
Composition metal_cn_composition(const Composition& material_pattern,const MetalCNVector&);
// The metal group velocity advects each actual CN isotope and the inert
// fraction. Individual metal concentration gradients within the group are
// not independent forces in this explicitly selected approximation.
MetalCNFaceResponse common_metal_cn_flux(const MetalSpeciesFaceResponse&,
    const Composition&,const Composition&,bool derivatives);
MetalCNTransportResult burn_metal_cn_and_diffuse(const Model&,const Model&,const PPCNNetwork&,
    const MixingRegions&,const MetalCNFlux&,double dt,const SpeciesTransportOptions& = {},
    std::span<const double> common_mixing_rates = {});
// The same physical inventories and implicit solver with macroscopic mixing
// alone. The common diffusivity transports initial D as well as H/He/metals.
// Microscopic D transport still requires an isotope-aware face provider.
std::vector<Composition> burn_metal_cn_and_transport(const Model&,const Model&,
    const PPCNNetwork&,const MixingRegions&,std::span<const double> diffusivity,
    double dt,double tolerance=1e-12);

struct MetalFluxReconstruction {
  std::vector<MetalSpeciesVector> face_rates,cell_balances,region_balances;
  MetalSpeciesVector integrated_balance{};
};
MetalFluxReconstruction reconstruct_metal_fluxes(const Model&,const Model&,const PPCNNetwork&,
    const MixingRegions&,std::span<const MetalCNBoundaryFlux>,double dt);
} // namespace ember
