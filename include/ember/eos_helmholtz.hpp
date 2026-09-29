#pragma once
#include "ember/eos.hpp"
#include <array>
#include <filesystem>
#include <iosfwd>
#include <span>
#include <string>
#include <vector>

namespace ember {

// Natural-log T/rho derivatives of material F/T; entries i+j<=3 are used.
using HelmholtzJet = std::array<std::array<double,4>,4>;
EosResponse helmholtz_response(double T, double rho, HelmholtzJet material);

// Fixed-composition material Helmholtz potential, phi=F/T, in ln(T) and
// ln[rho/(T/1e6)^1.5]. Biquintic Hermite interpolation is C2 across cells;
// all P/E/S and transport responses differentiate that SAME potential.
// Radiation is added analytically, once. No extrapolation or fallback.
class HelmholtzTableEos final : public Eos {
public:
  enum class Mixture { exact, allow_documented_proxy };
  explicit HelmholtzTableEos(const std::filesystem::path&, Mixture = Mixture::exact);
  EosState eval(double T, double rho, const Composition&) const override;
  EosResponse eval_with_derivatives(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  std::optional<DensityRange> density_range_near(double,const Composition&,double) const override;
  double rho_from_PT(double T, double P, const Composition&, double guess = 0) const override;
  const char* name() const override { return source_.c_str(); }
  const Composition& composition() const { return composition_; }
  // These evaluate the table's declared material, before adding radiation.
  HelmholtzJet material_jet(double T, double rho) const;
  DensityRange material_density_range(double T) const;
  DensityRange material_density_range_near(double T,double rho) const;
private:
  friend class MetalHelmholtzEos;
  friend class SmoothMetalHelmholtzEos;
  friend class VariableMetalHelmholtzEos;
  HelmholtzTableEos(std::istream&, Mixture);
  void write_binary(std::ostream&) const;
  void initialize_support();
  struct WeightedTable {const HelmholtzTableEos* table;double weight;};
  bool same_material_grid(const HelmholtzTableEos& other) const {return t_==other.t_ && q_==other.q_;}
  // Caller verifies the common grid once when constructing its family.
  static HelmholtzJet mixed_material_jet(double T,double rho,std::span<const WeightedTable>);
  struct WeightedCompositionTable {
    const HelmholtzTableEos* table;
    // Value, three first derivatives, six symmetric second derivatives.
    std::array<double,10> weight{};
  };
  static std::array<HelmholtzJet,10> mixed_composition_jets(double,double,
      std::span<const WeightedCompositionTable>,std::size_t channels);
  struct Node { bool valid{}; std::array<double,9> d{}; };
  void check_composition(const Composition&) const;
  std::pair<std::size_t,std::size_t> supported_q(std::size_t it) const;
  std::pair<std::size_t,std::size_t> supported_q(std::size_t it,std::size_t iq) const;
  bool supported_cell(std::size_t it,std::size_t iq) const {
    const auto n=q_.size();
    return nodes_[it*n+iq].valid && nodes_[it*n+iq+1].valid &&
      nodes_[(it+1)*n+iq].valid && nodes_[(it+1)*n+iq+1].valid;
  }
  std::string source_, proxy_;
  Composition composition_;
  std::vector<double> t_, q_; // natural logarithms; q=rho/(T/1e6)^1.5
  std::vector<Node> nodes_;
  std::vector<std::vector<std::pair<std::size_t,std::size_t>>> support_; // connected valid node intervals
};
} // namespace ember
