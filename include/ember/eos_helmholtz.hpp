#pragma once
#include "ember/eos.hpp"
#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace ember {

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
  double rho_from_PT(double T, double P, const Composition&, double guess = 0) const override;
  const char* name() const override { return source_.c_str(); }
private:
  struct Node { bool valid{}; std::array<double,9> d{}; };
  void check_composition(const Composition&) const;
  std::pair<std::size_t,std::size_t> supported_q(std::size_t it) const;
  std::string source_;
  Composition composition_;
  std::vector<double> t_, q_; // natural logarithms; q=rho/(T/1e6)^1.5
  std::vector<Node> nodes_;
};
} // namespace ember
