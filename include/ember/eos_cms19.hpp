#pragma once
#include "ember/eos.hpp"
#include <filesystem>
#include <memory>

namespace ember {
// CMS19 additive-volume H/He mixture, interpolating original TP density and
// entropy. Static structure only: the original energy columns have defects
// at model joins. Thermodynamic consistency errors must be audited before
// evolution; has_internal_energy() is false and E is unavailable (NaN).
// No electron abundance is supplied by these tables. No extrapolation.
class Cms19Eos final : public Eos {
public:
  enum class Metals { reject, helium_proxy };
  Cms19Eos(const std::filesystem::path& hydrogen, const std::filesystem::path& helium,
           Metals = Metals::reject);
  ~Cms19Eos() override;
  EosState eval(double T, double rho, const Composition&) const override;
  EosResponse eval_with_derivatives(double T, double rho, const Composition&) const override;
  double rho_from_PT(double T, double P, const Composition&, double guess = 0) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  bool has_internal_energy() const override { return false; }
  const char* name() const override { return "CMS19 H/He entropy EOS (static only)"; }
private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
} // namespace ember
