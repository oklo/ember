#pragma once
#include "ember/atmosphere.hpp"
#include <array>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>

namespace ember {

// A measured log-temperature/log-gas-pressure metal response applied to an
// existing atmosphere. Its H/He3 dependence is retained in the reference;
// the added metal response is measured at He3=0 and assumed separable over
// the small He3 interval declared by the response table. GS98 metal ratios
// remain fixed. This approximation requires explicit selection and checks.
class MetalResponseAtmosphere final : public Atmosphere {
public:
  enum class Approximation { separable_gs98_response };
  struct Options {
    // Optional, explicitly bounded linear continuation above reference Z.
    // Intended for a measured small nuclear-capture increase, never clipping.
    double maximum_reference_metal_excess{};
  };
  MetalResponseAtmosphere(const Eos &, const Atmosphere &reference,
                          const std::filesystem::path &, Approximation, Options);
  MetalResponseAtmosphere(const Eos &, const Atmosphere &reference,
                          std::istream &, Approximation, Options);
  AtmosphereState eval(double Teff, double gravity, const Composition &) const override;
  const char *name() const override { return source_.c_str(); }

private:
  void read(std::istream &);
  std::optional<std::array<std::size_t, 3>>
  stencil(const std::array<double, 3> &) const;
  // Natural-log response, and its derivatives in H, ln Teff and ln g.
  std::array<double, 4> interpolate(const std::vector<double> &,
                                  const std::array<double, 3> &) const;
  const Eos &eos_;
  const Atmosphere &reference_;
  Options options_;
  std::string source_, approximation_;
  double reference_z_{}, source_z_{}, maximum_helium3_{}, tau_{};
  std::array<std::vector<double>, 3> axes_;
  std::vector<double> delta_logT_, delta_logPg_;
  std::vector<bool> valid_;
};

} // namespace ember
