#pragma once
#include "ember/eos.hpp"
#include "ember/opacity.hpp"
#include <cstddef>

namespace ember {

struct AtmosphereState {
  double T{}, P{}, Pgas{}, rho{}, tau{};  // CGS; P includes a_rad*T^4/3
  // Logarithmic partials at fixed composition and matching optical depth.
  double dlnT_dlnTeff{}, dlnT_dlng{};
  double dlnP_dlnTeff{}, dlnP_dlng{};
};

// Pressure and temperature at a specified optical depth. The caller chooses
// an atmosphere explicitly: table errors must not trigger a silent fallback.
class Atmosphere {
public:
  virtual ~Atmosphere() = default;
  virtual AtmosphereState eval(double Teff, double gravity, const Composition&) const = 0;
  virtual const char* name() const = 0;
};

struct GreyAtmosphereOptions {
  double tau_match{2.0 / 3.0};
  double tau_top{1e-6};
  double tolerance{1e-8};           // integration error in ln Pgas and sensitivities
  std::size_t max_steps{10000};     // attempted steps, including rejected steps
};

// Plane-parallel, radiative Eddington atmosphere, with local EOS and radiative opacity.
// T^4 = (3/4) Teff^4 (tau+2/3), dP/dtau = g/kappa.
// Integrates positive gas pressure, retaining radiation's pressure gradient.
// The finite top uses Pgas=(g/kappa-F/c)*tau_top with self-consistent kappa;
// its truncation error is separate from the integration tolerance. Check
// convergence in tau_top when using a new opacity regime. No convective
// correction is applied to the T(tau) law; use a model atmosphere for that.
class GreyAtmosphere final : public Atmosphere {
public:
  GreyAtmosphere(const Eos& eos, const Opacity& opacity, GreyAtmosphereOptions options = {});
  AtmosphereState eval(double Teff, double gravity, const Composition&) const override;
  const char* name() const override { return "Eddington grey (varying opacity)"; }
private:
  const Eos& eos_;
  const Opacity& opacity_;
  GreyAtmosphereOptions options_;
};

} // namespace ember
