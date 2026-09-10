#pragma once
#include "ember/conduction.hpp"
#include <filesystem>
#include <string>
#include <vector>

namespace ember {
// Electron-fraction weighted resistivities at common electron density.
// This approximates the mixture collision sum and counts electron-electron
// scattering once. Pure-ion table correlations and isotope masses are a
// documented approximation, particularly at partial degeneracy.
class TabulatedConduction final : public Conduction {
public:
  explicit TabulatedConduction(const std::filesystem::path&);
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<Opacity::DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override { return label_.c_str(); }

private:
  struct Conductivity {
    double logK, dT, drho;
  };
  Conductivity pure(std::size_t, double logT, double logElectronMassDensity) const;
  Conductivity ion(double charge, double logT, double logElectronMassDensity) const;
  std::vector<double> logT_, logRho_, logZ_, mass_, logK_;
  std::string label_;
};

// Full-ionization conduction is suppressed in the neutral envelope. Smooth
// activation in log T between specified temperatures, with the derivative
// included. This is an explicit domain join; test its negligible influence.
class HotConduction final : public Conduction {
public:
  HotConduction(const Conduction&, double start = 3e5, double end = 1e6);
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<Opacity::DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override { return "ionized conduction with explicit envelope temperature join"; }

private:
  const Conduction& source_;
  double start_, end_;
};
} // namespace ember
