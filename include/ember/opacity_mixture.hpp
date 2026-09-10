#pragma once
#include "ember/opacity_blend.hpp"
#include "ember/opacity_table.hpp"
#include <memory>

namespace ember {
// Source (X,Z) family. Each plane uses the source's native T/density grid;
// ln opacity is interpolated linearly in Z. No metallicity extrapolation.
class MixtureOpacity final : public Opacity {
public:
  explicit MixtureOpacity(const std::filesystem::path& manifest);
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override { return name_.c_str(); }

private:
  std::pair<std::size_t, double> interval(const Composition&) const;
  Composition source(const Composition&, std::size_t) const;
  std::vector<double> z_;
  std::vector<std::unique_ptr<TabulatedOpacity>> tables_;
  std::string name_;
};

// Map H/He isotope number densities into atomic-mixture opacity tables.
// rho_source=s*rho and kappa=s*kappa_source preserve extinction per length.
// This omits isotope-dependent cross sections, CIA reduced-mass corrections
// and collision broadening isotope shifts; source metal patterns stay fixed.
class ElementalOpacity final : public Opacity {
public:
  explicit ElementalOpacity(const Opacity& table) : table_(table) {}
  OpacityState eval(double T, double rho, const Composition&) const override;
  std::optional<DensityRange> density_range(double T, const Composition&) const override;
  const char* name() const override {
    return "number-density mapped elemental opacity; He isotope cross sections approximated";
  }

private:
  struct Mapping {
    Composition c;
    double scale;
    std::array<double, 2> ds, dx, dz;
  };
  static Mapping map(const Composition&);
  const Opacity& table_;
};

// Versioned AESOPUS/TOPS source families with the same physical T joins as
// the historical stellar calculation, followed by the isotope mapping.
class StellarMixtureOpacity final : public Opacity {
public:
  explicit StellarMixtureOpacity(const std::filesystem::path& directory)
      : low_(directory / "aesopus21_gs98_mixture.dat"), hot_low_(directory / "tops_gs98_mixture_low.dat"),
        hot_high_(directory / "tops_gs98_mixture_high.dat"), hot_(hot_low_, hot_high_, 5.6, 5.7),
        radiative_(low_, hot_, 4.4, 4.5), mapped_(radiative_) {}
  OpacityState eval(double T, double rho, const Composition& c) const override {
    return mapped_.eval(T, rho, c);
  }
  std::optional<DensityRange> density_range(double T, const Composition& c) const override {
    return mapped_.density_range(T, c);
  }
  const char* name() const override {
    return "AESOPUS/TOPS X/Z family, number-density isotope mapping; GS98 metal proxy";
  }

private:
  MixtureOpacity low_, hot_low_, hot_high_;
  BlendedOpacity hot_, radiative_;
  ElementalOpacity mapped_;
};
} // namespace ember
