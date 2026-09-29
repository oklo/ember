#pragma once
// Hydrostatic outer layer between the tau=100 atmosphere and enclosed mass M-dM.
// The boundary arguments refer to the mesh radius r_b: Teff_b^4=L/(4 pi sigma r_b^2)
// and g_b=G M/r_b^2. photosphere() supplies the physical surface radius and Teff.
//
// Composition and luminosity are constant through the layer; mixing-length
// transport uses the supplied opacity. Thermodynamics can use the interior EOS
// or an explicit source table. The layer's mass, composition and thermal energy
// enter the evolution through the outer cell's reservoir; its internal energy
// distribution is approximated by the base state. Mass above tau=100 is neglected.
#include "ember/atmosphere.hpp"
#include "ember/opacity.hpp"
#include <array>
#include <atomic>
#include <mutex>
#include <string>
#include <vector>

namespace ember {

enum class EnvelopeMetals { reject, neutral, ionized };

class EnvelopeSource {
public:
  explicit EnvelopeSource(const std::string& path);
  struct State { double lnP, grad_ad, cp, delta, chiRho; };
  // Bicubic (Catmull-Rom) in (ln T, ln rho) per composition plane, bilinear between planes in (X, Y3).
  State eval(double lnT, double lnrho, double X, double Y3) const;
  double lnrho_from(double lnT, double lnP, double X, double Y3, double guess) const;
  struct PressureState { double rho, cp, delta, grad_ad, chiRho; };
  PressureState at_pressure(double lnT,double lnP,const Composition&,double& guess,EnvelopeMetals) const;
  std::array<double,2> X_interval() const { return {X_.front(), X_.back()}; }
  std::array<double,2> Y3_interval() const { return {Y3_.front(), Y3_.back()}; }
private:
  std::vector<double> X_, Y3_, lt_, lr_;
  std::vector<std::array<std::vector<double>,5>> planes_;   // [plane][quantity][iT*nr + ir]
};

struct EnvelopeSurface {
  double R{}, Teff{}, L{}, T100{}, P100{};        // true photospheric radius/temperature and luminosity
  double T_base{}, P_base{}, rho_base{}, r_base{};
  std::size_t steps{};
};

class EnvelopeAtmosphere final : public Atmosphere {
public:
  EnvelopeAtmosphere(const Atmosphere& top, const Opacity& opacity, const EnvelopeSource& source,
                     double alpha_mlt, double total_mass, double envelope_mass, double steps_per_unit_lnP = 40,
                     EnvelopeMetals metals = EnvelopeMetals::reject);
  EnvelopeAtmosphere(const Atmosphere& top, const Opacity& opacity, const Eos& eos,
                     double alpha_mlt, double total_mass, double envelope_mass, double steps_per_unit_lnP = 40);
  AtmosphereState eval(double Teff_b, double g_b, const Composition&) const override;
  AtmosphereState eval_value(double Teff_b, double g_b, const Composition&) const override;
  void evaluation_threads(std::size_t threads);
  // Modified Newton: reuse nearby boundary derivatives, never boundary values.
  // Refresh after at most seven uses or a composition/structure displacement.
  void jacobian_reuse(double radius);
  std::size_t jacobian_reused() const { return jacobian_reused_; }
  std::size_t jacobian_computed() const { return jacobian_computed_; }
  const char* name() const override { return "integrated outer envelope"; }
  // True surface for diagnostics (same solve as eval, no derivatives).
  EnvelopeSurface photosphere(double Teff_b, double g_b, const Composition&) const;
  // Initialize from physical surface quantities; the structure mesh ends at r_base.
  EnvelopeSurface from_photosphere(double Teff, double R, const Composition& c) const {
    auto result=integrate(Teff,R,c);
    ratio_.store(R/result.r_base);
    return result;
  }
  double envelope_mass() const { return dM_; }
  // Applies to the H/He source-table constructor, not the interior-EOS constructor.
  static constexpr double max_Z = 1e-6;
  // AtmosphereState::tau for this boundary: the match is at a mass depth, not an optical depth. This atmosphere is
  // the outermost selection and must not be wrapped by composition/gravity interval blends.
  static constexpr double tau_sentinel = -1;
private:
  EnvelopeSurface integrate(double Teff, double R, const Composition&) const;   // from tau100 down to M - dM
  EnvelopeSurface solve(double L, double r_b, const Composition&) const;         // root in R
  EnvelopeSource::PressureState at_pressure(double lnT,double lnP,const Composition&,double& guess) const;
  const Atmosphere& top_; const Opacity& opacity_;
  const EnvelopeSource* source_{};
  EnvelopeMetals metals_{EnvelopeMetals::reject};
  const Eos* eos_{};
  double alpha_, M_, dM_, per_unit_;
  std::size_t threads_{1};
  double jacobian_radius_{};
  struct JacobianCache {
    bool valid{};
    double logTeff{}, logg{};
    Composition composition{};
    std::array<double,4> derivatives{};
    unsigned uses{};
  };
  mutable JacobianCache jacobian_;
  mutable std::mutex jacobian_mutex_;
  mutable std::atomic<std::size_t> jacobian_reused_{0}, jacobian_computed_{0};
  mutable std::atomic<double> ratio_{0};   // last R/r_b, used as a radius-root starting guess
  mutable std::atomic<std::size_t> integrations_{0};
public:
  std::size_t integrations() const { return integrations_; }
};

} // namespace ember

namespace ember {

// Tabulated EnvelopeAtmosphere. Nodes on a uniform (ln Teff_b, ln g_b) grid for each explicit (X, Y3) plane of
// the layer source; each node stores ln T_b, ln P_b, ln R (true photosphere) and their first derivatives from the
// native reference. Bicubic Hermite in (ln Teff_b, ln g_b) (cross derivatives by centred differences of the stored
// slopes), bilinear in (X, Y3). Outside the grid or the composition interval it throws.
class EnvelopeMapAtmosphere final : public Atmosphere {
public:
  explicit EnvelopeMapAtmosphere(const std::string& path);
  AtmosphereState eval(double Teff_b, double g_b, const Composition&) const override;
  const char* name() const override { return "tabulated outer envelope"; }
  // True photospheric radius and Teff for the same base state.
  std::array<double,2> photosphere_R_Teff(double Teff_b, double g_b, const Composition&) const;
  double envelope_mass() const { return dM_; }
  double total_mass() const { return M_; }
  struct Values { double lnT, lnP, lnR, dT[2], dP[2], dR[2]; };
  Values values(double Teff_b, double g_b, const Composition&) const;
private:
  double M_{}, dM_{};
  std::vector<double> X_, Y3_, lte_, lg_;
  // [plane][q][it*ng+ig], q: 0..2 value (lnT, lnP, lnR), 3..8 slopes (dT/dlte, dT/dlg, dP/dlte, dP/dlg, dR/dlte, dR/dlg)
  std::vector<std::array<std::vector<double>,9>> planes_;
  std::array<double,9> plane_hermite(std::size_t plane, double x, double y) const;
};

}  // namespace ember
