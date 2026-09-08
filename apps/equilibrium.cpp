#include "ember/relaxation.hpp"
#include "ember/eos_composite.hpp"
#include "ember/eos_cms19.hpp"
#include "ember/eos_helmholtz.hpp"
#include "ember/atmosphere_table.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_opal.hpp"
#include "ember/opacity_tops.hpp"
#include "ember/opacity_blend.hpp"
#include "../examples/stellar_seed.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <string_view>

namespace {
// stdout is always JSON, including invalid inputs and unconverged models.
void json_string(std::string_view value) {
  std::putchar('"');
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') { std::putchar('\\'); std::putchar(c); }
    else if (c < 0x20) std::printf("\\u%04x", static_cast<unsigned int>(c));
    else std::putchar(c);
  }
  std::putchar('"');
}
void number(double value) {
  if (std::isfinite(value)) std::printf("%.16g", value);
  else std::printf("null");
}
template<class T> T parse(const char* text) {
  T value{};
  const char* end = text + std::strlen(text);
  const auto result = std::from_chars(text, end, value);
  if (result.ec != std::errc{} || result.ptr != end)
    throw std::invalid_argument("arguments must be complete numeric values");
  return value;
}
}

int main(int argc, char** argv) {
  using namespace ember;
  if (argc == 2 && std::string_view(argv[1]) == "--help") {
    std::puts("usage: ember-equilibrium [mesh_points [mass_solar [seed_radius_solar]]] [--eos ionized|cms19|freeeos]\n"
              "       [--hot-opacity opal|tops] [--tau-top value] [--seed-index 3|1.5]\n"
              "       [--atmosphere grey|cond-solar-proxy]\n"
              "Defaults: 128, 0.1, 0.2. Experimental static trial; convergence is not guaranteed.\n"
              "Default EOS: ionized. Default tau_top: 1e-6 (ionized), 1e-3 (cms19/freeeos).\n"
              "JSON on stdout, diagnostics on stderr; exit 0 only on convergence.");
    return 0;
  }
  try {
    std::vector<const char*> positional;
    std::string eos_name = "ionized";
    std::string hot_opacity = "opal";
    std::string atmosphere_name = "grey";
    double seed_index = 3;
    double tau_top = 0;
    bool tau_given = false;
    for (int i = 1; i < argc; ++i) {
      const std::string_view arg(argv[i]);
      if (arg == "--eos" || arg == "--tau-top" || arg == "--hot-opacity" || arg == "--seed-index"
          || arg == "--atmosphere") {
        if (++i == argc) throw std::invalid_argument("missing option value");
        if (arg == "--eos") eos_name = argv[i];
        else if (arg == "--hot-opacity") hot_opacity = argv[i];
        else if (arg == "--seed-index") seed_index = parse<double>(argv[i]);
        else if (arg == "--atmosphere") atmosphere_name = argv[i];
        else { tau_top = parse<double>(argv[i]); tau_given = true; }
      } else if (arg.starts_with("--")) throw std::invalid_argument("unknown option");
      else positional.push_back(argv[i]);
    }
    if (positional.size() > 3 || (eos_name != "ionized" && eos_name != "cms19" && eos_name != "freeeos")
        || (hot_opacity != "opal" && hot_opacity != "tops")
        || (atmosphere_name != "grey" && atmosphere_name != "cond-solar-proxy"))
      throw std::invalid_argument("invalid arguments; see --help");
    const bool cond = atmosphere_name == "cond-solar-proxy";
    if (cond && tau_given) throw std::invalid_argument("--tau-top applies only to the grey atmosphere");
    const bool cms = eos_name == "cms19";
    const bool free = eos_name == "freeeos";
    if (!tau_given) tau_top = cms || free ? 1e-3 : 1e-6;
    const std::size_t points = positional.size() > 0 ? parse<std::size_t>(positional[0]) : 128;
    const double mass = (positional.size() > 1 ? parse<double>(positional[1]) : 0.1) * constants::Msun;
    const double radius = (positional.size() > 2 ? parse<double>(positional[2]) : 0.2) * constants::Rsun;
    const std::string data = EMBER_DATA_DIR;
    std::unique_ptr<Eos> eos;
    if (cms) eos = std::make_unique<Cms19Eos>(data + "/eos/cms19_h_tp.dat", data + "/eos/cms19_he_tp.dat",
                                             Cms19Eos::Metals::helium_proxy);
    else if (free) eos = std::make_unique<HelmholtzTableEos>(data + "/eos/freeeos300_hhe_x070_potential.dat",
                                             HelmholtzTableEos::Mixture::allow_documented_proxy);
    else eos = std::make_unique<CompositeEos>();
    AesopusOpacity low(data + "/opacity/aesopus21_gs98_z020.dat");
    std::unique_ptr<Opacity> high;
    if (hot_opacity == "tops") high = std::make_unique<TopsOpacity>(data + "/opacity");
    else high = std::make_unique<OpalOpacity>(data + "/opacity/opal_gs98_z020.dat");
    // Retain the dense AESOPUS gas domain through the ionization layers;
    // require both tables only in their upper-temperature overlap.
    BlendedOpacity opacity(low, *high, 4.4, 4.5);
    PPChains nuclear;
    std::unique_ptr<Atmosphere> atmosphere;
    if (cond) atmosphere = std::make_unique<TabulatedAtmosphere>(*eos,
        data + "/atmosphere/cond_gn93_tau100_solar_proxy.dat",
        TabulatedAtmosphere::Mixture::allow_documented_proxy);
    else atmosphere = std::make_unique<GreyAtmosphere>(*eos, opacity, GreyAtmosphereOptions{.tau_top = tau_top});
    Physics physics{eos.get(), &opacity, &nuclear, 1.9};
    const auto seed = example::stellar_seed(points, mass, radius, solar_scaled(0.7, 0.02), nuclear, *atmosphere, seed_index);
    RelaxationOptions options; options.max_iterations = 100;
    const auto result = relax(seed, physics, *atmosphere, options);
    const auto& m = result.model;
    const double Teff = std::pow(m.y.back().L / (4.0 * M_PI * constants::sigma_SB
                                              * std::pow(m.r(points - 1), 2)), 0.25);
    std::fprintf(stderr, "%s\naccepted updates=%zu residual=%.8g correction=%.8g\n"
                 "R/Rsun=%.8g L/Lsun=%.8g Teff=%.8g Tc=%.8g rhoc=%.8g\n",
                 result.message.c_str(), result.iterations, result.residual, result.correction,
                 m.r(points - 1) / constants::Rsun, m.y.back().L / constants::Lsun,
                 Teff, m.T(0), m.rho(0));
    double L_nuclear = m.m[0] * nuclear.eval(m.T(0), m.rho(0), m.comp[0]).eps;
    bool increasing_radius = true, decreasing_temperature = true;
    double max_logR = -1000.0;
    double max_maxwell_defect = 0;
    double previous_pressure = 0, previous_defect = 0, defect_integral = 0;
    double pressure_integral = 0, gravity_integral = 0;
    std::size_t max_defect_point = 0;
    for (std::size_t i = 0; i < points; ++i) {
      max_logR = std::max(max_logR, std::log10(m.rho(i)) - 3.0 * (std::log10(m.T(i)) - 6.0));
      const auto e = eos->eval(m.T(i), m.rho(i), m.comp[i]);
      const double defect = e.P * e.delta / (m.rho(i) * m.T(i) * e.cp * e.grad_ad) - 1;
      if (std::abs(defect) > max_maxwell_defect) {
        max_maxwell_defect = std::abs(defect); max_defect_point = i;
      }
      if (i == 0) {
        pressure_integral = 3 * m.m[0] * e.P / m.rho(0);
        gravity_integral = .6 * constants::G * m.m[0] * m.m[0] / m.r(0);
        defect_integral = m.m[0] * defect * defect;
      } else {
        const double dm = m.m[i] - m.m[i - 1];
        pressure_integral += 1.5 * dm * (previous_pressure / m.rho(i - 1) + e.P / m.rho(i));
        gravity_integral += .5 * constants::G * dm * (m.m[i - 1] / m.r(i - 1) + m.m[i] / m.r(i));
        defect_integral += .5 * dm * (previous_defect * previous_defect + defect * defect);
      }
      previous_pressure = e.P; previous_defect = defect;
      if (i == 0) continue;
      increasing_radius &= m.y[i].lnr > m.y[i - 1].lnr;
      decreasing_temperature &= m.y[i].lnT < m.y[i - 1].lnT;
      L_nuclear += 0.5 * (m.m[i] - m.m[i - 1])
          * (nuclear.eval(m.T(i), m.rho(i), m.comp[i]).eps
             + nuclear.eval(m.T(i - 1), m.rho(i - 1), m.comp[i - 1]).eps);
    }
    const double surface_term = 4 * M_PI * std::pow(m.r(points - 1), 3) * previous_pressure;
    const double virial_error = (pressure_integral - surface_term) / gravity_integral - 1;
    char tau_literal[32]; std::snprintf(tau_literal, sizeof tau_literal, "%.16g", tau_top);
    std::printf("{\n  \"calculation\": \"experimental static stellar equilibrium trial\",\n"
                "  \"physics\": {\"eos\": \"%s\",\n"
                "    \"opacity\": \"AESOPUS 2.1 gas + %s GS98, Z=0.020, no extrapolation\",\n"
                "    \"opacity_blend_logT\": [4.4, 4.5],\n"
                "    \"atmosphere\": \"%s\",\n"
                "    \"tau_top\": %s, \"tau_match\": %.16g,\n"
                "    \"atmosphere_composition_approximation\": \"%s\",\n"
                "    \"nuclear\": \"pp chains, fixed composition, X_He3=0\",\n"
                "    \"composition\": \"X=0.7, Z=0.02, AAG21 resolved metals (GS98 opacity approximation)\",\n"
                "    \"alpha_mlt\": 1.9},\n"
                "  \"limitations\": \"%s; thin atmosphere; interior has no grain opacity, conduction, mixing, or evolution; alpha uncalibrated\",\n"
                "  \"converged\": %s, \"message\": ",
                cms ? "CMS19 pressure/entropy + radiation; metals represented by helium; static only"
                    : free ? "FreeEOS 3.0 EOS1 material Helmholtz potential + radiation; metals represented by helium; fixed composition"
                    : "ionized ions + radiation + FD electrons", hot_opacity == "tops" ? "TOPS ATOMIC" : "OPAL",
                cond ? "AMES-COND-2000 non-grey via MESA; untouched cool-dwarf cells" : "radiative Eddington grey",
                cond ? "null" : tau_literal, cond ? 100.0 : 2.0 / 3.0,
                cond ? "GN93 solar proxy; atmosphere and interior mixtures are not matched" : "none beyond EOS and opacity choices",
                cms ? "CMS19 thermodynamic consistency is approximate; internal energy disabled"
                    : free ? "fixed-composition EOS; source-fit joins regularized by C2 interpolation (see docs/FREEEOS.md); He3 and composition evolution unavailable"
                    : "EOS lacks partial ionization and molecules",
                result.converged ? "true" : "false");
    json_string(result.message);
    std::printf(",\n  \"points\": %zu, \"iterations\": %zu, \"residual\": ", points, result.iterations);
    number(result.residual); std::printf(", \"correction\": "); number(result.correction);
    std::printf(",\n  \"relative_virial_error\": %.16g, \"mass_rms_eos_maxwell_defect\": %.16g,\n"
                "  \"max_eos_defect_location\": {\"mass_fraction\": %.16g, \"temperature_K\": %.16g, \"density_g_cm3\": %.16g}",
                virial_error, std::sqrt(defect_integral / mass), m.m[max_defect_point] / mass,
                m.T(max_defect_point), m.rho(max_defect_point));
    std::printf(",\n  \"mass_g\": %.16g, \"seed_radius_cm\": %.16g,\n"
                "  \"seed_polytropic_index\": %.16g, \"max_eos_maxwell_defect\": %.16g,\n"
                "  \"radius_cm\": %.16g, \"luminosity_erg_s\": %.16g, \"Teff_K\": %.16g,\n"
                "  \"integrated_nuclear_luminosity_erg_s\": %.16g,\n"
                "  \"inner_mass_fraction\": %.16g, \"max_logR\": %.16g,\n"
                "  \"increasing_radius\": %s, \"decreasing_temperature\": %s,\n"
                "  \"history\": [\n", mass, radius, seed_index, max_maxwell_defect, m.r(points - 1), m.y.back().L, Teff,
                L_nuclear, m.m[0] / mass, max_logR, increasing_radius ? "true" : "false",
                decreasing_temperature ? "true" : "false");
    for (std::size_t i = 0; i < result.history.size(); ++i) {
      const auto& h = result.history[i];
      std::printf("    {\"residual\": %.16g, \"correction\": %.16g, \"damping\": %.16g, \"linear_error\": %.16g}%s\n",
                  h.residual, h.correction, h.damping, h.linear_error, i + 1 == result.history.size() ? "" : ",");
    }
    std::printf("  ],\n  \"columns\": [\"mass_g\", \"radius_cm\", \"density_g_cm3\", \"temperature_K\", \"luminosity_erg_s\"],\n  \"profile\": [\n");
    for (std::size_t i = 0; i < points; ++i)
      std::printf("    [%.16g, %.16g, %.16g, %.16g, %.16g]%s\n", m.m[i], m.r(i), m.rho(i),
                  m.T(i), m.y[i].L, i + 1 == points ? "" : ",");
    std::printf("  ]\n}\n");
    return result.converged ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    std::printf("{\"converged\": false, \"message\": "); json_string(e.what()); std::printf("}\n");
    return 1;
  }
}
