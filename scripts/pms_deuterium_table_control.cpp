// Numerical thermal coupling check using a relaxed late-PMS structure.
// A declared small D inventory is injected for this test only. This is NOT
// the initial condition or trajectory of the continuous physical PMS model.
#include "ember/stellar_seed.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/evolution.hpp"
#include "ember/deuterium.hpp"
#include "ember/eos_deuterium.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace ember;

struct ContractingSeed final : Nuclear {
  const Nuclear& nuclear;
  const double entropy_loss;
  ContractingSeed(const Nuclear& n, double rate) : nuclear(n), entropy_loss(rate) {}
  NuclearState eval(double T, double rho, const Composition& c) const override {
    auto s = nuclear.eval(T, rho, c);
    const double cooling = T * entropy_loss, total = s.eps + cooling;
    s.dlneps_dlnT = (s.eps * s.dlneps_dlnT + cooling) / total;
    s.dlneps_dlnRho *= s.eps / total;
    s.eps = total;
    return s;
  }
  const char* name() const override { return "initial-model entropy loss only"; }
};

int main(int argc, char** argv) {
  if (argc != 11) return 2;
  std::cout << std::setprecision(17);
  try {
    const double entropy_loss = std::stod(argv[8]), years = std::stod(argv[9]);
    if (!std::isfinite(entropy_loss) || !std::isfinite(years) ||
        entropy_loss <= 0 || years <= 0) throw std::invalid_argument("positive finite controls required");
    VariableMetalHelmholtzEos table_eos(argv[1], HelmholtzTableEos::Mixture::allow_documented_proxy);
    DeuteriumApproxEos eos(table_eos);
    MixtureOpacity molecular(argv[2]), warm(argv[3]), bridge(argv[4]), hot(argv[5]);
    BlendedOpacity mid(warm, bridge, 5.05, 5.10), upper(mid, hot, 5.6, 5.7);
    BlendedOpacity raw(molecular, upper, 4.4, 4.47);
    auto radiation = std::make_shared<ElementalOpacity>(raw);
    TabulatedConduction table(argv[6]);
    CombinedOpacity opacity(radiation, std::make_shared<HotConduction>(table));
    CompositionAtmosphereGrid table_atmosphere(eos, argv[7], CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    TraceDeuteriumAtmosphere atmosphere(eos,table_atmosphere);
    PPDeuterium nuclear;
    ContractingSeed seed_source(nuclear, entropy_loss);
    Physics physics{&eos, &opacity, &nuclear, 1.9, ConvectiveCriterion::ledoux};
    auto seed_physics = physics; seed_physics.nuclear = &seed_source;
    auto c = solar_scaled(.7, .02); c.basis = AbundanceBasis::baryon_mass;
    c.metal_inventory = MetalInventory::gs98;
    auto guess = example::stellar_seed(512, .1*constants::Msun, .165*constants::Rsun,
        c, seed_source, atmosphere, 1.5, 2800.);
    RelaxationOptions ro; ro.zone_threads = 2;
    const auto start = std::clock();
    auto seed = relax(guess, seed_physics, atmosphere, ro);
    std::cout << "{\"seed_converged\":" << seed.converged << ",\"seed_message\":" << std::quoted(seed.message);
    if (!seed.converged) { std::cout << "}\n"; return 1; }
    const double D=std::stod(argv[10]);
    if(!(D>0 && D<=1e-4))throw std::domain_error("invalid test deuterium inventory");
    for(auto& comp:seed.model.comp){comp.X[0]-=D;comp[Species::H2]=D;}
    const double dt = years*31557600.;
    EvolutionOptions options; options.relaxation = ro;options.abundance_tolerance=1e-14;
    const auto full = evolve_step(seed.model, physics, atmosphere, dt, options);
    const auto half1 = evolve_step(seed.model, physics, atmosphere, dt/2, options);
    const auto half2 = half1.converged ? evolve_step(half1.model, physics, atmosphere, dt/2, options) : half1;
    const auto record = [&](const char* name, const Model& m) {
      const auto w = nodal_mass_weights(m); double entropy=0, lnuc=0, u=0, deuterium=0, deuterium_power=0;
      for (std::size_t i=0; i<m.size(); ++i) {
        const auto e=eos.eval(m.T(i), m.rho(i), m.comp[i]);
        deuterium+=w[i]*m.comp[i][Species::H2];entropy+=w[i]*e.S; u+=w[i]*e.E;
        lnuc+=w[i]*nuclear.eval(m.T(i),m.rho(i),m.comp[i]).eps;
        auto other=m.comp[i];const double D=other[Species::H2];other[Species::H2]=0;
        deuterium_power+=w[i]*deuterium_capture(m.T(i),m.rho(i),other,D).source.eps;
      }
      const double r=m.r(m.size()-1), lum=m.y.back().L;
      std::cout << ',' << std::quoted(name) << ":{\"radius_Rsun\":" << r/constants::Rsun
        << ",\"luminosity_Lsun\":" << lum/constants::Lsun
        << ",\"Teff_K\":" << std::pow(lum/(4*M_PI*constants::sigma_SB*r*r), .25)
        << ",\"central_T_K\":" << m.T(0) << ",\"mean_entropy\":" << entropy/m.M
        << ",\"deuterium_g\":" << deuterium << ",\"deuterium_luminosity\":" << deuterium_power
        << ",\"internal_energy\":" << u << ",\"nuclear_luminosity\":" << lnuc
        << ",\"nuclear_fraction\":" << lnuc/lum << '}';
    };
    record("seed",seed.model);record("full",full.model);record("half1",half1.model); record("final",half2.model);
    double structure=0, abundance=0;
    for (std::size_t i=0;i<seed.model.size();++i) {
      structure=std::max({structure,std::abs(full.model.y[i].lnr-half2.model.y[i].lnr),
        std::abs(full.model.y[i].lnrho-half2.model.y[i].lnrho),std::abs(full.model.y[i].lnT-half2.model.y[i].lnT)});
      for(std::size_t j=0;j<NSPEC;++j) abundance=std::max(abundance,
        std::abs(full.model.comp[i].X[j]-half2.model.comp[i].X[j]));
    }
    std::cout << ",\"years\":" << years << ",\"entropy_loss_seed_only\":" << entropy_loss
      << ",\"maximum_log_structure_difference\":" << structure
      << ",\"maximum_abundance_difference\":" << abundance << ",\"steps\":[";
    bool first=true;
    for(const auto* s:{&full,&half1,&half2}) {
      if(!first)std::cout<<','; first=false;
      std::cout << "{\"converged\":" << s->converged << ",\"message\":" << std::quoted(s->message)
        << ",\"first_law_error\":" << s->luminosity_balance
        << ",\"nuclear_mass_error\":" << s->nuclear_mass_balance
        << ",\"gravitational_luminosity\":" << s->gravitational_luminosity
        << ",\"convective_mass_fraction\":" << s->convective_mass_fraction << '}';
    }
    std::cout << "],\"cpu_seconds\":" << double(std::clock()-start)/CLOCKS_PER_SEC << "}\n";
    return full.converged && half1.converged && half2.converged ? 0 : 1;
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}
