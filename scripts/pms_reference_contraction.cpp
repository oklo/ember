// Contracting PMS seed and first coupled interval with specified initial D/He3.
// No post-seed fuel injection. Requires a separately accepted composition-
// dependent atmosphere covering the actual model; there is no grey fallback.
// CLI: six material inputs, atmosphere, entropy-loss rate, years,
//      nine-species composition file, radius/Rsun, trial Teff, threads, work,
//      optional --seed-only (no evolutionary interval is attempted).
#include "ember/stellar_seed.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/evolution.hpp"
#include "ember/deuterium.hpp"
#include "ember/eos_deuterium.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include "ember/evolution_checkpoint.hpp"
#include <ctime>
#include <fstream>
#include <sstream>
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

struct ReportedAtmosphere final : Atmosphere {
  const Atmosphere& source;
  explicit ReportedAtmosphere(const Atmosphere& a):source(a){}
  AtmosphereState eval(double T,double g,const Composition& c) const override {
    try {return source.eval(T,g,c);}
    catch(const std::domain_error& e) {
      std::ostringstream message;message<<std::setprecision(17)<<e.what()
        << "; surface query Teff=" << T << " logg=" << std::log10(g)
        << " XH=" << c[Species::H1] << " XHe3=" << c[Species::He3]
        << " XD=" << c[Species::H2];
      throw std::domain_error(message.str());
    }
  }
  const char* name() const override {return source.name();}
};

int main(int argc, char** argv) {
  if (argc != 15 && argc != 16) return 2;
  const bool seed_only=argc==16 && std::string(argv[15])=="--seed-only";
  if(argc==16 && !seed_only)return 2;
  std::ostringstream report; report << std::setprecision(17);
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
    TraceDeuteriumAtmosphere isotope_atmosphere(eos,table_atmosphere);
    ReportedAtmosphere atmosphere(isotope_atmosphere);
    PPDeuterium nuclear;
    ContractingSeed seed_source(nuclear, entropy_loss);
    Physics physics{&eos, &opacity, &nuclear, 1.9, ConvectiveCriterion::ledoux};
    auto seed_physics = physics; seed_physics.nuclear = &seed_source;
    Composition c{}; c.basis = AbundanceBasis::baryon_mass;
    c.metal_inventory = MetalInventory::gs98;
    std::ifstream composition_input(argv[10]);
    for (auto& x:c.X) {
      composition_input >> x;
      if(!composition_input || !std::isfinite(x) || x<0)throw std::invalid_argument("invalid initial isotope input");
    }
    std::string extra;
    if(composition_input>>extra || std::abs(c.sum()-1)>2e-12 || c[Species::H2]<=0)
      throw std::invalid_argument("initial isotope input is not normalized or lacks D");
    const double radius=std::stod(argv[11]),teff=std::stod(argv[12]);
    const auto threads=std::stoul(argv[13]);
    if(!std::isfinite(radius+teff) || radius<=0 || teff<=0 || threads<1 || threads>4)
      throw std::invalid_argument("invalid physical starting controls");
    const std::filesystem::path work(argv[14]);
    if(!std::filesystem::is_directory(work))throw std::invalid_argument("missing output directory");
    for(const auto* name:{"seed","full","half1","final"})
      if(std::filesystem::exists(work/(std::string(name)+".checkpoint")))throw std::invalid_argument("checkpoint output already exists");
    driver::Identities identities{{"executable",driver::file_identity(argv[0])}};
    for(int i=1;i<=7;++i)identities[std::filesystem::canonical(argv[i]).string()]=driver::file_identity(argv[i]);
    identities[std::filesystem::canonical(argv[10]).string()]=driver::file_identity(argv[10]);
    const driver::Selections selections{"PMS contraction with initial D","PPDeuterium","GS98 variable EOS", "composition atmosphere", "Ledoux instantaneous convection"};
    const auto checkpoint=[&](const char* name,const Model& m,std::size_t accepted) {
      driver::Checkpoint state{m,years*31557600.,accepted,0};
      driver::write_checkpoint(work/(std::string(name)+".checkpoint"),state,selections,1e-14,identities);
    };
    auto guess = example::stellar_seed(512, .1*constants::Msun, radius*constants::Rsun,
        c, seed_source, atmosphere, 1.5, teff);
    RelaxationOptions ro; ro.zone_threads = threads;
    const auto start = std::clock();
    auto seed = relax(guess, seed_physics, atmosphere, ro);
    report << "{\"seed_converged\":" << seed.converged << ",\"seed_message\":" << std::quoted(seed.message);
    report << ",\"seed_only\":" << seed_only << ",\"seed_iterations\":" << seed.iterations
      << ",\"seed_residual\":" << seed.residual << ",\"seed_history\":[";
    bool first_iteration=true;
    for(const auto& h:seed.history) {
      if(!first_iteration)report<<',';first_iteration=false;
      report << "{\"residual\":" << h.residual << ",\"correction\":" << h.correction
        << ",\"damping\":" << h.damping << ",\"linear_error\":" << h.linear_error << '}';
    }
    report << ']';
    if (!seed.converged) {
      const auto& m=seed.model;
      const double r=m.r(m.size()-1);
      report << ",\"last_radius_Rsun\":" << r/constants::Rsun
        << ",\"last_Teff_K\":" << std::pow(m.y.back().L/(4*M_PI*constants::sigma_SB*r*r),.25)
        << ",\"last_central_T_K\":" << m.T(0) << "}\n";
      std::cout << report.str(); return 1;
    }
    // Relaxation holds the actual initial isotope abundances fixed. No
    // deuterium or helium is inserted after constructing the seed.
    for(const auto& comp:seed.model.comp)if(comp!=c)throw std::logic_error("seed changed initial fuel");
    checkpoint("seed",seed.model,0);
    if(seed_only) {
      const double r=seed.model.r(seed.model.size()-1),lum=seed.model.y.back().L;
      report << ",\"radius_Rsun\":" << r/constants::Rsun
        << ",\"luminosity_Lsun\":" << lum/constants::Lsun
        << ",\"Teff_K\":" << std::pow(lum/(4*M_PI*constants::sigma_SB*r*r),.25)
        << ",\"central_T_K\":" << seed.model.T(0)
        << ",\"entropy_loss_seed_only\":" << entropy_loss
        << ",\"accepted_evolutionary_intervals\":0,\"cpu_seconds\":"
        << double(std::clock()-start)/CLOCKS_PER_SEC << "}\n";
      std::cout << report.str();return 0;
    }
    const double dt = years*31557600.;
    EvolutionOptions options; options.relaxation = ro;options.abundance_tolerance=1e-14;
    const auto full = evolve_step(seed.model, physics, atmosphere, dt, options);
    const auto half1 = evolve_step(seed.model, physics, atmosphere, dt/2, options);
    const auto half2 = half1.converged ? evolve_step(half1.model, physics, atmosphere, dt/2, options) : half1;
    if(full.converged)checkpoint("full",full.model,1);
    if(half1.converged)checkpoint("half1",half1.model,1);
    if(half1.converged && half2.converged)checkpoint("final",half2.model,2);
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
      report << ',' << std::quoted(name) << ":{\"radius_Rsun\":" << r/constants::Rsun
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
    report << ",\"years\":" << years << ",\"entropy_loss_seed_only\":" << entropy_loss
      << ",\"maximum_log_structure_difference\":" << structure
      << ",\"maximum_abundance_difference\":" << abundance << ",\"steps\":[";
    bool first=true;
    for(const auto* s:{&full,&half1,&half2}) {
      if(!first)report<<','; first=false;
      report << "{\"converged\":" << s->converged << ",\"message\":" << std::quoted(s->message)
        << ",\"first_law_error\":" << s->luminosity_balance
        << ",\"nuclear_mass_error\":" << s->nuclear_mass_balance
        << ",\"gravitational_luminosity\":" << s->gravitational_luminosity
        << ",\"convective_mass_fraction\":" << s->convective_mass_fraction << '}';
    }
    report << "],\"cpu_seconds\":" << double(std::clock()-start)/CLOCKS_PER_SEC << "}\n";
    std::cout << report.str();
    return full.converged && half1.converged && half2.converged ? 0 : 1;
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cout << "{\"failed\":true,\"message\":" << std::quoted(e.what()) << "}\n"; return 1;
  }
}
