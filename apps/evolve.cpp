#include "../examples/stellar_seed.hpp"
#include "evolution_checkpoint.hpp"
#include "ember/eos_composition.hpp"
#include "ember/eos_mixture.hpp"
#include "ember/evolution.hpp"
#include "ember/evolution_proxies.hpp"
#include "ember/atmosphere_composition.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <future>
#include <string>
#include <stdexcept>

namespace {
// Each worker owns its cache for repeated residual and Jacobian boundary
// queries. Every composition component is part of the key.
class CachedAtmosphere final : public ember::Atmosphere {
public:
  explicit CachedAtmosphere(const ember::Atmosphere& source):source_(source) {}
  ember::AtmosphereState eval(double T,double g,const ember::Composition& c) const override {
    for(const auto& e:entries_)if(e.T==T && e.g==g && e.c.basis==c.basis && e.c.metal_inventory==c.metal_inventory && e.c.X==c.X)return e.state;
    const auto state=source_.eval(T,g,c);
    if(entries_.size()==16)entries_.erase(entries_.begin());
    entries_.push_back({T,g,c,state});return state;
  }
  const char* name() const override {return source_.name();}
private:
  struct Entry {double T,g;ember::Composition c;ember::AtmosphereState state;};
  const ember::Atmosphere& source_;
  mutable std::vector<Entry> entries_;
};
// Each worker owns this instance. Adjacent zones share an identical endpoint,
// and burning residual/Jacobian calls frequently repeat the same state.
// Retain complete PP responses keyed by every physical input; no approximate
// lookup, rounded key or reuse across a changed composition is allowed.
class CachedNuclear final : public ember::Nuclear {
public:
  explicit CachedNuclear(const ember::PPChains& source):source_(source) {}
  ember::NuclearResponse composition_response(double T,double rho,const ember::Composition& c) const override {
    for(const auto& e:entries_)if(e.T==T && e.rho==rho && e.c.basis==c.basis
        && e.c.metal_inventory==c.metal_inventory && e.c.X==c.X)return e.response;
    const auto response=source_.composition_response(T,rho,c);
    if(entries_.size()==16)entries_.erase(entries_.begin());
    entries_.push_back({T,rho,c,response});return response;
  }
  ember::NuclearState eval(double T,double rho,const ember::Composition& c) const override {
    return composition_response(T,rho,c).state;
  }
  const char* name() const override { return source_.name(); }
private:
  struct Entry {double T,rho;ember::Composition c;ember::NuclearResponse response;};
  const ember::PPChains& source_;
  mutable std::vector<Entry> entries_;
};
void json_string(const std::string& value) {
  std::putchar('"');
  for(unsigned char ch:value) {
    if(ch=='"' || ch=='\\') {std::putchar('\\');std::putchar(ch);}
    else if(ch<32) std::printf("\\u%04x",ch);
    else std::putchar(ch);
  }
  std::putchar('"');
}
template<class T> T number(const char* text) {
  const std::string value(text);T out{};
  const auto [end,error]=std::from_chars(value.data(),value.data()+value.size(),out);
  if(error!=std::errc{} || end!=value.data()+value.size()) throw std::invalid_argument("invalid number");
  return out;
}
}
int main(int argc,char** argv) {
  using namespace ember;
  if(argc==2 && std::string(argv[1])=="--help") {
    std::puts("usage: ember-evolve [points [duration_years [initial_step_years [tolerance_scale [nuclear_model [transport [atmosphere [eos [criterion]]]]]]]]]\nDefaults: 512, 1e8, 1e7, 1, sfii-svh, wd, cond-corrected. Fixed baryonic 0.1 Msun, X=.7 Z=.02 initially.\nNuclear: sfii-svh, sfii-debye, sfii-legacy-screening, legacy.\nTransport: wd (weakly damped 2021 conduction), classic, undamped, none, early (historical bounded physics).\nAtmosphere: cond-corrected, grey-convective, cond-y076, cond-top002, cond-alpha15, nongrey:/path/to/atmosphere.dat.\nEOS: proxy (historical H/He), metal:/path/to/family.dat (GS98 H/He3). Convection criterion: schwarzschild (default), ledoux, ledoux-diffusive (Langer alpha=.1, Kippenhahn alpha=1).\nImplicit pp burning and instantaneous convective mixing, step-doubling control.\nExtended X/Z opacity with elemental isotope mapping; grey convective composition correction anchored to solar COND.\nJSON on stdout, progress on stderr; exit zero only at requested duration.\nTrailing options: --checkpoint FILE writes the latest accepted state atomically; use a new path. --restart FILE restores internal variables and next step; duration is the target age. Physics, tolerances, executable bytes and original tables must match. Invoke with an executable file path.\n--checkpoint-after N writes once after N accepted steps in this invocation while continuing the calculation, for exact restart comparisons. Restart JSON history contains only the continued segment.\n--opacity-directory DIR selects a separately validated AESOPUS/TOPS family for extended transport; default is the original data/opacity. Selection and every source plane are bound to restart identity.\n--thermal-neutrinos none|plasma-hrw selects thermal losses (default none); plasma-hrw is plasma decay only, not all thermal channels.\n--step-workers 1|2 uses one or two CPU workers for the independent full-step and half-step estimates (default 1). Physics, accuracy tests and accepted-state ordering are unchanged.");return 0;
  }
  try {
    std::string restart_path,checkpoint_path,opacity_directory,thermal_neutrinos;
    unsigned step_workers=1;
    std::size_t checkpoint_after=0;bool after_requested=false,flags_started=false;
    std::vector<char*> positional{argv[0]};
    for(int i=1;i<argc;++i) {
      const std::string option=argv[i];
      if(option.starts_with("--")) {
        flags_started=true;
        if(i+1==argc || std::string(argv[i+1]).empty())throw std::invalid_argument("driver option requires a value");
        if(option=="--restart" && restart_path.empty())restart_path=argv[++i];
        else if(option=="--step-workers") {
          step_workers=number<unsigned>(argv[++i]);
          if(step_workers<1 || step_workers>2)throw std::invalid_argument("step-workers must be 1 or 2");
        }
        else if(option=="--checkpoint" && checkpoint_path.empty())checkpoint_path=argv[++i];
        else if(option=="--opacity-directory" && opacity_directory.empty())opacity_directory=argv[++i];
        else if(option=="--thermal-neutrinos" && thermal_neutrinos.empty())thermal_neutrinos=argv[++i];
        else if(option=="--checkpoint-after" && !after_requested) {
          checkpoint_after=number<std::size_t>(argv[++i]);after_requested=true;
          if(checkpoint_after==0 || checkpoint_after>10000)throw std::invalid_argument("invalid checkpoint step");
        } else throw std::invalid_argument("unknown or repeated driver option");
      } else {
        if(flags_started)throw std::invalid_argument("positional arguments must precede driver options");
        positional.push_back(argv[i]);
      }
    }
    if(after_requested && checkpoint_path.empty())throw std::invalid_argument("checkpoint-after requires checkpoint output");
    if(thermal_neutrinos.empty())thermal_neutrinos="none";
    if(thermal_neutrinos!="none" && thermal_neutrinos!="plasma-hrw")
      throw std::invalid_argument("unknown thermal-neutrino model");
    if(!checkpoint_path.empty() && std::filesystem::exists(checkpoint_path))
      throw std::invalid_argument("use a new checkpoint output path");
    argc=static_cast<int>(positional.size());argv=positional.data();
    if(argc>10) throw std::invalid_argument("too many arguments");
    const std::size_t points=argc>1?number<std::size_t>(argv[1]):512;
    constexpr double year=365.25*86400.;
    const double duration=(argc>2?number<double>(argv[2]):1e8)*year;
    double dt=(argc>3?number<double>(argv[3]):1e7)*year;
    const double tolerance_scale=argc>4?number<double>(argv[4]):1.;
    const std::string nuclear_model=argc>5?argv[5]:"sfii-svh";
    const std::string transport_model=argc>6?argv[6]:"wd";
    const std::string atmosphere_model=argc>7?argv[7]:"cond-corrected";
    const std::string eos_model=argc>8?argv[8]:"proxy";
    const std::string criterion=argc>9?argv[9]:"schwarzschild";
    if(eos_model!="proxy" && !eos_model.starts_with("metal:"))throw std::invalid_argument("unknown EOS model");
    if(criterion!="schwarzschild" && criterion!="ledoux" && criterion!="ledoux-diffusive")throw std::invalid_argument("unknown convection criterion");
    if(transport_model!="wd" && transport_model!="classic" && transport_model!="undamped" && transport_model!="none" && transport_model!="early")
      throw std::invalid_argument("unknown transport model");
    const bool nongrey=atmosphere_model.starts_with("nongrey:");
    if(!nongrey && atmosphere_model!="cond-corrected" && atmosphere_model!="grey-convective" && atmosphere_model!="cond-y076"
        && atmosphere_model!="cond-top002" && atmosphere_model!="cond-alpha15")throw std::invalid_argument("unknown atmosphere model");
    if(transport_model=="early" && atmosphere_model!="cond-corrected")throw std::invalid_argument("early transport uses its historical frozen atmosphere");
    if(eos_model.starts_with("metal:") && (!nongrey && atmosphere_model!="grey-convective"))
      throw std::invalid_argument("metal EOS requires a nongrey grid or grey-convective atmosphere");
    if(nuclear_model!="sfii-svh" && nuclear_model!="sfii-debye" && nuclear_model!="sfii-legacy-screening" && nuclear_model!="legacy")
      throw std::invalid_argument("unknown nuclear model");
    if(points<128 || points>8192 || !std::isfinite(duration) || duration<=0 || !std::isfinite(dt) || dt<=0
        || !std::isfinite(tolerance_scale) || tolerance_scale<.001 || tolerance_scale>100)
      throw std::invalid_argument("invalid resolution or duration");
    const std::string data=EMBER_DATA_DIR;
    const bool early=transport_model=="early";
    if(early && !opacity_directory.empty())
      throw std::invalid_argument("early transport requires its historical opacity inputs");
    if(opacity_directory.empty())opacity_directory=data+"/opacity";
    opacity_directory=std::filesystem::canonical(opacity_directory).string();
    std::unique_ptr<Eos> selected_eos;
    if(eos_model.starts_with("metal:"))
      selected_eos=std::make_unique<MetalHelmholtzEos>(eos_model.substr(6),HelmholtzTableEos::Mixture::allow_documented_proxy);
    else selected_eos=std::make_unique<CompositionHelmholtzEos>(data+(early?"/eos/freeeos300_hhe_composition.dat":"/eos/freeeos300_hhe_extended.dat"),HelmholtzTableEos::Mixture::allow_documented_proxy);
    const Eos& eos=*selected_eos;
    AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat");
    TopsOpacity high(data+"/opacity",TopsOpacity::Grid::composition);
    BlendedOpacity old_radiative(low,high,4.4,4.5);
    std::shared_ptr<Opacity> radiative=early?std::shared_ptr<Opacity>(std::make_shared<NominalAbundanceOpacity>(old_radiative)):
      std::shared_ptr<Opacity>(std::make_shared<StellarMixtureOpacity>(opacity_directory));
    TabulatedConduction conduction(data+"/conduction/"+(transport_model=="classic"?"condtab21_I":
      transport_model=="undamped"?"condtab21nd":"condtab21wd")+(eos_model.starts_with("metal:")?"_metals.dat":".dat"));
    const bool conductive=!early && transport_model!="none";
    CombinedOpacity opacity(radiative,conductive?std::make_shared<HotConduction>(conduction):nullptr);
    TabulatedAtmosphere cond(eos,data+"/atmosphere/cond_gn93_tau100_solar_proxy.dat",TabulatedAtmosphere::Mixture::allow_documented_proxy);
    FrozenCompositionAtmosphere frozen(eos,cond);
    auto composition=solar_scaled(.7,.02);composition.basis=AbundanceBasis::baryon_mass;
    if(eos_model.starts_with("metal:"))composition.metal_inventory=MetalInventory::gs98;
    ConvectiveAtmosphereOptions atmosphere_options;
    if(atmosphere_model=="cond-y076")atmosphere_options.henyey_y=.076;
    if(atmosphere_model=="cond-top002")atmosphere_options.tau_top=.002;
    if(atmosphere_model=="cond-alpha15")atmosphere_options.alpha=1.5;
    ConvectiveAtmosphere column(eos,*radiative,atmosphere_options);
    CompositionCorrectedAtmosphere corrected(eos,cond,column,composition);
    std::unique_ptr<CompositionAtmosphereGrid> grid;
    if(nongrey)grid=std::make_unique<CompositionAtmosphereGrid>(eos,atmosphere_model.substr(8),
      CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    const Atmosphere& raw_atmosphere=nongrey?static_cast<const Atmosphere&>(*grid):early?static_cast<const Atmosphere&>(frozen):atmosphere_model=="grey-convective"?
      static_cast<const Atmosphere&>(column):static_cast<const Atmosphere&>(corrected);
    CachedAtmosphere atmosphere(raw_atmosphere);
    PPChains nuclear_source(nuclear_model=="legacy"?PPRates::legacy:PPRates::solar_fusion_ii,
      nuclear_model=="sfii-svh"?PPScreening::salpeter_van_horn:nuclear_model=="sfii-debye"?
      PPScreening::debye_fermi:PPScreening::legacy_weak);
    CachedNuclear nuclear(nuclear_source);
    Physics physics{&eos,&opacity,&nuclear,1.9,
      criterion=="schwarzschild"?ConvectiveCriterion::schwarzschild:ConvectiveCriterion::ledoux,
      criterion=="ledoux-diffusive"?.1:0,criterion=="ledoux-diffusive"?1.:0};
    PlasmaNeutrinoLosses plasma_losses;
    if(thermal_neutrinos=="plasma-hrw")physics.neutrino_losses=&plasma_losses;
    // The full-step estimate is independent of both half steps. Its worker
    // owns separate caches; immutable EOS, opacity and rate sources are shared.
    CachedNuclear full_nuclear(nuclear_source);
    CachedAtmosphere full_atmosphere(raw_atmosphere);
    Physics full_physics=physics;full_physics.nuclear=&full_nuclear;
    double seed_radius=.15*constants::Rsun, seed_teff=0;
    if(grid) {
      const auto support=grid->support();
      seed_teff=std::exp(.5*(std::log(support.teff[0])+std::log(support.teff[1])));
      const double seed_g=std::exp(.5*(std::log(support.gravity[0])+std::log(support.gravity[1])));
      seed_radius=std::sqrt(constants::G*.1*constants::Msun/seed_g);
    }
    const driver::Selections selections{nuclear_model,transport_model,atmosphere_model,eos_model,criterion};
    driver::Identities identities;
    if(!restart_path.empty() || !checkpoint_path.empty())
      identities=driver::input_identities(argv[0],data,atmosphere_model,eos_model,opacity_directory);
    identities["thermal_neutrinos"]=thermal_neutrinos;
    Model model;std::size_t accepted=0,rejected=0;
    if(restart_path.empty()) {
      const auto seed=example::stellar_seed(points,.1*constants::Msun,seed_radius,composition,nuclear,atmosphere,1.5,seed_teff);
      const auto initial=relax(seed,physics,atmosphere);
      if(!initial.converged) throw std::runtime_error("initial equilibrium: "+initial.message);
      model=initial.model;
    } else {
      auto saved=driver::read_checkpoint(restart_path,points,.1*constants::Msun,composition,selections,tolerance_scale,identities);
      model=std::move(saved.model);dt=saved.next_dt;accepted=saved.accepted;rejected=saved.rejected;
      if(model.age>duration)throw std::invalid_argument("requested age precedes checkpoint");
    }
    const auto accepted_before=accepted,rejected_before=rejected;
    const double starting_age=model.age;
    bool checkpoint_written=false;
    auto save=[&]() {
      driver::Checkpoint state{model,dt,accepted,rejected};
      driver::check_state(state,points,composition);
      driver::write_checkpoint(checkpoint_path,state,selections,tolerance_scale,identities);
      checkpoint_written=true;
    };
    if(!checkpoint_path.empty() && !after_requested)save();
    struct Record {double age,step,R,L,Teff,X,Y3,error,energy,mass;std::size_t iterations;double convective,Tc,rhoc,Xs,Y3s,Lnuc,Lnucnu,Lgrav,MH,Lthermalnu;};
    std::vector<Record> history;
    const auto diagnostic_weights=nodal_mass_weights(model);
    auto record=[&](double step,double error,const EvolutionStep& result) {
      const double radius=model.r(points-1),lum=model.y.back().L;
      double convective=result.convective_mass_fraction;
      double nuclear_lum=result.nuclear_luminosity,nuclear_neutrino_lum=result.neutrino_luminosity,hydrogen_mass=0;
      double thermal_neutrino_lum=result.thermal_neutrino_luminosity;
      for(std::size_t i=0;i<points;++i) {
        hydrogen_mass+=diagnostic_weights[i]*model.comp[i].h1();
        if(step==0) {
          const auto n=nuclear.eval(model.T(i),model.rho(i),model.comp[i]);
          nuclear_lum+=diagnostic_weights[i]*n.eps;
          nuclear_neutrino_lum+=diagnostic_weights[i]*n.eps_neutrino;
          thermal_neutrino_lum+=diagnostic_weights[i]*evaluate_losses(
              physics.neutrino_losses,model.T(i),model.rho(i),model.comp[i]).eps;
        }
      }
      if(step==0){for(auto [begin,end]:convective_mixing_regions(model,physics))if(end>begin+1)for(std::size_t i=begin;i<end;++i)convective+=diagnostic_weights[i]/model.M;}
      history.push_back({model.age/year,step/year,radius/constants::Rsun,lum/constants::Lsun,
        std::pow(lum/(4*M_PI*constants::sigma_SB*radius*radius),.25),model.comp[0].X[0],model.comp[0].X[1],
        error,result.luminosity_balance,result.nuclear_mass_balance,result.coupling_iterations,
        convective,model.T(0),model.rho(0),model.comp.back().X[0],model.comp.back().X[1],
        nuclear_lum/constants::Lsun,nuclear_neutrino_lum/constants::Lsun,
        result.gravitational_luminosity/constants::Lsun,hydrogen_mass/constants::Msun,
        thermal_neutrino_lum/constants::Lsun});
    };
    record(0,0,{});
    EvolutionOptions options;std::string last_failure;bool success=true;
    // Checkpoints are written after accepted steps, so a restored state has
    // no unresolved consecutive rejection history. Keep the total rejected
    // count for diagnostics without treating normal adaptive retries over a
    // trillion-year trajectory as evidence that the solver has stalled.
    std::size_t consecutive_rejected=0;
    while(model.age<duration) {
      dt=std::min(dt,duration-model.age);
      if(accepted-accepted_before>=10000) {success=false;last_failure="accepted-step limit reached in this invocation";break;}
      if(accepted==std::numeric_limits<std::size_t>::max() || rejected==std::numeric_limits<std::size_t>::max())
        {success=false;last_failure="lifetime step counter overflow";break;}
      if(consecutive_rejected>100) {success=false;last_failure="consecutive rejection limit reached: "+last_failure;break;}
      if(dt<year) {success=false;last_failure="minimum timestep reached: "+last_failure;break;}
      EvolutionStep full,first,second;
      if(step_workers==2) {
        auto pending=std::async(std::launch::async,[&] {
          return evolve_step(model,full_physics,full_atmosphere,dt,options);
        });
        first=evolve_step(model,physics,atmosphere,.5*dt,options);
        if(first.converged)second=evolve_step(first.model,physics,atmosphere,.5*dt,options);
        full=pending.get();
      } else {
        full=evolve_step(model,physics,atmosphere,dt,options);
        if(full.converged)first=evolve_step(model,physics,atmosphere,.5*dt,options);
        if(first.converged)second=evolve_step(first.model,physics,atmosphere,.5*dt,options);
      }
      double error=0;
      if(full.converged && first.converged && second.converged) {
        for(std::size_t i=0;i<points;++i) {
          error=std::max({error,std::abs(full.model.y[i].lnT-second.model.y[i].lnT)/1e-5,
            std::abs(full.model.y[i].lnrho-second.model.y[i].lnrho)/1e-5,
            std::abs(full.model.y[i].lnr-second.model.y[i].lnr)/1e-5});
          for(std::size_t j=0;j<3;++j)
            error=std::max(error,std::abs(full.model.comp[i].X[j]-second.model.comp[i].X[j])/1e-8);
        }
        error=std::max(error,std::abs(full.model.y.back().L/second.model.y.back().L-1)/1e-4);
        error/=tolerance_scale;
      } else {
        error=2;last_failure=!full.converged?full.message:!first.converged?first.message:second.message;
      }
      if(!std::isfinite(error) || error>1) {++rejected;++consecutive_rejected;dt*=.5;std::fprintf(stderr,"retry dt=%.6g yr (%s; error %.3g)\n",dt/year,last_failure.c_str(),error);continue;}
      model=second.model;++accepted;consecutive_rejected=0;record(dt,error,second);last_failure.clear();
      std::fprintf(stderr,"age=%.8g yr, X=%.10g Y3=%.10g, R=%.9g L=%.9g, convective=%.8g, step error=%.3g\n",
        model.age/year,model.comp[0].X[0],model.comp[0].X[1],history.back().R,history.back().L,second.convective_mass_fraction,error);
      dt*=std::clamp(.9/std::sqrt(std::max(error,1e-6)),.5,2.);
      if(!checkpoint_path.empty() && (!after_requested || accepted-accepted_before==checkpoint_after))save();
    }
    if(success && !checkpoint_path.empty() && !checkpoint_written)
      throw std::runtime_error("requested checkpoint step was not reached");
    std::printf("{\n\"calculation\":\"composition-dependent pp evolution with instantaneous convective mixing\",\n\"converged\":%s,\n\"message\":",success?"true":"false");
    json_string(success?"requested duration reached":last_failure);
    if(!restart_path.empty()) {
      std::printf(",\n\"restart\":{\"source\":");json_string(restart_path);
      std::printf(",\"history_start_age_yr\":%.17g,\"accepted_steps_before_restart\":%zu,\"rejected_steps_before_restart\":%zu,\"scope\":\"history contains this continuation segment only; saved internal state and next timestep restored\"}",
        starting_age/year,accepted_before,rejected_before);
    }
    std::printf(",\n\"nuclear_model\":");json_string(nuclear_model);
    std::printf(",\n\"nuclear_physics\":");json_string(nuclear.name());
    std::printf(",\n\"thermal_neutrino_model\":");json_string(thermal_neutrinos);
    std::printf(",\n\"thermal_neutrino_physics\":");json_string(physics.neutrino_losses?physics.neutrino_losses->name():"thermal neutrinos omitted");
    std::printf(",\n\"thermal_neutrino_limitations\":");json_string(physics.neutrino_losses?"plasma decay only; fully ionized approximation; pair, photo, bremsstrahlung and recombination omitted":"all thermal-neutrino channels omitted");
    std::printf(",\n\"transport_model\":");json_string(transport_model);
    std::printf(",\n\"eos_model\":");json_string(eos_model);
    std::printf(",\n\"metal_inventory\":");json_string(composition.metal_inventory==MetalInventory::gs98?
      "GS98 elemental number pattern; inert carrier labels; representative metal isotopes":"carried isotope labels with Ne-like Zrest");
    std::printf(",\n\"eos\":");json_string(eos.name());
    std::printf(",\n\"convection_criterion\":");json_string(criterion);
    std::printf(",\n\"secular_mixing\":{\"semiconvection\":\"Langer mixing-only\",\"alpha_semiconvection\":%.17g,\"thermohaline\":\"Kippenhahn\",\"alpha_thermohaline\":%.17g}",
      physics.alpha_semiconvection,physics.alpha_thermohaline);
    std::printf(",\n\"atmosphere_model\":");json_string(early?"frozen":atmosphere_model);
    std::printf(",\n\"opacity\":");json_string(radiative->name());
    std::printf(",\n\"opacity_directory\":");json_string(opacity_directory);
    std::printf(",\n\"atmosphere\":");json_string(atmosphere.name());
    std::printf(",\n\"conduction\":");json_string(conductive?conduction.name():"omitted");
    if(!nongrey)std::printf(",\n\"atmosphere_integration\":{\"log_state_tolerance\":%.17g,\"sensitivity_tolerance\":%.17g,\"tau_top\":%.17g,\"alpha\":%.17g,\"henyey_y\":%.17g}",atmosphere_options.tolerance,atmosphere_options.sensitivity_tolerance,atmosphere_options.tau_top,atmosphere_options.alpha,atmosphere_options.henyey_y);
    if(nongrey) {
      std::printf(",\n\"atmosphere_grid_approximation\":");json_string(grid->approximation());
      std::printf(",\n\"atmosphere_tau_match\":%.17g",grid->tau_match());
      const auto support=grid->support();
      std::printf(",\n\"atmosphere_grid_support\":{\"hydrogen\":[%.17g,%.17g],\"helium3\":[%.17g,%.17g],\"teff_K\":[%.17g,%.17g],\"gravity_cm_s2\":[%.17g,%.17g]}",
        support.hydrogen[0],support.hydrogen[1],support.helium3[0],support.helium3[1],
        support.teff[0],support.teff[1],support.gravity[0],support.gravity[1]);
      if(grid->has_missing_states())
        std::printf(",\n\"atmosphere_grid_support_scope\":\"Outer bounds only; cells require all sixteen validated corner states, including the derivative stencil. Missing source states are rejected.\"");
    }
    std::printf(",\n\"nuclear_limitations\":\"reduced pp network: no pep, hep, ppIII or CNO; SFII pp curvature omitted; SVH approximate, not derived from FreeEOS; modern screening requires zeta<=.2\"");
    std::printf(",\n\"history_energy_scope\":\"Nuclear deposited and nuclear-neutrino luminosities use nodal baryonic mass weights. Thermal-neutrino luminosity is the positive sink from the explicitly selected thermal model and is subtracted from deposited nuclear plus gravothermal heating. Gravothermal luminosity describes the last accepted implicit half step; null at the initial or restored state. Hydrogen mass includes all remaining core and envelope fuel.\"");
    const auto screening=pp_screening(model.T(0),model.rho(0),model.comp[0],PPReaction::pp,
      nuclear_model=="sfii-svh"?PPScreening::salpeter_van_horn:nuclear_model=="sfii-debye"?
      PPScreening::debye_fermi:PPScreening::legacy_weak);
    std::printf(",\n\"central_screening\":{\"pp_log_factor\":%.17g,\"electron_susceptibility\":%.17g,\"gamma_e\":%.17g,\"pp_zeta\":%.17g}",
      screening.log_factor,screening.electron_susceptibility,screening.gamma_e,screening.zeta);
    std::printf(",\n\"points\":%zu,\n\"mass_basis\":\"conserved baryonic mass; nuclear rest mass release accounted in energy, Newtonian gravity\",\n\"mass_Msun\":0.1,\n\"age_origin\":\"specified static composition, not age since formation\",\n\"limitations\":\"source domains enforced; He isotope cross sections approximated, GS98 metals; atmosphere/EOS/opacity source mixture approximations; selected convection criterion; semiconvective heat flux and microscopic settling omitted; fixed mass mesh; conduction mixture and envelope join approximations\",\n\"step_error_tolerances\":{\"log_structure\":%.17g,\"absolute_abundance\":%.17g,\"relative_surface_luminosity\":%.17g},\n\"rejected_steps\":%zu,\n\"columns\":[\"age_yr\",\"step_yr\",\"R_Rsun\",\"L_Lsun\",\"Teff_K\",\"central_X\",\"central_Y3\",\"step_error_norm\",\"discrete_luminosity_balance\",\"nuclear_rest_mass_balance\",\"last_halfstep_coupling_iterations\",\"convective_mass_fraction\",\"central_T_K\",\"central_rho\",\"surface_X\",\"surface_Y3\",\"nuclear_deposited_Lsun\",\"nuclear_neutrino_Lsun\",\"last_halfstep_gravothermal_Lsun\",\"hydrogen_mass_Msun\",\"thermal_neutrino_Lsun\"],\n\"history\":[\n",points,1e-5*tolerance_scale,1e-8*tolerance_scale,1e-4*tolerance_scale,rejected);
    for(std::size_t i=0;i<history.size();++i) {
      const auto& r=history[i];
      std::printf("%s[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,",
        i?",\n":"",r.age,r.step,r.R,r.L,r.Teff,r.X,r.Y3,r.error,r.energy,r.mass,r.iterations,r.convective,r.Tc,r.rhoc,r.Xs,r.Y3s,r.Lnuc,r.Lnucnu);
      if(r.step==0)std::printf("null");else std::printf("%.17g",r.Lgrav);
      std::printf(",%.17g,%.17g]",r.MH,r.Lthermalnu);
    }
    std::printf("\n],\n\"profile_columns\":[\"mass_g\",\"radius_cm\",\"density_g_cm3\",\"temperature_K\",\"luminosity_erg_s\",\"X\",\"Y3\",\"Y4\"],\n\"profile\":[\n");
    for(std::size_t i=0;i<points;++i)
      std::printf("%s[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",i?",\n":"",model.m[i],model.r(i),model.rho(i),model.T(i),model.y[i].L,model.comp[i].X[0],model.comp[i].X[1],model.comp[i].X[2]);
    std::puts("\n]}\n");return success?0:1;
  } catch(const std::exception& e) {
    std::printf("{\"converged\":false,\"message\":");json_string(e.what());std::puts("}");return 1;
  }
}
