#include "lifetime_driver.hpp"
#include "contracting_seed.hpp"
#include "lifetime_checks.hpp"
#include "convective_material_heat.hpp"
#include "evolution_checkpoint.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/eos_deuterium.hpp"
#include "ember/metal_microscopic_transport.hpp"
#include "ember/opacity_mixture.hpp"
#include <charconv>
#include <chrono>
#include <ctime>
#include <iostream>
#include <set>

namespace ember::driver {
namespace {
constexpr double year=31557600.;
struct Settings {
  std::map<std::string,std::string> values;
  explicit Settings(const fs::path& path) {
    std::ifstream in(path);std::string line;
    if(!in)throw std::invalid_argument("missing lifetime configuration");
    while(std::getline(in,line)) {
      if(line.empty() || line.front()=='#')continue;
      std::istringstream row(line);std::string key,value,extra;
      if(!(row>>key>>std::quoted(value)) || row>>extra || !values.emplace(key,value).second)
        throw std::invalid_argument("invalid or duplicate lifetime setting: "+line);
    }
  }
  std::string get(const std::string& key) {
    auto found=values.find(key);if(found==values.end())throw std::invalid_argument("missing lifetime setting: "+key);
    auto value=found->second;values.erase(found);return value;
  }
  double number(const std::string& key) {
    const auto text=get(key);double x{};auto r=std::from_chars(text.data(),text.data()+text.size(),x);
    if(r.ec!=std::errc{} || r.ptr!=text.data()+text.size() || !std::isfinite(x))
      throw std::invalid_argument("invalid numeric lifetime setting: "+key);
    return x;
  }
};
void add_identity(Identities& ids,const fs::path& path) {
  const auto p=fs::canonical(path);ids[p.string()]=file_identity(p);
}
void add_family_identities(Identities& ids,const fs::path& path,bool eos) {
  add_identity(ids,path);std::ifstream in(path);std::string label,file;std::size_t count{};
  if(eos) {
    std::getline(in,label);count=1;
    for(int k=0;k<3;++k) {
      std::size_t n{};in>>label>>n;if(n<2 || n>100)throw std::runtime_error("invalid EOS identity axis");
      count*=n;double coordinate{};for(std::size_t i=0;i<n;++i)in>>coordinate;
    }
  }else {
    int version{};in>>label>>version>>count>>label>>label;
    if(version!=1 || count<2 || count>100)throw std::runtime_error("invalid opacity identity family");
  }
  for(std::size_t i=0;i<count;++i) {
    if(!eos){double Z{};in>>Z;}
    in>>std::quoted(file);if(!in)throw std::runtime_error("missing family identity source");
    add_identity(ids,path.parent_path()/file);
  }
}
double positive(const char* text) {
  std::string s(text);double x{};auto r=std::from_chars(s.data(),s.data()+s.size(),x);
  if(r.ec!=std::errc{} || r.ptr!=s.data()+s.size() || !std::isfinite(x) || x<=0)
    throw std::invalid_argument("positive finite lifetime argument required");
  return x;
}
bool has_D(const Model& m) {return std::any_of(m.comp.begin(),m.comp.end(),[](const auto& c){return c[Species::H2]!=0;});}
}

int lifetime_main(int argc,char** argv) {
  if(argc<5) {
    std::cerr<<"usage: ember-evolve --lifetime CONFIG TARGET_YEARS NEW_OUTPUT_DIRECTORY"
      " [--restart CHECKPOINT] [--max-steps N] [--cpu-seconds S]\n"
      "Common Hayashi-to-remnant driver under integration. Source-domain failures stop; no atmosphere fallback.\n";
    return 2;
  }
  try {
    const auto wall_start=std::chrono::steady_clock::now();const auto cpu_start=std::clock();
    const fs::path config=fs::canonical(argv[2]),work=argv[4];const double target=positive(argv[3])*year;
    if(fs::exists(work))throw std::invalid_argument("lifetime output directory already exists");
    std::string restart;std::size_t maximum_steps=500;double maximum_cpu=240;
    for(int i=5;i<argc;i+=2) {
      if(i+1>=argc)throw std::invalid_argument("lifetime option needs a value");
      const std::string key=argv[i];
      if(key=="--restart" && restart.empty())restart=argv[i+1];
      else if(key=="--max-steps") {
        const double value=positive(argv[i+1]);if(value>100000 || std::floor(value)!=value)throw std::invalid_argument("invalid step budget");
        maximum_steps=static_cast<std::size_t>(value);
      }else if(key=="--cpu-seconds")maximum_cpu=positive(argv[i+1]);
      else throw std::invalid_argument("unknown lifetime argument: "+key);
    }
    Settings cfg(config);
    if(cfg.get("version")!="1")throw std::invalid_argument("unsupported lifetime configuration version");
    const auto path=[&](const char* key){return fs::canonical(config.parent_path()/cfg.get(key));};
    const auto eos_path=path("eos"),low_path=path("opacity_low"),warm_path=path("opacity_warm"),bridge_path=path("opacity_bridge"),hot_path=path("opacity_hot");
    const auto conduction_path=path("conduction"),atmosphere_path=path("atmosphere"),collision_path=path("collisions"),composition_path=path("composition");
    const double mass=cfg.number("mass_Msun")*constants::Msun,radius=cfg.number("initial_radius_Rsun")*constants::Rsun;
    const double teff=cfg.number("initial_Teff_K"),entropy_loss=cfg.number("initial_entropy_loss");
    const double count=cfg.number("points"),threads=cfg.number("zone_threads"),minimum_temperature=cfg.number("screened_minimum_T_K");
    double dt=cfg.number("initial_step_years")*year;const double maximum_dt=cfg.number("maximum_step_years")*year;
    const double structure_tolerance=cfg.number("structure_tolerance"),species_tolerance=cfg.number("species_tolerance"),energy_tolerance=cfg.number("energy_tolerance");
    const double abundance_tolerance=cfg.number("coupling_abundance_tolerance");
    const double inventory_tolerance=cfg.number("inventory_abundance_tolerance");
    if(!cfg.values.empty())throw std::invalid_argument("unknown lifetime setting: "+cfg.values.begin()->first);
    if(!(mass>0 && radius>0 && teff>0 && entropy_loss>0 && dt>0 && maximum_dt>=dt && minimum_temperature>0
        && count>=128 && count<=8192 && std::floor(count)==count && threads>=1 && threads<=4 && std::floor(threads)==threads
        && structure_tolerance>0 && structure_tolerance<=1e-3 && species_tolerance>0 && species_tolerance<=1e-6
        && energy_tolerance>0 && energy_tolerance<=.01 && abundance_tolerance>0 && abundance_tolerance<=1e-12
        && inventory_tolerance>=abundance_tolerance && inventory_tolerance<=std::min(1e-12,.01*species_tolerance)))
      throw std::invalid_argument("lifetime physical or accuracy setting out of range");
    const auto points=static_cast<std::size_t>(count);
    Composition initial;initial.basis=AbundanceBasis::baryon_mass;initial.metal_inventory=MetalInventory::gs98;
    std::ifstream input(composition_path);for(auto& x:initial.X)x=read_representable_double(input);
    std::string extra;if(input>>extra || std::abs(initial.sum()-1)>2e-12 || initial[Species::H2]<=0)
      throw std::invalid_argument("initial composition must be normalized and include primordial D");
    for(double x:initial.X)if(x<0)throw std::invalid_argument("negative initial abundance");
    initial.cn_molality=initial_gs98_cn(initial);initial=explicit_cn_material(initial);

    Identities identities{{"executable",file_identity(argv[0])}};add_identity(identities,config);
    add_family_identities(identities,eos_path,true);
    for(const auto& p:{low_path,warm_path,bridge_path,hot_path})add_family_identities(identities,p,false);
    for(const auto& p:{conduction_path,atmosphere_path,collision_path,composition_path})add_identity(identities,p);
    const Selections selections{"continuous Hayashi start; face energy","initial D + SFIII/SVH pp/CN; physical metal masses",
      "plasma neutrino losses","assessed whole-star instantaneous convection; total composition heat after initial D",
      "whole-convective transport limit; radiative species boundary requires microscopic transport"};
    VariableMetalHelmholtzEos table_eos(eos_path,HelmholtzTableEos::Mixture::allow_documented_proxy);DeuteriumApproxEos eos(table_eos);
    MixtureOpacity low(low_path),warm(warm_path),bridge(bridge_path),hot(hot_path);
    BlendedOpacity mid(warm,bridge,5.05,5.10),upper(mid,hot,5.6,5.7),raw(low,upper,4.4,4.47);
    auto radiation=std::make_shared<ElementalOpacity>(raw);TabulatedConduction table_conduction(conduction_path);
    auto conduction=std::make_shared<HotConduction>(table_conduction);CombinedOpacity combined(radiation,conduction);
    CompositionAtmosphereGrid table_atmosphere(eos,atmosphere_path,CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    TraceDeuteriumAtmosphere atmosphere(eos,table_atmosphere);
    PPCNNetwork nuclear(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,PPRates::solar_fusion_iii);
    PlasmaNeutrinoLosses losses;ScreenedCollisionTransport collisions(collision_path.string());
    ScreenedMetalMicroscopicTransport microscopic(table_eos,collisions,true,minimum_temperature,{true,true,true},true);
    ConvectiveMaterialHeat convective_heat(table_eos,*conduction);
    Physics early{&eos,&combined,&nuclear,1.9,ConvectiveCriterion::ledoux};early.neutrino_losses=&losses;early.explicit_metal_mixing_only=true;
    auto later=early;later.opacity=radiation.get();later.microscopic=&convective_heat;later.explicit_metal_mixing_only=false;
    EvolutionOptions options;options.relaxation.zone_threads=static_cast<std::size_t>(threads);options.abundance_tolerance=abundance_tolerance;
    Checkpoint state;
    if(restart.empty()) {
      ContractingSource seed_source(nuclear,entropy_loss);auto seed_physics=early;seed_physics.nuclear=&seed_source;
      const auto guess=contracting_guess(points,mass,radius,teff,initial,seed_physics,atmosphere);
      auto solved=relax(guess,seed_physics,atmosphere,options.relaxation);
      if(!solved.converged)throw std::runtime_error("Hayashi initial model: "+solved.message);
      for(const auto& c:solved.model.comp)if(c!=initial)throw std::runtime_error("initial relaxation changed isotope inventory");
      state={std::move(solved.model),dt,0,0};
    }else {
      std::ifstream checkpoint_header(restart);std::string magic;int version{};
      if(!(checkpoint_header>>magic>>version) || magic!="EMBER_EVOLUTION_CHECKPOINT" || (version!=5 && version!=6))
        throw std::invalid_argument("lifetime restart requires a volume-face checkpoint");
      state=read_checkpoint(restart,points,mass,initial,selections,abundance_tolerance,identities,
                            LuminosityGrid::volume_faces,version==6);
    }
    if(state.model.age>=target)throw std::invalid_argument("target must exceed the saved age");
    fs::create_directories(work);write_checkpoint(work/"seed.checkpoint",state,selections,abundance_tolerance,identities);
    std::ofstream history(work/"history.jsonl"),attempts(work/"attempts.jsonl");history<<std::setprecision(17);attempts<<std::setprecision(17);
    const auto record=[&](double step,double error,const HomogeneousCheck& guard) {
      const auto& m=state.model;const auto w=nodal_mass_weights(m);double H=0,Y3=0,D=0,Lnuc=0;
      for(std::size_t i=0;i<m.size();++i){H+=w[i]*m.comp[i].X[0];Y3+=w[i]*m.comp[i].X[1];D+=w[i]*m.comp[i][Species::H2];Lnuc+=w[i]*nuclear.eval(m.T(i),m.rho(i),m.comp[i]).eps;}
      const double R=m.r(m.size()-1),L=m.y.back().L;
      history<<"{\"years\":"<<m.age/year<<",\"step_years\":"<<step/year<<",\"accepted\":"<<state.accepted<<",\"rejected\":"<<state.rejected
        <<",\"radius_Rsun\":"<<R/constants::Rsun<<",\"Teff_K\":"<<std::pow(L/(4*M_PI*constants::sigma_SB*R*R),.25)
        <<",\"luminosity_Lsun\":"<<L/constants::Lsun<<",\"nuclear_fraction\":"<<Lnuc/L<<",\"central_T_K\":"<<m.T(0)
        <<",\"H_mass_g\":"<<H<<",\"He3_mass_g\":"<<Y3<<",\"D_mass_g\":"<<D<<",\"surface_X\":"<<m.comp.back().X[0]
        <<",\"surface_Z\":"<<m.comp.back().Z()<<",\"error_norm\":"<<error<<",\"initial_D_approximation\":"<<has_D(m)
        <<",\"omitted_mixing_heat_fraction\":"<<(has_D(m)?guard.maximum_heat_fraction:0)<<",\"gross_mixing_heat_fraction\":"<<guard.maximum_gross_heat_fraction
        <<",\"burn_gradient_estimate\":"<<guard.burn_gradient_estimate<<",\"drift_gradient_proxy\":"<<guard.drift_gradient_proxy
        <<",\"kinetic_heat_proxy\":"<<guard.kinetic_heat_proxy
        <<",\"convective_travel_years\":"<<guard.travel_years<<",\"cpu_seconds\":"<<double(std::clock()-cpu_start)/CLOCKS_PER_SEC<<"}\n";history.flush();
    };
    std::size_t accepted_here=0,failed=0;std::string stop;bool reached=false;double error=0;
    try {
      const auto assess=[&](const Model& m,std::span<const MetalSpeciesVector> rates) {
        const bool initial_D=has_D(m);convective_heat.diagnostic_rates=rates;
        try {
          const auto result=check_initial_convection(m,initial_D?early:later,table_eos,nuclear,initial_D);
          convective_heat.diagnostic_rates={};return result;
        }catch(...) {convective_heat.diagnostic_rates={};throw;}
      };
      HomogeneousCheck guard=assess(state.model,state.metal_heat_rates);
      record(0,0,guard);
      for(;;) {
        if(state.model.age>=target || target-state.model.age<=8*(std::nextafter(target,INFINITY)-target)){stop="requested age";reached=true;break;}
        if(accepted_here>=maximum_steps){stop="accepted-step budget";break;}
        if(double(std::clock()-cpu_start)/CLOCKS_PER_SEC>=maximum_cpu){stop="CPU budget";break;}
        const double ds=std::min({state.next_dt,maximum_dt,target-state.model.age});
        if(state.model.age+ds==state.model.age || ds<=0)throw std::runtime_error("timestep is below clock resolution");
        const bool initial_D=has_D(state.model);auto selected=options;
        selected.previous_metal_heat_rates=state.metal_heat_rates;
        selected.convective_mixing=ConvectiveMixing::instantaneous;
        const auto& physics=initial_D?early:later;
        const auto full=evolve_step(state.model,physics,atmosphere,ds,selected);
        const auto h1=full.converged?evolve_step(state.model,physics,atmosphere,ds/2,selected):full;
        auto second_options=selected;second_options.previous_metal_heat_rates=h1.total_metal_species_rates;
        const auto h2=h1.converged?evolve_step(h1.model,physics,atmosphere,ds/2,second_options):h1;
        const bool converged=full.converged && h1.converged && h2.converged;
        const auto af=check_interval(state.model,full,ds,nuclear,inventory_tolerance);
        const auto a1=check_interval(state.model,h1,ds/2,nuclear,inventory_tolerance);
        const auto a2=check_interval(h1.model,h2,ds/2,nuclear,inventory_tolerance);
        const bool audited=af.pass && a1.pass && a2.pass;error=1e30;
        HomogeneousCheck end_guard;
        if(converged) {
          error=0;
          for(std::size_t i=0;i<points;++i) {
            for(const auto v:{Var::lnr,Var::lnrho,Var::lnT})error=std::max(error,std::abs(full.model.y[i][v]-h2.model.y[i][v])/structure_tolerance);
            const auto a=metal_cn_abundances(full.model.comp[i]),b=metal_cn_abundances(h2.model.comp[i]);
            for(std::size_t k=0;k<METAL_CN_SIZE;++k)error=std::max(error,std::abs(a[k]-b[k])/species_tolerance);
          }
          const double half_power=.5*(h1.nuclear_luminosity+h2.nuclear_luminosity),power_scale=std::max(std::abs(half_power),.5*(std::abs(h1.model.y.back().L)+std::abs(h2.model.y.back().L)));
          error=std::max(error,std::abs(full.nuclear_luminosity-half_power)/(energy_tolerance*power_scale));
          error=std::max(error,std::abs(full.model.y.back().L/h2.model.y.back().L-1)/(4*structure_tolerance));
          (void)assess(full.model,full.total_metal_species_rates);
          (void)assess(h1.model,h1.total_metal_species_rates);
          end_guard=assess(h2.model,h2.total_metal_species_rates);
        }
        const bool take=converged && audited && error<=1;
        const auto reason=!full.converged?full.message:!h1.converged?h1.message:h2.message;
        attempts<<"{\"start_years\":"<<state.model.age/year<<",\"step_years\":"<<ds/year<<",\"converged\":"<<converged<<",\"audit_pass\":"<<audited
          <<",\"accepted\":"<<take<<",\"error_norm\":"<<error<<",\"message\":"<<std::quoted(reason)
          <<",\"species_error\":["<<af.maximum_species_error<<','<<a1.maximum_species_error<<','<<a2.maximum_species_error
          <<"],\"mass_error_surface\":["<<af.mass_error_surface<<','<<a1.mass_error_surface<<','<<a2.mass_error_surface
          <<"],\"first_law\":["<<af.first_law<<','<<a1.first_law<<','<<a2.first_law<<"]}\n";attempts.flush();
        if(take) {
          state.model=h2.model;state.metal_heat_rates=h2.total_metal_species_rates;++state.accepted;++accepted_here;failed=0;
          state.next_dt=std::min(maximum_dt,ds*std::clamp(.8/std::sqrt(std::max(error,1e-8)),.5,1.5));
          record(ds,error,end_guard);write_checkpoint(work/"latest.checkpoint",state,selections,abundance_tolerance,identities);
          std::cerr<<std::setprecision(4)<<"age="<<state.model.age/year<<" yr; step error="<<error<<"; accepted="<<state.accepted<<'\n';
        }else {
          ++state.rejected;state.next_dt=ds*.5;
          if(converged && !audited)throw std::runtime_error("physical inventory/energy audit failed");
          if(reason.find("radiative species boundary")!=std::string::npos)
            throw std::runtime_error(reason);
          if(++failed>=8)throw std::runtime_error("repeated step rejection: "+reason);
        }
      }
    }catch(const std::exception& e){stop=e.what();}
    write_checkpoint(work/"final.checkpoint",state,selections,abundance_tolerance,identities);
    std::ofstream report(work/"report.json");report<<std::setprecision(17)
      <<"{\"requested_age_reached\":"<<reached<<",\"age_years\":"<<state.model.age/year<<",\"accepted_this_invocation\":"<<accepted_here
      <<",\"accepted_total\":"<<state.accepted<<",\"rejected_total\":"<<state.rejected<<",\"stop_reason\":"<<std::quoted(stop)
      <<",\"cpu_seconds\":"<<double(std::clock()-cpu_start)/CLOCKS_PER_SEC
      <<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_start).count()
      <<",\"full_production_track_validated\":false,\"scope\":\"Common driver under integration; cool microscopic transport and atmosphere corridor unfinished\"}\n";
    std::cout<<std::setprecision(4)<<"Saved age "<<state.model.age/year<<" yr: "<<stop<<'\n';
    return reached || stop=="accepted-step budget" || stop=="CPU budget"?0:1;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
} // namespace ember::driver
