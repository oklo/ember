#include "lifetime_driver.hpp"
#include "lifetime_extrapolation.hpp"
#include "ember/convection.hpp"
#include "ember/contracting_seed.hpp"
#include "ember/convective_evolution_checks.hpp"
#include "ember/convective_material_heat.hpp"
#include "ember/envelope_transport.hpp"
#include "ember/envelope_atmosphere.hpp"
#include "ember/evolution_checkpoint.hpp"
#include "ember/runtime_identity.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include "ember/atmosphere_fixed_metal.hpp"
#include "ember/atmosphere_overlap.hpp"
#include "ember/atmosphere_metal_chain.hpp"
#include "ember/atmosphere_metal_interval.hpp"
#include "ember/atmosphere_hydrogen_interval.hpp"
#include "ember/atmosphere_hydrogen_envelope.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/eos_deuterium.hpp"
#include "ember/metal_microscopic_transport.hpp"
#include "ember/opacity_radiative.hpp"
#include "ember/opacity_conductive_interior.hpp"
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
      " [--restart CHECKPOINT] [--restart-step-years Y] [--max-steps N] [--cpu-seconds S]\n"
      "Common Hayashi-to-remnant driver under integration. Source-domain failures stop; no atmosphere fallback.\n";
    return 2;
  }
  try {
    const auto wall_start=std::chrono::steady_clock::now();const auto cpu_start=std::clock();
    const fs::path config=fs::canonical(argv[2]),work=argv[4];const double target=positive(argv[3])*year;
    if(fs::exists(work))throw std::invalid_argument("lifetime output directory already exists");
    std::string restart;std::size_t maximum_steps=500;double maximum_cpu=240,restart_step_years=0;
    for(int i=5;i<argc;i+=2) {
      if(i+1>=argc)throw std::invalid_argument("lifetime option needs a value");
      const std::string key=argv[i];
      if(key=="--restart" && restart.empty())restart=argv[i+1];
      else if(key=="--restart-step-years" && restart_step_years==0)restart_step_years=positive(argv[i+1]);
      else if(key=="--max-steps") {
        const double value=positive(argv[i+1]);if(value>100000 || std::floor(value)!=value)throw std::invalid_argument("invalid step budget");
        maximum_steps=static_cast<std::size_t>(value);
      }else if(key=="--cpu-seconds")maximum_cpu=positive(argv[i+1]);
      else throw std::invalid_argument("unknown lifetime argument: "+key);
    }
    if(restart_step_years>0 && restart.empty())
      throw std::invalid_argument("--restart-step-years requires --restart");
    Settings cfg(config);
    if(cfg.get("version")!="1")throw std::invalid_argument("unsupported lifetime configuration version");
    const auto path=[&](const char* key){return fs::canonical(config.parent_path()/cfg.get(key));};
    const auto eos_path=path("eos"),low_path=path("opacity_low"),warm_path=path("opacity_warm"),bridge_path=path("opacity_bridge"),hot_path=path("opacity_hot");
    const auto low_metal_interpolation=cfg.values.contains("eos_low_metal_interpolation")
        ?cfg.get("eos_low_metal_interpolation"):"cubic";
    if(low_metal_interpolation!="cubic" && low_metal_interpolation!="quadratic")
      throw std::invalid_argument("unknown low-metal EOS interpolation");
    const auto ion_quantum=cfg.values.contains("eos_ion_quantum")?cfg.get("eos_ion_quantum"):"none";
    if(ion_quantum!="none" && ion_quantum!="liquid_bc22")
      throw std::invalid_argument("unknown quantum-ion EOS selection");
    fs::path cold_eos_path;
    if(cfg.values.contains("eos_cold_potential"))cold_eos_path=path("eos_cold_potential");
    const auto conduction_path=path("conduction"),atmosphere_path=path("atmosphere"),collision_path=path("collisions"),composition_path=path("composition");
    const auto conduction_envelope=cfg.values.contains("conduction_envelope")
        ?cfg.get("conduction_envelope"):"temperature_join";
    if(conduction_envelope!="temperature_join" && conduction_envelope!="ionized")
      throw std::invalid_argument("unknown conduction envelope selection");
    fs::path cold_opacity_path;
    if(cfg.values.contains("opacity_cold_dense"))cold_opacity_path=path("opacity_cold_dense");
    const double cold_opacity_scale=cfg.values.contains("opacity_cold_dense_scale")
        ?cfg.number("opacity_cold_dense_scale"):1;
    if(!std::isfinite(cold_opacity_scale) || cold_opacity_scale<.1 || cold_opacity_scale>10
        || (cold_opacity_path.empty() && cold_opacity_scale!=1))
      throw std::invalid_argument("invalid cold dense-gas opacity selection");
    const double conductive_opacity=cfg.values.contains("opacity_conductive_interior")
        ?cfg.number("opacity_conductive_interior"):0;
    const double conductive_opacity_scale=cfg.values.contains("opacity_conductive_scale")
        ?cfg.number("opacity_conductive_scale"):1;
    if((conductive_opacity!=0 && conductive_opacity!=1) || !std::isfinite(conductive_opacity_scale)
        || conductive_opacity_scale<.1 || conductive_opacity_scale>10
        || (conductive_opacity==0 && conductive_opacity_scale!=1))
      throw std::invalid_argument("invalid conductive-interior radiative opacity selection");
    const double envelope_opacity=cfg.values.contains("opacity_conductive_envelope")
        ?cfg.number("opacity_conductive_envelope"):0;
    const double envelope_opacity_scale=cfg.values.contains("opacity_conductive_envelope_scale")
        ?cfg.number("opacity_conductive_envelope_scale"):1;
    const double envelope_opacity_uncertainty=cfg.values.contains("opacity_conductive_envelope_uncertainty")
        ?cfg.number("opacity_conductive_envelope_uncertainty"):.001;
    if((envelope_opacity!=0 && envelope_opacity!=1) || !std::isfinite(envelope_opacity_scale)
        || envelope_opacity_scale<.01 || envelope_opacity_scale>100 || (envelope_opacity==0 && envelope_opacity_scale!=1)
        || !std::isfinite(envelope_opacity_uncertainty) || envelope_opacity_uncertainty<=0 || envelope_opacity_uncertainty>.01
        || (envelope_opacity==0 && envelope_opacity_uncertainty!=.001))
      throw std::invalid_argument("invalid conductive-envelope opacity selection");
    std::optional<RadiativeOpacity::Extension> opacity_extension;
    if(cfg.values.contains("opacity_hydrogen_response")) {
      opacity_extension=RadiativeOpacity::Extension{path("opacity_hydrogen_response"),
        cfg.number("opacity_minimum_Z"),cfg.number("opacity_maximum_Z"),cfg.number("opacity_maximum_X")};
      if(cfg.values.contains("opacity_dense_hydrogen_maximum_logR"))
        opacity_extension->dense_hydrogen_maximum_logR=cfg.number("opacity_dense_hydrogen_maximum_logR");
      if(cfg.values.contains("opacity_dense_hydrogen_maximum_logT"))
        opacity_extension->dense_hydrogen_maximum_logT=cfg.number("opacity_dense_hydrogen_maximum_logT");
    }
    fs::path main_atmosphere_path;
    AtmosphereOverlap::Options atmosphere_join;
    if(cfg.values.contains("atmosphere_main_sequence")) {
      main_atmosphere_path=path("atmosphere_main_sequence");
      atmosphere_join={cfg.number("atmosphere_join_log_g_low"),cfg.number("atmosphere_join_log_g_high"),
        cfg.number("atmosphere_join_hydrogen_low"),cfg.number("atmosphere_join_hydrogen_high")};
    }
    const double mass=cfg.number("mass_Msun")*constants::Msun,radius=cfg.number("initial_radius_Rsun")*constants::Rsun;
    const double teff=cfg.number("initial_Teff_K"),entropy_loss=cfg.number("initial_entropy_loss");
    const double count=cfg.number("points"),threads=cfg.number("zone_threads"),minimum_temperature=cfg.number("screened_minimum_T_K");
    double dt=cfg.number("initial_step_years")*year;const double maximum_dt=cfg.number("maximum_step_years")*year;
    const double structure_tolerance=cfg.number("structure_tolerance"),species_tolerance=cfg.number("species_tolerance"),energy_tolerance=cfg.number("energy_tolerance");
    const double abundance_tolerance=cfg.number("coupling_abundance_tolerance");
    const double inventory_tolerance=cfg.number("inventory_abundance_tolerance");
    const auto structure_prediction=cfg.values.contains("structure_prediction")?cfg.get("structure_prediction"):"none";
    if(structure_prediction!="none" && structure_prediction!="linear")
      throw std::invalid_argument("unknown structure prediction method");
    const double collision_radius=cfg.values.contains("collision_taylor_radius")?cfg.number("collision_taylor_radius"):0;
    const double collision_verify=cfg.values.contains("collision_verify_reuse")?cfg.number("collision_verify_reuse"):0;
    if(!(collision_radius>=0 && collision_radius<=1e-4 && (collision_verify==0 || collision_verify==1))
        || (collision_verify!=0 && collision_radius==0))
      throw std::invalid_argument("invalid collision reuse settings");
    const auto optional_number=[&](const char* key) {return cfg.values.contains(key)?cfg.number(key):0.;};
    const double envelope_jacobian_radius=optional_number("envelope_jacobian_radius");
    if(!std::isfinite(envelope_jacobian_radius) || envelope_jacobian_radius<0 || envelope_jacobian_radius>.01)
      throw std::invalid_argument("envelope_jacobian_radius must lie in [0,0.01]");
    const double eos_radius=optional_number("eos_taylor_radius");
    const double buoyancy_spacing=optional_number("buoyancy_reuse_spacing");
    const double screening_spacing=optional_number("screening_reuse_spacing");
    const double quantum_screening=optional_number("quantum_screening_zeta_max");
    if(!std::isfinite(quantum_screening) || quantum_screening<0 || quantum_screening>1.6)
      throw std::invalid_argument("quantum_screening_zeta_max must lie in [0,1.6]");
    const double verify_responses=optional_number("verify_response_reuse");
    const double linearized_burning=optional_number("linearized_burning");
    if(linearized_burning!=0 && linearized_burning!=1)
      throw std::invalid_argument("linearized_burning must be zero or one");
    const double abundance_cap=cfg.values.contains("abundance_cap")?cfg.number("abundance_cap"):.001;
    if(!std::isfinite(abundance_cap) || abundance_cap<=0 || abundance_cap>.01)
      throw std::invalid_argument("abundance cap must lie in (0,0.01]");
    const double cap_after_mixing=optional_number("abundance_cap_after_mixing");
    if(cap_after_mixing!=0 && cap_after_mixing!=1)
      throw std::invalid_argument("abundance_cap_after_mixing must be 0 or 1");
    const double coupling_stop=optional_number("coupling_stop_tolerance");
    const double richardson=optional_number("richardson_extrapolation");
    if(richardson!=0 && richardson!=1)throw std::invalid_argument("richardson_extrapolation must be 0 or 1");
    const double material_heat=optional_number("material_heat_tolerance");
    const double verification_residual=optional_number("verification_residual_tolerance");
    const double verification_correction=optional_number("verification_correction_tolerance");
    for(double t:{coupling_stop,material_heat,verification_residual,verification_correction})
      if(!std::isfinite(t) || t<0)throw std::invalid_argument("invalid coupling or verification tolerance");
    if(coupling_stop>1e-3*species_tolerance || material_heat>1e-3*energy_tolerance
        || verification_residual>1e-3*structure_tolerance || verification_correction>1e-2*structure_tolerance)
      throw std::invalid_argument("coupling or verification tolerance must remain below time accuracy");
    for(double r:{eos_radius,buoyancy_spacing,screening_spacing})
      if(!std::isfinite(r) || r<0 || r>1e-4)throw std::invalid_argument("response reuse radius must lie in [0,1e-4]");
    if(verify_responses!=0 && verify_responses!=1)throw std::invalid_argument("verify_response_reuse must be zero or one");
    const auto transport_selection=cfg.values.contains("transport")?cfg.get("transport"):"whole_convective";
    if(transport_selection!="whole_convective" && transport_selection!="screened_core")
      throw std::invalid_argument("unknown lifetime transport selection");
    const bool screened_core=transport_selection=="screened_core";
    const auto mixing_selection=cfg.values.contains("convective_mixing")?cfg.get("convective_mixing"):"instantaneous";
    ConvectiveMixing mixing_mode=ConvectiveMixing::instantaneous;
    if(mixing_selection=="finite_implicit")mixing_mode=ConvectiveMixing::finite_implicit;
    else if(mixing_selection=="finite_lagged")mixing_mode=ConvectiveMixing::finite_lagged;
    else if(mixing_selection!="instantaneous")throw std::invalid_argument("unknown convective mixing selection");
    const double instantaneous_below=cfg.values.contains("instantaneous_mixing_below_T_K")
        ?cfg.number("instantaneous_mixing_below_T_K"):0;
    if(instantaneous_below<0 || (mixing_mode==ConvectiveMixing::instantaneous && instantaneous_below!=0)
        || (mixing_mode!=ConvectiveMixing::instantaneous && !screened_core))
      throw std::invalid_argument("finite convection requires screened-core transport and a nonnegative mixing temperature");
    const double heat_lower=screened_core && cfg.values.contains("screened_heat_lower_T_K")
        ?cfg.number("screened_heat_lower_T_K"):minimum_temperature;
    const double heat_upper=screened_core?cfg.number("screened_heat_upper_T_K"):minimum_temperature;
    const double mixing_gradient=screened_core?cfg.number("maximum_relative_mixing_gradient"):0;
    if(screened_core && (!(heat_upper>heat_lower && heat_lower>=minimum_temperature)
        || !(mixing_gradient>0 && mixing_gradient<=.01)))
      throw std::invalid_argument("invalid screened-core heat overlap or mixing approximation");
    const auto atmosphere_metals=cfg.values.contains("atmosphere_metals")?cfg.get("atmosphere_metals"):"strict";
    const double atmosphere_delta_Z=cfg.values.contains("atmosphere_maximum_delta_Z")
        ?cfg.number("atmosphere_maximum_delta_Z"):0;
    fs::path metal_atmosphere_path;
    double metal_join_low=0,metal_join_high=0;
    if(cfg.values.contains("atmosphere_metal_chain")) {
      metal_atmosphere_path=path("atmosphere_metal_chain");
      metal_join_low=cfg.number("atmosphere_metal_join_low");
      metal_join_high=cfg.number("atmosphere_metal_join_high");
      if(main_atmosphere_path.empty() || atmosphere_metals!="bounded_fixed_Z" ||
          !(metal_join_low>=0 && metal_join_high>metal_join_low && metal_join_high<1))
        throw std::invalid_argument("metal-dependent atmosphere requires a bounded reference and ordered overlap");
    }
    fs::path hydrogen_atmosphere_path;
    if(cfg.values.contains("atmosphere_hydrogen_interval"))hydrogen_atmosphere_path=path("atmosphere_hydrogen_interval");
    fs::path hydrogen_envelope_path;
    if(cfg.values.contains("atmosphere_hydrogen_envelope"))hydrogen_envelope_path=path("atmosphere_hydrogen_envelope");
    fs::path envelope_map_path;
    if(cfg.values.contains("envelope_map"))envelope_map_path=path("envelope_map");
    std::unique_ptr<EnvelopeMapAtmosphere> deep_envelope;
    if(!envelope_map_path.empty()) {
      deep_envelope=std::make_unique<EnvelopeMapAtmosphere>(envelope_map_path.string());
      if(deep_envelope->total_mass()!=mass)
        throw std::invalid_argument("envelope map requires a matching total stellar mass");
    }
    fs::path envelope_source_path;
    if(cfg.values.contains("envelope_source"))envelope_source_path=path("envelope_source");
    const auto envelope_eos=cfg.values.contains("envelope_eos")?cfg.get("envelope_eos"):"";
    if(!envelope_eos.empty() && (envelope_eos!="interior" || !envelope_source_path.empty()))
      throw std::invalid_argument("envelope EOS must be interior, without a separate layer source");
    const auto envelope_metals=cfg.values.contains("envelope_metals")?cfg.get("envelope_metals"):"reject";
    if((envelope_metals!="reject" && envelope_metals!="neutral" && envelope_metals!="ionized")
        || (envelope_metals!="reject" && envelope_source_path.empty()))
      throw std::invalid_argument("invalid envelope metal approximation");
    const bool direct_envelope=!envelope_source_path.empty() || envelope_eos=="interior";
    const double envelope_fraction=cfg.values.contains("envelope_mass_fraction")?cfg.number("envelope_mass_fraction"):0;
    const double envelope_thermal_limit=cfg.values.contains("envelope_thermal_fraction_limit")
        ?cfg.number("envelope_thermal_fraction_limit"):.001;
    if(!(std::isfinite(envelope_thermal_limit) && envelope_thermal_limit>0
        && envelope_thermal_limit<=.005 && (!direct_envelope || envelope_thermal_limit<=energy_tolerance))
        || (!direct_envelope && envelope_thermal_limit!=.001))
      throw std::invalid_argument("envelope thermal fraction limit requires a native envelope and must not exceed the energy tolerance or 0.005");
    if(direct_envelope && (deep_envelope
        || !(envelope_fraction>0 && envelope_fraction<=.01)))
      throw std::invalid_argument("native envelope requires mass fraction in (0,0.01] and no envelope map");
    if(!direct_envelope && envelope_fraction!=0)
      throw std::invalid_argument("envelope mass fraction requires a source or the interior EOS");
    const double selected_envelope_mass=deep_envelope?deep_envelope->envelope_mass():mass*envelope_fraction;
    if((atmosphere_metals!="strict" && atmosphere_metals!="bounded_fixed_Z") ||
        (atmosphere_metals=="strict"?atmosphere_delta_Z!=0:atmosphere_delta_Z<=0))
      throw std::invalid_argument("invalid atmosphere metal approximation selection");
    if(!cfg.values.empty())throw std::invalid_argument("unknown lifetime setting: "+cfg.values.begin()->first);
    if(!(mass>0 && radius>0 && teff>0 && entropy_loss>0 && dt>0 && maximum_dt>=dt && minimum_temperature>0
        && count>=128 && count<=8192 && std::floor(count)==count && threads>=1 && threads<=16 && std::floor(threads)==threads
        && structure_tolerance>0 && structure_tolerance<=1e-3 && species_tolerance>0 && species_tolerance<=1e-4
        && energy_tolerance>0 && energy_tolerance<=.01 && abundance_tolerance>0 && abundance_tolerance<=1e-12
        && inventory_tolerance>0 && inventory_tolerance<=std::min(1e-12,.01*species_tolerance)))
      throw std::invalid_argument("lifetime physical or accuracy setting out of range");
    const auto points=static_cast<std::size_t>(count);
    Composition initial;initial.basis=AbundanceBasis::baryon_mass;initial.metal_inventory=MetalInventory::gs98;
    std::ifstream input(composition_path);for(auto& x:initial.X)x=read_representable_double(input);
    std::string extra;if(input>>extra || std::abs(initial.sum()-1)>2e-12 || initial[Species::H2]<=0)
      throw std::invalid_argument("initial composition must be normalized and include primordial D");
    for(double x:initial.X)if(x<0)throw std::invalid_argument("negative initial abundance");
    initial.cn_molality=initial_gs98_cn(initial);initial=explicit_cn_material(initial);

    RuntimeIdentity identity;identity.file("executable",argv[0]);
    if(conduction_envelope=="ionized")
      identity.values["conduction.envelope"]="fully_ionized_source.v1";
    if(deep_envelope) {
      identity.file("atmosphere.envelope_map",envelope_map_path);
      identity.number("atmosphere.envelope_mass",deep_envelope->envelope_mass());
      identity.values["atmosphere.envelope_thermal"]="base_state_reservoir.v1";
    }
    if(direct_envelope) {
      if(envelope_metals!="reject")identity.values["atmosphere.envelope_metals"]=envelope_metals+".additive_volume.v1";
      if(!envelope_source_path.empty())identity.file("atmosphere.envelope_source",envelope_source_path);
      else identity.values["atmosphere.envelope_eos"]="interior.v1";
      identity.number("atmosphere.envelope_mass",selected_envelope_mass);
      identity.values["atmosphere.envelope_thermal"]="base_state_reservoir.v1";
      identity.number("atmosphere.envelope_lnP_steps",20);
    }
    identity.values["opacity.composition_extension"]=opacity_extension?"hydrogen_share.linear_Z.source_log_X.v1":"none";
    if(!cold_opacity_path.empty()) {
      identity.file("opacity.cold_dense.table",cold_opacity_path);
      identity.values["opacity.cold_dense"]="computed_baryonic_gas.C2_logR_T_trace_Z.v1";
      identity.number("opacity.cold_dense_scale",cold_opacity_scale);
    }
    if(envelope_opacity==1) {
      identity.values["opacity.conductive_envelope"]=conduction_envelope=="ionized"
          ?"hydrogen_density_continuation.logR_overlap.v2":"hydrogen_density_continuation.v6";
      identity.number("opacity.conductive_envelope_scale",envelope_opacity_scale);
      if(envelope_opacity_uncertainty!=.001)
        identity.number("opacity.conductive_envelope_uncertainty",envelope_opacity_uncertainty);
    }
    if(conductive_opacity==1) {
      identity.values["opacity.conductive_interior"]=conduction_envelope=="ionized"
          ?"source_slope.fixed_density_overlap.ionized_cold.v1":"source_slope.fixed_density_overlap.v4";
      identity.number("opacity.conductive_transport_uncertainty_limit",.001);
      identity.number("opacity.conductive_scale",conductive_opacity_scale);
    }
    if(opacity_extension) {
      identity.family("opacity_hydrogen_response",opacity_extension->hydrogen_response,false);
      identity.number("opacity.minimum_Z",opacity_extension->minimum_Z);
      identity.number("opacity.maximum_Z",opacity_extension->maximum_Z);
      identity.number("opacity.maximum_X",opacity_extension->maximum_X);
      if(opacity_extension->dense_hydrogen_maximum_logR!=1.8)
        identity.number("opacity.dense_hydrogen_maximum_logR",opacity_extension->dense_hydrogen_maximum_logR);
      if(opacity_extension->dense_hydrogen_maximum_logT!=6.1)
        identity.number("opacity.dense_hydrogen_maximum_logT",opacity_extension->dense_hydrogen_maximum_logT);
    }
    identity.values["transport.selection"]=transport_selection;
    if(mixing_mode!=ConvectiveMixing::instantaneous)
      identity.number("convection.instantaneous_mixing_below_T_K",instantaneous_below);
    if(screened_core) {
      if(heat_lower!=minimum_temperature)
        identity.number("transport.screened_heat_lower_T_K",heat_lower);
      identity.number("transport.screened_heat_upper_T_K",heat_upper);
      identity.number("transport.maximum_relative_mixing_gradient",mixing_gradient);
    }
    identity.values["atmosphere.metals"]=atmosphere_metals;
    identity.number("atmosphere.maximum_delta_Z",atmosphere_delta_Z);
    if(!metal_atmosphere_path.empty()) {
      identity.metal_atmosphere_chain("atmosphere.metal_chain",metal_atmosphere_path);
      identity.number("atmosphere.metal_join_low",metal_join_low);
      identity.number("atmosphere.metal_join_high",metal_join_high);
    }
    if(!hydrogen_atmosphere_path.empty())
      identity.hydrogen_atmosphere_interval("atmosphere.hydrogen_interval",hydrogen_atmosphere_path);
    if(!hydrogen_envelope_path.empty())
      identity.hydrogen_envelope("atmosphere.hydrogen_envelope",hydrogen_envelope_path);
    identity.values["atmosphere.overlap"]=main_atmosphere_path.empty()?"none":"gravity_hydrogen.v1";
    if(!main_atmosphere_path.empty()) {
      identity.file("atmosphere_main_sequence",main_atmosphere_path);
      identity.number("atmosphere.join.log_g_low",atmosphere_join.log_g_low);
      identity.number("atmosphere.join.log_g_high",atmosphere_join.log_g_high);
      identity.number("atmosphere.join.hydrogen_low",atmosphere_join.hydrogen_low);
      identity.number("atmosphere.join.hydrogen_high",atmosphere_join.hydrogen_high);
    }
    for(const auto& [key,value]:std::map<std::string,double>{{"mass_g",mass},{"initial_radius_cm",radius},
        {"initial_Teff_K",teff},{"initial_entropy_loss",entropy_loss},{"points",count},
        {"screened_minimum_T_K",minimum_temperature},{"structure_tolerance",structure_tolerance},
        {"species_tolerance",species_tolerance},{"energy_tolerance",energy_tolerance},
        {"coupling_abundance_tolerance",abundance_tolerance},{"inventory_abundance_tolerance",inventory_tolerance}})
      identity.number("configuration."+key,value);
    identity.number("solver.homogeneous_abundance_tolerance",std::min(1e-15,abundance_tolerance));
    if(envelope_thermal_limit!=.001)
      identity.number("envelope.thermal_fraction_limit",envelope_thermal_limit);
    if(structure_prediction=="linear")identity.number("solver.structure_prediction",1);
    if(linearized_burning==1)identity.number("solver.linearized_burning",1);
    if(abundance_cap!=.001)identity.number("solver.abundance_cap",abundance_cap);
    if(collision_radius>0)identity.number("solver.collision_taylor_radius",collision_radius);
    if(eos_radius>0)identity.number("solver.eos_taylor_radius",eos_radius);
    if(buoyancy_spacing>0)identity.number("solver.buoyancy_reuse_spacing",buoyancy_spacing);
    if(screening_spacing>0)identity.number("solver.screening_reuse_spacing",screening_spacing);
    if(quantum_screening>0) {
      identity.values["nuclear.quantum_screening"]="svh.cd09_finite_zeta.v1";
      identity.number("nuclear.quantum_screening_zeta_max",quantum_screening);
    }
    for(const auto& [key,value]:std::map<std::string,double>{{"coupling_stop_tolerance",coupling_stop},
        {"material_heat_tolerance",material_heat},{"verification_residual_tolerance",verification_residual},
        {"verification_correction_tolerance",verification_correction}})
      if(value>0)identity.number("solver."+key,value);
    identity.family("eos",eos_path,true);
    if(ion_quantum=="liquid_bc22")
      identity.values["eos.ion_quantum"]="bc22.liquid.common_ne.linear_mixture.full_expression.v7";
    if(low_metal_interpolation=="quadratic")
      identity.values["eos.low_metal_interpolation"]="quadratic.C2_to_cubic.v1";
    if(!cold_eos_path.empty()){
      identity.file("eos.cold_potential",cold_eos_path);
      identity.values["eos.cold_model"]=VariableMetalHelmholtzEos::cold_model_identifier;
    }
    for(const auto& [role,p]:std::map<std::string,fs::path>{{"opacity_low",low_path},{"opacity_warm",warm_path},
        {"opacity_bridge",bridge_path},{"opacity_hot",hot_path}})identity.family(role,p,false);
    for(const auto& [role,p]:std::map<std::string,fs::path>{{"conduction",conduction_path},{"atmosphere",atmosphere_path},
        {"collisions",collision_path},{"composition",composition_path}})identity.file(role,p);
    if(cap_after_mixing==1)
      identity.values["solver.abundance_cap_reference"]="instantaneously_mixed_previous.v1";
    if(richardson==1) {
      identity.values["integrator"]="richardson.full_two_half.assessed.v3";
      identity.number("integrator.residual_limit",1e-5);
      identity.number("integrator.energy_tolerance",.001);
    }
    const auto& identities=identity.values;
    const Selections selections{"lifetime.volume_faces.v1","ppcn.sfiii.svh.physical_metals.v1",
      "losses.plasma_neutrino.v1","convection."+mixing_selection+".material_heat_after_D.v1",
      screened_core?"transport.screened_core.material_envelope.v1":"transport.whole_convective_limit.v1"};
    Checkpoint state;
    // Reject incompatible or damaged restarts before allocating/parsing EOS
    // and opacity objects. The configuration's execution controls may change.
    if(!restart.empty()) {
      std::ifstream checkpoint_header(restart);std::string magic;int version{};
      if(!(checkpoint_header>>magic>>version) || magic!="EMBER_EVOLUTION_CHECKPOINT" || (version!=5 && version!=6 && version!=7))
        throw std::invalid_argument("lifetime restart requires a volume-face checkpoint");
      bool saved_heat=version==6;
      if(version==7) {
        int cn{},metal{},heat{};std::string grid;checkpoint_header>>cn>>metal>>grid>>heat;
        saved_heat=heat!=0;
      }
      state=read_checkpoint(restart,points,mass,initial,selections,abundance_tolerance,identities,
                            LuminosityGrid::volume_faces,saved_heat,selected_envelope_mass);
      if(state.model.age>=target)throw std::invalid_argument("target must exceed the saved age");
      if(restart_step_years>0)state.next_dt=std::min(restart_step_years,maximum_dt/year)*year;
    }
    VariableMetalHelmholtzEos table_eos(eos_path,HelmholtzTableEos::Mixture::allow_documented_proxy,
        low_metal_interpolation=="quadratic"?VariableMetalHelmholtzEos::LowMetalInterpolation::quadratic
                                            :VariableMetalHelmholtzEos::LowMetalInterpolation::cubic,
        cold_eos_path,ion_quantum=="liquid_bc22");
    DeuteriumApproxEos eos(table_eos);
    const RadiativeOpacity::Tables opacity_tables{low_path,warm_path,bridge_path,hot_path,
        cold_opacity_path,cold_opacity_scale};
    auto source_radiation=opacity_extension?std::make_shared<RadiativeOpacity>(opacity_tables,*opacity_extension)
        :std::make_shared<RadiativeOpacity>(opacity_tables);
    auto table_conduction=std::make_shared<TabulatedConduction>(conduction_path);
    std::shared_ptr<Conduction> conduction=table_conduction;
    if(conduction_envelope=="temperature_join")
      conduction=std::make_shared<HotConduction>(*table_conduction);
    std::shared_ptr<ConductiveInteriorOpacity> continued_radiation;
    std::shared_ptr<Opacity> radiation=source_radiation;
    if(conductive_opacity==1) {
      auto domain=ConductiveInteriorOpacity::Domain{};
      if(conduction_envelope=="ionized")domain.minimum_T=6e5;
      continued_radiation=std::make_shared<ConductiveInteriorOpacity>(*source_radiation,*conduction,.001,conductive_opacity_scale,domain);
      radiation=continued_radiation;
    }
    std::shared_ptr<ConductiveInteriorOpacity> envelope_radiation;
    if(envelope_opacity==1) {
      // Complete the density join before the hot hydrogen source reaches
      // log R = 3.5 at its lowest selected temperature, log T = 5.6.
      // That corner lies at rho = 199.5 g/cm^3.
      // Recover the supported radiative source across the conduction turn-on.
      const auto domain=conduction_envelope=="ionized"
          ?ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain()
          :ConductiveInteriorOpacity::hydrogen_envelope_domain();
      envelope_radiation=std::make_shared<ConductiveInteriorOpacity>(*radiation,*conduction,envelope_opacity_uncertainty,envelope_opacity_scale,domain);
      radiation=envelope_radiation;
    }
    CombinedOpacity combined(radiation,conduction);
    const auto selected_envelope_metals=envelope_metals=="neutral"?EnvelopeMetals::neutral
        :(envelope_metals=="ionized"?EnvelopeMetals::ionized:EnvelopeMetals::reject);
    std::unique_ptr<EnvelopeSource> envelope_source;
    std::unique_ptr<EnvelopeSourceDensity> envelope_density;
    if(!envelope_source_path.empty()) {
      envelope_source=std::make_unique<EnvelopeSource>(envelope_source_path.string());
      envelope_density=std::make_unique<EnvelopeSourceDensity>(*envelope_source,selected_envelope_metals);
    }
    const PressureDensity& atmosphere_density=envelope_density
        ?static_cast<const PressureDensity&>(*envelope_density):eos;
    CompositionAtmosphereGrid table_atmosphere(atmosphere_density,atmosphere_path,CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    std::unique_ptr<FixedMetalAtmosphere> fixed_metal_atmosphere;
    if(atmosphere_metals=="bounded_fixed_Z")fixed_metal_atmosphere=std::make_unique<FixedMetalAtmosphere>(
        atmosphere_density,table_atmosphere,atmosphere_delta_Z);
    const Atmosphere& contraction_atmosphere=fixed_metal_atmosphere
        ?static_cast<const Atmosphere&>(*fixed_metal_atmosphere):table_atmosphere;
    std::unique_ptr<CompositionAtmosphereGrid> main_atmosphere;
    std::unique_ptr<FixedMetalAtmosphere> fixed_main_atmosphere;
    std::unique_ptr<AtmosphereOverlap> atmosphere_overlap;
    if(!main_atmosphere_path.empty()) {
      main_atmosphere=std::make_unique<CompositionAtmosphereGrid>(atmosphere_density,main_atmosphere_path,
          CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
      if(atmosphere_metals=="bounded_fixed_Z")fixed_main_atmosphere=std::make_unique<FixedMetalAtmosphere>(
          atmosphere_density,*main_atmosphere,atmosphere_delta_Z);
      const Atmosphere& main=fixed_main_atmosphere?static_cast<const Atmosphere&>(*fixed_main_atmosphere):*main_atmosphere;
      atmosphere_overlap=std::make_unique<AtmosphereOverlap>(atmosphere_density,contraction_atmosphere,main,atmosphere_join);
    }
    const Atmosphere& reference_boundary=atmosphere_overlap
        ?static_cast<const Atmosphere&>(*atmosphere_overlap):contraction_atmosphere;
    std::unique_ptr<MetalAtmosphereChain> metal_atmosphere;
    std::unique_ptr<MetalIntervalAtmosphere> metal_overlap;
    if(!metal_atmosphere_path.empty()) {
      const double z=main_atmosphere->reference_metallicity();
      if(metal_join_high>z || metal_join_low<z-atmosphere_delta_Z)
        throw std::invalid_argument("metal atmosphere overlap exceeds the fixed-Z reference allowance");
      metal_atmosphere=std::make_unique<MetalAtmosphereChain>(atmosphere_density,*main_atmosphere,metal_atmosphere_path,z);
      metal_overlap=std::make_unique<MetalIntervalAtmosphere>(atmosphere_density,*metal_atmosphere,
          reference_boundary,metal_join_low,metal_join_high);
    }
    const Atmosphere& metal_boundary=metal_overlap
        ?static_cast<const Atmosphere&>(*metal_overlap):reference_boundary;
    std::unique_ptr<HydrogenIntervalAtmosphere> hydrogen_atmosphere;
    if(!hydrogen_atmosphere_path.empty())hydrogen_atmosphere=std::make_unique<HydrogenIntervalAtmosphere>(
        atmosphere_density,metal_boundary,hydrogen_atmosphere_path);
    const Atmosphere& hydrogen_boundary=hydrogen_atmosphere
        ?static_cast<const Atmosphere&>(*hydrogen_atmosphere):metal_boundary;
    std::unique_ptr<HydrogenEnvelopeAtmosphere> hydrogen_envelope;
    if(!hydrogen_envelope_path.empty())hydrogen_envelope=std::make_unique<HydrogenEnvelopeAtmosphere>(
        atmosphere_density,hydrogen_boundary,hydrogen_envelope_path);
    TraceDeuteriumAtmosphere thin_atmosphere(atmosphere_density,hydrogen_envelope
        ?static_cast<const Atmosphere&>(*hydrogen_envelope):hydrogen_boundary);
    std::unique_ptr<EnvelopeAtmosphere> native_envelope;
    if(envelope_source) {
      native_envelope=std::make_unique<EnvelopeAtmosphere>(thin_atmosphere,combined,*envelope_source,
          1.9,mass,selected_envelope_mass,20,selected_envelope_metals);
    }
    if(envelope_eos=="interior")native_envelope=std::make_unique<EnvelopeAtmosphere>(
        thin_atmosphere,combined,eos,1.9,mass,selected_envelope_mass,20);
    if(native_envelope) {
      native_envelope->evaluation_threads(static_cast<std::size_t>(threads));
      native_envelope->jacobian_reuse(envelope_jacobian_radius);
    }
    const Atmosphere& atmosphere=native_envelope?static_cast<const Atmosphere&>(*native_envelope)
        :(deep_envelope?static_cast<const Atmosphere&>(*deep_envelope):thin_atmosphere);
    PPCNNetwork nuclear(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,PPRates::solar_fusion_iii);
    PlasmaNeutrinoLosses losses;ScreenedCollisionTransport collisions(collision_path.string());
    ScreenedMetalMicroscopicTransport microscopic(table_eos,collisions,true,minimum_temperature,{true,true,true},true);
    set_composition_buoyancy_reuse(buoyancy_spacing,verify_responses==1);
    set_screening_reuse(screening_spacing);
    set_quantum_screening(quantum_screening);
    if(eos_radius>0)microscopic.use_eos_taylor(eos_radius,verify_responses==1);
    std::shared_ptr<CollisionTaylorCache> collision_reuse;
    if(collision_radius>0) {
      collision_reuse=std::make_shared<CollisionTaylorCache>(collision_radius,points,collision_verify==1);
      microscopic.use_collision_taylor(collision_reuse);
    }
    ConvectiveMaterialHeat convective_heat(table_eos,*conduction);
    EnvelopeTransport envelope_heat(convective_heat,microscopic,heat_lower,
        screened_core?heat_upper:1.5*minimum_temperature);
    Physics early{&eos,&combined,&nuclear,1.9,ConvectiveCriterion::ledoux};early.neutrino_losses=&losses;early.explicit_metal_mixing_only=true;
    auto later=early;later.opacity=radiation.get();later.microscopic=&convective_heat;later.explicit_metal_mixing_only=false;
    if(screened_core)later.microscopic=&envelope_heat;
    EvolutionOptions options;options.relaxation.zone_threads=static_cast<std::size_t>(threads);options.abundance_tolerance=abundance_tolerance;
    options.linearized_burning=linearized_burning==1;
    options.max_abundance_change=abundance_cap;
    options.abundance_cap_after_mixing=cap_after_mixing==1;
    options.homogeneous_abundance_tolerance=std::min(1e-15,abundance_tolerance);
    options.coupling_stop_tolerance=coupling_stop;
    if(material_heat>0)options.material_heat_tolerance=material_heat;
    options.verification_residual_tolerance=verification_residual;
    options.verification_correction_tolerance=verification_correction;
    if(restart.empty()) {
      ContractingSource seed_source(nuclear,entropy_loss);auto seed_physics=early;seed_physics.nuclear=&seed_source;
      double seed_radius=radius,seed_teff=teff;
      if(native_envelope) {
        const auto surface=native_envelope->from_photosphere(teff,radius,initial);
        seed_radius=surface.r_base;seed_teff=teff*std::sqrt(radius/seed_radius);
      }
      const auto guess=contracting_guess(points,mass,seed_radius,seed_teff,initial,seed_physics,atmosphere,selected_envelope_mass);
      auto solved=relax(guess,seed_physics,atmosphere,options.relaxation);
      if(!solved.converged)throw std::runtime_error("Hayashi initial model: "+solved.message);
      for(const auto& c:solved.model.comp)if(c!=initial)throw std::runtime_error("initial relaxation changed isotope inventory");
      state={std::move(solved.model),dt,0,0};
    }
    if(state.model.age>=target)throw std::invalid_argument("target must exceed the saved age");
    fs::create_directories(work);write_checkpoint(work/"seed.checkpoint",state,selections,abundance_tolerance,identities);
    std::ofstream execution(work/"execution.json");execution<<std::setprecision(17)
      <<"{\"configuration_path\":"<<std::quoted(config.string())<<",\"configuration_fnv1a\":"<<std::quoted(file_identity(config))
      <<",\"restart_path\":"<<std::quoted(restart)<<",\"zone_threads\":"<<threads
      <<",\"restart_step_years\":"<<restart_step_years
      <<",\"effective_initial_step_years\":"<<state.next_dt/year
      <<",\"initial_step_years\":"<<dt/year<<",\"maximum_step_years\":"<<maximum_dt/year
      <<",\"maximum_steps\":"<<maximum_steps<<",\"maximum_cpu_seconds\":"<<maximum_cpu<<"}\n";
    std::ofstream history(work/"history.jsonl"),attempts(work/"attempts.jsonl");history<<std::setprecision(17);attempts<<std::setprecision(17);
    const auto record=[&](double step,double error,const HomogeneousCheck& guard) {
      const auto& m=state.model;const auto w=nodal_mass_weights(m);double H=0,Y3=0,D=0,Lnuc=0;
      std::vector<double> eps(m.size());std::size_t burning_peak=0,he3_peak=0;
      for(std::size_t i=0;i<m.size();++i) {
        H+=w[i]*m.comp[i][Species::H1];Y3+=w[i]*m.comp[i][Species::He3];D+=w[i]*m.comp[i][Species::H2];
        eps[i]=nuclear.eval(m.T(i),m.rho(i),m.comp[i]).eps;Lnuc+=w[i]*eps[i];
        if(eps[i]>eps[burning_peak])burning_peak=i;
        if(m.comp[i][Species::He3]>m.comp[he3_peak][Species::He3])he3_peak=i;
      }
      // Contiguous cells above half the peak specific nuclear power. This
      // describes the resolved burning region; it is not a convergence test.
      std::size_t left=burning_peak,right=burning_peak;
      double inner_q=0,outer_q=0;std::size_t burning_cells=0;
      if(eps[burning_peak]>0) {
        const double half=.5*eps[burning_peak];
        while(left>0 && eps[left-1]>=half)--left;
        while(right+1<m.size() && eps[right+1]>=half)++right;
        inner_q=left==0?0:.5*(m.m[left-1]+m.m[left])/m.M;
        outer_q=right+1==m.size()?1:.5*(m.m[right]+m.m[right+1])/m.M;
        burning_cells=right-left+1;
      }
      double R=m.r(m.size()-1);const double L=m.y.back().L;
      double Teff=std::pow(L/(4*M_PI*constants::sigma_SB*R*R),.25);
      if(deep_envelope) {
        const auto surface=deep_envelope->photosphere_R_Teff(Teff,constants::G*m.M/(R*R),m.comp.back());
        R=surface[0];Teff=surface[1];
      } else if(native_envelope) {
        const auto surface=native_envelope->photosphere(Teff,constants::G*m.M/(R*R),m.comp.back());
        R=surface.R;Teff=surface.Teff;
      }
      history<<"{\"years\":"<<m.age/year<<",\"step_years\":"<<step/year<<",\"accepted\":"<<state.accepted<<",\"rejected\":"<<state.rejected
        <<",\"radius_Rsun\":"<<R/constants::Rsun<<",\"Teff_K\":"<<Teff
        <<",\"luminosity_Lsun\":"<<L/constants::Lsun<<",\"nuclear_fraction\":"<<Lnuc/L<<",\"central_T_K\":"<<m.T(0)
        <<",\"central_density_g_cm3\":"<<m.rho(0)<<",\"central_X\":"<<m.comp.front()[Species::H1]
        <<",\"central_He3\":"<<m.comp.front()[Species::He3]
        <<",\"eps_nuc_peak_erg_g_s\":"<<eps[burning_peak]
        <<",\"eps_nuc_peak_q\":"<<(eps[burning_peak]>0?m.m[burning_peak]/m.M:0)
        <<",\"burning_half_max_inner_q\":"<<inner_q<<",\"burning_half_max_outer_q\":"<<outer_q
        <<",\"burning_half_max_cells\":"<<burning_cells
        <<",\"He3_peak\":"<<m.comp[he3_peak][Species::He3]
        <<",\"He3_peak_q\":"<<(m.comp[he3_peak][Species::He3]>0?m.m[he3_peak]/m.M:0)
        <<",\"H_mass_g\":"<<H<<",\"He3_mass_g\":"<<Y3<<",\"D_mass_g\":"<<D<<",\"surface_X\":"<<m.comp.back().X[0]
        <<",\"surface_Z\":"<<m.comp.back().Z()<<",\"error_norm\":"<<error<<",\"initial_D_approximation\":"<<has_D(m)
        <<",\"omitted_mixing_heat_fraction\":"<<(has_D(m)?guard.maximum_heat_fraction:0)<<",\"gross_mixing_heat_fraction\":"<<guard.maximum_gross_heat_fraction
        <<",\"burn_gradient_estimate\":"<<guard.burn_gradient_estimate<<",\"drift_gradient_proxy\":"<<guard.drift_gradient_proxy
        <<",\"kinetic_heat_proxy\":"<<guard.kinetic_heat_proxy
        <<",\"convective_mass_fraction\":"<<guard.convective_mass_fraction
        <<",\"radiative_boundaries\":"<<guard.radiative_boundaries
        <<",\"relative_mixing_gradient\":"<<guard.maximum_relative_mixing_gradient
        <<",\"convective_travel_years\":"<<guard.travel_years<<",\"cpu_seconds\":"<<double(std::clock()-cpu_start)/CLOCKS_PER_SEC<<"}\n";history.flush();
    };
    EvolutionControlOptions control;
    control.step=options;control.step.convective_mixing=ConvectiveMixing::instantaneous;
    control.target_age=target;control.maximum_dt=maximum_dt;
    control.structure_tolerance=structure_tolerance;control.species_tolerance=species_tolerance;
    control.energy_tolerance=energy_tolerance;control.maximum_steps=maximum_steps;
    control.predict_structure=structure_prediction=="linear";
    control.richardson_extrapolation=richardson==1;
    control.maximum_cpu_seconds=maximum_cpu;
    // A failed trial must not replace the accepted model. Retry with a shorter
    // interval under the same audits; the controller bounds repeated rejection.
    control.audit_failure_is_fatal=false;
    HomogeneousCheck guard;
    double envelope_maximum_thermal_fraction=0,envelope_maximum_nuclear_fraction=0,envelope_maximum_base_T=0;
    EvolutionControlHooks hooks;
    hooks.physics=[&](const Model& m)->const Physics& {return has_D(m)?early:later;};
    hooks.configure_step=[&](const Model& m,EvolutionOptions& selected) {
      selected.convective_mixing=has_D(m)?ConvectiveMixing::instantaneous:mixing_mode;
      selected.instantaneous_mixing_below_T=has_D(m)?0:instantaneous_below;
    };
    hooks.species_difference=[](const Composition& x,const Composition& y) {
      const auto a=metal_cn_abundances(x),b=metal_cn_abundances(y);double difference=0;
      for(std::size_t k=0;k<METAL_CN_SIZE;++k)difference=std::max(difference,std::abs(a[k]-b[k]));
      return difference;
    };
    hooks.audit=[&](const Model& old,const EvolutionStep& step,double duration) {
      auto audit=check_interval(old,step,duration,nuclear,inventory_tolerance);
      if(step.converged && native_envelope && selected_envelope_mass>0) {
        const auto& m=step.model;const auto i=m.size()-1;
        const auto a=eos.eval(old.T(i),old.rho(i),old.comp[i]);
        const auto b=eos.eval(m.T(i),m.rho(i),m.comp[i]);
        const double light=std::min(old.y.back().L,m.y.back().L);
        // Conservative base-reservoir scale: absolute internal-energy change
        // plus pressure work, without cancellation. It is an approximation
        // check, separate from the full-star conservation audit.
        const double power=selected_envelope_mass*(std::abs(b.E-a.E)
          +std::max(a.P,b.P)*std::abs(1/m.rho(i)-1/old.rho(i)))/duration;
        const double fraction=power/light;
        if(!std::isfinite(fraction) || light<=0 || fraction>envelope_thermal_limit) {
          std::ostringstream message;
          message<<std::setprecision(4)<<"native envelope thermal reservoir fraction "<<fraction
            <<" exceeds assessed limit "<<envelope_thermal_limit;
          throw std::domain_error(message.str());
        }
        envelope_maximum_thermal_fraction=std::max(envelope_maximum_thermal_fraction,fraction);
      }
      return audit;
    };
    hooks.assess=[&](const Model& m,std::span<const std::array<double,3>> rates) {
      if(selected_envelope_mass>0) {
        const double minimum_envelope_T=native_envelope?1e4:(cold_eos_path.empty()?2e5:3e5);
        const double maximum_envelope_T=native_envelope?2e6:1e6;
        if(m.T(m.size()-1)<minimum_envelope_T || m.T(m.size()-1)>maximum_envelope_T)
          throw std::domain_error("envelope base outside assessed thermal regime");
        envelope_maximum_base_T=std::max(envelope_maximum_base_T,m.T(m.size()-1));
        if(native_envelope && m.T(m.size()-1)>1e6) {
          const double fraction=selected_envelope_mass*std::abs(nuclear.eval(m.T(m.size()-1),m.rho(m.size()-1),m.comp.back()).eps)/m.y.back().L;
          if(!std::isfinite(fraction) || fraction>.0001)
            throw std::domain_error("native envelope nuclear reservoir exceeds assessed heating fraction");
          envelope_maximum_nuclear_fraction=std::max(envelope_maximum_nuclear_fraction,fraction);
        }
        const auto regions=convective_mixing_regions(m,has_D(m)?early:later,options.relaxation.zone_threads);
        if(regions.back().second-regions.back().first<2)
          throw std::domain_error("envelope reservoir requires a connected convective base");
      }
      const bool initial_D=has_D(m);convective_heat.diagnostic_rates=rates;envelope_heat.diagnostic_rates=rates;
      try {
        auto selected=options;hooks.configure_step(m,selected);
        guard=screened_core && !initial_D && !rates.empty()
          ?check_envelope_transport(m,later,envelope_heat,rates,minimum_temperature,mixing_gradient,selected)
          :check_initial_convection(m,initial_D?early:later,table_eos,nuclear,initial_D);
        convective_heat.diagnostic_rates={};envelope_heat.diagnostic_rates={};
      }catch(...) {convective_heat.diagnostic_rates={};envelope_heat.diagnostic_rates={};throw;}
    };
    std::ofstream extrapolation_log;
    if(richardson==1)extrapolation_log.open(work/"richardson.jsonl");
    const LifetimeExtrapolation assessment{hooks,options,atmosphere,nuclear,convective_heat,envelope_heat,
        static_cast<std::size_t>(threads),1e-5,.001,species_tolerance,extrapolation_log};
    if(richardson==1)hooks.assess_extrapolated=assessment;
    hooks.cpu_seconds=[&]{return double(std::clock()-cpu_start)/CLOCKS_PER_SEC;};
    hooks.terminal_failure=[](std::string_view reason) {
      return reason.find("radiative species boundary")!=std::string_view::npos;
    };
    hooks.attempted=[&](const EvolutionAttempt& attempt) {
      const auto& af=attempt.audits[0];const auto& a1=attempt.audits[1];const auto& a2=attempt.audits[2];
      attempts<<"{\"start_years\":"<<attempt.start_age/year<<",\"step_years\":"<<attempt.dt/year<<",\"converged\":"<<attempt.converged<<",\"audit_pass\":"<<attempt.audit_pass
        <<",\"accepted\":"<<attempt.accepted<<",\"error_norm\":"<<attempt.error_norm<<",\"message\":"<<std::quoted(attempt.message)
        <<",\"species_error\":["<<af.maximum_species_error<<','<<a1.maximum_species_error<<','<<a2.maximum_species_error
        <<"],\"mass_error_surface\":["<<af.mass_error_surface<<','<<a1.mass_error_surface<<','<<a2.mass_error_surface
        <<"],\"first_law\":["<<af.first_law<<','<<a1.first_law<<','<<a2.first_law<<"]}\n";attempts.flush();
    };
    hooks.accepted=[&](const EvolutionState&,double duration,double error) {
      record(duration,error,guard);
      if(duration>0) {
        write_checkpoint(work/"latest.checkpoint",state,selections,abundance_tolerance,identities);
        std::cerr<<std::setprecision(4)<<"age="<<state.model.age/year<<" yr; step error="<<error<<"; accepted="<<state.accepted<<'\n';
      }
    };
    const auto outcome=ember::evolve(state,atmosphere,control,hooks);
    {
      std::ofstream out(work/"envelope_approximation.json");
      out<<std::setprecision(17)<<"{\"maximum_base_T_K\":"<<envelope_maximum_base_T
        <<",\"thermal_fraction_limit\":"<<envelope_thermal_limit
        <<",\"maximum_thermal_fraction\":"<<envelope_maximum_thermal_fraction
        <<",\"maximum_nuclear_fraction\":"<<envelope_maximum_nuclear_fraction<<"}\n";
    }
    if(native_envelope) {
      std::ofstream diagnostic(work/"envelope_evaluations.json");diagnostic<<std::setprecision(17)
        <<"{\"jacobian_radius\":"<<envelope_jacobian_radius
        <<",\"integrations\":"<<native_envelope->integrations()
        <<",\"jacobians_computed\":"<<native_envelope->jacobian_computed()
        <<",\"jacobians_reused\":"<<native_envelope->jacobian_reused()<<"}\n";
    }
    if(envelope_radiation) {
      std::ofstream diagnostic(work/"conductive_envelope_opacity.json");diagnostic<<std::setprecision(17)
        <<"{\"evaluations\":"<<envelope_radiation->continued_evaluations()
        <<",\"maximum_transport_uncertainty\":"<<envelope_radiation->maximum_transport_uncertainty()
        <<",\"opacity_scale\":"<<envelope_opacity_scale<<"}\n";
    }
    if(continued_radiation) {
      std::ofstream diagnostic(work/"conductive_opacity.json");diagnostic<<std::setprecision(17)
        <<"{\"evaluations\":"<<continued_radiation->continued_evaluations()
        <<",\"maximum_relative_transport_uncertainty\":"
        <<continued_radiation->maximum_transport_uncertainty()<<",\"opacity_scale\":"<<conductive_opacity_scale<<"}\n";
    }
    if(richardson==1)std::cerr<<"Richardson accepted "<<outcome.richardson_accepted
      <<", declined "<<outcome.richardson_declined<<'\n';
    if(eos_radius>0 || buoyancy_spacing>0 || screening_spacing>0) {
      const auto e=microscopic.eos_reuse_statistics();const auto b=composition_buoyancy_reuse_statistics();
      std::ofstream out(work/"response_reuse.json");out<<std::setprecision(17)
        <<"{\"eos_radius\":"<<eos_radius<<",\"eos_hits\":"<<e.hits<<",\"eos_exact\":"<<e.exact
        <<",\"eos_verified\":"<<e.verified<<",\"worst_potential\":"<<e.worst_potential<<",\"worst_enthalpy\":"<<e.worst_enthalpy
        <<",\"buoyancy_spacing\":"<<buoyancy_spacing<<",\"buoyancy_hits\":"<<b.hits<<",\"buoyancy_exact\":"<<b.exact
        <<",\"buoyancy_verified\":"<<b.verified<<",\"worst_absolute_B_error\":"<<b.worst_absolute_B_error
        <<",\"screening_spacing\":"<<screening_spacing<<"}\n";
    }
    if(collision_reuse) {
      const auto stats=collision_reuse->statistics();
      std::ofstream reuse_report(work/"collision_reuse.json");reuse_report<<std::setprecision(17)
        <<"{\"radius\":"<<collision_radius<<",\"hits\":"<<stats.hits<<",\"exact\":"<<stats.misses
        <<",\"verified\":"<<stats.verified<<",\"maximum_coefficient_group_relative_error\":"<<stats.worst_relative_error<<"}\n";
    }
    const auto& stop=outcome.stop_reason;const bool reached=outcome.requested_age_reached;
    const auto accepted_here=outcome.accepted_this_invocation;
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
