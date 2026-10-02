#include "ember/eos_variable_metal.hpp"
#include "ember/envelope_transport.hpp"
#include "ember/convective_material_heat.hpp"
#include "ember/convective_evolution_checks.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/atmosphere_metal_chain.hpp"
#include "ember/opacity_radiative.hpp"
#include "ember/conduction_table.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <map>
#include "ember/evolution_checkpoint.hpp"
#include <chrono>
using namespace ember;
int main(int argc,char**argv) {
 if(argc!=3 && argc!=4)return 2;
 try {
  const auto start=std::clock();const std::filesystem::path data=argv[1];std::ifstream input(data/"production/configuration.txt");std::map<std::string,std::string> cfg;std::string k,v;
  while(input>>k>>std::quoted(v))cfg[k]=v;
  for(const auto* role:{"eos","opacity_low","opacity_warm","opacity_bridge","opacity_hot","conduction","collisions"})
    cfg[role]=(data/"production"/cfg.at(role)).string();
  cfg["atmosphere_main_sequence"]=(data/"atmosphere/lifetime/main_sequence_reference.dat").string();
  cfg["atmosphere_metal_chain"]=(data/"atmosphere/lifetime/metal_chain_extended.dat").string();
  VariableMetalHelmholtzEos eos(cfg.at("eos"),HelmholtzTableEos::Mixture::allow_documented_proxy);
  RadiativeOpacity opacity({cfg.at("opacity_low"),cfg.at("opacity_warm"),cfg.at("opacity_bridge"),cfg.at("opacity_hot")},
      {data/"opacity/lifetime/hydrogen_response.dat",0,.16,1});
  CompositionAtmosphereGrid source(eos,cfg.at("atmosphere_main_sequence"),CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
  MetalAtmosphereChain atmosphere(eos,source,cfg.at("atmosphere_metal_chain"),.02);
  TabulatedConduction table(cfg.at("conduction"));HotConduction conduction(table);
  ScreenedCollisionTransport collisions(cfg.at("collisions"));ScreenedMetalMicroscopicTransport microscopic(eos,collisions,true,2e6,{true,true,true},true);
  driver::ConvectiveMaterialHeat material(eos,conduction);EnvelopeTransport transport(material,microscopic,2e6,3e6);
  PPCNNetwork nuclear(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,PPRates::solar_fusion_iii);PlasmaNeutrinoLosses losses;
  Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux};physics.microscopic=&transport;physics.neutrino_losses=&losses;
  Model original;std::size_t count;std::ifstream in(argv[2]);in>>original.age>>count;
  for(std::size_t i=0;i<count;++i) {
   double m,r,rho,T,L,h,y,z;CNAbundances cn;in>>m>>r>>rho>>T>>L>>h>>y>>cn[0]>>cn[1]>>cn[2]>>z;
   auto c=solar_scaled(h,z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
   c.X[1]=y;c.X[2]-=y;c.cn_molality=cn;c.cn_mass_convention=CNMassConvention::explicit_metal_mass;
   original.m.push_back(m);original.y.push_back({std::log(r),std::log(rho),std::log(T),L});original.comp.push_back(c);
  }
  original.M=original.m.back();std::vector<MetalSpeciesVector> rates(count-1);
  for(auto& rate:rates)for(auto& value:rate)in>>value;
  if(!in)throw std::runtime_error("invalid physical fixture");
  const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
  Model m=original;m.luminosity_grid=LuminosityGrid::volume_faces;
  for(std::size_t i=0;i+1<count;++i)m.y[i].L=.5*(original.y[i].L+original.y[i+1].L);
  transport.diagnostic_rates=rates;material.diagnostic_rates=rates;
  const auto regions=convective_mixing_regions(m,physics);
  require(regions.size()>1 && regions.size()<count,"fixture must contain radiative and mixed regions");
  // Finite abundance gradients require pressure inversions even inside the
  // cool convective envelope. Check actual nearby states without requiring
  // thermodynamic stability at an unrelated high-density table endpoint.
  std::size_t pressure_inversions=0;
  for(std::size_t i=0;i+1<count;++i)if(std::max(m.T(i),m.T(i+1))<2e6) {
    const double T=.5*(m.T(i)+m.T(i+1));
    const double P=.5*(eos.eval(m.T(i),m.rho(i),m.comp[i]).P+
        eos.eval(m.T(i+1),m.rho(i+1),m.comp[i+1]).P);
    const double density=eos.rho_from_PT(T,P,m.comp[i],.5*(m.rho(i)+m.rho(i+1)));
    require(std::abs(eos.eval(T,density,m.comp[i]).P/P-1)<2e-12,
            "cool-envelope density inversion must satisfy the pressure constraint");
    ++pressure_inversions;
  }
  require(pressure_inversions>0,"fixture must test cool-envelope pressure inversion");
  double maximum_flux=0;
  for(auto [a,b]:regions)if(b<count) {
    const auto f=transport.metal_eval(b-1,m.m[b-1],m.m[b],m.y[b-1],m.comp[b-1],m.y[b],m.comp[b],false);
    maximum_flux=std::max(maximum_flux,std::abs(f.species.rate[2]));
  }
  require(maximum_flux>1e8,"real screened species law must produce nonzero metal exchange");
  transport.diagnostic_rates={};material.diagnostic_rates={};
  EvolutionOptions options;options.relaxation.zone_threads=2;options.abundance_tolerance=1e-12;
  options.homogeneous_abundance_tolerance=1e-15;
  options.previous_metal_heat_rates=rates;
  const bool finite=argc==4;
  if(finite) {
    options.convective_mixing=ConvectiveMixing::finite_implicit;
    options.instantaneous_mixing_below_T=2e6;
    const auto blocks=instantaneous_mixing_regions(m,physics,options);
    require(blocks.size()>regions.size() && blocks.size()<count,
            "fixture must exercise finite convection and cool instantaneous blocks");
    auto all=options;all.instantaneous_mixing_below_T=1e99;
    require(instantaneous_mixing_regions(m,physics,all)==regions,
            "a temperature choice must never collapse a radiative face");
    auto unsupported=options;unsupported.instantaneous_mixing_below_T=0;
    const auto rejected=evolve_step(m,physics,atmosphere,31557600.,unsupported);
    require(!rejected.converged && rejected.model.comp==m.comp,
            "finite convection must not omit unsupported cool microscopic exchange");
  }
  const double dt=1e6*31557600.;const auto first=evolve_step(m,physics,atmosphere,dt,options);
  if(!first.converged)throw std::runtime_error(first.message);
  const auto audit=driver::check_interval(m,first,dt,nuclear,1e-14);
  require(audit.pass,"real screened evolution must pass the unchanged isotope and energy audits");
  if(finite) {
    transport.diagnostic_rates=first.total_metal_species_rates;
    material.diagnostic_rates=first.total_metal_species_rates;
    const auto guard=driver::check_envelope_transport(first.model,physics,transport,
        first.total_metal_species_rates,2e6,.01,options);
    require(std::abs(guard.convective_mass_fraction-first.convective_mass_fraction)<1e-14,
            "finite faces must retain physical convection diagnostics");
    const auto blocks=instantaneous_mixing_regions(first.model,physics,options);
    for(auto [a,b]:blocks)for(auto i=a;i<b;++i)
      require(first.model.comp[i]==first.model.comp[a],"selected instantaneous block is not homogeneous");
    transport.diagnostic_rates={};material.diagnostic_rates={};
    const auto half=evolve_step(m,physics,atmosphere,dt/2,options);
    require(half.converged,"finite first half must converge");
    auto next=options;next.previous_metal_heat_rates=half.total_metal_species_rates;
    const auto end=evolve_step(half.model,physics,atmosphere,dt/2,next);
    require(end.converged && driver::check_interval(m,half,dt/2,nuclear,1e-14).pass
        && driver::check_interval(half.model,end,dt/2,nuclear,1e-14).pass,
        "finite half steps must each conserve species and energy");
  }
  const auto preserved=first.model;
  auto wrong=first;
  auto abundance=metal_cn_abundances(wrong.model.comp[256]);abundance[0]+=1e-6;
  wrong.model.comp[256]=metal_cn_composition(wrong.model.comp[256],abundance);
  require(!driver::check_interval(m,wrong,dt,nuclear,1e-14).pass,
          "audit must reject an injected hydrogen inventory error");
  const auto stamp=std::chrono::high_resolution_clock::now().time_since_epoch().count();
  const auto checkpoint=std::filesystem::temp_directory_path()/("ember-screened-"+std::to_string(stamp)+".checkpoint");
  const driver::Selections selections{"volume faces","ppCN","plasma neutrinos","instantaneous","screened core"};
  const driver::Identities identities{{"executable",driver::file_identity(argv[0])},
      {"purpose","real screened transport integration test"}};
  driver::write_checkpoint(checkpoint,{first.model,dt,1,0,first.total_metal_species_rates},selections,1e-12,identities);
  const auto restored=driver::read_checkpoint(checkpoint,count,m.M,m.comp[0],selections,1e-12,identities,LuminosityGrid::volume_faces,true);
  std::filesystem::remove(checkpoint);
  options.previous_metal_heat_rates=first.total_metal_species_rates;
  const auto follow=evolve_step(first.model,physics,atmosphere,dt,options);
  options.previous_metal_heat_rates=restored.metal_heat_rates;
  const auto replay=evolve_step(restored.model,physics,atmosphere,dt,options);
  require(follow.converged && replay.converged,"screened transport restart must converge");
  require(follow.model.comp==replay.model.comp && follow.total_metal_species_rates==replay.total_metal_species_rates,
          "screened restart must reproduce species and carried heat exactly");
  for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<NVAR;++j)
    require(follow.model.y[i][static_cast<Var>(j)]==replay.model.y[i][static_cast<Var>(j)],
            "screened restart must reproduce structure exactly");
  require(driver::check_interval(first.model,follow,dt,nuclear,1e-14).pass,
          "screened continuation must retain its conservation accuracy");
  require(first.model.comp==preserved.comp,"fault injection altered the retained control");
  std::cout<<std::setprecision(4)<<"PASS real screened diffusion across "<<regions.size()-1
      <<" radiative faces; maximum metal rate "<<maximum_flux<<" g/s; inventory residual "
      <<audit.maximum_species_error<<"; exact restart; CPU "<<double(std::clock()-start)/CLOCKS_PER_SEC<<" s\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
