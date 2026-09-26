// Sequential conditional whole-star steps; retain ages and log-state checkpoints.
// Neither the cooler opacity sampling nor the omitted cool kinetic heat
// is accepted by this experiment. No output is a selected trajectory point.
#include "conditional_envelope_heat.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/opacity_blend.hpp"
#include "ember/conduction_table.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>

using namespace ember;
namespace {
std::string scalar(double x) {
  if (!std::isfinite(x)) return "null";
  std::ostringstream s; s << std::setprecision(17) << x; return s.str();
}
struct ScaledOpacity final : Opacity {
  const Opacity& source; double factor;
  ScaledOpacity(const Opacity& s, double f):source(s),factor(f) {
    if (!(f>0) || !std::isfinite(f)) throw std::invalid_argument("invalid opacity factor");
  }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    auto s=source.eval(T,rho,c); s.kappa*=factor; return s;
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    return source.density_range(T,c);
  }
  const char* name() const override {return "conditional coherent bridge-opacity perturbation";}
};
double teff(const Model& m) {
  return std::pow(m.y.back().L/(4*std::numbers::pi*constants::sigma_SB*
    m.r(m.y.size()-1)*m.r(m.y.size()-1)),.25);
}
bool identical(const Model& a,const Model& b) {
  if (a.m!=b.m || a.age!=b.age || a.M!=b.M) return false;
  for (std::size_t i=0;i<a.y.size();++i)
    if (a.comp[i].X!=b.comp[i].X || a.y[i].lnr!=b.y[i].lnr ||
        a.y[i].lnrho!=b.y[i].lnrho || a.y[i].lnT!=b.y[i].lnT || a.y[i].L!=b.y[i].L) return false;
  return true;
}
}
int main(int argc,char** argv) {
  if (argc!=9 && argc!=10) return 2;
  const bool select_radiation=argc==10 && std::string(argv[9])=="select-radiation";
  const bool select_cn=argc==10 && std::string(argv[9])=="select-cn";
  const bool default_cn=argc==10 && std::string(argv[9])=="cn";
  const bool default_radiation=select_cn || default_cn || (argc==10 && std::string(argv[9])=="radiation");
  if(argc==10 && !select_radiation && !default_radiation)return 2;
  try {
    const auto begin=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    MixtureOpacity molecular(argv[2]),warm(argv[3]),bridge(argv[4]),hot(argv[5]);
    TabulatedConduction conduction(argv[6]);
    CompositionAtmosphereGrid atmosphere(eos,argv[7],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
    PPCNO cn,nitrogen(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn,PPRates::solar_fusion_iii,0.);
    PlasmaNeutrinoLosses losses;
    ScreenedCollisionTransport collisions(argv[8]);
    std::cerr << "load_seconds " << std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count() << '\n';
    std::cout << std::setprecision(17); std::string line;
    while (std::getline(std::cin,line)) {
      try {
        std::istringstream input(line); double factor,dt,abundance_tolerance,heat_tolerance,age; std::size_t n; int ions,logarithmic;
        int parcel_radiation=default_radiation;
        int cn_carbon=default_cn?1:-1;
        if(select_cn && (!(input>>cn_carbon) || cn_carbon < -1 || cn_carbon > 1))return 2;
        const Nuclear& nuclear=cn_carbon<0?static_cast<const Nuclear&>(pp):
          cn_carbon==0?static_cast<const Nuclear&>(nitrogen):static_cast<const Nuclear&>(cn);
        if(select_radiation && (!(input>>parcel_radiation) || (parcel_radiation!=0 && parcel_radiation!=1)))return 2;
        if (!(input>>ions>>factor>>dt>>abundance_tolerance>>heat_tolerance>>age>>logarithmic>>n) || n!=512 || !(dt>0) || !(abundance_tolerance>0) || !std::isfinite(age) || age<0 || (logarithmic!=0 && logarithmic!=1)) return 2;
        Model old; old.age=age;
        for (std::size_t i=0;i<n;++i) {
          double m,r,rho,T,L,x,y3;
          if (!(input>>m>>r>>rho>>T>>L>>x>>y3)) return 2;
          auto c=solar_scaled(x,.02); c.basis=AbundanceBasis::baryon_mass;
          c.metal_inventory=MetalInventory::gs98; c.X[1]=y3; c.X[2]-=y3;
          old.m.push_back(m); old.y.push_back({logarithmic?r:std::log(r),logarithmic?rho:std::log(rho),logarithmic?T:std::log(T),L}); old.comp.push_back(c);
        }
        // Physical ages are retained; evolved checkpoints carry native logarithms.
        old.M=old.m.back(); const Model preserved=old;
        ScaledOpacity shifted(bridge,factor);
        BlendedOpacity warm_bridge(warm,shifted,5.05,5.10);
        BlendedOpacity joined_hot(warm_bridge,hot,5.6,5.7);
        BlendedOpacity atomic(molecular,joined_hot,4.4,4.47);
        auto radiative=std::make_shared<ElementalOpacity>(atomic);
        auto material=std::make_shared<HotConduction>(conduction);
        CombinedOpacity opacity(radiative,material);
        ConditionalEnvelopeHeat transport(eos,collisions,*material,ions,2e6,3e6,parcel_radiation);
        Physics physics{&eos,radiative.get(),&nuclear,1.9,ConvectiveCriterion::ledoux,.1,1.,&losses};
        physics.microscopic=&transport;
        std::size_t good=0; std::vector<std::string> errors;
        for (std::size_t i=0;i<n;++i) {
          try {
            const auto e=eos.eval(old.T(i),old.rho(i),old.comp[i]);
            const auto k=opacity.eval(old.T(i),old.rho(i),old.comp[i]);
            if (!(e.P>0 && k.kappa>0) || !std::isfinite(k.dlnk_dlnT) ||
                !std::isfinite(k.dlnk_dlnRho) || !std::isfinite(k.dlnk_dX) || !std::isfinite(k.dlnk_dY3))
              throw std::runtime_error("invalid EOS pressure or opacity/derivative");
            ++good;
          } catch (const std::exception& e) {errors.push_back("node "+std::to_string(i)+": "+e.what());}
        }
        // The unchanged conditional heat provider has retained independent
        // derivatives. Do not repeat that batch at every sequential step.
        const double gravity=constants::G*old.M/std::pow(old.r(n-1),2);
        try { (void)atmosphere.eval(teff(old),gravity,old.comp.back()); }
        catch (const std::exception& e) {errors.push_back(std::string("atmosphere: ")+e.what());}
        std::cout << "{\"kind\":\"preflight\",\"factor\":" << factor << ",\"dt_seconds\":" << dt
          << ",\"abundance_tolerance\":" << abundance_tolerance
          << ",\"supported_nodes\":" << good << ",\"initial_Teff\":" << teff(old) << ",\"errors\":[";
        for (std::size_t i=0;i<errors.size();++i) {if(i) std::cout<<',';std::cout<<std::quoted(errors[i]);}
        std::cout << "]}\n" << std::flush;
        if (!errors.empty()) continue;
        EvolutionOptions options; options.max_coupling_iterations=30; options.max_abundance_change=.001;
        options.abundance_tolerance=abundance_tolerance; options.material_heat_tolerance=heat_tolerance;
        const auto start=std::chrono::steady_clock::now();
        const auto result=evolve_step(old,physics,atmosphere,dt,options);
        const double step_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        MixingRegions final_regions;
        if(result.converged) {
          transport.prescribed=result.total_species_rates;
          final_regions=convective_mixing_regions(result.model,physics);
        }
        double maxT=0,maxrho=0,maxX=0;
        for (std::size_t i=0;i<n;++i) {
          maxT=std::max(maxT,std::abs(result.model.y[i].lnT-old.y[i].lnT));
          maxrho=std::max(maxrho,std::abs(result.model.y[i].lnrho-old.y[i].lnrho));
          for (std::size_t k=0;k<NSPEC;++k) maxX=std::max(maxX,std::abs(result.model.comp[i].X[k]-old.comp[i].X[k]));
        }
        std::cout << "{\"kind\":\"step\",\"factor\":" << factor << ",\"dt_seconds\":" << dt
          << ",\"cn_carbon\":" << cn_carbon
          << ",\"parcel_radiation\":" << parcel_radiation
          << ",\"abundance_tolerance\":" << abundance_tolerance
          << ",\"converged\":" << result.converged << ",\"message\":" << std::quoted(result.message)
          << ",\"iterations\":" << result.coupling_iterations << ",\"residual\":" << scalar(result.residual)
          << ",\"material_heat_residual\":" << scalar(result.material_heat_residual)
          << ",\"maximum_species_face\":" << transport.maximum_species_face
          << ",\"screening_ions\":" << ions
          << ",\"correction\":" << scalar(result.correction) << ",\"abundance_residual\":" << scalar(result.abundance_residual)
          << ",\"luminosity_balance\":" << scalar(result.luminosity_balance)
          << ",\"nuclear_mass_balance\":" << scalar(result.nuclear_mass_balance)
          << ",\"convective_mass_fraction\":" << scalar(result.convective_mass_fraction)
          << ",\"elapsed_age_seconds\":" << scalar(result.model.age-old.age) << ",\"input_preserved\":" << identical(old,preserved)
          << ",\"initial_age_seconds\":" << old.age << ",\"age_seconds\":" << result.model.age
          << ",\"seconds\":" << step_seconds
          << ",\"max_dlnT\":" << maxT << ",\"max_dlnrho\":" << maxrho << ",\"max_dX\":" << maxX
          << ",\"Teff\":" << teff(result.model) << ",\"model\":[";
        for (std::size_t i=0;i<n;++i) {
          if(i) std::cout<<','; const auto& m=result.model;
          std::cout<<'['<<m.m[i]<<','<<m.r(i)<<','<<m.rho(i)<<','<<m.T(i)<<','<<m.y[i].L<<','<<m.comp[i].X[0]<<','<<m.comp[i].X[1]<<']';
        }
        std::cout << "],\"model_log\":[";
        for(std::size_t i=0;i<n;++i) {
          if(i)std::cout<<',';const auto& m=result.model;
          std::cout<<'['<<m.m[i]<<','<<m.y[i].lnr<<','<<m.y[i].lnrho<<','<<m.y[i].lnT<<','<<m.y[i].L<<','<<m.comp[i].X[0]<<','<<m.comp[i].X[1]<<']';
        }
        std::cout << "],\"nuclear_dXdt\":[";
        for(std::size_t i=0;i<n;++i) {
          if(i)std::cout<<',';const auto& m=result.model;
          const auto source=nuclear.eval(m.T(i),m.rho(i),m.comp[i]);
          std::cout<<'['<<source.dXdt[0]<<','<<source.dXdt[1]<<']';
        }
        std::cout << "],\"mixing_regions\":[";
        for(std::size_t i=0;i<final_regions.size();++i) {
          if(i)std::cout<<',';
          std::cout<<'['<<final_regions[i].first<<','<<final_regions[i].second<<']';
        }
        std::cout << "],\"total_species_rates\":[";
        for(std::size_t i=0;i<result.total_species_rates.size();++i) {
          if(i)std::cout<<',';
          std::cout<<'['<<result.total_species_rates[i][0]<<','<<result.total_species_rates[i][1]<<']';
        }
        std::cout << ']';
        if(select_radiation && result.converged) {
          // Compare the two prescriptions at identical final states and total
          // rates. These diagnostics do not feed back into the step.
          ConditionalEnvelopeHeat material_only(eos,collisions,*material,ions),
            with_radiation(eos,collisions,*material,ions,2e6,3e6,true);
          double conductivity_difference=0;
          std::cout<<",\"extra_radiation_luminosity\":[";
          const auto& m=result.model;
          for(std::size_t i=0;i+1<n;++i) {
            const auto a=material_only.heat_with_total_species_rate(i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],result.total_species_rates[i],false);
            const auto b=with_radiation.heat_with_total_species_rate(i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],result.total_species_rates[i],false);
            if(i)std::cout<<',';
            std::cout<<b.carried_luminosity-a.carried_luminosity;
            conductivity_difference=std::max(conductivity_difference,std::abs(b.conductivity-a.conductivity));
          }
          std::cout<<"],\"conductivity_difference\":"<<conductivity_difference;
        }
        std::cout << "}\n" << std::flush;
      } catch (const std::exception& e) {std::cout << "{\"error\":" << std::quoted(e.what()) << "}\n" << std::flush;}
    }
  } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
