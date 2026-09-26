// Isolated whole-star solver control. The candidate opacity is not accepted
// for a trajectory: its cooler TOPS density sampling fails the stated target.
// No microscopic diffusion is enabled in this integration control.
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
  if (argc!=8) return 2;
  try {
    const auto begin=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    MixtureOpacity molecular(argv[2]),warm(argv[3]),bridge(argv[4]),hot(argv[5]);
    TabulatedConduction conduction(argv[6]);
    CompositionAtmosphereGrid atmosphere(eos,argv[7],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
    PlasmaNeutrinoLosses losses;
    std::cerr << "load_seconds " << std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count() << '\n';
    std::cout << std::setprecision(17); std::string line;
    while (std::getline(std::cin,line)) {
      try {
        std::istringstream input(line); double factor,dt,abundance_tolerance; std::size_t n;
        if (!(input>>factor>>dt>>abundance_tolerance>>n) || n!=512 || !(dt>0) || !(abundance_tolerance>0)) return 2;
        Model old;
        for (std::size_t i=0;i<n;++i) {
          double m,r,rho,T,L,x,y3;
          if (!(input>>m>>r>>rho>>T>>L>>x>>y3)) return 2;
          auto c=solar_scaled(x,.02); c.basis=AbundanceBasis::baryon_mass;
          c.metal_inventory=MetalInventory::gs98; c.X[1]=y3; c.X[2]-=y3;
          old.m.push_back(m); old.y.push_back({std::log(r),std::log(rho),std::log(T),L}); old.comp.push_back(c);
        }
        // Age zero is the origin of this diagnostic interval, not formation.
        old.M=old.m.back(); const Model preserved=old;
        ScaledOpacity shifted(bridge,factor);
        BlendedOpacity warm_bridge(warm,shifted,5.05,5.10);
        BlendedOpacity joined_hot(warm_bridge,hot,5.6,5.7);
        BlendedOpacity atomic(molecular,joined_hot,4.4,4.47);
        auto radiative=std::make_shared<ElementalOpacity>(atomic);
        CombinedOpacity opacity(radiative,std::make_shared<HotConduction>(conduction));
        Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux,.1,1.,&losses};
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
        options.abundance_tolerance=abundance_tolerance;
        const auto start=std::chrono::steady_clock::now();
        const auto result=evolve_step(old,physics,atmosphere,dt,options);
        double maxT=0,maxrho=0,maxX=0;
        for (std::size_t i=0;i<n;++i) {
          maxT=std::max(maxT,std::abs(result.model.y[i].lnT-old.y[i].lnT));
          maxrho=std::max(maxrho,std::abs(result.model.y[i].lnrho-old.y[i].lnrho));
          for (std::size_t k=0;k<NSPEC;++k) maxX=std::max(maxX,std::abs(result.model.comp[i].X[k]-old.comp[i].X[k]));
        }
        std::cout << "{\"kind\":\"step\",\"factor\":" << factor << ",\"dt_seconds\":" << dt
          << ",\"abundance_tolerance\":" << abundance_tolerance
          << ",\"converged\":" << result.converged << ",\"message\":" << std::quoted(result.message)
          << ",\"iterations\":" << result.coupling_iterations << ",\"residual\":" << scalar(result.residual)
          << ",\"correction\":" << scalar(result.correction) << ",\"abundance_residual\":" << scalar(result.abundance_residual)
          << ",\"luminosity_balance\":" << scalar(result.luminosity_balance)
          << ",\"nuclear_mass_balance\":" << scalar(result.nuclear_mass_balance)
          << ",\"convective_mass_fraction\":" << scalar(result.convective_mass_fraction)
          << ",\"elapsed_age_seconds\":" << scalar(result.model.age) << ",\"input_preserved\":" << identical(old,preserved)
          << ",\"seconds\":" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
          << ",\"max_dlnT\":" << maxT << ",\"max_dlnrho\":" << maxrho << ",\"max_dX\":" << maxX
          << ",\"Teff\":" << teff(result.model) << ",\"model\":[";
        for (std::size_t i=0;i<n;++i) {
          if(i) std::cout<<','; const auto& m=result.model;
          std::cout<<'['<<m.m[i]<<','<<m.r(i)<<','<<m.rho(i)<<','<<m.T(i)<<','<<m.y[i].L<<','<<m.comp[i].X[0]<<','<<m.comp[i].X[1]<<']';
        }
        std::cout << "]}\n" << std::flush;
      } catch (const std::exception& e) {std::cout << "{\"error\":" << std::quoted(e.what()) << "}\n" << std::flush;}
    }
  } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
