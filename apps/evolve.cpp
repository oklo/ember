#include "../examples/stellar_seed.hpp"
#include "ember/eos_composition.hpp"
#include "ember/evolution.hpp"
#include "ember/evolution_proxies.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <string>
#include <stdexcept>

namespace {
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
    std::puts("usage: ember-evolve [points [duration_years [initial_step_years [tolerance_scale]]]]\nDefaults: 512, 1e8, 1e7, 1. Fixed baryonic 0.1 Msun, X=.7 Z=.02 initially.\nImplicit pp burning and instantaneous Schwarzschild mixing, step-doubling control.\nInitial-evolution experiment: nominal isotope opacity; frozen COND boundary, |delta X| and Y3 <= .005.\nJSON on stdout, progress on stderr; exit zero only at requested duration.");return 0;
  }
  try {
    if(argc>5) throw std::invalid_argument("too many arguments");
    const std::size_t points=argc>1?number<std::size_t>(argv[1]):512;
    constexpr double year=365.25*86400.;
    const double duration=(argc>2?number<double>(argv[2]):1e8)*year;
    double dt=(argc>3?number<double>(argv[3]):1e7)*year;
    const double tolerance_scale=argc>4?number<double>(argv[4]):1.;
    if(points<128 || points>16384 || !std::isfinite(duration) || duration<=0 || !std::isfinite(dt) || dt<=0
        || !std::isfinite(tolerance_scale) || tolerance_scale<.001 || tolerance_scale>100)
      throw std::invalid_argument("invalid resolution or duration");
    const std::string data=EMBER_DATA_DIR;
    CompositionHelmholtzEos eos(data+"/eos/freeeos300_hhe_composition.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
    AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat");
    TopsOpacity high(data+"/opacity",TopsOpacity::Grid::composition);
    BlendedOpacity radiative(low,high,4.4,4.5);NominalAbundanceOpacity opacity(radiative);
    TabulatedAtmosphere cond(eos,data+"/atmosphere/cond_gn93_tau100_solar_proxy.dat",TabulatedAtmosphere::Mixture::allow_documented_proxy);
    FrozenCompositionAtmosphere atmosphere(eos,cond);PPChains nuclear;
    Physics physics{&eos,&opacity,&nuclear,1.9};
    auto composition=solar_scaled(.7,.02);composition.basis=AbundanceBasis::baryon_mass;
    const auto seed=example::stellar_seed(points,.1*constants::Msun,.15*constants::Rsun,composition,nuclear,atmosphere,1.5);
    const auto initial=relax(seed,physics,atmosphere);
    if(!initial.converged) throw std::runtime_error("initial equilibrium: "+initial.message);
    Model model=initial.model;
    struct Record {double age,step,R,L,Teff,X,Y3,error,energy,mass;std::size_t iterations;};
    std::vector<Record> history;
    auto record=[&](double step,double error,const EvolutionStep& result) {
      const double radius=model.r(points-1),lum=model.y.back().L;
      history.push_back({model.age/year,step/year,radius/constants::Rsun,lum/constants::Lsun,
        std::pow(lum/(4*M_PI*constants::sigma_SB*radius*radius),.25),model.comp[0].X[0],model.comp[0].X[1],
        error,result.luminosity_balance,result.nuclear_mass_balance,result.coupling_iterations});
    };
    record(0,0,{});
    EvolutionOptions options;std::size_t rejected=0;std::string last_failure;bool success=true;
    while(model.age<duration) {
      dt=std::min(dt,duration-model.age);
      if(history.size()>10000 || rejected>100 || dt<year) {success=false;last_failure="step limit or minimum timestep reached: "+last_failure;break;}
      const auto full=evolve_step(model,physics,atmosphere,dt,options);
      EvolutionStep first,second;
      if(full.converged) first=evolve_step(model,physics,atmosphere,.5*dt,options);
      if(first.converged) second=evolve_step(first.model,physics,atmosphere,.5*dt,options);
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
      if(!std::isfinite(error) || error>1) {++rejected;dt*=.5;std::fprintf(stderr,"retry dt=%.6g yr (%s; error %.3g)\n",dt/year,last_failure.c_str(),error);continue;}
      model=second.model;record(dt,error,second);
      std::fprintf(stderr,"age=%.8g yr, X=%.10g Y3=%.10g, R=%.9g L=%.9g, step error=%.3g\n",
        model.age/year,model.comp[0].X[0],model.comp[0].X[1],history.back().R,history.back().L,error);
      dt*=std::clamp(.9/std::sqrt(std::max(error,1e-6)),.5,2.);
    }
    std::printf("{\n\"calculation\":\"bounded initial evolution with coupled pp burning and instantaneous convective mixing\",\n\"converged\":%s,\n\"message\":",success?"true":"false");
    json_string(success?"requested duration reached":last_failure);
    std::printf(",\n\"points\":%zu,\n\"mass_basis\":\"conserved baryonic mass; nuclear rest mass release accounted in energy, Newtonian gravity\",\n\"mass_Msun\":0.1,\n\"age_origin\":\"specified static composition, not age since formation\",\n\"eos\":\"FreeEOS potential family; He3 number-density mapping; metals as He4\",\n\"opacity\":\"AESOPUS/TOPS varying nominal X, fixed Z=.02; He3 as He4\",\n\"atmosphere\":\"frozen GN93 COND tau100 T/P boundary, composition-dependent EOS inversion\",\n\"limitations\":\"early evolution only: |Xsurface-.7|<=.005, Y3<=.005; isotope opacity and atmosphere-mixture approximations; Schwarzschild criterion, no Ledoux/diffusion/conduction/remeshing\",\n\"step_error_tolerances\":{\"log_structure\":%.17g,\"absolute_abundance\":%.17g,\"relative_surface_luminosity\":%.17g},\n\"rejected_steps\":%zu,\n\"columns\":[\"age_yr\",\"step_yr\",\"R_Rsun\",\"L_Lsun\",\"Teff_K\",\"central_X\",\"central_Y3\",\"step_error_norm\",\"discrete_luminosity_balance\",\"nuclear_rest_mass_balance\",\"last_halfstep_coupling_iterations\"],\n\"history\":[\n",points,1e-5*tolerance_scale,1e-8*tolerance_scale,1e-4*tolerance_scale,rejected);
    for(std::size_t i=0;i<history.size();++i) {
      const auto& r=history[i];
      std::printf("%s[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu]",
        i?",\n":"",r.age,r.step,r.R,r.L,r.Teff,r.X,r.Y3,r.error,r.energy,r.mass,r.iterations);
    }
    std::printf("\n],\n\"profile_columns\":[\"mass_g\",\"radius_cm\",\"density_g_cm3\",\"temperature_K\",\"luminosity_erg_s\",\"X\",\"Y3\",\"Y4\"],\n\"profile\":[\n");
    for(std::size_t i=0;i<points;++i)
      std::printf("%s[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",i?",\n":"",model.m[i],model.r(i),model.rho(i),model.T(i),model.y[i].L,model.comp[i].X[0],model.comp[i].X[1],model.comp[i].X[2]);
    std::puts("\n]}\n");return success?0:1;
  } catch(const std::exception& e) {
    std::printf("{\"converged\":false,\"message\":");json_string(e.what());std::puts("}");return 1;
  }
}
