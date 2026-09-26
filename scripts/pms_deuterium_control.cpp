// Numerical coupling control: explicit D, analytic fully ionized EOS and
// constant opacity. This is NOT a physical low-mass PMS boundary or track.
#include "ember/stellar_seed.hpp"
#include "ember/eos_composite.hpp"
#include "ember/deuterium_burning.hpp"
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace ember;
class TestOpacity final:public Opacity {
 public:
  OpacityState eval(double,double,const Composition&) const override{return {1.,0.,0.,0.,0.,0.};}
  const char* name()const override{return "constant opacity, numerical PMS coupling control";}
};
class ContractingSeed final:public Nuclear {
 public:
  ContractingSeed(const Nuclear& source,double loss):source_(source),loss_(loss){}
  NuclearState eval(double T,double rho,const Composition& c)const override {
    auto s=source_.eval(T,rho,c);const double added=T*loss_,old=s.eps;
    s.eps+=added;s.dlneps_dlnT=(old*s.dlneps_dlnT+added)/s.eps;s.dlneps_dlnRho*=old/s.eps;return s;
  }
  const char* name()const override{return "specified seed entropy loss only";}
 private:const Nuclear& source_;double loss_;
};
int main(int argc,char** argv) {
  if(argc!=4 && argc!=5)return 2;
  const auto start=std::clock();std::cout<<std::setprecision(17);
  try {
    const double radius=std::stod(argv[1]),loss=std::stod(argv[2]),years=std::stod(argv[3]);
    const double seed_index=argc==5?std::stod(argv[4]):1.5;
    if(!(radius>0 && loss>0 && years>0))return 2;
    CompositeEos eos;TestOpacity opacity;GreyAtmosphere atmosphere(eos,opacity);
    PPDeuterium nuclear;ContractingSeed seed_source(nuclear,loss);
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
    Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux};
    auto construction=physics;construction.nuclear=&seed_source;
    auto guess=example::stellar_seed(192,.1*constants::Msun,radius*constants::Rsun,c,seed_source,atmosphere,seed_index,3000.);
    RelaxationOptions ro;ro.zone_threads=1;ro.max_iterations=200;
    auto initial=relax(guess,construction,atmosphere,ro);
    std::cout<<"{\"seed_converged\":"<<initial.converged<<",\"seed_message\":"<<std::quoted(initial.message);
    std::cout<<",\"seed_iterations\":"<<initial.iterations<<",\"seed_residual\":"<<initial.residual
        <<",\"seed_correction\":"<<initial.correction;
    if(!initial.converged){std::cout<<"}\n";return 1;}
    // Coupling diagnostic only: add a declared test inventory to an already
    // relaxed D-free structure. This is not the initial condition of the
    // physical PMS sequence. A thermal step must account for its consumption.
    for(auto& mixture:initial.model.comp){mixture.X[0]-=2e-5;mixture[Species::H2]=2e-5;}
    const auto record=[&](const char* label,const Model& m) {
      const auto w=nodal_mass_weights(m);double D=0,he3=0,entropy=0;
      for(std::size_t i=0;i<m.size();++i){D+=w[i]*m.comp[i][Species::H2];he3+=w[i]*m.comp[i].X[1];entropy+=w[i]*eos.eval(m.T(i),m.rho(i),m.comp[i]).S;}
      const double r=m.r(m.size()-1),L=m.y.back().L;
      std::cout<<','<<std::quoted(label)<<":{\"radius_Rsun\":"<<r/constants::Rsun<<",\"Teff_K\":"
          <<std::pow(L/(4*M_PI*constants::sigma_SB*r*r),.25)<<",\"central_T_K\":"<<m.T(0)
          <<",\"deuterium_g\":"<<D<<",\"helium3_g\":"<<he3<<",\"mean_entropy\":"<<entropy/m.M<<'}';
    };
    record("initial",initial.model);
    EvolutionOptions options;options.relaxation=ro;
    const double dt=years*31557600.;
    const auto full=evolve_step(initial.model,physics,atmosphere,dt,options);
    const auto half1=evolve_step(initial.model,physics,atmosphere,dt/2,options);
    const auto half2=half1.converged?evolve_step(half1.model,physics,atmosphere,dt/2,options):half1;
    std::cout<<",\"steps\":[";bool first=true;
    for(const auto* s:{&full,&half1,&half2}) {
      if(!first)std::cout<<',';first=false;
      std::cout<<"{\"converged\":"<<s->converged<<",\"message\":"<<std::quoted(s->message)
          <<",\"first_law_error\":"<<s->luminosity_balance<<",\"nuclear_mass_error\":"<<s->nuclear_mass_balance
          <<",\"convective_mass_fraction\":"<<s->convective_mass_fraction
          <<",\"nuclear_luminosity\":"<<s->nuclear_luminosity<<",\"contraction_luminosity\":"<<s->gravitational_luminosity<<'}';
    }
    std::cout<<']';
    if(full.converged && half2.converged){record("full",full.model);record("two_half",half2.model);}
    std::cout<<",\"cpu_seconds\":"<<double(std::clock()-start)/CLOCKS_PER_SEC<<"}\n";
    return full.converged && half1.converged && half2.converged?0:1;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
