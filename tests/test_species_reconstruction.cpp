#include "ember/species_transport.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace ember;
namespace {
int checks=0;
void require(bool ok,const char* message) {++checks;if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double tolerance=2e-13) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tolerance,"flux reconstruction comparison failed");
}
template<class F> void rejects(F f) {
  bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}
  require(rejected,"invalid reconstruction input accepted");
}
struct Source final:Nuclear {
  NuclearState eval(double T,double rho,const Composition&)const override {
    NuclearState r;r.dXdt[0]=-.003*T;r.dXdt[1]=.002*rho;
    r.dXdt[2]=-r.dXdt[0]-r.dXdt[1];return r;
  }
  NuclearResponse composition_response(double T,double rho,const Composition& c)const override {
    NuclearResponse r;r.state=eval(T,rho,c);return r;
  }
  const char* name()const override{return "prescribed spatially varying source";}
};
Model model() {
  Model m;m.M=10;m.m={.2,.4,1.,2.,5.,10.};
  for(std::size_t i=0;i<m.m.size();++i) {
    auto c=solar_scaled(.3,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.02;c.X[2]-=c.X[1];m.comp.push_back(c);
    const double index=static_cast<double>(i);
    m.y.push_back({std::log(1.+index),std::log(1.+.2*index),std::log(1.+.1*index),1.});
  }
  return m;
}
void manufactured() {
  const Source source;
  // Construct OLD abundances from prescribed new abundances, local reactions,
  // and a known all-face flux. Neither a regional average nor reconstruction
  // is used to create the reference flux.
  for(const MixingRegions& regions:std::vector<MixingRegions>{
      {{0,6}},{{0,2},{2,3},{3,6}},{{0,1},{1,2},{2,3},{3,4},{4,5},{5,6}}})
    for(double dt:{.01,.2,1.}) {
      auto current=model(),previous=current;const auto w=nodal_mass_weights(current);
      std::vector<SpeciesVector> rates(5);
      for(std::size_t i=0;i<5;++i)rates[i]={.004*std::sin(1.+static_cast<double>(i)),-.002*std::cos(1.+static_cast<double>(i))};
      std::vector<SpeciesBoundaryFlux> boundaries;
      for(std::size_t r=0;r<regions.size();++r) {
        for(std::size_t i=regions[r].first;i<regions[r].second;++i) {
          current.comp[i].X[0]+=.01*static_cast<double>(r);current.comp[i].X[2]-=.01*static_cast<double>(r);
        }
        if(r+1<regions.size())boundaries.push_back({regions[r].second-1,rates[regions[r].second-1]});
      }
      for(std::size_t i=0;i<6;++i) {
        const auto s=source.eval(current.T(i),current.rho(i),current.comp[i]);
        const SpeciesVector left=i?rates[i-1]:SpeciesVector{},right=i<5?rates[i]:SpeciesVector{};
        previous.comp[i]=current.comp[i];
        for(std::size_t k=0;k<2;++k)
          previous.comp[i].X[k]-=dt*(s.dXdt[k]-(right[k]-left[k])/w[i]);
        previous.comp[i].X[2]=1-previous.comp[i].Z()-previous.comp[i].X[0]-previous.comp[i].X[1];
      }
      const auto result=reconstruct_species_fluxes(current,previous,source,regions,boundaries,dt);
      require(result.face_rates.size()==5,"missing internal face");
      for(std::size_t i=0;i<5;++i)for(std::size_t k=0;k<2;++k)close(result.face_rates[i][k],rates[i][k]);
      for(const auto& f:boundaries)require(result.face_rates[f.face]==f.rate,"supplied boundary was changed");
      for(const auto& v:result.cell_balances)for(double x:v)close(x,0,2e-16);
      for(const auto& v:result.region_balances)for(double x:v)close(x,0,2e-16);
      for(double x:result.integrated_balance)close(x,0,2e-16);
    }
}
void solver_and_unhidden_error() {
  const auto previous=model();const Source source;const MixingRegions regions{{0,2},{2,3},{3,6}};
  SpeciesFlux drift=[](std::size_t i,const Composition& a,const Composition& b,bool d) {
    SpeciesFaceResponse r;
    for(std::size_t k=0;k<2;++k) {
      const double strength=.01*(1+static_cast<double>(i))*(1+static_cast<double>(k));
      r.rate[k]=strength*(a.X[k]-b.X[k])+.0001*(1+static_cast<double>(i));
      if(d){r.dleft[k][k]=strength;r.dright[k][k]=-strength;}
    }
    return r;
  };
  const double dt=.3;const auto solution=burn_and_diffuse(previous,previous,source,regions,drift,dt);
  auto current=previous;current.comp=solution.composition;
  const auto result=reconstruct_species_fluxes(current,previous,source,regions,solution.boundary_fluxes,dt);
  for(double x:result.integrated_balance)close(x,0,2e-15);
  const auto w=nodal_mass_weights(current);
  // Check EACH local continuity equation directly from the original source,
  // rather than only a region/global sum which discards internal fluxes.
  for(std::size_t i=0;i<6;++i) {
    const auto s=source.eval(current.T(i),current.rho(i),current.comp[i]);
    const auto left=i?result.face_rates[i-1]:SpeciesVector{},right=i<5?result.face_rates[i]:SpeciesVector{};
    for(std::size_t k=0;k<2;++k)
      close(w[i]*(current.comp[i].X[k]-previous.comp[i].X[k])/dt+(right[k]-left[k]),w[i]*s.dXdt[k],2e-14);
  }
  // A finite numerical imbalance must be visible. Do not force both ends
  // closed by quietly changing the reconstructed source or total inventory.
  for(std::size_t i=3;i<6;++i){current.comp[i].X[0]+=.0001;current.comp[i].X[2]-=.0001;}
  const auto bad=reconstruct_species_fluxes(current,previous,source,regions,solution.boundary_fluxes,dt);
  const double error=.0001*(w[3]+w[4]+w[5])/current.M;
  close(bad.region_balances[2][0],error,2e-16);close(bad.cell_balances[5][0],error,2e-16);
  close(bad.integrated_balance[0],error,2e-16);
  for(const auto& f:solution.boundary_fluxes)require(bad.face_rates[f.face]==f.rate,"imbalance changed boundary flux");
}
void invalid() {
  const auto m=model();const Source source;const MixingRegions region{{0,6}};
  rejects([&]{reconstruct_species_fluxes(m,m,source,region,{},0);});
  rejects([&]{reconstruct_species_fluxes(m,m,source,{{1,6}},{},1);});
  rejects([&]{reconstruct_species_fluxes(m,m,source,{{0,5}},{},1);});
  rejects([&]{reconstruct_species_fluxes(m,m,source,{{0,2},{2,6}},{},1);});
  auto bad=m;bad.comp[0].X[0]+=.001;bad.comp[0].X[2]-=.001;
  rejects([&]{reconstruct_species_fluxes(bad,m,source,region,{},1);});
  bad=m;bad.comp[0].X[3]+=.001;bad.comp[0].X[2]-=.001;
  rejects([&]{reconstruct_species_fluxes(bad,m,source,region,{},1);});
  bad=m;bad.m[0]=bad.m[1];
  rejects([&]{reconstruct_species_fluxes(bad,bad,source,region,{},1);});
  std::vector<SpeciesBoundaryFlux> boundary{{2,{0,0}}};
  rejects([&]{reconstruct_species_fluxes(m,m,source,{{0,2},{2,6}},boundary,1);});
  boundary={{1,{std::numeric_limits<double>::quiet_NaN(),0}}};
  rejects([&]{reconstruct_species_fluxes(m,m,source,{{0,2},{2,6}},boundary,1);});
}
}
int main() {
  try {manufactured();solver_and_unhidden_error();invalid();std::cout<<checks<<" species reconstruction checks passed\n";return 0;}
  catch(const std::exception& e){std::cerr<<"after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}
