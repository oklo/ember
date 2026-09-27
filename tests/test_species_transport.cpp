#include "ember/species_transport.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace ember;
namespace {
int checks=0;
void require(bool condition,const char* message) {
  ++checks;if(!condition)throw std::runtime_error(message);
}
void close(double a,double b,double tol=2e-12) {
  require(std::isfinite(a+b) && std::abs(a-b)<=tol,"species comparison failed");
}
template<class F> void rejects(F f) {
  bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}
  require(rejected,"invalid species input accepted");
}
struct Reaction final:Nuclear {
  // Independent first-order H -> He3 -> He4 test network, in mass fractions.
  double rate;
  explicit Reaction(double value):rate(value){}
  NuclearResponse composition_response(double T,double rho,const Composition& c) const override {
    NuclearResponse r;const double a=rate*T,b=rate*rho;
    r.state.dXdt[0]=-a*c.X[0];r.state.dXdt[1]=a*c.X[0]-b*c.X[1];r.state.dXdt[2]=b*c.X[1];
    r.d_dXdt_dX[0][0]=-a;r.d_dXdt_dX[1][0]=a;r.d_dXdt_dX[1][1]=-b;r.d_dXdt_dX[2][1]=b;
    return r;
  }
  NuclearState eval(double T,double rho,const Composition& c)const override{return composition_response(T,rho,c).state;}
  const char* name()const override{return "analytic sequential reaction";}
};
Model model() {
  Model m;m.M=10;m.m={.2,.4,1.,2.,5.,10.};
  for(std::size_t i=0;i<m.m.size();++i) {
    auto c=solar_scaled(.12+.02*static_cast<double>(i),.02);
    c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=i%2?.03:0.;c.X[2]-=c.X[1];m.comp.push_back(c);
    m.y.push_back({std::log(1.+static_cast<double>(i)),std::log(1.+.2*static_cast<double>(i)),
                  std::log(1.+.1*static_cast<double>(i)),1.});
  }
  return m;
}
// Independent dense Gaussian elimination on the complete region equations.
std::vector<double> dense(std::vector<std::vector<double>> a,std::vector<double> b) {
  const auto n=b.size();
  for(std::size_t i=0;i<n;++i) {
    auto pivot=i;for(std::size_t j=i+1;j<n;++j)if(std::abs(a[j][i])>std::abs(a[pivot][i]))pivot=j;
    std::swap(a[i],a[pivot]);std::swap(b[i],b[pivot]);require(a[i][i]!=0,"dense singularity");
    for(std::size_t j=i+1;j<n;++j) {
      const double f=a[j][i]/a[i][i];b[j]-=f*b[i];
      for(std::size_t k=i+1;k<n;++k)a[j][k]-=f*a[i][k];
    }
  }
  std::vector<double> x(n);
  for(std::size_t i=n;i-->0;) {
    double v=b[i];for(std::size_t j=i+1;j<n;++j)v-=a[i][j]*x[j];x[i]=v/a[i][i];
  }
  return x;
}
SpeciesFlux fick(double strength) {
  return [=](std::size_t face,const Composition& a,const Composition& b,bool derivatives) {
    SpeciesFaceResponse f;
    for(std::size_t j=0;j<2;++j) {
      const double g=strength*(1.+.2*static_cast<double>(face))*(j? .7:1.);
      f.rate[j]=g*(a.X[j]-b.X[j]);
      if(derivatives){f.dleft[j][j]=g;f.dright[j][j]=-g;}
    }
    return f;
  };
}
void reaction_diffusion() {
  const auto m=model();const auto w=nodal_mass_weights(m);
  const std::vector<MixingRegions> partitions{{{0,1},{1,2},{2,3},{3,4},{4,5},{5,6}},{{0,2},{2,3},{3,6}},{{0,6}}};
  for(const auto& regions:partitions)for(double rate:{0.,.03,3.})for(double strength:{0.,.01,10.})for(double dt:{.01,1.,100.}) {
    Reaction nuclear(rate);const auto f=fick(strength);
    const auto result=burn_and_diffuse(m,m,nuclear,regions,f,dt);
    const auto n=regions.size();std::vector<std::vector<double>> a(2*n,std::vector<double>(2*n));std::vector<double>b(2*n);
    for(std::size_t r=0;r<n;++r)for(std::size_t cell=regions[r].first;cell<regions[r].second;++cell) {
      a[2*r][2*r]+=w[cell]*(1+dt*rate*m.T(cell));
      a[2*r+1][2*r+1]+=w[cell]*(1+dt*rate*m.rho(cell));
      a[2*r+1][2*r]-=dt*w[cell]*rate*m.T(cell);
      for(std::size_t j=0;j<2;++j)b[2*r+j]+=w[cell]*m.comp[cell].X[j];
    }
    for(std::size_t r=0;r+1<n;++r)for(std::size_t j=0;j<2;++j) {
      const double g=dt*strength*(1.+.2*static_cast<double>(regions[r].second-1))*(j?.7:1.);
      a[2*r+j][2*r+j]+=g;a[2*r+2+j][2*r+2+j]+=g;
      a[2*r+j][2*r+2+j]-=g;a[2*r+2+j][2*r+j]-=g;
    }
    const auto exact=dense(a,b);
    for(std::size_t r=0;r<n;++r)for(std::size_t cell=regions[r].first;cell<regions[r].second;++cell) {
      for(std::size_t j=0;j<2;++j)close(result.composition[cell].X[j],exact[2*r+j]);
      close(result.composition[cell].sum(),1.);
      for(std::size_t j=3;j<NSPEC;++j)close(result.composition[cell].X[j],m.comp[cell].X[j],1e-16);
    }
    for(double v:result.integrated_balance)close(v,0.,2e-13);
    require(result.boundary_fluxes.size()==n-1,"wrong region flux count");
  }
}
void stiff_chain() {
  // A uniform mode must survive even when formation of the conventional
  // diagonals would subtract diffusive terms twelve orders larger than mass.
  constexpr std::size_t n=40;
  for(double strength:{0.,1.,1e6,1e12}) {
    std::vector<SpeciesMatrix> e(n),a(n-1),b(n-1);std::vector<SpeciesVector> rhs(n),d(n-1);
    for(std::size_t i=0;i<n;++i) {
      const double mass=.01+static_cast<double>(i)/n;
      e[i]={{{{mass,.1*mass}},{{-.2*mass,1.4*mass}}}};
      rhs[i]={mass*.17+.1*mass*.03,-.2*mass*.17+1.4*mass*.03};
    }
    for(std::size_t i=0;i+1<n;++i) {
      const double g=strength*(1.+static_cast<double>(i)/n);
      a[i]=b[i]={{{{g,.2*g}},{{-.1*g,.8*g}}}};
    }
    const auto result=solve_species_flux_chain(e,a,b,rhs,d);
    for(const auto& x:result){close(x[0],.17);close(x[1],.03);}
  }
}
void pp_and_inputs() {
  auto m=model();for(auto& y:m.y){y.lnT=std::log(1e7);y.lnrho=std::log(1000.);}
  const PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
  const MixingRegions regions{{0,2},{2,3},{3,6}};
  const auto reference=burn_and_mix(m,m,pp,regions,1e14,1e-14);
  const auto result=burn_and_diffuse(m,m,pp,regions,fick(0),1e14);
  for(std::size_t i=0;i<m.size();++i)for(std::size_t j=0;j<NSPEC;++j)close(result.composition[i].X[j],reference[i].X[j],2e-14);
  SpeciesTransportOptions independent;independent.abundance_tolerance=1e-6;
  independent.integrated_balance_tolerance=1e-14;
  const auto bounded=burn_and_diffuse(m,m,pp,regions,fick(0),1e14,independent);
  require(bounded.abundance_correction<=independent.abundance_tolerance,"local Newton correction bound");
  for(double v:bounded.integrated_balance)require(std::abs(v)<=1e-14,"independent integrated PP balance");
  independent.integrated_balance_tolerance=-1;
  rejects([&]{burn_and_diffuse(m,m,pp,regions,fick(0),1e14,independent);});
  const Reaction zero(0);
  rejects([&]{burn_and_diffuse(m,m,zero,{{0,2},{1,6}},fick(1),1);});
  rejects([&]{burn_and_diffuse(m,m,zero,regions,fick(1),0);});
  auto bad=m;bad.comp[0].X[3]+=.001;bad.comp[0].X[2]-=.001;
  rejects([&]{burn_and_diffuse(m,bad,zero,regions,fick(1),1);});
  bad=m;bad.M*=2;
  rejects([&]{burn_and_diffuse(m,bad,zero,regions,fick(1),1);});
  struct MetalReaction final:Nuclear {
    NuclearState eval(double,double,const Composition&)const override {
      NuclearState s;s.dXdt[2]=-.001;s.dXdt[3]=.001;return s;
    }
    NuclearResponse composition_response(double T,double rho,const Composition& c)const override {
      NuclearResponse r;r.state=eval(T,rho,c);return r;
    }
    const char* name()const override{return "unsupported metal reaction";}
  } metal_reaction;
  rejects([&]{burn_and_diffuse(m,m,metal_reaction,regions,fick(1),1);});
  // A fully mixed region does not call internal microscopic faces: their
  // species flux cancels. Their heat still belongs in the structure solve.
  const auto mixed=burn_and_diffuse(m,m,zero,{{0,6}},[](auto,const auto&,const auto&,bool)->SpeciesFaceResponse {
    throw std::runtime_error("internal region face incorrectly evaluated");
  },1);
  require(mixed.boundary_fluxes.empty(),"unexpected internal boundary");
}
void interior_guess() {
  const auto m=model();const Reaction zero(0),reaction(.03);
  const MixingRegions regions{{0,2},{2,3},{3,6}};
  for(const Nuclear* nuclear:{static_cast<const Nuclear*>(&zero),static_cast<const Nuclear*>(&reaction)}) {
    const auto reference=burn_and_diffuse(m,m,*nuclear,regions,fick(1),2);
    for(double x:{.03,.4,.7}) {
      auto guess=m.comp;
      for(auto& c:guess){c.X[0]=x;c.X[1]=.07;c.X[2]=1-c.Z()-x-.07;}
      SpeciesTransportOptions opts;opts.initial_guess=guess;
      const auto result=burn_and_diffuse(m,m,*nuclear,regions,fick(1),2,opts);
      for(std::size_t i=0;i<m.size();++i)for(std::size_t k=0;k<3;++k)
        close(result.composition[i].X[k],reference.composition[i].X[k],2e-13);
      for(double balance:result.integrated_balance)close(balance,0,2e-13);
    }
  }
  SpeciesTransportOptions opts;auto guess=m.comp;guess.pop_back();opts.initial_guess=guess;
  rejects([&]{burn_and_diffuse(m,m,zero,regions,fick(1),2,opts);});
  guess=m.comp;guess[0].X[0]=-.01;opts.initial_guess=guess;
  rejects([&]{burn_and_diffuse(m,m,zero,regions,fick(1),2,opts);});
  // An interior-only force law can begin from old local zeros. The seed is
  // not a source term: compare with the independently checked Fick solution.
  const auto ordinary=fick(1);
  SpeciesFlux interior=[&](std::size_t i,const Composition& a,const Composition& b,bool d) {
    for(std::size_t k=0;k<2;++k)if(a.X[k]<=0 || b.X[k]<=0)
      throw std::domain_error("test force requires an interior candidate");
    return ordinary(i,a,b,d);
  };
  const MixingRegions cells{{0,1},{1,2},{2,3},{3,4},{4,5},{5,6}};
  const auto reference=burn_and_diffuse(m,m,zero,cells,ordinary,2);
  SpeciesTransportOptions seeded;seeded.seed_present_species=true;
  const auto result=burn_and_diffuse(m,m,zero,cells,interior,2,seeded);
  require(result.interior_guess_adjusted,"old zero did not receive an interior guess");
  for(std::size_t i=0;i<m.size();++i)for(std::size_t k=0;k<3;++k)
    close(result.composition[i].X[k],reference.composition[i].X[k],2e-13);
}
}
int main() {
  try {stiff_chain();reaction_diffusion();pp_and_inputs();interior_guess();std::cout<<checks<<" species transport checks passed\n";return 0;}
  catch(const std::exception& e){std::cerr<<"after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}
