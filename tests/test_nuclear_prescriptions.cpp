#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <limits>

using namespace ember;
int main() {
  int failures=0;
  auto check=[&](bool ok,const char* what,double value=0) {
    failures+=!ok;std::printf("[%s] %s (%.10g)\n",ok?"PASS":"FAIL",what,value);
  };
  auto rejects=[](auto f) {try{f();}catch(const std::exception&){return true;}return false;};
  std::ifstream input(std::string(EMBER_TEST_DATA_DIR)+"/nuclear_reference.dat");
  double rate_error=0,screen_error=0;int nr=0,ns=0;
  for(std::string line;std::getline(input,line);) {
    if(line.empty() || line[0]=='#')continue;
    std::istringstream row(line);char type;row>>type;
    if(type=='R') {
      double T,value,slope;int r;row>>T>>r>>value>>slope;
      const auto result=pp_bare_rate(T,static_cast<PPReaction>(r),PPRates::solar_fusion_ii);
      rate_error=std::max({rate_error,std::abs(result.molar_rate/value-1),std::abs(result.dlnrate_dlnT/slope-1)});++nr;
    } else {
      double T,rho,X,X3,eta,theta,ge,weak,svh;row>>T>>rho>>X>>X3>>eta>>theta>>ge>>weak>>svh;
      Composition comp;comp.basis=AbundanceBasis::baryon_mass;comp.X[0]=X;comp.X[1]=X3;comp.X[2]=1-X-X3;
      const auto a=pp_screening(T,rho,comp,PPReaction::pp,PPScreening::salpeter_van_horn);
      const auto b=pp_screening(T,rho,comp,PPReaction::pp,PPScreening::debye_fermi);
      screen_error=std::max({screen_error,std::abs(a.electron_eta-eta)/std::max(1.,std::abs(eta)),
        std::abs(a.electron_susceptibility/theta-1),std::abs(a.gamma_e/ge-1),
        std::abs(a.log_factor/svh-1),std::abs(b.log_factor/weak-1)});++ns;
    }
    if(!row) {check(false,"reference file parses");return 1;}
  }
  check(nr==18 && rate_error<2e-9,"SFII rates and slopes match independent adaptive energy integration",rate_error);
  check(ns==10 && screen_error<2e-9,"screening and Fermi susceptibility match independent bracketed quadrature",screen_error);
  double thermal=0,composition=0,baryons=0,mass=0;
  for(auto basis:{AbundanceBasis::baryon_mass,AbundanceBasis::atomic_mass})
  for(auto screening:{PPScreening::salpeter_van_horn,PPScreening::debye_fermi})
  for(auto [T,rho]:{std::pair{4.5e6,360.},std::pair{1.55e7,150.},std::pair{1e6,1000.}})
  for(double X3:{0.,.004,.05}) {
    PPChains nuc(PPRates::solar_fusion_ii,screening);
    auto comp=solar_scaled(.5,.02);comp.basis=basis;comp.X[1]=X3;comp.X[2]-=X3;
    const auto response=nuc.composition_response(T,rho,comp);const auto& s=response.state;
    const double step=1e-5;
    const double dt=(std::log(nuc.eval(T*std::exp(step),rho,comp).eps)-std::log(nuc.eval(T*std::exp(-step),rho,comp).eps))/(2*step);
    const double dr=(std::log(nuc.eval(T,rho*std::exp(step),comp).eps)-std::log(nuc.eval(T,rho*std::exp(-step),comp).eps))/(2*step);
    thermal=std::max({thermal,std::abs(dt-s.dlneps_dlnT),std::abs(dr-s.dlneps_dlnRho)});
    double release=0,sum=0,scale=0;
    for(std::size_t j=0;j<NSPEC;++j) {
      release-=s.dXdt[j]*nuclides[j].A/comp.abundance_weight(j)*constants::c*constants::c;
      sum+=s.dXdt[j];scale+=std::abs(s.dXdt[j]);
    }
    if(basis==AbundanceBasis::baryon_mass) baryons=std::max(baryons,std::abs(sum)/scale);
    mass=std::max(mass,std::abs(release/(s.eps+s.eps_neutrino)-1));
    for(std::size_t j=0;j<NSPEC;++j) {
      const double dx=1e-7;auto plus=comp,minus=comp;plus.X[j]+=dx;
      minus.X[j]+=comp.X[j]==0?2*dx:-dx;
      const auto p=nuc.eval(T,rho,plus),m=nuc.eval(T,rho,minus);
      auto diff=[&](double a,double b,double base){return comp.X[j]==0?(-3*base+4*a-b)/(2*dx):(a-b)/(2*dx);};
      composition=std::max(composition,std::abs(diff(p.eps,m.eps,s.eps)-response.deps_dX[j])/std::max(s.eps,std::abs(response.deps_dX[j])));
      for(std::size_t i=0;i<3;++i)
        composition=std::max(composition,std::abs(diff(p.dXdt[i],m.dXdt[i],s.dXdt[i])-response.d_dXdt_dX[i][j])/std::max(scale,std::abs(response.d_dXdt_dX[i][j])));
    }
  }
  check(thermal<2e-7,"screened heating derivatives follow actual rates and electron degeneracy",thermal);
  check(composition<1e-5,"all species screening/burning derivatives, both abundance bases and zero He3",composition);
  check(baryons<1e-15 && mass<1e-12,"modern rates conserve baryons and release the atomic mass defect",mass);
  auto comp=solar_scaled(.7,.02);comp.basis=AbundanceBasis::baryon_mass;
  const auto dilute=pp_screening(1e7,1e-12,comp,PPReaction::pp,PPScreening::salpeter_van_horn);
  const auto classical=pp_screening(1e7,1e-12,comp,PPReaction::pp,PPScreening::legacy_weak);
  check(std::abs(dilute.log_factor/classical.log_factor-1)<1e-4 && std::abs(dilute.electron_susceptibility-1)<1e-10,
    "dilute SVH screening recovers classical Debye charge sum",dilute.log_factor/classical.log_factor);
  const auto dense=pp_screening(1e6,1000,comp,PPReaction::he3_he3,PPScreening::salpeter_van_horn);
  const auto denser=pp_screening(1e6,1001,comp,PPReaction::he3_he3,PPScreening::salpeter_van_horn);
  check(dense.log_factor>2 && denser.log_factor>dense.log_factor && dense.dlog_dlnRho>0,
    "intermediate screening grows smoothly beyond the old exp(2) cap",dense.log_factor);
  PPChains modern(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
  check(rejects([&]{modern.eval(2.001e7,100,comp);}) && rejects([&]{modern.eval(1e6,1e6,comp);}),
    "unsupported rate temperatures and quantum-ion conditions are rejected");
  check(rejects([&]{modern.eval(0,1,comp);}) && rejects([&]{modern.eval(1e6,-1,comp);})
    && rejects([&]{modern.eval(std::numeric_limits<double>::quiet_NaN(),1,comp);}),"invalid thermal states are rejected");
  check(modern.eval(1e4,100,comp).eps==0,"retained cold cutoff has no burning");
  // Same particle inventory in the two mass conventions, including density.
  auto atomic=comp;atomic.basis=AbundanceBasis::atomic_mass;double scale=0;
  for(std::size_t j=0;j<NSPEC;++j){atomic.X[j]*=nuclides[j].A/mass_numbers[j];scale+=atomic.X[j];}
  for(double& x:atomic.X)x/=scale;
  const auto b=modern.eval(4.5e6,360,comp),a=modern.eval(4.5e6,360*scale,atomic);
  check(std::abs(a.eps*scale/b.eps-1)<1e-12,"changing mass convention preserves heat per physical volume",a.eps*scale/b.eps);
  return failures?1:0;
}
