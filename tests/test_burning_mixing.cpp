#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

using namespace ember;
static int failures=0;
static void check(bool ok,const char* label,double value=0) {
  failures+=!ok;std::printf("[%s] %s (%.9g)\n",ok?"PASS":"FAIL",label,value);
}
template<class F> bool throws(F f) {try{f();}catch(const std::exception&){return true;}return false;}
int main() {
  PPChains nuclear;
  double baryons=0,mass=0,jacobian=0;
  for(double y3:{0.,.001,.05}) for(auto [T,rho]:{std::pair{6e6,100.},std::pair{1.2e7,150.},std::pair{1e6,1e4}}) {
    auto comp=solar_scaled(.5,.02);comp.basis=AbundanceBasis::baryon_mass;comp.X[1]=y3;comp.X[2]-=y3;
    const auto response=nuclear.composition_response(T,rho,comp);const auto& s=response.state;
    double sum=0,release=0,scale=0;
    for(std::size_t j=0;j<NSPEC;++j) {sum+=s.dXdt[j];scale+=std::abs(s.dXdt[j]);
      release-=s.dXdt[j]*nuclides[j].A/mass_numbers[j]*constants::c*constants::c;}
    baryons=std::max(baryons,std::abs(sum)/scale);
    mass=std::max(mass,std::abs(release/(s.eps+s.eps_neutrino)-1));
    for(std::size_t j=0;j<NSPEC;++j) {
      const double h=1e-7;auto p=comp,m=comp;p.X[j]+=h;m.X[j]-=h;
      // At zero He3 use a second-order forward difference: the quadratic
      // He3+He3 term can dominate a first-order secant even at small h.
      if(comp.X[j]==0)m.X[j]=2*h;
      const auto plus=nuclear.eval(T,rho,p),minus=nuclear.eval(T,rho,m);
      auto difference=[&](double a,double b,double base) {
        return comp.X[j]==0?(-3*base+4*a-b)/(2*h):(a-b)/(2*h);
      };
      for(std::size_t row=0;row<3;++row)
        jacobian=std::max(jacobian,std::abs(difference(plus.dXdt[row],minus.dXdt[row],s.dXdt[row])-response.d_dXdt_dX[row][j])
          /std::max(scale,std::abs(response.d_dXdt_dX[row][j])));
      jacobian=std::max(jacobian,std::abs(difference(plus.eps,minus.eps,s.eps)-response.deps_dX[j])
        /std::max(s.eps,std::abs(response.deps_dX[j])));
    }
  }
  check(baryons<1e-15,"pp stoichiometry conserves baryons without renormalization",baryons);
  check(mass<1e-12,"atomic mass defects equal deposited heat plus escaping neutrinos",mass);
  check(jacobian<1e-5,"screened nuclear composition Jacobian matches independent differences including zero He3",jacobian);
  // Pure initial pp burning isolates the Debye screening factor. Compare
  // two dilute/denser states against the independent classical charge sum.
  auto comp=solar_scaled(.5,.02);comp.basis=AbundanceBasis::baryon_mass;
  const double T=1e7,lo=1e-4,hi=100,e2=std::pow(4.803204673e-10,2);
  double charges=0;
  for(std::size_t j=0;j<NSPEC;++j) charges+=comp.X[j]/mass_numbers[j]*(nuclides[j].Z*nuclides[j].Z+nuclides[j].Z);
  const double H=e2/(constants::kB*T)*std::sqrt(4*M_PI*e2*constants::NA*charges/(constants::kB*T));
  const double measured=std::log((nuclear.eval(T,hi,comp).dXdt[1]/hi)/(nuclear.eval(T,lo,comp).dXdt[1]/lo));
  check(std::abs(measured-H*(std::sqrt(hi)-std::sqrt(lo)))<1e-13,"screening uses ne + sum(ni Zi^2), with no extra mean charge",measured);
  Model model;model.M=10;model.m={1,2,4,10};
  for(int i=0;i<4;++i) {
    auto c=comp;c.X[0]+=.01*i;c.X[2]-=.01*i;model.comp.push_back(c);
    model.y.push_back({std::log(1.+i),std::log(100.),std::log(8e6),1.});
  }
  const auto w=nodal_mass_weights(model);
  check(std::abs(std::accumulate(w.begin(),w.end(),0.)-model.M)<1e-14,"nodal quadrature includes unresolved central mass");
  Model cold=model;for(auto& y:cold.y)y.lnT=std::log(1e4);
  const MixingRegions split{{0,2},{2,4}};
  const auto mixed=burn_and_mix(cold,model,nuclear,split,1e15);
  bool conservative=true;
  for(std::size_t j=0;j<NSPEC;++j) {
    double old=0,now=0;for(std::size_t i=0;i<4;++i){old+=w[i]*model.comp[i].X[j];now+=w[i]*mixed[i].X[j];}
    conservative&=std::abs(now-old)<1e-14;
  }
  check(conservative && mixed[0].X==mixed[1].X && mixed[2].X==mixed[3].X && mixed[1].X!=mixed[2].X,
        "mixing conserves every species and does not cross disconnected regions");
  const double dt=1e17;const auto burnt=burn_and_mix(model,model,nuclear,{{0,4}},dt,1e-14);
  double residual=0;
  for(std::size_t j=0;j<3;++j) {
    double old=0,rate=0;for(std::size_t i=0;i<4;++i){old+=w[i]*model.comp[i].X[j]/model.M;
      rate+=w[i]*nuclear.eval(model.T(i),model.rho(i),burnt[i]).dXdt[j]/model.M;}
    residual=std::max(residual,std::abs(burnt[0].X[j]-old-dt*rate));
  }
  double old_h=0;for(std::size_t i=0;i<4;++i) old_h+=w[i]*model.comp[i].X[0]/model.M;
  check(residual<1e-12 && burnt[0].X[0]<old_h && burnt[0].X[1]>0 && burnt[0].X[2]>0,
        "stiff implicit burn and mixing solve integrated species equations with positive abundances",residual);
  check(throws([&]{burn_and_mix(model,model,nuclear,{{0,2},{1,4}},dt);}),"overlapping mixing regions are rejected");
  auto invalid=model;invalid.comp[0].basis=AbundanceBasis::atomic_mass;
  check(throws([&]{burn_and_mix(model,invalid,nuclear,{{0,4}},dt);}),"evolution rejects legacy atomic mass fractions");
  return failures?1:0;
}
