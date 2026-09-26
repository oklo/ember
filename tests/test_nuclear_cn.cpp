#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace ember;
int main() {
  int checks=0,failures=0;
  auto check=[&](bool pass,const char* name) {
    ++checks;if(!pass){++failures;std::printf("FAIL %s\n",name);}
  };
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  double worst_thermal=0,worst_composition=0,worst_energy=0;
  for(auto basis:{AbundanceBasis::baryon_mass,AbundanceBasis::atomic_mass})
  for(auto rates:{PPRates::solar_fusion_ii,PPRates::solar_fusion_iii})
  for(auto screen:{PPScreening::salpeter_van_horn,PPScreening::debye_fermi})
  for(auto [T,rho]:{std::pair{5e6,300.},std::pair{1.32e7,1e4},std::pair{1.9e7,300.}})
  for(double carbon:{0.,.45,1.}) {
    auto c=solar_scaled(.01,.02);c.basis=basis;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.003;c.X[2]=.967;
    CNCycle cn(rates,screen,carbon);
    const auto response=cn.composition_response(T,rho,c);const auto& s=response.state;
    const double h=1e-5;
    const auto a=cn.eval(T*std::exp(h),rho,c),b=cn.eval(T*std::exp(-h),rho,c);
    const auto d=cn.eval(T,rho*std::exp(h),c),e=cn.eval(T,rho*std::exp(-h),c);
    const double t=(std::log(a.eps)-std::log(b.eps))/(2*h);
    const double r=(std::log(d.eps)-std::log(e.eps))/(2*h);
    worst_thermal=std::max({worst_thermal,std::abs(t-s.dlneps_dlnT),std::abs(r-s.dlneps_dlnRho)});
    check(std::abs(t-s.dlneps_dlnT)<1e-6 && std::abs(r-s.dlneps_dlnRho)<1e-6,"thermal derivatives");
    double released=0,sum=0;
    for(std::size_t i=0;i<NSPEC;++i) {
      released-=s.dXdt[i]*nuclides[i].A/c.abundance_weight(i)*constants::c*constants::c;
      sum+=s.dXdt[i];
      if(i!=0 && i!=2)check(s.dXdt[i]==0,"closed cycle preserves other species");
      auto up=c,dn=c;const double step=std::max(c.X[i]*1e-4,1e-8);
      const bool edge=c.X[i]<step;up.X[i]+=step;dn.X[i]+=edge?2*step:-step;
      const auto su=cn.eval(T,rho,up),sd=cn.eval(T,rho,dn);
      auto difference=[&](double a,double b,double base){return edge?(-3*base+4*a-b)/(2*step):(a-b)/(2*step);};
      double err=std::abs(difference(su.eps,sd.eps,s.eps)-response.deps_dX[i])/(s.eps+std::abs(response.deps_dX[i]));
      for(std::size_t j=0;j<NSPEC;++j)
        err=std::max(err,std::abs(difference(su.dXdt[j],sd.dXdt[j],s.dXdt[j])-response.d_dXdt_dX[j][i])/
          std::max(1e-200,std::abs(s.dXdt[0])+std::abs(response.d_dXdt_dX[j][i])));
      worst_composition=std::max(worst_composition,err);check(err<2e-6,"all abundance derivatives including screening");
    }
    worst_energy=std::max(worst_energy,std::abs(released/(s.eps+s.eps_neutrino)-1));
    check(std::abs(released/(s.eps+s.eps_neutrino)-1)<1e-11,"mass defect equals deposited plus neutrino energy");
    check(basis!=AbundanceBasis::baryon_mass || sum==0,"integer baryon conservation");
    const auto pp=PPChains(PPRates::solar_fusion_ii,screen).composition_response(T,rho,c);
    const auto both=PPCNO(PPRates::solar_fusion_ii,screen,rates,carbon).composition_response(T,rho,c);
    check(both.state.eps==pp.state.eps+s.eps,"combined energy");
    check(both.state.eps_neutrino==pp.state.eps_neutrino+s.eps_neutrino,"combined neutrinos");
    check(std::abs(both.state.dlneps_dlnT*both.state.eps-pp.state.dlneps_dlnT*pp.state.eps-s.dlneps_dlnT*s.eps)
      <1e-12*both.state.eps,"combined thermal slope");
    for(std::size_t i=0;i<NSPEC;++i) {
      check(both.state.dXdt[i]==pp.state.dXdt[i]+s.dXdt[i],"combined abundance source");
      check(both.deps_dX[i]==pp.deps_dX[i]+response.deps_dX[i],"combined energy derivative");
      for(std::size_t j=0;j<NSPEC;++j)
        check(both.d_dXdt_dX[i][j]==pp.d_dXdt_dX[i][j]+response.d_dXdt_dX[i][j],"combined abundance derivative");
    }
  }
  CNCycle cn;
  auto c=solar_scaled(.02,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  auto zero=c;zero.X[0]=0;
  auto z=cn.composition_response(1.3e7,1e4,zero);
  check(z.state.eps==0 && z.deps_dX[0]>0 && z.d_dXdt_dX[0][0]<0,"initially absent hydrogen response");
  zero=c;for(std::size_t j=3;j<NSPEC;++j)zero.X[j]=0;
  z=cn.composition_response(1.3e7,1e4,zero);
  check(z.state.eps==0 && z.deps_dX[3]>0,"initially absent catalysts response");
  const auto original=cn.eval(1.3e7,1e4,c);zero=c;
  for(std::size_t j=3;j<NSPEC;++j)zero.X[j]=0;zero.X[7]=c.Z();
  check(std::abs(cn.eval(1.3e7,1e4,zero).eps/original.eps-1)<1e-14,"GS98 proxy labels do not define isotopes");
  check(cn.eval(9e4,1.,c).eps==0,"explicit negligible-burning cutoff");
  check(rejects([&]{cn.eval(1e6,1e6,c);}),"quantum-ion domain rejected");
  check(rejects([&]{cn.eval(2.01e7,1e4,c);}),"S-factor temperature domain rejected");
  check(rejects([&]{cn.eval(std::numeric_limits<double>::quiet_NaN(),1e4,c);}),"nonfinite state rejected");
  zero=c;zero.metal_inventory=MetalInventory::carried_isotopes;
  check(rejects([&]{cn.eval(1.3e7,1e4,zero);}),"literal isotope inventory not silently treated as GS98");
  check(rejects([]{CNCycle n(PPRates::legacy);}),"legacy CN rate unavailable");
  check(rejects([]{CNCycle n(PPRates::solar_fusion_iii,PPScreening::legacy_weak);}),"capped weak screening unavailable");
  check(rejects([]{CNCycle n(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,1.01);}),"invalid converted fraction rejected");
  std::printf("checks %d failures %d thermal %.17g composition %.17g energy %.17g\n",
    checks,failures,worst_thermal,worst_composition,worst_energy);
  return failures?1:0;
}
