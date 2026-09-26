#include "ember/nuclear_cn.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace ember;
int main() {
  int checks=0,failures=0;
  auto check=[&](bool ok,const char* label){++checks;if(!ok){++failures;std::printf("FAIL %s\n",label);}};
  auto rejects=[](auto fn){try{fn();}catch(const std::exception&){return true;}return false;};
  auto close=[](double a,double b,double tol=1e-10){return std::abs(a-b)<=tol*std::max({1e-200,std::abs(a),std::abs(b)});};
  double worst_thermal=0,worst_composition=0,worst_mass=0;
  for(auto prescription:{PPRates::solar_fusion_ii,PPRates::solar_fusion_iii})
  for(auto screening:{PPScreening::debye_fermi,PPScreening::salpeter_van_horn})
  for(auto [T,rho]:{std::pair{5e6,300.},std::pair{8e6,300.},std::pair{1.59e7,8500.}}) {
    auto c=solar_scaled(.2,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.003;c.X[2]=.777;
    const auto initial=initial_gs98_cn(c);
    CNAbundances y{.4*initial[0],.15*initial[0],initial[2]+.45*initial[0]};
    CNNetwork network(prescription,screening);
    const auto response=network.response(T,rho,c,y);const auto& n=response.physical;const auto& s=n.state;
    double baryons=0,mass=0;
    for(std::size_t i=0;i<NSPEC;++i) {
      baryons+=s.dXdt[i];
      mass-=s.dXdt[i]*(nuclides[i].A/mass_numbers[i]-1)*constants::c*constants::c;
    }
    worst_mass=std::max(worst_mass,std::abs(mass/(s.eps+s.eps_neutrino)-1));
    check(close(mass,s.eps+s.eps_neutrino,1e-11),"physical rest-mass release");
    check(std::abs(baryons)<1e-13*std::abs(s.dXdt[0]),"physical baryon conservation");
    check(std::abs(s.dXdt[3]/12+s.dXdt[4]/13+s.dXdt[5]/14)<1e-13*std::abs(s.dXdt[0]),"catalyst number conservation");
    check(s.dXdt[1]==0 && s.dXdt[6]==0 && s.dXdt[7]==0,"unaffected species");
    const double h=1e-5;
    const double t=(std::log(network.response(T*std::exp(h),rho,c,y).physical.state.eps)
      -std::log(network.response(T*std::exp(-h),rho,c,y).physical.state.eps))/(2*h);
    const double r=(std::log(network.response(T,rho*std::exp(h),c,y).physical.state.eps)
      -std::log(network.response(T,rho*std::exp(-h),c,y).physical.state.eps))/(2*h);
    worst_thermal=std::max({worst_thermal,std::abs(t-s.dlneps_dlnT),std::abs(r-s.dlneps_dlnRho)});
    check(std::abs(t-s.dlneps_dlnT)<2e-6 && std::abs(r-s.dlneps_dlnRho)<2e-6,"thermal derivatives");
    for(std::size_t j=0;j<NSPEC;++j) {
      auto up=c,dn=c;const double step=std::max(1e-8,c.X[j]*1e-4);
      const bool edge=c.X[j]<step;up.X[j]+=step;dn.X[j]+=edge?2*step:-step;
      const auto u=network.response(T,rho,up,y).physical.state,d=network.response(T,rho,dn,y).physical.state;
      auto difference=[&](double a,double b,double base){return edge?(-3*base+4*a-b)/(2*step):(a-b)/(2*step);};
      double error=std::abs(difference(u.eps,d.eps,s.eps)-n.deps_dX[j])/(s.eps+std::abs(n.deps_dX[j]));
      for(std::size_t i=0;i<NSPEC;++i)
        error=std::max(error,std::abs(difference(u.dXdt[i],d.dXdt[i],s.dXdt[i])-n.d_dXdt_dX[i][j])/
          (std::abs(s.dXdt[0])+std::abs(n.d_dXdt_dX[i][j])));
      worst_composition=std::max(worst_composition,error);check(error<2e-6,"bulk composition derivatives");
    }
    for(std::size_t k=0;k<3;++k) {
      CNAbundances pure{};pure[k]=y[k];const auto source=network.response(T,rho,c,pure);
      check(close(source.physical.state.eps/y[k],response.deps_dY[k]),"independent catalyst energy derivative");
      for(std::size_t i=0;i<NSPEC;++i)
        check(close(source.physical.state.dXdt[i]/y[k],response.d_dXdt_dY[i][k]),"independent catalyst source derivative");
      if(k<2)check(source.physical.state.dXdt[2]==0,"carbon conversion makes no helium");
      else check(source.physical.state.dXdt[3]>0 && source.physical.state.dXdt[2]>0,"cycle returns carbon and produces helium");
    }
    const double flow=1e-25;
    CNAbundances equilibrium{};for(std::size_t k=0;k<3;++k)equilibrium[k]=flow/response.frequency[k];
    const auto eq=network.response(T,rho,c,equilibrium).physical.state;
    check(close(eq.dXdt[0],-4*flow) && close(eq.dXdt[2],4*flow),"stationary complete cycle fuel");
    check(close(eq.eps+eq.eps_neutrino,flow*(4*nuclides[0].A-nuclides[2].A)*constants::c*constants::c,1e-11),"stationary complete cycle energy");
    auto zero=c;zero.X[0]=0;const auto z=network.response(T,rho,zero,y);
    check(z.physical.state.eps==0 && z.physical.deps_dX[0]>0,"zero-hydrogen onset response");
  }
  auto c=solar_scaled(.2,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  const auto initial=initial_gs98_cn(c);const double total=initial[0]+initial[2];
  for(double fraction:{0.,.1,.5,1.}) {
    CNAbundances y{initial[0]*(1-fraction),0,initial[2]+fraction*initial[0]};
    const auto ledger=cn_physical_ledger(c,y);
    check(close(ledger.extra_metal_mass,2*initial[0]*fraction,1e-12),"two captured protons stored in metals");
    check(std::abs(ledger.electron_molality_difference)<1e-17,"full conversion electron-count match");
    check(close(ledger.helium4+ledger.metal_fraction+c.X[0]+c.X[1],1.),"physical composition normalization");
    auto original=c;original.X[0]+=ledger.extra_metal_mass;original.X[2]-=ledger.extra_metal_mass;
    const double proxy=(nuclides[0].A-nuclides[2].A/4)*ledger.extra_metal_mass*constants::c*constants::c;
    const double physical=(nuclides[3].A+2*nuclides[0].A-nuclides[5].A)*initial[0]*fraction*constants::c*constants::c;
    check(close(proxy-ledger.rest_energy_difference,physical,1e-11),"lookup plus binding correction restores physical Q");
  }
  const auto intermediate=cn_physical_ledger(c,{0,initial[0],initial[2]});
  check(close(intermediate.electron_molality_difference,-initial[0]/2),"C13 intermediate electron-count mismatch retained");
  for(auto rate:{std::array{1.,4.,.001},std::array{1.,0.,0.},std::array{0.,0.,1.},std::array{0.,0.,0.}})
  for(double dt:{0.,1e-10,1.,1e5,1e100,1e300}) {
    const auto y=cn_backward_euler(initial,rate,dt);
    check(std::all_of(y.begin(),y.end(),[](double v){return std::isfinite(v) && v>=0;}),"implicit positivity at arbitrary stiffness");
    check(close(y[0]+y[1]+y[2],total,1e-14),"implicit catalyst conservation at arbitrary stiffness");
    if(dt<1e6) {
      const std::array<double,3> r{rate[0]*y[0],rate[1]*y[1],rate[2]*y[2]};
      for(std::size_t k=0;k<3;++k) {
        const double residual=y[k]-initial[k]-dt*(r[(k+2)%3]-r[k]);
        check(std::abs(residual)<2e-11*total,"backward-Euler equation");
      }
    }
  }
  CNNetwork n;
  check(n.response(9e4,1,c,initial).physical.state.eps==0,"cold negligible-burning limit");
  check(rejects([&]{n.response(2.01e7,1,c,initial);}),"hot extrapolation rejected");
  check(rejects([&]{n.response(1e6,1e6,c,initial);}),"quantum screening domain rejected");
  auto atomic=c;atomic.basis=AbundanceBasis::atomic_mass;
  check(rejects([&]{n.response(1e7,1e3,atomic,initial);}),"atomic and baryon bases not silently mixed");
  check(cn_physical_ledger(c,{total*2,0,0}).metal_fraction>c.Z(),"local catalyst enrichment represented");
  check(rejects([&]{cn_backward_euler(initial,{-1.,1.,1.},1.);}),"negative rate rejected");
  std::printf("checks %d failures %d thermal %.17g composition %.17g mass %.17g\n",checks,failures,worst_thermal,worst_composition,worst_mass);
  return failures?1:0;
}
