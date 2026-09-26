#include "ember/cn_burning.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
int main() {
  int checks=0,failures=0;
  auto check=[&](bool ok,const char* label){++checks;if(!ok){++failures;std::printf("FAIL %s\n",label);}};
  PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);CNNetwork cn;
  double maxmass=0,maxres=0;std::size_t maxit=0;
  for(bool late:{false,true})for(bool mixed:{false,true})for(double years:{1e4,1e6,1e8}) {
    Model model;model.M=1;model.m={.1,.4,.7,1.};std::vector<CNAbundances> old_cn;
    for(std::size_t i=0;i<4;++i) {
      const double X=late?.0008:.1+.1*i;
      auto c=solar_scaled(X,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=late?1e-9:.001;c.X[2]=1-X-c.X[1]-.02;model.comp.push_back(c);
      model.y.push_back({std::log(1.+i),std::log(late?8500.:300.),std::log(late?1.3e7:5e6+i*1e6),1.});
      const auto y=initial_gs98_cn(c);const double f=.2*i;
      old_cn.push_back({y[0]*(1-f),0,y[2]+f*y[0]});
    }
    const MixingRegions regions=mixed?MixingRegions{{0,4}}:MixingRegions{{0,1},{1,2},{2,3},{3,4}};
    const double dt=years*3.15576e7;
    try {
      const auto result=burn_cn_and_mix(model,model,old_cn,pp,cn,regions,dt,1e-14);
      maxmass=std::max(maxmass,std::abs(result.nuclear_mass_balance));
      maxres=std::max(maxres,result.maximum_equation_residual);maxit=std::max(maxit,result.maximum_iterations);
      check(result.maximum_equation_residual<1e-14,"coupled abundance equations");
      // The shortest cold burn changes H by only about 1e-9; the energy
      // audit differences that finite change and has corresponding roundoff.
      check(std::abs(result.nuclear_mass_balance)<2e-5,"physical energy equals composition mass defect");
      const std::array<double,4> weights{.25,.3,.3,.15};
      for(auto [begin,end]:regions) {
        std::array<double,5> residual{};
        for(std::size_t i=begin;i<end;++i) {
          const auto& c=result.lookup[i];const auto& y=result.catalysts[i];
          const auto p=pp.eval(model.T(i),model.rho(i),c);
          const auto n=cn.response(model.T(i),model.rho(i),c,y).physical.state;
          for(std::size_t k=0;k<2;++k)
            residual[k]+=weights[i]*(c.X[k]-model.comp[i].X[k]-dt*(p.dXdt[k]+n.dXdt[k]));
          for(std::size_t k=0;k<3;++k)
            residual[k+2]+=weights[i]*(mass_numbers[k+3]*(y[k]-old_cn[i][k])-dt*n.dXdt[k+3]);
          const auto physical=cn_physical_ledger(c,y);
          check(c.X[0]>=0 && physical.helium4>=0 && *std::min_element(y.begin(),y.end())>=0,"positive actual fuel and catalysts");
          if(mixed)check(c.X==result.lookup[0].X && y==result.catalysts[0],"simultaneous homogeneous convective burning");
        }
        check(std::abs(*std::max_element(residual.begin(),residual.end(),[](double a,double b){return std::abs(a)<std::abs(b);}))<1e-12,"independent integrated H/He3/C/N equations");
      }
      std::printf("late %d mixed %d years %.4g iterations %zu residual %.4g mass %.4g\n",late,mixed,years,result.maximum_iterations,result.maximum_equation_residual,result.nuclear_mass_balance);
    }catch(const std::exception& e){++failures;std::printf("FAILED CASE late %d mixed %d years %.4g: %s\n",late,mixed,years,e.what());}
  }
  std::printf("checks %d failures %d max_iterations %zu max_residual %.17g max_mass %.17g\n",checks,failures,maxit,maxres,maxmass);
  return failures?1:0;
}
