#include "ember/dense_nuclear.hpp"
#include "ember/nuclear.hpp"
#include "ember/nuclear_cn.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
int main() {
  int failures=0;
  auto check=[&](bool ok,const char* msg,double x=0.) {failures+=!ok;std::printf("[%s] %s %.4g\n",ok?"PASS":"FAIL",msg,x);};
  auto comp=solar_scaled(.01,.13);comp.basis=AbundanceBasis::baryon_mass;
  comp.metal_inventory=MetalInventory::gs98;comp.X[1]=.04;comp.X[2]-=.04;
  comp.cn_molality=initial_gs98_cn(comp);comp.cn_mass_convention=CNMassConvention::explicit_metal_mass;
  const NuclearPair pair{1,1,1,1,4.09e-25};
  // Independent Python evaluation of the volume rates in Y06 Eqs.33,36--39,
  // divided by the pair number densities. Integer reacting masses here.
  Composition reference_mix;reference_mix.basis=AbundanceBasis::baryon_mass;
  reference_mix.X[0]=.1;reference_mix.X[1]=.1;reference_mix.X[2]=.8;
  struct Reference { double T,rho;NuclearPair pair;double optimal,high; };
  const Reference references[]{
    {10,5e4,{1,1,1,1,4.09e-25},4.8319342600615378e-24,1.6678906088082759e-21},
    {2e5,5e4,{2,2,3,4,.000561},1.7096809115359542e-50,6.7427719230545118e-44},
    {1e6,1e5,{1,7,1,14,.00168},8.9145063715164426e-41,7.9535316441627214e-38}
  };
  for(auto model:{DenseNuclearModel::uniform_optimal,DenseNuclearModel::uniform_high}) {
    for(const auto& ref:references) {
      const double expected=model==DenseNuclearModel::uniform_optimal?ref.optimal:ref.high;
      const double error=dense_nuclear_rate(ref.T,ref.rho,reference_mix,ref.pair,model).molar_rate/expected-1;
      check(std::abs(error)<1e-10,"independent volume-rate reference",error);
    }
    const auto cold=dense_nuclear_rate(10,5e4,comp,pair,model);
    const auto warmer=dense_nuclear_rate(1000,5e4,comp,pair,model);
    check(cold.molar_rate>0&&std::abs(cold.molar_rate/warmer.molar_rate-1)<1e-12
          &&std::abs(cold.dlnrate_dlnT)<1e-12,"finite temperature-independent quantum limit",cold.molar_rate);
    double worst=0;
    for(double T:{1e4,1e5,2e5,1e6})for(double rho:{100.,5e4,1e6}) {
      const auto out=dense_nuclear_rate(T,rho,comp,pair,model);const double h=1e-5;
      auto compare=[&](double num,double exact){worst=std::max(worst,std::abs(num-exact)/std::max(1.,std::abs(exact)));};
      compare((std::log(dense_nuclear_rate(T*std::exp(h),rho,comp,pair,model).molar_rate)
              -std::log(dense_nuclear_rate(T*std::exp(-h),rho,comp,pair,model).molar_rate))/(2*h),out.dlnrate_dlnT);
      compare((std::log(dense_nuclear_rate(T,rho*std::exp(h),comp,pair,model).molar_rate)
              -std::log(dense_nuclear_rate(T,rho*std::exp(-h),comp,pair,model).molar_rate))/(2*h),out.dlnrate_dlnRho);
      for(std::size_t j=0;j<NSPEC;++j) {
        auto plus=comp,minus=comp;const double dx=std::min(h,comp.X[j]/4);
        if(dx==0)continue;plus.X[j]+=dx;minus.X[j]-=dx;
        compare((std::log(dense_nuclear_rate(T,rho,plus,pair,model).molar_rate)
                -std::log(dense_nuclear_rate(T,rho,minus,pair,model).molar_rate))/(2*dx),out.dlnrate_dX[j]);
      }
    }
    check(worst<2e-5,"dense coefficient derivatives match finite differences",worst);
  }
  // Eq.13: independent classical Gamow-peak coefficient for constant S.
  using namespace constants;
  const double T=1e7,mu=amu/2,hb=h/(2*M_PI),e2=4.803204673e-10*4.803204673e-10;
  const double tau=std::cbrt(27*M_PI*M_PI*mu*e2*e2/(2*kB*T*hb*hb));
  const double reference=NA*4*std::sqrt(2*kB*T*tau/(9*mu))*pair.s0*1.602176634e-30/(kB*T)*std::exp(-tau);
  const auto dilute=dense_nuclear_rate(T,1e-20,comp,pair,DenseNuclearModel::uniform_optimal);
  check(std::abs(dilute.molar_rate/reference-1)<1e-7,"classical dilute limit matches Gamow expression",dilute.molar_rate/reference-1);

  set_quantum_screening(1.6);
  const PPCNNetwork network(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn,PPRates::solar_fusion_iii);
  const auto old=network.composition_response(4e6,400,comp);
  for(auto model:{DenseNuclearModel::uniform_optimal,DenseNuclearModel::uniform_high}) {
    set_dense_nuclear_model(model);
    const auto warm=network.composition_response(4e6,400,comp);
    check(warm.state.eps==old.state.eps&&warm.state.dXdt==old.state.dXdt,"warm network unchanged below the quantum join");
    double worst=0,baryons=0,energy=0;
    for(double t:{3e4,8e4,1.2e5,2e5,3e5,5e5,1e6}) {
      const auto out=network.composition_response(t,5e4,comp);const double step=1e-5;
      const auto plus=network.eval(t*std::exp(step),5e4,comp),minus=network.eval(t*std::exp(-step),5e4,comp);
      const double thermal=(std::log(plus.eps)-std::log(minus.eps))/(2*step);
      worst=std::max(worst,std::abs(thermal-out.state.dlneps_dlnT)/std::max(1.,std::abs(thermal)));
      const double density=(std::log(network.eval(t,5e4*std::exp(step),comp).eps)
                          -std::log(network.eval(t,5e4*std::exp(-step),comp).eps))/(2*step);
      worst=std::max(worst,std::abs(density-out.state.dlneps_dlnRho)/std::max(1.,std::abs(density)));
      double sum=0,rest=0,norm=0;
      for(std::size_t i=0;i<NSPEC;++i) {sum+=out.state.dXdt[i];norm+=std::abs(out.state.dXdt[i]);}
      baryons=std::max(baryons,std::abs(sum)/norm);
      const auto pp=PPChains(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn).composition_response(t,5e4,comp);
      for(std::size_t i=0;i<NSPEC;++i)rest-=pp.state.dXdt[i]*(nuclides[i].A/mass_numbers[i]-1)*c_light*c_light;
      energy=std::max(energy,std::abs(rest/(pp.state.eps+pp.state.eps_neutrino)-1));
      for(std::size_t j:{std::size_t(0),std::size_t(1),std::size_t(2)}) {
        auto a=comp,b=comp;const double dx=std::min(1e-6,comp.X[j]/4);a.X[j]+=dx;b.X[j]-=dx;
        const double num=(network.eval(t,5e4,a).eps-network.eval(t,5e4,b).eps)/(2*dx);
        worst=std::max(worst,std::abs(num-out.deps_dX[j])/std::max(out.state.eps,std::abs(num)));
      }
    }
    check(worst<3e-5,"joined network derivatives match finite differences",worst);
    check(baryons<2e-12,"cold pp and CN networks conserve baryons",baryons);
    check(energy<2e-12,"cold pp heat plus neutrinos matches nuclear mass loss",energy);
  }
  set_dense_nuclear_model(DenseNuclearModel::none);set_quantum_screening(0);
  return failures?1:0;
}
