// Optional profile audit; compile as documented in docs/NUCLEAR.md.
#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <cstdio>
#include <cmath>
#include <iostream>
using namespace ember;
int main(int argc,char** argv) {
  const bool metal=argc==2 && std::string(argv[1])=="--gs98";
  if(argc!=1 && !metal)return 1;
  double T,rho,X,Y3;
  while(std::cin>>T>>rho>>X>>Y3) {
    auto comp=solar_scaled(X,.02);comp.basis=AbundanceBasis::baryon_mass;comp.X[1]=Y3;comp.X[2]-=Y3;
    if(metal)comp.metal_inventory=MetalInventory::gs98;
    try {
      std::printf("%.17g %.17g",T,rho);
      for(auto r:{PPReaction::pp,PPReaction::he3_he3,PPReaction::he3_he4}) {
        const auto modern=pp_bare_rate(T,r,PPRates::solar_fusion_ii),old=pp_bare_rate(T,r,PPRates::legacy);
        std::printf(" %.17g",old.molar_rate>0?modern.molar_rate/old.molar_rate:1.);
      }
      for(auto screening:{PPScreening::legacy_weak,PPScreening::debye_fermi,PPScreening::salpeter_van_horn}) {
        const auto pp=pp_screening(T,rho,comp,PPReaction::pp,screening);
        const auto he=pp_screening(T,rho,comp,PPReaction::he3_he3,screening);
        const PPChains nuclear(PPRates::solar_fusion_ii,screening);
        std::printf(" %.17g %.17g %.17g",pp.log_factor,he.log_factor,nuclear.eval(T,rho,comp).eps);
      }
      const auto p=pp_screening(T,rho,comp,PPReaction::pp,PPScreening::salpeter_van_horn);
      const auto n=PPChains(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn).eval(T,rho,comp);
      const double npp=.5*rho*X*X*pp_bare_rate(T,PPReaction::pp,PPRates::solar_fusion_ii).molar_rate*std::exp(p.log_factor);
      const auto he=pp_screening(T,rho,comp,PPReaction::he3_he4,PPScreening::salpeter_van_horn);
      const double n34=rho*(Y3/3)*(comp.X[2]/4)*pp_bare_rate(T,PPReaction::he3_he4,PPRates::solar_fusion_ii).molar_rate*std::exp(he.log_factor);
      std::printf(" %.17g %.17g %.17g %.17g %.17g %.17g\n",p.electron_eta,p.electron_susceptibility,p.gamma_e,p.zeta,n.eps_neutrino,npp>0?n34/npp:0.);
    } catch(const std::exception& e) {std::fprintf(stderr,"%s\n",e.what());return 1;}
  }
  return std::cin.eof()?0:1;
}
