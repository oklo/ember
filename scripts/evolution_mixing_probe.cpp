// Read the eight-column 0.1-Msun GS98 profile and evaluate the actual Ledoux
// mixing-region routine. Printed physical profiles are reconstructed with
// log/exp roundoff; this is an offline diagnostic, never a restart converter.
#include "ember/conduction_table.hpp"
#include "ember/eos_composition.hpp"
#include "ember/eos_mixture.hpp"
#include "ember/evolution.hpp"
#include "ember/opacity_mixture.hpp"
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

int main(int argc,char** argv) {
  using namespace ember;
  try {
    if(argc!=3)throw std::invalid_argument("usage: mixing-probe EOS_FAMILY OPACITY_DIRECTORY < profile.txt");
    MetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    auto radiative=std::make_shared<StellarMixtureOpacity>(argv[2]);
    TabulatedConduction conduction("data/conduction/condtab21wd_metals.dat");
    CombinedOpacity opacity(radiative,std::make_shared<HotConduction>(conduction));
    Physics physics{&eos,&opacity,nullptr,1.9,ConvectiveCriterion::ledoux};
    Model model;
    double m,r,rho,T,L,X,Y3,Y4;
    int fields;
    while((fields=std::scanf("%lf %lf %lf %lf %lf %lf %lf %lf",&m,&r,&rho,&T,&L,&X,&Y3,&Y4))==8) {
      if(!(m>0 && r>0 && rho>0 && T>0) || !std::isfinite(L))
        throw std::invalid_argument("invalid physical profile state");
      auto c=solar_scaled(X,.02);
      c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=Y3;c.X[2]=Y4;
      if(std::abs(c.sum()-1)>1e-10)throw std::invalid_argument("unnormalized composition");
      for(double v:c.X)if(!std::isfinite(v) || v<0)throw std::invalid_argument("invalid abundance");
      model.m.push_back(m);model.y.push_back({std::log(r),std::log(rho),std::log(T),L});
      model.comp.push_back(c);
    }
    if(fields!=EOF || model.size()<2)throw std::invalid_argument("incomplete profile");
    model.M=model.m.back();
    const auto weights=nodal_mass_weights(model);
    const auto regions=convective_mixing_regions(model,physics);
    double total=0;
    for(auto [begin,end]:regions) {
      double mass=0;for(std::size_t i=begin;i<end;++i)mass+=weights[i];
      if(end>begin+1)total+=mass/model.M;
      std::printf("{\"begin\":%zu,\"end_exclusive\":%zu,\"convective\":%s,\"mass_fraction\":%.17g,\"inner_enclosed_fraction\":%.17g,\"outer_enclosed_fraction\":%.17g}\n",
          begin,end,end>begin+1?"true":"false",mass/model.M,model.m[begin]/model.M,model.m[end-1]/model.M);
    }
    std::fprintf(stderr,"convective mass fraction %.4g; central region %s\n",total,
        regions.front().second>1?"convective":"nonconvective");
  } catch(const std::exception& e) {
    std::fprintf(stderr,"%s\n",e.what());return 1;
  }
}
