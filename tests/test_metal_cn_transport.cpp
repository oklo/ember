#include "ember/metal_cn_transport.hpp"
#include "ember/constants.hpp"
#include "../src/metal_cn_source.hpp"
#include "ember/evolution_checkpoint.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace ember;
namespace {
Composition composition(double h,double conversion,double enrichment,double inert_factor=1) {
  auto c=solar_scaled(h,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=.001;c.X[2]-=.001;const auto y=initial_gs98_cn(c);
  c.cn_molality=CNAbundances{enrichment*y[0]*(1-conversion),enrichment*.03*y[0],enrichment*(y[2]+(conversion-.03)*y[0])};
  c=explicit_cn_material(c);auto v=metal_cn_abundances(c);v[5]*=inert_factor;
  return metal_cn_composition(c,v);
}
MetalSpeciesFaceResponse material_flux(const Composition& a,const Composition& b,double factor) {
  MetalSpeciesFaceResponse f;f.rate={factor*(.1+a.X[0]-b.X[0]),factor*.3*(a.X[1]-b.X[1]),factor*(.002+a.Z()-b.Z())};
  f.dleft[0][0]=factor;f.dright[0][0]=-factor;
  f.dleft[1][1]=.3*factor;f.dright[1][1]=-.3*factor;
  f.dleft[2][2]=factor;f.dright[2][2]=-factor;return f;
}
}
int main() {
  int checks=0,failures=0;double max_flux=0,max_source=0,max_balance=0,max_energy=0,max_stiff_mismatch=0;
  auto check=[&](bool ok,const char* label,double error=0){++checks;if(!ok){++failures;std::printf("FAIL %s %.4g\n",label,error);}};
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  auto old=solar_scaled(.4,.02);old.basis=AbundanceBasis::baryon_mass;old.metal_inventory=MetalInventory::gs98;
  old.X[1]=.01;old.X[2]-=.01;const auto initial=initial_gs98_cn(old);
  old.cn_molality=CNAbundances{.3*initial[0],.05*initial[0],initial[2]+.65*initial[0]};
  const auto before=cn_physical_ledger(old,*old.cn_molality);const auto converted=explicit_cn_material(old);
  const auto after=cn_physical_ledger(converted,*converted.cn_molality);
  check(before.hydrogen==after.hydrogen && before.helium3==after.helium3 && before.helium4==after.helium4,"material conversion preserves hydrogen and helium");
  check(before.molality==after.molality && std::abs(before.metal_fraction-after.metal_fraction)<1e-16,"material conversion preserves actual CN and total metals");
  check(explicit_cn_material(converted)==converted,"explicit composition conversion is idempotent");
  for(double z:{1e-50,1e-32,1e-12,.1}) {
    const CNAbundances y{z/12,0,0};
    const double rounded_low=std::nextafter(12*y[0],0.);
    check(cn_inert_metal_fraction(rounded_low,y)==0,"CN zero remainder tolerates one rounding step");
    check(rejects([&]{cn_inert_metal_fraction(z*(1-1e-12),y);}),"actual CN excess remains rejected");
    check(rejects([&]{cn_inert_metal_fraction(0,y);}),"nonzero CN cannot inhabit zero total metals");
    const double inert=cn_inert_metal_fraction(z,CNAbundances{z/24,0,0});
    check(std::abs(inert/z-.5)<1e-15,"positive inert-metal remainder unchanged");
  }
  // Trace He4 remains representable even when sum(independent) rounds to one.
  {
    MetalCNVector v{};v[0]=std::nextafter(1.,0.);v[1]=1-v[0]-1e-20;
    const auto c=metal_cn_composition(converted,v);
    check(c.X[2]>0 && c.X[2]<2e-20,"positive trace helium is not lost in the unit sum",c.X[2]);
    v[1]=1-v[0]+1e-20;
    check(rejects([&]{metal_cn_composition(converted,v);}),"negative trace helium is rejected without a floor");
  }
  auto left=composition(.2,.1,.7,.6),right=composition(.55,.8,1.3,1.4);
  const auto midpoint=mean_composition(left,right);
  const auto physical_left=metal_cn_abundances(left),physical_right=metal_cn_abundances(right);
  const auto physical_midpoint=metal_cn_abundances(midpoint);
  for(std::size_t k=0;k<METAL_CN_SIZE;++k)
    check(std::abs(physical_midpoint[k]-.5*(physical_left[k]+physical_right[k]))<2e-16,
        "face average preserves the actual isotope and inert-metal inventories");
  check(midpoint==mean_composition(right,left),"face composition is independent of endpoint order");
  auto incompatible=right;incompatible.cn_molality.reset();
  check(rejects([&]{mean_composition(left,incompatible);}),"face average rejects inconsistent CN inventory");
  incompatible=right;incompatible.cn_mass_convention=CNMassConvention::fixed_metal_proxy;
  check(rejects([&]{mean_composition(left,incompatible);}),"face average rejects mixed metal conventions");
  for(double sign:{-1.,1.}) {
    const auto base=material_flux(left,right,sign);const auto f=common_metal_cn_flux(base,left,right,true);
    double metals=0;for(std::size_t row=2;row<6;++row)metals+=f.rate[row];
    check(std::abs(metals-base.rate[2])<1e-15,"metal group flux partitions exactly among CN and inert mass");
    for(std::size_t side=0;side<2;++side)for(std::size_t col=0;col<6;++col) {
      const auto at=metal_cn_abundances(side?right:left);const double h=1e-5*std::max(at[col],1e-5);
      auto up=at,dn=at;up[col]+=h;dn[col]-=h;auto a=left,b=left,c=right,d=right;
      if(side){c=metal_cn_composition(right,up);d=metal_cn_composition(right,dn);}
      else{a=metal_cn_composition(left,up);b=metal_cn_composition(left,dn);}
      const auto plus=common_metal_cn_flux(material_flux(a,c,sign),a,c,false),minus=common_metal_cn_flux(material_flux(b,d,sign),b,d,false);
      for(std::size_t row=0;row<6;++row) {
        const double exact=side?f.dright[row][col]:f.dleft[row][col];
        const double e=std::abs((plus.rate[row]-minus.rate[row])/(2*h)-exact)/std::max(std::abs(exact),1e-5);
        max_flux=std::max(max_flux,e);check(e<2e-5,"six-species flux derivative",e);
      }
    }
  }
  PPCNNetwork nuclear;const double T=1.4e7,rho=300;
  const auto source=detail::metal_cn_source(T,rho,left,nuclear.pp(),nuclear.cn());
  for(std::size_t col=0;col<6;++col) {
    const auto at=metal_cn_abundances(left);const double h=1e-5*std::max(at[col],1e-5);auto up=at,dn=at;up[col]+=h;dn[col]-=h;
    const auto a=detail::metal_cn_source(T,rho,metal_cn_composition(left,up),nuclear.pp(),nuclear.cn());
    const auto b=detail::metal_cn_source(T,rho,metal_cn_composition(left,dn),nuclear.pp(),nuclear.cn());
    for(std::size_t row=0;row<6;++row) {
      const double exact=source.jacobian[row][col];const double e=std::abs((a.source[row]-b.source[row])/(2*h)-exact)/std::max(1e-24,std::abs(exact));
      max_source=std::max(max_source,e);check(e<2e-5,"physical CN source derivative",e);
    }
  }
  const auto mapped=nuclear.eval(T,rho,left);double zrate=0,sum=0;
  for(std::size_t k=0;k<NSPEC;++k){sum+=mapped.dXdt[k];if(k>=3)zrate+=mapped.dXdt[k];}
  check(std::abs(sum)<1e-12*std::abs(mapped.dXdt[0]),"mapped physical nuclear source conserves baryons");
  check(std::abs(zrate-(source.source[2]+source.source[3]+source.source[4]))<1e-28,"nuclear captures update actual total metal mass");
  Model model;model.M=10;model.m={1,4,10};model.comp={left,composition(.35,.5,1),right};
  model.y={{0,0,0,1},{std::log(2.),0,0,1},{std::log(3.),0,0,1}};const auto weights=nodal_mass_weights(model);
  for(bool warm:{false,true})for(int mode=0;mode<4;++mode) {
    for(std::size_t i=0;i<3;++i){model.y[i].lnT=std::log(warm?1.2e7+1e6*static_cast<double>(i):1e4);model.y[i].lnrho=std::log(warm?300.:1.);}
    const double dt=warm?1e13:2.;const MixingRegions regions=mode==2?MixingRegions{{0,2},{2,3}}:MixingRegions{{0,1},{1,2},{2,3}};
    const MetalCNFlux flux=[&](std::size_t,const Composition& a,const Composition& b,bool derivatives) {
      return mode?common_metal_cn_flux(material_flux(a,b,1/dt),a,b,derivatives):MetalCNFaceResponse{};
    };
    SpeciesTransportOptions options;options.abundance_tolerance=1e-14;
    const std::vector<double> mixing(2,mode==1?.3/dt:(mode==3?1e12/dt:0.));
    const auto result=burn_metal_cn_and_diffuse(model,model,nuclear,regions,flux,dt,options,mixing);
    if(warm && mode==0) {
      auto separate=options;separate.abundance_tolerance=1e-6;
      const auto loose=burn_metal_cn_and_diffuse(model,model,nuclear,regions,flux,dt,separate,mixing);
      double loose_balance=0;for(double x:loose.integrated_balance)loose_balance=std::max(loose_balance,std::abs(x));
      check(loose_balance>1e-14,"independent-balance control must exercise a loose conservation residual",loose_balance);
      separate.integrated_balance_tolerance=1e-14;
      const auto bounded=burn_metal_cn_and_diffuse(model,model,nuclear,regions,flux,dt,separate,mixing);
      for(double x:bounded.integrated_balance)check(std::abs(x)<=1e-14,
          "looser local Newton accuracy must preserve the selected integrated balance",x);
      check(bounded.abundance_correction<=separate.abundance_tolerance,"independent local correction bound");
    }
    auto current=model;current.comp=result.composition;const auto reconstruction=reconstruct_metal_fluxes(current,model,nuclear,regions,result.boundary_fluxes,dt);
    MetalCNVector total{};double rest=0,heat=0;
    for(std::size_t i=0;i<3;++i) {
      const auto a=metal_cn_abundances(model.comp[i]),b=metal_cn_abundances(current.comp[i]);
      const auto reaction=detail::metal_cn_source(current.T(i),current.rho(i),current.comp[i],nuclear.pp(),nuclear.cn());
      for(std::size_t k=0;k<6;++k) {check(b[k]>=0,"positive physical abundance");total[k]+=weights[i]/model.M*(b[k]-a[k]-dt*reaction.source[k]);}
      check(current.comp[i].X[2]>=0,"positive physical helium");
      double delta=nuclear.rest_energy_correction(current.comp[i])-nuclear.rest_energy_correction(model.comp[i]);
      for(std::size_t k=0;k<NSPEC;++k)delta+=(current.comp[i].X[k]-model.comp[i].X[k])*(nuclides[k].A/mass_numbers[k]-1)*constants::c*constants::c;
      rest-=weights[i]*delta/dt;const auto n=nuclear.eval(current.T(i),current.rho(i),current.comp[i]);heat+=weights[i]*(n.eps+n.eps_neutrino);
    }
    for(double e:total){max_balance=std::max(max_balance,std::abs(e));check(std::abs(e)<3e-14,"global physical isotope and inert-metal balance",e);}
    for(const auto& cell:reconstruction.cell_balances)for(double e:cell) {
      check(std::abs(e)<3e-14,"reconstructed material continuity including extreme mixing",e);
      if(mode==3)max_stiff_mismatch=std::max(max_stiff_mismatch,std::abs(e));
    }
    if(warm){const double e=std::abs(rest/heat-1);max_energy=std::max(max_energy,e);check(e<2e-6,"nuclear heat agrees with changing physical rest mass",e);}
    if(mode)check(current.comp.front().Z()!=model.comp.front().Z(),"total metal mass actually redistributed");
    for(auto [a,b]:regions)for(std::size_t i=a+1;i<b;++i)check(current.comp[i]==current.comp[a],"mixed composition homogeneous");
  }
  // Seeding absent metals in a hydrogen-rich cell must not create negative He4.
  {auto seed_model=model;
   auto depleted=metal_cn_abundances(seed_model.comp.front());depleted[0]=1e-8;depleted[1]=1e-14;
   seed_model.comp.front()=metal_cn_composition(left,depleted);
   const MetalCNVector surface{.99,.0099999,0,0,0,0,0};
   seed_model.comp.back()=metal_cn_composition(right,surface);
   for(auto& y:seed_model.y){y.lnT=std::log(1e4);y.lnrho=0;}
   const MixingRegions regions{{0,1},{1,2},{2,3}};
   double maximum_core_H=0;
   const MetalCNFlux zero=[&](std::size_t face,const Composition& a,const Composition&,bool){
     if(face==0)maximum_core_H=std::max(maximum_core_H,a.X[0]);return MetalCNFaceResponse{};
   };
   SpeciesTransportOptions options;options.seed_present_species=true;
   try {
     const auto result=burn_metal_cn_and_diffuse(seed_model,seed_model,nuclear,regions,zero,1.,options);
     double error=0;
     for(std::size_t i=0;i<seed_model.size();++i) {
       const auto start=metal_cn_abundances(seed_model.comp[i]),next=metal_cn_abundances(result.composition[i]);
       for(std::size_t k=0;k<METAL_CN_SIZE;++k)error=std::max(error,std::abs(next[k]-start[k]));
       check(result.composition[i].X[2]>=0,"seeded solve preserves nonnegative helium");
     }
     check(error<1e-12,"positive starting guess adds no material to accepted state",error);
     check(maximum_core_H<=1.011e-8,"starting guess does not refill depleted core hydrogen",maximum_core_H);
   } catch(const std::exception& e) {
     std::printf("seed regression: %s\n",e.what());
     check(false,"positive starting guess reaches the unchanged state");
   }
  }
  // A uniform mixed region must remain uniform in the initial iterate,
  // including trace reference helium. Unequal nodal weights exercise rounding.
  {
    auto uniform=model;MetalCNVector v{};v[0]=std::nextafter(1.,0.);v[1]=(1-v[0])*.9;
    const auto c=metal_cn_composition(left,v);uniform.comp.assign(uniform.size(),c);
    for(auto& y:uniform.y){y.lnT=std::log(1e4);y.lnrho=0;}
    SpeciesTransportOptions options;options.initial_guess=uniform.comp;
    const MetalCNFlux zero=[&](std::size_t,const Composition& a,const Composition& b,bool){
      if(a.X[2]<=0 || b.X[2]<=0)throw std::domain_error("test: lost trace helium");
      return MetalCNFaceResponse{};
    };
    try {
      const auto result=burn_metal_cn_and_diffuse(uniform,uniform,nuclear,{{0,2},{2,3}},zero,1.,options);
      bool same=true;for(const auto& next:result.composition)same&=next.X==c.X;
      check(same,"uniform warm start preserves trace composition");
    }catch(const std::exception& e){std::printf("trace warm start: %s\n",e.what());check(false,"uniform trace warm start succeeds");}
  }
  // Steep gradients in vanishing metals make the chemical-potential
  // derivative large even though the actual mass flux is negligible.
  // Keep every trace species, including its isotope ratios, without a floor.
  {
    auto trace=model;
    for(std::size_t i=0;i<trace.size();++i) {
      auto v=metal_cn_abundances(left);const double z=i==1?1e-48:1e-30;
      double old_z=0;for(std::size_t k=2;k<6;++k)old_z+=v[k];
      for(std::size_t k=2;k<6;++k)v[k]*=z/old_z;
      v[0]=.99;v[1]=1e-4;trace.comp[i]=metal_cn_composition(left,v);
      trace.y[i].lnT=std::log(1e4);trace.y[i].lnrho=0;
    }
    const MetalCNFlux flux=[](std::size_t,const Composition& a,const Composition& b,bool derivatives) {
      const double za=a.Z(),zb=b.Z(),gradient=std::log(za/zb),mean=.5*(za+zb);
      MetalSpeciesFaceResponse f;f.rate[2]=mean*gradient;f.rate[0]=-f.rate[2];
      f.dleft[2][2]=.5*gradient+mean/za;f.dright[2][2]=.5*gradient-mean/zb;
      f.dleft[0][2]=-f.dleft[2][2];f.dright[0][2]=-f.dright[2][2];
      return common_metal_cn_flux(f,a,b,derivatives);
    };
    try {
      const auto result=burn_metal_cn_and_diffuse(trace,trace,nuclear,{{0,1},{1,2},{2,3}},flux,1.);
      for(std::size_t i=0;i<trace.size();++i) {
        const auto x=metal_cn_abundances(result.composition[i]);
        for(std::size_t k=2;k<6;++k)check(x[k]>0 && x[k]<1e-28,"trace metals retained without a floor");
      }
      for(double b:result.integrated_balance)check(std::abs(b)<1e-14,"trace chemical-potential flux conserves material",b);
    } catch(const std::exception& e) {std::printf("trace gradient: %s\n",e.what());check(false,"extreme trace-metal gradient has a nonsingular solve");}
  }
  const auto path=std::filesystem::temp_directory_path()/"ember-metal-cn-checkpoint-test.restart";
  check(!std::filesystem::exists(path),"checkpoint test uses a new file");
  driver::Checkpoint state{model,1.,1,0};driver::Selections selections{"metal-cn","a","b","c","d"};driver::Identities identities{{"executable","test-only"}};
  driver::write_checkpoint(path,state,selections,1.,identities);
  const auto restored=driver::read_checkpoint(path,3,model.M,model.comp.front(),selections,1.,identities);
  check(restored.model.comp==model.comp,"variable-metal physical checkpoint round trip");
  check(rejects([&]{driver::read_checkpoint(path,3,model.M,old,selections,1.,identities);}),"physical checkpoint rejects proxy convention");
  std::filesystem::remove(path);
  check(rejects([&]{burn_cn_and_diffuse(model,model,nuclear,{{0,3}},[](auto,const auto&,const auto&,bool){return CNSpeciesFaceResponse{};},1.);}),"old solver rejects physical metal convention");
  std::printf("%d checks, %d failures; flux %.4g source %.4g balance %.4g mass-energy %.4g extreme-mixing-local-mismatch %.4g\n",checks,failures,max_flux,max_source,max_balance,max_energy,max_stiff_mismatch);
  return failures?1:0;
}
