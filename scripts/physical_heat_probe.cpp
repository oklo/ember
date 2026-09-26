#include "ember/screened_microscopic_transport.hpp"
#include "ember/evolution.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace ember;
namespace {
template<std::size_t N> void array(const std::array<double,N>& a) {
  std::cout<<'[';for(std::size_t i=0;i<N;++i)std::cout<<(i?",":"")<<a[i];std::cout<<']';
}
template<std::size_t N> void matrix(const std::array<std::array<double,N>,N>& a) {
  std::cout<<'[';for(std::size_t i=0;i<N;++i){if(i)std::cout<<',';array(a[i]);}std::cout<<']';
}
Composition composition(double x,double y) {
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;
  c.metal_inventory=MetalInventory::gs98;c.X[1]=y;c.X[2]-=y;return c;
}
struct AnalyticRadiation final:Opacity {
  OpacityState eval(double T,double rho,const Composition&)const override {
    return {20*std::pow(T/1e7,-.8)*std::pow(rho/1e3,.3),-.8,.3};
  }
  const char* name()const override{return "analytic radiation for transport integration controls";}
};
void heat_output(const MicroscopicHeatResponse& r) {
  std::cout<<"\"carried\":"<<r.carried_luminosity<<",\"conductivity\":"<<r.conductivity;
  std::cout<<",\"dcarried_lo\":";array(r.dcarried_lo);std::cout<<",\"dcarried_hi\":";array(r.dcarried_hi);
  std::cout<<",\"dconductivity_lo\":";array(r.dconductivity_lo);
  std::cout<<",\"dconductivity_hi\":";array(r.dconductivity_hi);
}
}
int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);
    const PPChains nuclear(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
    const AnalyticRadiation opacity;
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line)) {
      try {
        std::istringstream in(line);std::string mode;int ions,derivatives,mask;
        double ma,mb,xa,ya,xb,yb;Point a,b;
        if(!(in>>mode>>ions>>derivatives>>mask>>ma>>mb>>a.lnr>>a.lnrho>>a.lnT>>a.L>>xa>>ya
            >>b.lnr>>b.lnrho>>b.lnT>>b.L>>xb>>yb))return 2;
        const auto ca=composition(xa,ya),cb=composition(xb,yb);
        ScreenedMicroscopicTransport provider(eos,collisions,ions,2e6,{bool(mask&1),bool(mask&2)});
        const auto begin=std::chrono::steady_clock::now();
        if(mode=="heat" || mode=="face") {
          MicroscopicHeatResponse r;
          if(mode=="heat")r=microscopic_heat(provider,0,ma,mb,a,ca,b,cb,derivatives);
          else r=microscopic_face(provider,0,ma,mb,a,ca,b,cb,derivatives);
          const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
          std::cout<<'{';heat_output(r);std::cout<<",\"seconds\":"<<seconds<<"}\n";
        }else if(mode=="zone") {
          Model m;m.M=mb;m.m={ma,mb};m.y={a,b};m.comp={ca,cb};
          Physics p{&eos,&opacity,&nuclear,1.9};p.microscopic=&provider;
          p.criterion=ConvectiveCriterion::ledoux;p.alpha_semiconvection=.01;p.alpha_thermohaline=.01;
          ZoneResidual r;
          if(derivatives)r=zone_residual(m,0,p,0);else r.f=zone_equations(m,0,p,0);
          const auto regions=convective_mixing_regions(m,p);
          const auto mixing=secular_mixing_diffusivities(m,p);
          auto baseline=p;baseline.microscopic=nullptr;
          const auto plain=zone_residual(m,0,baseline,0);
          bool unchanged=true;
          for(std::size_t j=0;j<3;++j)unchanged &= r.f[j]==plain.f[j];
          const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
          std::cout<<"{\"residual\":";array(r.f);
          std::cout<<",\"left\":";matrix(r.dfdy_lo);std::cout<<",\"right\":";matrix(r.dfdy_hi);
          std::cout<<",\"nonthermal_equations_unchanged\":"<<unchanged<<",\"regions\":"<<regions.size()
            <<",\"secular_diffusivity\":"<<mixing[0]<<",\"seconds\":"<<seconds<<"}\n";
        }else return 2;
        std::cout<<std::flush;
      }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
