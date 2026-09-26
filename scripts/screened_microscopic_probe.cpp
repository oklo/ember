#include "ember/screened_microscopic_transport.hpp"
#include "ember/material_flux.hpp"
#include "ember/eos_component.hpp"
#include "ember/constants.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {
template<std::size_t N> void array(const std::array<double,N>& a) {
  std::cout<<'[';for(std::size_t i=0;i<N;++i) {if(i)std::cout<<',';std::cout<<a[i];}std::cout<<']';
}
template<std::size_t N,std::size_t M> void matrix(const std::array<std::array<double,M>,N>& a) {
  std::cout<<'[';for(std::size_t i=0;i<N;++i) {if(i)std::cout<<',';array(a[i]);}std::cout<<']';
}
ember::Composition composition(double x,double y) {
  auto c=ember::solar_scaled(x,.02);c.basis=ember::AbundanceBasis::baryon_mass;
  c.metal_inventory=ember::MetalInventory::gs98;c.X[1]=y;c.X[2]-=y;return c;
}
}
int main(int argc,char** argv) {
  if(argc!=3)return 2;
  using namespace ember;
  try {
    const auto start=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);
    ScreenedMicroscopicTransport electron_screen(eos,collisions,false,2e6),ion_screen(eos,collisions,true,2e6);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line)) {
      std::istringstream in(line);std::string mode;in>>mode;
      try {
        if(mode=="regular") {
          double x,y,T,rho;if(!(in>>x>>y>>T>>rho))return 2;
          const auto r=eos.regular_composition_potential(T,rho,composition(x,y));
          std::cout<<"{\"phi\":"<<r.phi<<",\"gradient\":";array(r.gradient);
          std::cout<<",\"gradient_T\":";array(r.dgradient_dlnT);
          std::cout<<",\"gradient_rho\":";array(r.dgradient_dlnRho);
          std::cout<<",\"hessian\":";matrix(r.hessian);std::cout<<"}\n";
        } else if(mode=="face") {
          int ions,derivatives,reference;double ma,mb,xa,ya,xb,yb;Point a,b;
          if(!(in>>ions>>derivatives>>reference>>ma>>mb>>a.lnr>>a.lnrho>>a.lnT>>a.L>>xa>>ya
              >>b.lnr>>b.lnrho>>b.lnT>>b.L>>xb>>yb))return 2;
          const auto ca=composition(xa,ya),cb=composition(xb,yb);
          const auto begin=std::chrono::steady_clock::now();
          const auto r=microscopic_face(ions?ion_screen:electron_screen,0,ma,mb,a,ca,b,cb,derivatives);
          const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
          // Construct the independent full 3x3 material flux before writing a
          // reply, so a failed reference cannot leave malformed partial JSON.
          MaterialFaceFlux f;CollisionTransportResponse k;double stiffness=0,length=0;
          if(reference) {
            const double T=std::exp(.5*(a.lnT+b.lnT)),rho=.5*(std::exp(a.lnrho)+std::exp(b.lnrho));
            const auto c=composition(.5*(xa+xb),.5*(ya+yb));
            const auto e=ElectronGas{}.eval(T,rho,c);stiffness=e.dP_dlnRho/(constants::NA*rho*c.mu_elec_inv());
            length=ScreenedCollisionTransport::screening_length(T,rho,c.h1(),c.X[1],c.Z(),stiffness,ions);
            k=collisions.eval(T,rho,c.h1(),c.X[1],c.Z(),length);
            const auto h=eos.composition_heat(T,rho,c,{true,true},false);
            const auto full=full_material_mobility(k.mobility,h.exchange_enthalpy,k.energy_scale,1e15,1e8);
            const auto pa=eos.composition_potential(std::exp(a.lnT),std::exp(a.lnrho),ca);
            const auto pb=eos.composition_potential(std::exp(b.lnT),std::exp(b.lnrho),cb);
            const double radius=.5*(std::exp(a.lnr)+std::exp(b.lnr));
            f=material_face_flux(full,pa.gradient,pb.gradient,std::exp(a.lnT),std::exp(b.lnT),rho,4*M_PI*radius*radius,mb-ma,1e15,1e8);
          }
          std::cout<<"{\"rate\":";array(r.species.rate);std::cout<<",\"dleft\":";matrix(r.species.dleft);
          std::cout<<",\"dright\":";matrix(r.species.dright);
          std::cout<<",\"carried\":"<<r.carried_luminosity<<",\"conductivity\":"<<r.conductivity;
          std::cout<<",\"dcarried_lo\":";array(r.dcarried_lo);std::cout<<",\"dcarried_hi\":";array(r.dcarried_hi);
          std::cout<<",\"dconductivity_lo\":";array(r.dconductivity_lo);std::cout<<",\"dconductivity_hi\":";array(r.dconductivity_hi);
          std::cout<<",\"seconds\":"<<seconds;
          if(reference) {
            std::cout<<",\"reference_rate\":";array(f.species_rate);
            std::cout<<",\"reference_carried\":"<<f.carried_luminosity<<",\"conductive\":"<<f.conductive_luminosity
              <<",\"material\":"<<f.material_luminosity<<",\"entropy\":"<<f.entropy_production
              <<",\"stiffness\":"<<stiffness<<",\"length\":"<<length<<",\"eta\":"<<k.eta<<",\"B\":"<<k.b_thermal;
          }
          std::cout<<"}\n";
        }else return 2;
      }catch(const std::exception& e) {std::cout<<"{\"error\":\""<<e.what()<<"\"}\n";}
    }
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
