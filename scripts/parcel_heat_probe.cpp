#include "conditional_envelope_heat.hpp"
#include "ember/conduction_table.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace ember;
namespace {
void scalar(double x){if(std::isfinite(x))std::cout<<x;else std::cout<<"null";}
template<std::size_t N> void array(const std::array<double,N>& a){
  std::cout<<'[';for(std::size_t i=0;i<N;++i){if(i)std::cout<<',';scalar(a[i]);}std::cout<<']';
}
Composition composition(double x,double y){
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=y;c.X[2]-=y;return c;
}
void heat(const MicroscopicHeatResponse& r){
  std::cout<<"{\"Q\":";scalar(r.carried_luminosity);std::cout<<",\"K\":";scalar(r.conductivity);
  std::cout<<",\"h_total\":";array(r.total_rate_enthalpy);
  std::cout<<",\"dQlo\":";array(r.dcarried_lo);std::cout<<",\"dQhi\":";array(r.dcarried_hi);
  std::cout<<",\"dKlo\":";array(r.dconductivity_lo);std::cout<<",\"dKhi\":";array(r.dconductivity_hi);std::cout<<'}';
}
}
int main(int argc,char** argv){
  if(argc!=4)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);TabulatedConduction table(argv[3]);HotConduction conduction(table);
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line)){
      try {
        std::istringstream in(line);std::string mode;in>>mode;
        if(mode=="eos"){
          double x,y,T,rho;int a,b;if(!(in>>x>>y>>T>>rho>>a>>b))return 2;
          const auto c=composition(x,y);const auto h=eos.composition_heat(T,rho,c,{bool(a),bool(b)},true);
          std::cout<<"{\"material\":";array(h.exchange_enthalpy);std::cout<<",\"radiation\":";array(h.radiation_enthalpy);
          std::cout<<",\"d_material\":[";array(h.enthalpy_partials[0]);std::cout<<',';array(h.enthalpy_partials[1]);
          std::cout<<"],\"d_radiation\":[";array(h.radiation_enthalpy_partials[0]);std::cout<<',';array(h.radiation_enthalpy_partials[1]);
          const auto e=eos.eval(T,rho,c);std::cout<<"],\"P\":"<<e.P<<",\"rho\":"<<rho<<",\"h_bulk\":"<<e.E+e.P/rho<<'}';
        }else if(mode=="bulk"){
          double x,y,T,P;if(!(in>>x>>y>>T>>P))return 2;
          const auto c=composition(x,y);const double rho=eos.rho_from_PT(T,P,c);const auto e=eos.eval(T,rho,c);
          std::cout<<"{\"P\":"<<e.P<<",\"rho\":"<<rho<<",\"h_rad\":"<<4*constants::a_rad*T*T*T*T/(3*rho)
            <<",\"h_bulk\":"<<e.E+e.P/rho<<'}';
        }else if(mode=="face"){
          int combined,derivatives,mask;double ma,mb,x0,y0,x1,y1;Point a,b;SpeciesVector total;
          if(!(in>>combined>>derivatives>>mask>>ma>>mb>>a.lnr>>a.lnrho>>a.lnT>>a.L>>x0>>y0
               >>b.lnr>>b.lnrho>>b.lnT>>b.L>>x1>>y1>>total[0]>>total[1]))return 2;
          const auto ca=composition(x0,y0),cb=composition(x1,y1);
          MicroscopicHeatResponse old,corrected;SpeciesVector micro{};
          if(combined){
            ConditionalEnvelopeHeat p0(eos,collisions,conduction,true),p1(eos,collisions,conduction,true,2e6,3e6,true);
            old=p0.heat_with_total_species_rate(0,ma,mb,a,ca,b,cb,total,derivatives);
            corrected=p1.heat_with_total_species_rate(0,ma,mb,a,ca,b,cb,total,derivatives);
          }else{
            ScreenedMicroscopicTransport p0(eos,collisions,true,2e6,{bool(mask&1),bool(mask&2)}),
                p1(eos,collisions,true,2e6,{bool(mask&1),bool(mask&2)},true);
            old=p0.heat_with_total_species_rate(0,ma,mb,a,ca,b,cb,total,derivatives);
            corrected=p1.heat_with_total_species_rate(0,ma,mb,a,ca,b,cb,total,derivatives);
            micro=p0.eval(0,ma,mb,a,ca,b,cb,false).species.rate;
          }
          std::cout<<"{\"old\":";heat(old);std::cout<<",\"corrected\":";heat(corrected);
          std::cout<<",\"micro\":";array(micro);std::cout<<'}';
        }else return 2;
        std::cout<<'\n'<<std::flush;
      }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
