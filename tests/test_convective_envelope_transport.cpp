#include "ember/convective_envelope_transport.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
using namespace ember;
namespace {
int checks=0;
void check(bool b,const char* s){++checks;if(!b)throw std::runtime_error(s);}
struct Source final:MetalMicroscopicTransport {
  double offset;bool require_hot;mutable std::size_t queries{};
  Source(double a,bool b):offset(a),require_hot(b){}
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,const Point& lo,
      const Composition&,const Point& hi,const Composition&,bool)const override {
    ++queries;if(require_hot && std::min(lo.lnT,hi.lnT)<std::log(2e6))throw std::domain_error("hot domain");
    MetalMicroscopicFaceResponse s;s.species.rate={3,-1,-.5};s.species.dleft[0][0]=4;return s;
  }
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,const Point& lo,
      const Composition&,const Point& hi,const Composition&,const MetalSpeciesVector& rate,bool d)const override {
    ++queries;if(require_hot && std::min(lo.lnT,hi.lnT)<std::log(2e6))throw std::domain_error("hot domain");
    const double K=offset+std::exp(.1*lo.lnrho+.2*hi.lnT),h=offset+.7*lo.lnT+.3*hi.lnrho;
    MicroscopicHeatResponse s;s.conductivity=K;s.carried_luminosity=h*rate[0]+2*rate[1]+3*rate[2];
    s.total_rate_enthalpy={h,2};s.total_metal_rate_enthalpy=3;
    if(d){s.dconductivity_lo[1]=.1*(K-offset);s.dconductivity_hi[2]=.2*(K-offset);s.dcarried_lo[2]=.7*rate[0];s.dcarried_hi[1]=.3*rate[0];}
    return s;
  }
  const char* name()const override{return "analytic transport fixture";}
};
}
int main(){try {
  Source hot(7,true),cool(2,false);ConvectiveEnvelopeTransport provider(hot,cool,2e6,3e6);
  auto c=solar_scaled(.5,.02);c.basis=AbundanceBasis::baryon_mass;
  Point a{1,2,0,1e30},b{2,1,0,2e30};MetalSpeciesVector rate{2,3,4};double worst=0;
  for(auto pair:{std::pair{1e6,4e6},std::pair{2.2e6,2.7e6},std::pair{2.6e6,4e6},std::pair{4e6,5e6}}){
    a.lnT=std::log(pair.first);b.lnT=std::log(pair.second);
    const auto h=provider.heat_with_total_metal_rate(0,1,2,a,c,b,c,rate,true);
    for(int side=0;side<2;++side)for(std::size_t k=0;k<3;++k){
      const auto v=static_cast<Var>(k);auto ap=a,am=a,bp=b,bm=b;constexpr double step=2e-5;
      (side?bp:ap)[v]+=step;(side?bm:am)[v]-=step;
      const auto plus=provider.heat_with_total_metal_rate(0,1,2,ap,c,bp,c,rate,false),minus=provider.heat_with_total_metal_rate(0,1,2,am,c,bm,c,rate,false);
      const double dK=(plus.conductivity-minus.conductivity)/(2*step),dQ=(plus.carried_luminosity-minus.carried_luminosity)/(2*step);
      worst=std::max({worst,std::abs(dK-(side?h.dconductivity_hi:h.dconductivity_lo)[k])/std::max(1.,std::abs(dK)),std::abs(dQ-(side?h.dcarried_hi:h.dcarried_lo)[k])/std::max(1.,std::abs(dQ))});
    }
    const double linear=h.total_rate_enthalpy[0]*rate[0]+h.total_rate_enthalpy[1]*rate[1]+h.total_metal_rate_enthalpy*rate[2];
    check(std::abs(linear-h.carried_luminosity)<1e-12,"lost linear total-composition enthalpy");
  }
  check(worst<2e-7,"transition derivatives differ from finite differences");
  for(double T:{1e6,2.5e6}){
    a.lnT=b.lnT=std::log(T);bool rejected=false;
    try{(void)provider.metal_eval(0,1,2,a,c,b,c,true);}catch(const std::domain_error&){rejected=true;}
    check(rejected,"cold radiative exchange silently accepted");
  }
  a.lnT=b.lnT=std::log(4e6);const auto flux=provider.metal_eval(0,1,2,a,c,b,c,true),original=hot.metal_eval(0,1,2,a,c,b,c,true);
  check(flux.species.rate==original.species.rate && flux.species.dleft==original.species.dleft,"hot species exchange was tapered or changed");
  a.lnT=std::log(1e6);const auto before=hot.queries;
  (void)provider.heat_with_total_metal_rate(0,1,2,a,c,b,c,rate,true);
  check(hot.queries==before,"cold heat queries hot source outside its domain");
  // At each endpoint, the added derivative vanishes with the square of
  // distance. Both sides approach the unchanged physical source response.
  for(double T:{2e6,3e6}) {
    a.lnT=b.lnT=std::log(T);const auto exact=provider.heat_with_total_metal_rate(0,1,2,a,c,b,c,rate,true);
    const auto source=(T==2e6?cool:hot).heat_with_total_metal_rate(0,1,2,a,c,b,c,rate,true);
    check(exact.conductivity==source.conductivity && exact.carried_luminosity==source.carried_luminosity,"endpoint differs from source");
    check(exact.dconductivity_lo==source.dconductivity_lo && exact.dcarried_hi==source.dcarried_hi,"endpoint derivative differs from source");
  }
  std::cout<<"checks="<<checks<<" maximum_derivative_error="<<worst<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
