#include "ember/cold_helium_base.hpp"
#include <cmath>
#include <iostream>
using namespace ember;
Composition comp(double x,double y3,double z) { auto c=solar_scaled(x,z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;c.X[1]=y3;c.X[2]-=y3;return c; }
int main() {
  ColdHeliumOptions o;o.mixture_phase=true;o.liquid_continuation_gamma=200;
  auto c=comp(1e-4,1e-4,.136);
  auto at=[&](double T,double d,const Composition& x){return cold_helium_solid_response(T,d,x,o);};
  int failures=0;double worst=0;
  for(double rho:{5.06e4,1.5e3}) {
  c=comp(1e-4,1e-4,rho<3e3?.01:.136);
  const double scale=std::cbrt(rho/5.06e4);
  double lo=1.65e5*scale,hi=2.1e5*scale;
  const double target=.5*cold_helium_join_weight(lo,rho,c,o);
  while(at(lo,rho,c)[0]<target)lo/=1.2;
  while(at(hi,rho,c)[0]>target)hi*=1.2;
  for(int i=0;i<60;++i){double m=(lo+hi)/2;if(at(m,rho,c)[0]>target)lo=m;else hi=m;}
  const double mid=(lo+hi)/2;
  for(double T:{mid*.96,mid,mid*1.04}) {
    const auto w=at(T,rho,c);
    if(w[0]<0||w[0]>1)++failures;
    for(int v=0;v<5;++v) {
      const double h=v<2?1e-6:1e-7;
      auto p=c,m=c;double tp=T,tm=T,rp=rho,rm=rho;
      if(v==0){tp*=std::exp(h);tm*=std::exp(-h);}
      else if(v==1){rp*=std::exp(h);rm*=std::exp(-h);}
      else if(v<4){p.X[v-2]+=h;p.X[2]-=h;m.X[v-2]-=h;m.X[2]+=h;}
      else {p=comp(c.X[0],c.X[1],c.Z()+h);m=comp(c.X[0],c.X[1],c.Z()-h);}
      const double fd=(at(tp,rp,p)[0]-at(tm,rm,m)[0])/(2*h);
      const double e=std::abs(fd-w[v+1])/std::max(1e-6,std::abs(w[v+1]));worst=std::max(worst,e);
    }
  }
  if(worst>3e-5)++failures;
  if(at(1e6,rho,c)!=std::array<double,6>{})++failures;
  auto off=o;off.mixture_phase=false;
  if(cold_helium_solid_response(mid,rho,c,off)!=std::array<double,6>{})++failures;
  const auto hydrogen=comp(.7,0,.02);
  if(at(2e5,rho,hydrogen)!=std::array<double,6>{})++failures;
  }
  std::cout<<"phase mobility: weight derivative error "<<worst<<", failures "<<failures<<'\n';
  return failures?1:0;
}
