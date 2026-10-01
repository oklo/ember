#include "ember/collision_transport.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>
using namespace ember;
namespace {
void check(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
template<class F> void rejects(F f) {
  bool rejected=false;try {f();}catch(const std::domain_error&){rejected=true;}
  check(rejected,"invalid collision request was reused");
}
}
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  try {
    ScreenedCollisionTransport table(argv[1]),other(argv[1]);
    CollisionTaylorCache cache(1e-4,4,true);
    constexpr double T=5e6,rho=100,X=.4,Y3=.004,Z=.02,length=1e-9;
    const auto first=cache.derivatives(table,0,T,rho,X,Y3,Z,length);
    const auto nearby=cache.derivatives(table,0,T*std::exp(2e-5),rho,X,Y3,Z,length);
    check(cache.statistics().hits==1 && cache.statistics().verified==1,"nearby state was not reused");
    check(cache.statistics().worst_relative_error<1e-6,"Taylor response failed nearby exact comparison");
    check(nearby.value.conductivity>0 && first.value.active_species==nearby.value.active_species,"invalid response");
    cache.value(other,0,T,rho,X,Y3,Z,length);
    check(cache.statistics().misses==2,"different source table reused the old anchor");
    // A disappearing trace species must change the active subspace, even
    // though the absolute abundance change is smaller than the reuse radius.
    cache.value(table,1,T,rho,X,1e-6,Z,length);
    const auto absent=cache.derivatives(table,1,T,rho,X,0,Z,length);
    check(!absent.value.active_species[1] && !absent.defined[3],"vanished species remained active");
    const auto present=cache.derivatives(table,1,T,rho,X,1e-6,Z,length);
    check(present.value.active_species[1] && present.defined[3],"new species remained absent");
    const auto misses=cache.statistics().misses;
    cache.value(table,1,T,rho,X,5e-7,Z,length);
    check(cache.statistics().misses==misses+1,"large fractional trace change was reused");
    rejects([&]{cache.value(table,1,T,rho,X,-1e-6,Z,length);});
    // Cross the pair-table b boundary by less than the reuse radius.
    std::ifstream source(argv[1]);std::string magic;std::size_t nx,ny;
    source>>magic>>nx>>ny;std::vector<double> tx(nx),ty(ny);
    for(auto& x:tx)source>>x;for(auto& y:ty)source>>y;
    constexpr double me=9.1093837015e-28,kb=1.380649e-16,hbar=1.054571817e-27;
    const double edge=std::sqrt(std::exp(ty.back()-2e-5)*hbar*hbar/(8*me*kb*T));
    cache.value(table,2,T,rho,X,Y3,Z,edge);
    rejects([&]{cache.value(table,2,T,rho,X,Y3,Z,edge*std::exp(2e-5));});
    if(tx.back()>128.1) {
      // Independent Sommerfeld density and finite differences across the
      // change of quadrature, including a more degenerate state when covered.
      constexpr double pi=std::numbers::pi,amu=1.66053906660e-24;
      constexpr double temperature=1e6,hydrogen=.5,he3=.01;
      constexpr double Ye=hydrogen+2*he3/3+(1-hydrogen-he3)/2;
      const double length=std::sqrt(std::exp(.5*(ty.front()+ty.back()))*hbar*hbar/(8*me*kb*temperature));
      for(double eta:{127.99999,128.00001,std::min(512.,.9*tx.back())}) {
        const double F=(2./3)*std::pow(eta,1.5)+pi*pi/12/std::sqrt(eta)
            +7*std::pow(pi,4)/960/std::pow(eta,2.5)
            +31*std::pow(pi,6)/4608/std::pow(eta,4.5);
        const double ne=std::pow(2*me*kb*temperature,1.5)/(2*pi*pi*std::pow(hbar,3))*F;
        const double density=ne*amu/Ye;
        const auto r=table.bulk_metal_derivatives(temperature,density,hydrogen,he3,0,length);
        check(std::abs(r.value.eta/eta-1)<1e-10,"degenerate density integral disagrees with Sommerfeld expansion");
        constexpr double h=2e-5;
        const auto lo=table.bulk_metal_eval(temperature*std::exp(-h),density,hydrogen,he3,0,length);
        const auto hi=table.bulk_metal_eval(temperature*std::exp(h),density,hydrogen,he3,0,length);
        const double derivative=(hi.conductivity-lo.conductivity)/(2*h);
        check(std::abs(derivative-r.partials[0].conductivity)/r.value.conductivity<2e-5,
            "degenerate conductivity derivative failed across quadrature switch");
      }
    }
    std::cout<<"collision reuse: exact comparison, species changes, table identity and domain passed\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
