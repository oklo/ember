#include "ember/ion_classical.hpp"
#include "ember/detail/ion_free_energy.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
using namespace ember;
int main(){
 try {
  int failed=0;auto check=[&](bool ok,const char* why,double error=0.){
    std::printf("%s %s %.4g\n",ok?"PASS":"FAIL",why,error);failed+=!ok;
  };
  auto mixture=[](double X,double Y3,double Z){auto c=solar_scaled(X,Z);
    c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=Y3;c.X[2]-=Y3;return c;};
  // Independent MESA 26.04.1 Skye OCP function: f and its first three
  // derivatives with respect to ln Gamma. No stellar tables are required.
  const double references[][5]={
    {0.01, -0.00057591941205951719, -0.00086191316453028269, -0.0012876876839579153, -0.0019180294867539225},
    {0.10000000000000001, -0.017582449199645335, -0.025682144831328097, -0.036886478358032139, -0.051631218113886539},
    {1, -0.43729308096885899, -0.5720432858463973, -0.70917038376755892, -0.8232972009108015},
    {10, -7.1047770053801118, -7.9987356884429639, -8.5830733871656903, -8.8608720227368103},
    {100, -84.332878961573826, -87.525993240956282, -89.053978258375722, -89.634192336087011},
    {175, -150.28885731032946, -154.4342243963259, -156.32212596585623, -157.03892784389595},
    {200, -172.36879400505595, -176.77277314780315, -178.75928515984617, -179.52042837538275}
  };
  constexpr double charge=4.80320471257e-10,rho=1e4;
  double reference_error=0;
  for(auto c:{mixture(1,0,0),mixture(0,1,0),mixture(0,0,0)}) {
    const double A=c.X[0]>0?1:(c.X[1]>0?3:4),Z=c.X[0]>0?1:2;
    const double ne=rho*(Z/A)/constants::amu,ae=std::cbrt(3/(4*std::numbers::pi*ne));
    for(const auto& r:references) {
      const double T=charge*charge*std::pow(Z,5./3)/(ae*constants::kB*r[0]);
      const auto f=ion_classical_liquid_jets(T,rho,c,1)[0];
      for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j) {
        // Gamma=200 is an intentional one-sided caloric boundary.
        if(r[0]==200 && i+j>1)continue;
        const double expected=r[1+i+j]*constants::R_gas/A*(i%2?-1:1)*std::pow(1./3,j);
        reference_error=std::max(reference_error,std::abs(f[i][j]/expected-1));
      }
    }
  }
  check(reference_error<1e-9,"pure-ion potential and thermal derivatives match independent source",reference_error);
  using J=detail::Taylor3<1>;
  const auto weak=detail::classical_ocp_free_energy(detail::exp(J::variable(std::log(1e-10),0)));
  const double dh=-std::pow(1e-10,1.5)/std::sqrt(3.);
  check(std::abs(weak.value()/dh-1)<1e-5 && std::abs(weak.derivative({1})/(1.5*dh)-1)<1e-5,
        "dilute Debye-Huckel limit is retained without cancellation");
  const auto left=detail::classical_ocp_free_energy(detail::exp(J::variable(std::log(1e-3*(1-1e-8)),0)));
  const auto right=detail::classical_ocp_free_energy(detail::exp(J::variable(std::log(1e-3*(1+1e-8)),0)));
  check(std::abs(left.value()/right.value()-1)<4e-8,"small-Gamma series joins the full fit");
  auto c=mixture(.1,.08,.13);const double T=1e6,h=1e-5;
  const auto f=ion_classical_liquid_jets(T,rho,c);
  double thermal_error=0,composition_error=0;
  for(unsigned k=0;k<10;++k){
    const unsigned order=k==0?0:(k<4?1:2);
    for(unsigned d=0;d<2;++d){
      const auto p=ion_classical_liquid_jets(T*std::exp(d==0?h:0),rho*std::exp(d==1?h:0),c);
      const auto m=ion_classical_liquid_jets(T*std::exp(d==0?-h:0),rho*std::exp(d==1?-h:0),c);
      for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j+order<3;++j){
        const double fd=(p[k][i][j]-m[k][i][j])/(2*h),exact=f[k][i+(d==0)][j+(d==1)];
        thermal_error=std::max(thermal_error,std::abs(fd-exact)/std::max(std::abs(f[0][0][0]),std::abs(exact)));
      }
    }
  }
  auto perturb=[](Composition v,unsigned k,double dx){
    if(k<2)v.X[k]+=dx;
    else {const double z=v.Z();for(unsigned i=METAL_BEGIN;i<METAL_END;++i)v.X[i]*=(z+dx)/z;}
    v.X[2]-=dx;return v;
  };
  for(unsigned a=0;a<3;++a){
    const auto p=ion_classical_liquid_jets(T,rho,perturb(c,a,h)),m=ion_classical_liquid_jets(T,rho,perturb(c,a,-h));
    for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j<=2;++j)
      composition_error=std::max(composition_error,std::abs((p[0][i][j]-m[0][i][j])/(2*h)-f[1+a][i][j])/std::abs(f[0][0][0]));
    for(unsigned b=0;b<3;++b){
      const auto lo=std::min(a,b),hi=std::max(a,b);const unsigned k=4+lo*3-lo*(lo-1)/2+hi-lo;
      for(unsigned i=0;i<2;++i)for(unsigned j=0;i+j<=1;++j)
        composition_error=std::max(composition_error,std::abs((p[1+b][i][j]-m[1+b][i][j])/(2*h)-f[k][i][j])/std::abs(f[0][0][0]));
    }
  }
  check(thermal_error<2e-7,"thermal and chemical responses differentiate the same potential",thermal_error);
  check(composition_error<2e-7,"composition derivatives include common electron-density response",composition_error);
  const auto helium=mixture(0,0,0);const double ne=rho*.5/constants::amu;
  const double Tb=charge*charge*std::pow(2.,5./3)/(std::cbrt(3/(4*std::numbers::pi*ne))*constants::kB*200);
  const auto cold=ion_classical_liquid_jets(.5*Tb,rho,helium,1)[0];
  check(std::abs(cold[1][0]+cold[2][0])<1e-12*std::abs(cold[0][0]),"outside fit range constant-entropy continuation has zero ion-correlation Cv");
  const auto warm=ion_classical_liquid_jets(Tb*(1+1e-8),rho,helium,1)[0];
  const auto cool=ion_classical_liquid_jets(Tb*(1-1e-8),rho,helium,1)[0];
  check(std::abs(warm[0][0]/cool[0][0]-1)<3e-8 && std::abs(warm[1][0]/cool[1][0]-1)<3e-8,
        "potential and first derivatives are continuous at fit boundary");
  bool refused=false;try{(void)ion_classical_liquid_jets(-T,rho,c);}catch(const std::domain_error&){refused=true;}
  check(refused,"invalid thermodynamic state refused");
  return failed?1:0;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}
