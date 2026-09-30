// Checks of the cold dense-He base and its optional join to the production EOS.
// References: Ember's Fermi-Dirac IdealEos (electrons) and the pinned eos22.f EXCOR7 (exchange-correlation,
// values below from EOSFI22 minus its separately probed ion and screening parts). The join checks need the
// production family: set EMBER_COLD_HELIUM_FAMILY to its path, otherwise they are skipped.
#include "ember/cold_helium_base.hpp"
#include "ember/constants.hpp"
#include "ember/eos.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/ion_quantum.hpp"
#include "ember/detail/ion_ocp_components.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* what,double v){if(!ok){++failures;std::cerr<<"FAIL "<<what<<' '<<v<<'\n';}}
Composition comp(double X,double Y3,double Z){auto c=solar_scaled(X,Z);c.basis=AbundanceBasis::baryon_mass;
  c.metal_inventory=MetalInventory::gs98;c.X[1]=Y3;c.X[2]-=Y3;return c;}
}
int main(){
  {ColdHeliumOptions o;
   for(int kind=0;kind<3;++kind){auto c=comp(.01,1e-4,.02);
     if(kind==0){c.X[0]=-1e-4;c.X[2]+=.0101;}
     if(kind==1)c.X[2]=std::numeric_limits<double>::quiet_NaN();
     if(kind==2)c.X[2]-=.01;
     bool rejected=false;
     try{cold_helium_validate(6e5,2e4,c,o);}catch(const std::domain_error&){rejected=true;}
     check(rejected,"invalid composition rejected",kind);}}
  const double per=4*constants::amu/constants::kB;
  // 1. Electrons against exact Fermi-Dirac electrons (IdealEos minus ideal ions and radiation).
  {IdealEos ideal;double wp=0,wc=0;
   for(double rho:{3e3,1e4,1e5})for(double T:{3e5,4.5e5,6e5}){
     auto c=comp(0,0,0);const auto j=cold_helium_component_jet(0,T,rho,c);const auto s=ideal.eval(T,rho,c);
     const double ni=rho/(4*constants::amu),ar=constants::a_rad;
     const double Pe=s.P-ni*constants::kB*T-ar*std::pow(T,4)/3,Ee=s.E-1.5*ni*constants::kB*T/rho-ar*std::pow(T,4)/rho;
     const double cve=s.cv-1.5*ni*constants::kB/rho-4*ar*std::pow(T,3)/rho;
     wp=std::max({wp,std::abs(rho*T*j[0][1]/Pe-1),std::abs(-T*j[1][0]/Ee-1)});
     wc=std::max(wc,std::abs((-j[1][0]-j[2][0])-cve)*per);}
   check(wp<1e-6,"electron P and E vs Fermi-Dirac",wp);
   check(wc<1e-3,"electron heat capacity vs Fermi-Dirac (k per ion)",wc);}
  // 2. Exchange-correlation per ion (k units) vs EOSFI22: T, rho, f, u, p.
  {const double ref[][5]={{1.5e5,3e3,-16.439430,-16.440649,-5.367333},{3e5,1e4,-12.189379,-12.190182,-3.998982},
     {6e5,1e5,-13.005405,-13.005587,-4.293405}};double w=0;
   for(const auto& r:ref){const auto j=cold_helium_component_jet(1,r[0],r[1],comp(0,0,0));const double s=4/constants::R_gas;
     w=std::max({w,std::abs(j[0][0]*s-r[2]),std::abs(-j[1][0]*s-r[3]),std::abs(j[0][1]*s-r[4])});}
   // Published eq. A.1 constants and exact theta = T/T_F vs eos22.f single-precision literals (0.543 r_s/Gamma_e): ~1e-5 k.
   check(w<3e-5,"exchange-correlation vs EXCOR7",w);}
  // 3. Third-order and composition consistency of the base.
  {const auto c=comp(2e-4,2e-6,1e-5);double w3=0,wx=0;
   for(double T:{3e5,4e5,6e5}){const double rho=1.5e4,h=1e-4;
     const auto m=cold_helium_material_jets(T,rho,c,10),tp=cold_helium_material_jets(T*std::exp(h),rho,c,10),
         tm=cold_helium_material_jets(T*std::exp(-h),rho,c,10),rp=cold_helium_material_jets(T,rho*std::exp(h),c,10),
         rm=cold_helium_material_jets(T,rho*std::exp(-h),c,10);
     const double sc=std::abs(m[0][0][0]);
     for(int i=0;i<3;++i)for(int j=0;i+j<=2;++j)if(i+j==2){
       w3=std::max(w3,std::abs((tp[0][i][j]-tm[0][i][j])/(2*h)-m[0][i+1][j])/sc);
       w3=std::max(w3,std::abs((rp[0][i][j]-rm[0][i][j])/(2*h)-m[0][i][j+1])/sc);}
     const double hx=1e-6;const double fd=(cold_helium_material_jets(T,rho,comp(2e-4+hx,2e-6,1e-5),1)[0][0][0]
         -cold_helium_material_jets(T,rho,comp(2e-4-hx,2e-6,1e-5),1)[0][0][0])/(2*hx);
     wx=std::max(wx,std::abs(fd-m[1][0][0])/std::abs(m[1][0][0]));}
   check(w3<1e-7,"base third-order entries vs differences",w3);
   check(wx<1e-6,"base H composition gradient vs difference",wx);}
  // 5. Guards.
  {auto throws=[](auto f){try{f();return false;}catch(const std::domain_error&){return true;}};
   check(throws([&]{cold_helium_material_jets(3e5,1.5e4,comp(.06,0,0),1);}),"composition gate",0);
   check(throws([&]{cold_helium_material_jets(3e5,500,comp(0,0,0),1);}),"pressure-ionization density floor",0);
   // The Sommerfeld electron approximation is never used silently outside kT/eps_F <= 0.05.
   check(throws([&]{cold_helium_material_jets(2.5e6,1e3,comp(0,0,0),1);}),"electron degeneracy guard",0);}
  // 8. Join (needs the production family).
  if(const char* path=std::getenv("EMBER_COLD_HELIUM_FAMILY")){
    VariableMetalHelmholtzEos production(path,HelmholtzTableEos::Mixture::allow_documented_proxy,
        VariableMetalHelmholtzEos::LowMetalInterpolation::quadratic,{},true);
    const ColdHeliumTable table=[&](double T,double r,const Composition& x,std::size_t n){return production.material_jets(T,r,x,n);};
    const auto c=comp(6.5e-5,3.3e-7,1.2e-4);const double rho=1.4e4;ColdHeliumOptions o;
    auto b=cold_helium_material_jets(o.join_cold,rho,c);const auto g=cold_helium_alignment_jets(table,o.join_cold,rho,c);
    const auto jc=cold_helium_joined_jets(table,o.join_cold,rho,c),t2=table(o.join_hot,rho,c,10),
        jh=cold_helium_joined_jets(table,o.join_hot*(1-1e-12),rho,c);
    double wedge=0,wanchor=0;
    for(int k=0;k<10;++k)for(int i=0;i<4;++i)for(int j=0;i+j+(k==0?0:k<4?1:2)<=3;++j){
      wedge=std::max(wedge,std::abs(jc[k][i][j]-b[k][i][j]-g[k][i][j])/(1+std::abs(jc[k][i][j])));
      wedge=std::max(wedge,std::abs(jh[k][i][j]-t2[k][i][j])/(1+std::abs(t2[k][i][j])));}
    for(double T:{o.join_cold,o.join_hot}){auto bb=cold_helium_material_jets(T,rho,c);
      const auto gg=cold_helium_alignment_jets(table,T,rho,c);const auto tt=table(T,rho,c,10);
      for(int k=0;k<10;++k)for(int j=0;j+(k==0?0:k<4?1:2)<=3;++j)
        wanchor=std::max(wanchor,std::abs(bb[k][0][j]+gg[k][0][j]-tt[k][0][j])/(1+std::abs(tt[k][0][j])));}
    check(wedge<1e-6,"join continuity at window edges",wedge);
    check(wanchor<1e-9,"alignment reproduces the table at both anchors",wanchor);
    // Default (production quantum term on trace H): supported down to T_p,H/T = 4, about 2.3e5 K at 2.8e4 g/cc.
    int bad=0;for(double r:{7e3,1.4e4,2.8e4})for(int i=0;i<=600;++i){const double T=2.5e5*std::pow(2e6/2.5e5,i/600.);
      try{(void)helmholtz_response(T,r,cold_helium_joined_jets(table,T,r,c,1)[0]);}catch(const std::exception&){++bad;}}
    check(bad==0,"joined responses stable from 2.5e5 to 2e6 K",bad);
    // Composition join in X: jets (all channels) vs central differences of lower-order entries, inside both blends.
    ColdHeliumOptions mix;
    {double wfd=0;const double h=1e-6;
     for(double T:{4.3e5,6e5})for(double X:{.015,.03,.045}){
       auto at=[&](double x,double y3,double z){return cold_helium_joined_jets(table,T,3e3,comp(x,y3,z),10,mix);};
       const double y=2e-4,z=1e-4,hz=1e-6;const auto m=at(X,y,z),xp=at(X+h,y,z),xm=at(X-h,y,z),zp=at(X,y,z+hz),zm=at(X,y,z-hz);
       auto rel=[&](double fd,double v){return std::abs(fd-v)/(1e-2*std::abs(m[0][1][0])+std::abs(v));};
       for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j<=2;++j){
         wfd=std::max(wfd,rel((xp[0][i][j]-xm[0][i][j])/(2*h),m[1][i][j]));
         if(i+j<=1){wfd=std::max(wfd,rel((xp[1][i][j]-xm[1][i][j])/(2*h),m[4][i][j]));
                    wfd=std::max(wfd,rel((zp[1][i][j]-zm[1][i][j])/(2*hz),m[6][i][j]));}}}
     check(wfd<1e-5,"X-joined composition channels vs differences",wfd);}
    // The optional EOS join: one potential for eval, forces and heat; production unchanged at w=0.
    {const ColdHeliumOptions jo=mix;VariableMetalHelmholtzEos joined(path,HelmholtzTableEos::Mixture::allow_documented_proxy,
        VariableMetalHelmholtzEos::LowMetalInterpolation::quadratic,{},true,jo);
     const auto hot=comp(.06,1e-3,1e-6);double same=0;
     for(auto [T,r]:{std::pair{4.2e5,2.4e3},std::pair{9e5,2e4}}){
       const auto a=joined.eval_with_derivatives(T,r,hot).state,reference=production.eval_with_derivatives(T,r,hot).state;
       same=std::max({same,std::abs(a.P-reference.P),std::abs(a.cv-reference.cv)});}
     check(same==0,"joined EOS identical to production at zero join weight",same);
     const auto core=comp(1e-5,1e-8,.01);const auto s=joined.eval_with_derivatives(3.5e5,2e4,core).state;
     check(s.P>0&&s.cv>0,"joined EOS evaluates below the table window",s.cv);
     const auto pot=joined.composition_potential(3.5e5,2e4,core);const auto heat=joined.composition_heat(3.5e5,2e4,core);
     check(std::isfinite(pot.gradient[0]+pot.hessian[0][0]+heat.exchange_enthalpy[0]+heat.enthalpy_partials[2][4]),
         "joined forces and heat finite",0);}
  } else std::cout<<"join checks skipped (EMBER_COLD_HELIUM_FAMILY not set)\n";
  std::cout<<"cold helium base: failures "<<failures<<'\n';
  return failures?1:0;
}
