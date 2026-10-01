// Checks of the opt-in same-composition crystallization of the cold dense-He mixture
// (ion_mixture_phase_difference_jets and ColdHeliumOptions::mixture_phase). Not selected by default.
#include "ember/cold_helium_base.hpp"
#include "ember/constants.hpp"
#include "ember/detail/ion_ocp_components.hpp"
#include "ember/ion_phase.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* what,double v){if(!ok){++failures;std::cerr<<"FAIL "<<what<<' '<<v<<'\n';}}
Composition comp(double X,double Y3,double Z){auto c=solar_scaled(X,Z);c.basis=AbundanceBasis::baryon_mass;
  c.metal_inventory=MetalInventory::gs98;c.X[1]=Y3;c.X[2]-=Y3;return c;}
ColdHeliumOptions phase(double width=.005){ColdHeliumOptions o;o.mixture_phase=true;o.phase_width=width;return o;}
// Energy and heat capacity per mass from a material F/T jet in (ln T, ln rho).
double energy(const HelmholtzJet& j,double T){return -T*j[1][0];}
double heat(const HelmholtzJet& j){return -j[1][0]-j[2][0];}
double weight(double T,double rho,const Composition& c,double width=.005){
  MixturePhaseOptions p;p.width=width;return ion_mixture_phase_weight(T,rho,c,p)[0];}
// Temperature where the solid-solution weight is one half.
double onset(double rho,const Composition& c,double width=.005){
  double lo=4e4,hi=4e5;for(int i=0;i<200;++i){const double m=std::sqrt(lo*hi);(weight(m,rho,c,width)>.5?lo:hi)=m;}
  return std::sqrt(lo*hi);}
}
int main(){
  // 1. Pure-He4 limit: identical to the reviewed pure-He soft minimum.
  {double worst=0;
   for(double rho:{1e4,5e4,1e5}){const auto c=comp(0,0,0);const double Tm=onset(rho,c);
     for(double s=.8;s<=1.25;s+=.01){const double T=Tm*s;
       const auto a=ion_mixture_phase_difference_jets(T,rho,c,1)[0],b=ion_phase_difference_jets(T,rho,c,1)[0];
       for(int i=0;i<4;++i)for(int j=0;i+j<=3;++j)
         worst=std::max(worst,std::abs(a[i][j]-b[i][j])/std::max(constants::R_gas/4,std::abs(b[i][j])));}}
   check(worst<1e-10,"pure-He4 limit equals the pure-He soft minimum (per ion, kT)",worst);
   std::printf("pure-He limit: max difference %.2e\n",worst);}
  // 2. Actual core states (932 K model, q 0, 0.2, 0.49, 0.87): heat capacity, energy integral, latent heat, stiffness.
  struct State{double rho,X,Y3,Z;};
  const State states[]={{50556.6,2.3e-8,1e-14,.13466},{33030.3,2e-7,1.6e-12,.042403},{19676.3,4.2e-6,4.1e-9,.0012127},
                        {15461.5,1.5e-5,7.1e-8,1.34e-4}};
  for(const auto& s:states){
    // Lowest temperature the base admits for the actual trace hydrogen (hydrogen quantum limit), reported;
    // the window checks below use the same metals and He3 without hydrogen.
    {const auto ch=comp(s.X,s.Y3,s.Z);double lo=2e4,hi=4e5;
     for(int i=0;i<100;++i){const double m=std::sqrt(lo*hi);bool ok=true;
       try{cold_helium_validate(m,s.rho,ch,phase());}catch(const std::domain_error&){ok=false;}(ok?hi:lo)=m;}
     std::printf("rho %.4g X %.2g: base admits T >= %.4g K; ",s.rho,s.X,hi);}
    const auto c=comp(0,s.Y3,s.Z);const double Tm=onset(s.rho,c);
    auto lopt=phase();lopt.minimum_solid_gamma=1e9;   // same limits, no phase term
    const auto liquid=[&](double T){return cold_helium_material_jets(T,s.rho,c,1,lopt)[0];};
    const auto full=[&](double T){return cold_helium_material_jets(T,s.rho,c,1,phase())[0];};
    const double lo=std::log(Tm)-.25,hi=std::log(Tm)+.25;const int n=8001;
    double cvmin=1e300,integral=0,prevT=0,prevCv=0,stiff=1e300,press=1e300;
    for(int i=0;i<n;++i){const double T=std::exp(lo+(hi-lo)*i/(n-1));const auto j=full(T);const double cv=heat(j);
      cvmin=std::min(cvmin,cv);press=std::min(press,j[0][1]);stiff=std::min(stiff,j[0][1]+j[0][2]);
      if(i)integral+=.5*(cv+prevCv)*(T-prevT);prevT=T;prevCv=cv;}
    const double Tl=std::exp(lo),Th=std::exp(hi),dU=energy(full(Th),Th)-energy(full(Tl),Tl);
    // Latent heat: jump between the liquid extrapolated from above and the solid solution extrapolated from below.
    MixturePhaseOptions p;const auto dd=[&](double T){return ion_mixture_phase_difference_jets(T,s.rho,c,1,p)[0];};
    // Solid-branch energy = liquid energy + energy of the full difference far below (weight 1): D itself.
    const double uS_lo=energy(full(Tl),Tl),uL_hi=energy(liquid(Th),Th);
    const double uL_lo=energy(liquid(Tl),Tl);
    const double latent_lo=uL_lo-uS_lo;   // liquid minus solid energy at the low edge (solid-side latent heat)
    const double lat_mid=[&]{const auto dj=dd(Tm);return -Tm*dj[1][0];}();
    (void)uL_hi;(void)lat_mid;
    const double far_hi=std::abs(dd(Th)[0][0])/std::abs(liquid(Th)[0][0]);
    check(cvmin>0,"positive heat capacity through the window",cvmin);
    check(std::abs(dU-integral)/std::abs(dU)<1e-6,"energy change equals the heat-capacity integral",std::abs(dU-integral)/std::abs(dU));
    check(press>0&&stiff>0,"positive material pressure and isothermal stiffness",std::min(press,stiff));
    check(far_hi<1e-12,"liquid recovered above the window",far_hi);
    check(latent_lo>0,"solid solution releases latent heat",latent_lo);
    std::printf("rho %.4g Z %.3g X %.2g: T_half %.5g K, min Cv %.4g, |dU-intCv|/dU %.2e, latent (kT per ion at T-edge) %.4f\n",
      s.rho,s.Z,s.X,Tm,cvmin,std::abs(dU-integral)/std::abs(dU),
      latent_lo/(constants::R_gas*Tl*(c.X[0]+c.X[1]/3+c.X[2]/4+c.Z()*0.0551)));
    // Width sensitivity (c/2, 2c): the half-weight temperature is width independent (D=0 there).
    for(double w:{.0025,.01}){const double Tw=onset(s.rho,c,w);
      check(std::abs(Tw/Tm-1)<1e-9,"half-weight temperature independent of width",Tw/Tm-1);}
  }
  // 3. Composition channels against central differences inside the window (solid weight near one half).
  {const double rho=19676.3;const auto base=comp(4.2e-6,1e-4,.0012127);  /* He3 raised so central differences stay non-negative */const double T=onset(rho,base)*1.002;
   // The phase term itself (the unchanged liquid base is roundoff limited near 1e-5 at these steps).
   const auto j=ion_mixture_phase_difference_jets(T,rho,base,10);
   double worst1=0,worst2=0;const double h=1e-6;
   for(int a=0;a<3;++a){
     auto shift=[&](int k,double d){auto c=base;c.X[k==0?0:1]+=d;c.X[2]-=d;return c;};
     // Metals: move the whole GS98 pattern (solar_scaled keeps fractions in X[3..]); rebuild by Z instead.
     auto at=[&](double d){if(a<2)return shift(a,d);return comp(4.2e-6,1e-4,.0012127+d);};
     const auto p=ion_mixture_phase_difference_jets(T,rho,at(h),4),m=ion_mixture_phase_difference_jets(T,rho,at(-h),4);
     for(int i=0;i<3;++i)for(int k=0;i+k<=2;++k){
       const double fd=(p[0][i][k]-m[0][i][k])/(2*h),an=j[1+a][i][k];
       worst1=std::max(worst1,std::abs(fd-an)/std::max(1e-3*std::abs(j[1+a][0][0]),std::abs(an)));}
     for(int b=0;b<3;++b){const int lo=std::min(a,b),hi=std::max(a,b);const int ch=4+(lo==0?hi:(lo==1?2+hi:5));
       for(int i=0;i<2;++i)for(int k=0;i+k<=1;++k){
         const double fd=(p[1+b][i][k]-m[1+b][i][k])/(2*h),an=j[ch][i][k];
         worst2=std::max(worst2,std::abs(fd-an)/std::max(1e-3*std::abs(j[ch][0][0]),std::abs(an)));}}
   }
   check(worst1<1e-6,"composition gradients vs central differences",worst1);
   check(worst2<1e-5,"composition Hessian vs central differences",worst2);
   std::printf("composition channels at weight %.3f: gradient %.2e, Hessian %.2e\n",weight(T,rho,base),worst1,worst2);}
  // 4. Guards.
  {auto throws=[](auto f){try{f();return false;}catch(const std::domain_error&){return true;}};
   {const auto h=comp(.003,0,0);const double rho=2e4,T=1.11e5;   // Gamma_He ~ 140
    bool refused=false;
    try{(void)ion_mixture_phase_difference_jets(T,rho,h,1);}catch(const std::domain_error& e){
      refused=std::string(e.what()).find("hydrogen")!=std::string::npos;}
    check(refused,"hydrogen beyond trace with solid weight refused",T);}
   // Trace hydrogen (X=1e-3) shifts the half-weight temperature only slightly (its lattice branch is unassessed).
   {const double t0=onset(2e4,comp(0,0,0)),t1=onset(2e4,comp(1e-3,0,0));
    check(std::abs(t1/t0-1)<.02,"trace hydrogen onset shift below 2%",t1/t0-1);
    std::printf("trace hydrogen X=1e-3: half-weight temperature shift %.4f\n",t1/t0-1);}
   check(throws([&]{(void)cold_helium_material_jets(2e4,5e4,comp(0,0,.13),1,phase());}),"supercooled helium limit",0);
   // The Gamma_He cut at 90 is invisible over the core composition range and density.
   int cut=0;
   // Metal-rich material occurs only at rho >= 2e4 in the model; lower densities are near-pure helium.
   for(double rho=3e3;rho<=1.2e5;rho*=1.3)for(double Z:{0.,.001,.005,.02,.08,.16}){if(rho<2e4&&Z>.005)continue;const auto c=comp(0,0,Z);
     const double ne=rho*(.5*(1-Z)+gs98_ion_moment(1)*Z)/constants::amu;
     const double g1=std::pow(2.,5./3)*315775.02480407/(std::cbrt(3/(4*M_PI*ne))/5.29177210903e-9);
     for(double g=90;g<=100;g+=.5)try{(void)ion_mixture_phase_difference_jets(g1/g,rho,c,1);}catch(const std::domain_error& e){
       if(std::string(e.what()).find("Gamma cut")!=std::string::npos)++cut;}}
   check(cut==0,"solid weight negligible at the helium Gamma cut",cut);
   // Default path untouched: no phase term without the option.
   ColdHeliumOptions off;const auto a=cold_helium_material_jets(2.2e5,5e4,comp(0,0,.1),10,off);
   ColdHeliumOptions on=phase();on.minimum_solid_gamma=1e9;const auto b=cold_helium_material_jets(2.2e5,5e4,comp(0,0,.1),10,on);
   double diff=0;for(int k=0;k<10;++k)for(int i=0;i<4;++i)for(int j=0;i+j<=3;++j)diff=std::max(diff,std::abs(a[k][i][j]-b[k][i][j]));
   check(diff==0,"phase option adds nothing below the solid-coupling cut",diff);}
  // Continued-liquid scenarios retain the value and first two derivatives at
  // their join, remain above the bcc ground energy at large coupling, and
  // conserve latent/thermal energy with the actual trace-H core composition.
  for(double gstar:{175.,200.,300.}) {
    using J=detail::Taylor3<1>;
    const auto at=[&](double g){return detail::ioffe::fition9_continued(J::variable(g,0),gstar);};
    const auto a=at(gstar*(1-1e-7)),b=at(gstar*(1+1e-7));
    double join=0;
    for(unsigned d=0;d<=2;++d)join=std::max(join,std::abs(a.derivative({d})-b.derivative({d}))/std::max(1.,std::abs(a.derivative({d}))));
    check(join<1e-6,"continued liquid C2 at its join",join);
    const auto high=at(1e4);check(1e4*high.derivative({1})>-.895929256*1e4,
      "continued liquid energy above bcc ground energy",1e4*high.derivative({1}));
  }
  {
    auto o=phase();o.liquid_continuation_gamma=200.;
    const auto c=comp(2.06e-8,9.45e-15,.1364);const double rho=50633.,lo=1.65e5,hi=2.1e5;
    auto state=[&](double T){return cold_helium_material_jets(T,rho,c,1,o)[0];};
    double integral=0,previous=heat(state(lo)),worst=1e300;const unsigned n=1600;
    for(unsigned i=1;i<=n;++i){const double t=lo+(hi-lo)*i/n;const double cv=heat(state(t));
      integral+=(previous+cv)*.5*(hi-lo)/n;previous=cv;worst=std::min(worst,cv);}
    const double delta=energy(state(hi),hi)-energy(state(lo),lo);
    check(worst>0,"continued mixture positive heat capacity",worst);
    check(std::abs(integral/delta-1)<1e-5,"continued mixture heat integral",integral/delta-1);
    o.liquid_continuation_gamma=100.;bool invalid=false;
    try{(void)cold_helium_material_jets(2e5,rho,c,1,o);}catch(const std::domain_error&){invalid=true;}
    check(invalid,"unsupported continuation choice refused",0.);
  }
  std::cout<<"mixture phase: failures "<<failures<<'\n';
  return failures?1:0;
}
