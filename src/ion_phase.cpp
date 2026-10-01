#include "ember/ion_phase.hpp"
#include "ember/detail/ion_mixture.hpp"
#include "ember/detail/taylor3.hpp"
#include "ember/detail/ion_ocp_components.hpp"
#include "ember/constants.hpp"
#include "ember/gs98_mixture.hpp"
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace ember {
namespace {
using J=detail::Taylor3<2>;   // variables: 0 = ln T, 1 = ln n_e

template<class T> T pow(const T& x,double p){return detail::exp(p*detail::log(x));}
// eos22.f writes many screening constants as default-real (single precision) literals; reproduce them exactly.
constexpr double sp(float x){return x;}

template<class J> using Plasma=detail::ioffe::Plasma<J>;
using detail::ioffe::plasma;
Plasma<J> plasma(double logT,double logne,double A,double Z) {
  return detail::ioffe::plasma(J::variable(logT,0),J::variable(logne,1),A,Z);
}

// LIQUBC quantum-liquid free energy per ion (kT), Baiko & Chugunov (2022) eqs. 33-34.
template<class J> J liqubc(const J& rsi,const J& tpt) {
  constexpr double P1=.351,P2=.294,P3=90.,EPS3=2e-3,BIG=33.;
  const auto ci1=P1*rsi/(P3+rsi);
  const std::array<J,3> ci{ci1,J(P2),detail::sqrt(1-ci1*ci1-P2*P2)};
  J f(0);
  for(const auto& c:ci) {
    const auto x=c*tpt;
    if(x.value()<EPS3)f=f+x*x/24.;
    else if(x.value()>BIG)f=f+x*.5-detail::log(x);
    else f=f+detail::log(1-detail::exp(-x))+x*.5-detail::log(x);
  }
  return f;
}
// HLfit12 (bcc) thermal free energy of the harmonic lattice per ion (kT), zero point excluded.
template<class J> J hlfit12(const J& eta) {
  if(!(eta.value()>1e-5&&eta.value()<1e5))throw std::domain_error("ion phase: HLfit12 asymptotic regime not ported");
  constexpr double ALPHA=.265764,BETA=.334547,GAMMA=.932446,A1=.1839,A2=.593586,A3=.0054814,A4=5.01813e-4,
      A6=3.9247e-7,A8=5.8356e-11,B0=261.66,B2=7.07997,B4=.0409484,B5=.000397355,B6=5.11148e-5,B7=2.19749e-6,
      C9=.004757014,C11=.0047770935;
  const double B9=A6*C9,B11=A8*C11;
  const auto e2=eta*eta,e3=e2*eta,e4=e3*eta,e5=e4*eta,e6=e5*eta,e7=e6*eta,e8=e7*eta;
  const auto up=1+A1*eta+A2*e2+A3*e3+A4*e4+A6*e6+A8*e8;
  const auto dn=B0+B2*e2+B4*e4+B5*e5+B6*e6+B7*e7+e8*(B9*eta+B11*e3);
  return detail::log(1-detail::exp(-ALPHA*eta))+detail::log(1-detail::exp(-BETA*eta))
      +detail::log(1-detail::exp(-GAMMA*eta))-up/dn;
}
// ANHBC quantum anharmonic crystal free energy per ion (kT), Baiko & Chugunov (2022) eq. 15.
template<class J> J anhbc(const J& g,const J& tpt) {
  constexpr double A1CL=10.2,A2CL=248.,A3CL=2.03e5,A1Q=-.62/6.,A2Q=-.56,A3Q=2.35,A11=-10.,A12=6e-3;
  const double A13=-A1CL-A11,A14=(A1Q-A11*A12)/A13,A21=std::pow(-2*A2Q/A2CL,4./3);
  const auto q=g/tpt,t2=tpt*tpt,t4=t2*t2;
  const auto a1=A11/t2/(1+A12*t2)+A13/t2/(1+A14*t2)+A1Q;
  const auto root=detail::sqrt(detail::sqrt(1+A21*t4));
  const auto a2=-A2CL*.5*(root/tpt)*(root/tpt)*(root/tpt);
  const auto a3=-A3CL/3./t4+A3Q;
  const auto q2=q*q;
  return g*(a1/q2+a2/(q2*q)+a3/(q2*q2));
}
// FSCRsol8 electron-ion screening in the bcc crystal, per ion (kT).
template<class J> J fscr_solid(const Plasma<J>& p) {
  const double Z=std::max(p.Z,1.),ZLN=std::log(Z),Z13=std::cbrt(Z);
  constexpr double AP1=sp(1.1866f),AP2=sp(.684f),AP3=sp(17.9f),AP4=sp(41.5f),PX=sp(.205f),ENAT=sp(2.7182818285f);
  const auto xsr=p.x,x2=xsr*xsr;
  const double P1=sp(.00352f)*(1-AP1/std::pow(Z,sp(.267f))+sp(.27f)/Z);
  const double P2=1+sp(2.25f)/Z13*(1+AP2*std::pow(Z,5)+sp(.222f)*std::pow(Z,6))/(1+sp(.222f)*std::pow(Z,6));
  const auto finf=detail::sqrt(P2/x2+1)*(Z13*Z13*P1);
  const double R1=AP4/(1+ZLN),R2=sp(.395f)*ZLN+sp(.347f)/Z/std::sqrt(Z),R3=1/(1+ZLN*std::sqrt(ZLN)*sp(.01f)+sp(.097f)/(Z*Z));
  const auto q1=(R1+AP3*x2)/(1+R2*x2);
  J sup;
  if(p.tpt.value()<6/PX) {
    const auto y1=detail::exp(PX*p.tpt*PX*p.tpt);
    sup=detail::sqrt(detail::log(1+y1)/detail::log(ENAT-(ENAT-2)/y1));
  } else sup=PX*p.tpt;
  const auto gr3=pow(p.gami/sup,R3);
  return -p.gami*finf*(1+q1/gr3);
}
template<class J> J liquid_f(const Plasma<J>& p,IonScreening s,double gstar=0) {
  const auto fid=detail::ioffe::ideal_ion(p);
  auto f=detail::ioffe::fition9_continued(p.gami,gstar)+fid+liqubc(p.rsi,p.tpt);
  if(s==IonScreening::fitted_screening)f=f+detail::ioffe::fscr_liquid(p);
  return f;
}
template<class J> J solid_f(const Plasma<J>& p,IonScreening s) {
  constexpr double CM=.895929256,U1=.5113875;
  auto f=-CM*p.gami+hlfit12(p.tpt)+1.5*U1*p.tpt+anhbc(p.gami,p.tpt);
  if(s==IonScreening::fitted_screening)f=f+fscr_solid(p);
  return f;
}
HelmholtzJet to_jet(const J& f,double A) {
  HelmholtzJet out{};
  for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)out[i][j]=constants::R_gas/A*f.derivative({i,j});
  return out;
}
} // namespace

HelmholtzJet ion_ocp_liquid_jet(double logT,double logne,double A,double Z,IonScreening s) {
  return to_jet(liquid_f(plasma(logT,logne,A,Z),s),A);
}
HelmholtzJet ion_ocp_solid_jet(double logT,double logne,double A,double Z,IonScreening s) {
  return to_jet(solid_f(plasma(logT,logne,A,Z),s),A);
}

std::array<HelmholtzJet,10> ion_phase_difference_jets(
    double T,double rho,const Composition& c,std::size_t channels,const IonPhaseOptions& o) {
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0)||(channels!=1&&channels!=4&&channels!=10))
    throw std::domain_error("ion phase: invalid state or derivative request");
  if(!(o.width>0&&o.width<=.05)||!(o.max_non_helium>=0&&o.max_non_helium<.01))
    throw std::domain_error("ion phase: invalid options");
  if(c.basis!=AbundanceBasis::baryon_mass||c.metal_inventory!=MetalInventory::gs98||c[Species::H2]!=0)
    throw std::domain_error("ion phase: requires baryonic GS98 material after D mapping");
  const double non_helium=c.X[0]+c.X[1]+c.Z();
  if(!(non_helium<=o.max_non_helium)) {
    std::ostringstream m;m<<std::setprecision(4)<<"ion phase: pure-He reference outside composition domain (non-He4 "
      <<non_helium<<" > "<<o.max_non_helium<<"); a mixture phase treatment is required";
    throw std::domain_error(m.str());
  }
  auto evaluate=[&](double logT,double logne,double A,double Z)->HelmholtzJet {
    if(!(A==4&&Z==2))return HelmholtzJet{};   // no phase term for other species
    const auto p=plasma(logT,logne,A,Z);
    const double rsi=p.rsi.value(),tpt=p.tpt.value(),g=p.gami.value();
    if(rsi<500||rsi>1.2e5||tpt>30) {
      std::ostringstream m;m<<std::setprecision(4)<<"ion phase: outside Baiko-Chugunov fitted domain (R_S="<<rsi
        <<", Tp/T="<<tpt<<')';
      throw std::domain_error(m.str());
    }
    if(g<o.minimum_solid_gamma)return HelmholtzJet{};  // liquid only; spurious solid root excluded
    const auto fl=liquid_f(p,o.screening),fs=solid_f(p,o.screening);
    const auto d=(fs-fl)/o.width;   // soft minimum written relative to the liquid
    J delta;
    if(d.value()>=0)delta=-o.width*detail::log1p(detail::exp(-d));
    else delta=(fs-fl)-o.width*detail::log1p(detail::exp(d));
    // The cut at minimum_solid_gamma must be invisible: require negligible solid weight near it.
    if(g<o.minimum_solid_gamma+10&&d.value()<20)
      throw std::domain_error("ion phase: solid branch not negligible near the Gamma cut");
    return to_jet(delta,A);
  };
  return detail::common_density_ion_jets(T,rho,c,channels,evaluate);
}
namespace {
using J5=detail::Taylor3<5>;   // variables: ln T, ln rho, X_H, X_He3, Z (He4 is the remainder)
struct MixtureTerms {J5 D,W;double gamma_he;};
void check_mixture(double T,double rho,const Composition& c,std::size_t channels,const MixturePhaseOptions& o) {
  if(!(o.liquid_continuation_gamma==0. || (std::isfinite(o.liquid_continuation_gamma)
      && o.liquid_continuation_gamma>=175. && o.liquid_continuation_gamma<=300.)))
    throw std::domain_error("mixture phase: invalid metal-liquid continuation");
  for(double x:c.X)if(!std::isfinite(x)||x<0)
    throw std::domain_error("mixture phase: invalid abundance");
  if(std::abs(c.sum()-1)>1e-8)
    throw std::domain_error("mixture phase: abundances do not sum to one");
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0)||(channels!=1&&channels!=4&&channels!=10))
    throw std::domain_error("mixture phase: invalid state or derivative request");
  if(!(o.width>0&&o.width<=.05)||!(o.minimum_solid_gamma>=60)||!(o.trace_hydrogen>=0&&o.trace_hydrogen<=.01))
    throw std::domain_error("mixture phase: invalid options");
  if(c.basis!=AbundanceBasis::baryon_mass||c.metal_inventory!=MetalInventory::gs98||c[Species::H2]!=0)
    throw std::domain_error("mixture phase: requires baryonic GS98 material after D mapping");
}
MixtureTerms mixture_terms(double T,double rho,const Composition& c,const MixturePhaseOptions& o) {
  const J5 xh=J5::variable(c.X[0],2),x3=J5::variable(c.X[1],3),xz=J5::variable(c.Z(),4);
  const J5 x4=1-xh-x3-xz;
  const double Ye=c.X[0]+(2./3)*c.X[1]+.5*c.X[2]+gs98_ion_moment(1)*c.Z();
  const double logT=std::log(T),logne=std::log(rho*Ye/constants::amu);
  const double gamma_he=detail::ioffe::plasma(logT,logne,4.,2.).gami;
  J5 D(0);
  if(gamma_he>=o.minimum_solid_gamma) {
    for(const double A:{3.,4.}) {
      const auto p=detail::ioffe::plasma(logT,logne,A,2.);
      if(p.rsi<500||p.rsi>1.2e5||p.tpt>30) {
        std::ostringstream m;m<<std::setprecision(4)<<"mixture phase: He isotope outside Baiko-Chugunov fitted domain (R_S="
          <<p.rsi<<", Tp/T="<<p.tpt<<')';
        throw std::domain_error(m.str());
      }
    }
    // Each ion depends only on T and electron density. Apply the shared
    // mixture chain rule once, instead of differentiating every ion through
    // all five independent variables. The nonlinear minimum below still
    // differentiates the complete dimensional D/W, including its width.
    const auto d=detail::common_density_ion_jets(T,rho,c,10,[&](double lt,double ln,double A,double Z){
      const auto p=plasma(lt,ln,A,Z);
      return to_jet(solid_f(p,o.screening)-liquid_f(p,o.screening,o.liquid_continuation_gamma),A);
    });
    for(std::size_t k=0;k<J5::count;++k) {
      const auto p=J5::layout.powers[k];
      const unsigned order=p[2]+p[3]+p[4];
      if(order>2)continue;  // no requested response uses three composition derivatives
      unsigned channel=0,factor=1;
      for(auto power:p)factor*=power==3?6:(power==2?2:1);
      if(order==1)for(unsigned a=0;a<3;++a)if(p[2+a])channel=1+a;
      if(order==2) {
        unsigned ch=4;
        for(unsigned a=0;a<3;++a)for(unsigned b=a;b<3;++b,++ch) {
          J5::Powers q{};++q[2+a];++q[2+b];
          if(q[2]==p[2]&&q[3]==p[3]&&q[4]==p[4])channel=ch;
        }
      }
      D.c[k]=d[channel][p[0]][p[1]]/factor;
    }
  }
  double metal_moment=0;
  for(const auto& m:gs98_metals)metal_moment+=m.fraction/m.mass_number;
  const J5 W=(o.width*constants::R_gas)*(xh+x3/3.+x4/4.+metal_moment*xz);
  return {D,W,gamma_he};
}
}  // namespace

std::array<double,2> ion_mixture_phase_weight(double T,double rho,const Composition& c,const MixturePhaseOptions& o) {
  check_mixture(T,rho,c,1,o);
  const auto m=mixture_terms(T,rho,c,o);
  if(m.gamma_he<o.minimum_solid_gamma)return {0.,INFINITY};
  const double y=m.D.value()/m.W.value();
  return {y>=0?std::exp(-y)/(1+std::exp(-y)):1/(1+std::exp(y)),y};
}

std::array<double,6> ion_mixture_phase_weight_response(double T,double rho,const Composition& c,
    const MixturePhaseOptions& o) {
  check_mixture(T,rho,c,10,o);
  const auto m=mixture_terms(T,rho,c,o);
  std::array<double,6> out{};
  if(m.gamma_he<o.minimum_solid_gamma)return out;
  const J5 y=m.D/m.W;
  const double e=std::exp(-std::abs(y.value()));
  out[0]=y.value()>=0?e/(1+e):1/(1+e);
  const double slope=-e/((1+e)*(1+e));
  for(unsigned k=0;k<5;++k){J5::Powers p{};p[k]=1;out[1+k]=slope*y.derivative(p);}
  return out;
}

std::array<HelmholtzJet,10> ion_mixture_phase_difference_jets(
    double T,double rho,const Composition& c,std::size_t channels,const MixturePhaseOptions& o) {
  check_mixture(T,rho,c,channels,o);
  const auto m=mixture_terms(T,rho,c,o);
  std::array<HelmholtzJet,10> out{};
  if(m.gamma_he<o.minimum_solid_gamma)return out;   // liquid; the forced-solid formula re-crosses spuriously below
  const J5 y=m.D/m.W;
  // The coupling cut and the hydrogen limit must be invisible: negligible solid weight near/beyond them.
  if(m.gamma_he<o.minimum_solid_gamma+10&&y.value()<20)
    throw std::domain_error("mixture phase: solid weight not negligible near the helium Gamma cut");
  if(c.X[0]>o.trace_hydrogen&&y.value()<20) {
    std::ostringstream s;s<<std::setprecision(4)<<"mixture phase: hydrogen "<<c.X[0]
      <<" beyond the trace limit where the solid solution is not negligible (D/W="<<y.value()<<')';
    throw std::domain_error(s.str());
  }
  const J5 g=y.value()>=0?-detail::log1p(detail::exp(-y)):y-detail::log1p(detail::exp(y));
  const J5 f=m.W*g;
  using P=J5::Powers;
  for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)out[0][i][j]=f.derivative(P{i,j,0,0,0});
  if(channels>1)for(unsigned a=0;a<3;++a)for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j<=2;++j) {
    P p{i,j,0,0,0};p[2+a]=1;out[1+a][i][j]=f.derivative(p);
  }
  if(channels>4){unsigned ch=4;
    for(unsigned a=0;a<3;++a)for(unsigned b=a;b<3;++b,++ch)for(unsigned i=0;i<2;++i)for(unsigned j=0;i+j<=1;++j) {
      P p{i,j,0,0,0};++p[2+a];++p[2+b];out[ch][i][j]=f.derivative(p);
    }}
  return out;
}
} // namespace ember
