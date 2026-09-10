#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include "differential.hpp"
#include "fermi.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <unordered_map>

namespace ember {
using namespace constants;
namespace {
constexpr double e2 = 4.803204673e-10 * 4.803204673e-10;
constexpr double mev = 1.602176634e-6;
struct Rate { double v, dlnv_dlnT; };
Rate pp(double T9) {                       // p(p,e+ nu)d  - the bottleneck
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 3.381 / t913;
  // Historical fit; not a direct implementation of Solar Fusion II.
  const double f = 4.01e-15 / t923 * std::exp(-tau)
                 * (1.0 + 0.123 * t913 + 1.09 * t923 + 0.938 * T9);
  const double poly = 1.0 + 0.123 * t913 + 1.09 * t923 + 0.938 * T9;
  const double dpoly = (0.123 * t913 / 3.0 + 1.09 * t923 * 2.0 / 3.0 + 0.938 * T9) / poly;
  return {f, -2.0 / 3.0 + tau / 3.0 + dpoly};
}

Rate he3he3(double T9) {                   // He3(He3,2p)He4  - the ppI branch
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 12.276 / t913;
  const double poly = 1.0 - 0.034 * t913 + 0.213 * t923 - 0.026 * T9;
  const double f = 6.04e10 / t923 * std::exp(-tau) * poly;
  const double dpoly = (-0.034 * t913 / 3.0 + 0.213 * t923 * 2.0 / 3.0 - 0.026 * T9) / poly;
  return {f, -2.0 / 3.0 + tau / 3.0 + dpoly};
}

Rate he3he4(double T9) {                   // He3(alpha,gamma)Be7 - ppII/ppIII
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 12.826 / t913;
  const double f = 5.46e6 / t923 * std::exp(-tau);
  return {f, -2.0 / 3.0 + tau / 3.0};
}


struct Reaction {
  double z1,z2,m1,m2,s0,s1,s2; // S derivatives in MeV barn, barn, barn/MeV
};
Reaction reaction(PPReaction r) {
  // Kinematic nuclear masses exclude electrons (electronic binding neglected).
  const double m1=nuclides[0].A*amu-me, m3=nuclides[1].A*amu-2*me;
  const double m4=nuclides[2].A*amu-2*me;
  switch(r) {
  case PPReaction::pp: return {1,1,m1,m1,4.01e-25,4.01e-25*11.2,0};
  case PPReaction::he3_he3: return {2,2,m3,m3,5.21,-4.9,22.};
  case PPReaction::he3_he4: return {2,2,m3,m4,.00056,-.00036,.000151};
  }
  throw std::invalid_argument("PPReaction: unknown reaction");
}
double gamow_energy(const Reaction& r) {
  const double hbar=h/(2*M_PI), mu=r.m1*r.m2/(r.m1+r.m2);
  return 2*mu*std::pow(M_PI*r.z1*r.z2*e2/hbar,2);
}
// Gauss--Legendre integration in ln(E/E0), split at the Gamow peak.
struct Quadrature {
  static constexpr int n=48;
  std::array<double,n> x{},w{};
  Quadrature() {
    for(int i=0;i<n;++i) {
      double z=std::cos(M_PI*(i+.75)/(n+.5));
      for(int it=0;it<30;++it) {
        double p=1,previous=0;
        for(int j=1;j<=n;++j) {double old=previous;previous=p;p=((2*j-1)*z*previous-(j-1)*old)/j;}
        const double dp=n*(z*p-previous)/(z*z-1),dz=p/dp;
        z-=dz;
        if(std::abs(dz)<1e-15) break;
      }
      double p=1,previous=0;
      for(int j=1;j<=n;++j) {double old=previous;previous=p;p=((2*j-1)*z*previous-(j-1)*old)/j;}
      const double dp=n*(z*p-previous)/(z*z-1);
      x[i]=z;w[i]=2/((1-z*z)*dp*dp);
    }
  }
};
void validate(double T,double rho,const Composition& comp) {
  if(!std::isfinite(T) || !std::isfinite(rho) || T<=0 || rho<=0)
    throw std::domain_error("PPChains: positive finite temperature and density required");
  for(double x:comp.X) if(!std::isfinite(x) || x<0)
    throw std::domain_error("PPChains: finite nonnegative abundances required");
  // No normalization here: composition_response exposes unconstrained partials.
  if(!(comp.mu_elec_inv()>0)) throw std::domain_error("PPChains: empty charged mixture");
}
struct Susceptibility { double eta,theta,dtheta_dlnT,dtheta_dlnne; };
Susceptibility electrons(double T,double ne) {
  const double beta=kB*T/(me*c_light*c_light);
  const double norm=8*M_PI*std::pow(me*c_light/h,3);
  const double lambda=h/std::sqrt(2*M_PI*me*kB*T);
  const double nd=std::log(ne*lambda*lambda*lambda/2);
  const double xf=std::cbrt(3*ne/norm);
  double eta=nd<0?nd:xf*xf/((std::sqrt(1+xf*xf)+1)*beta);
  for(int it=0;it<100;++it) {
    const auto f=fermi::evaluate(eta,beta);
    if(!(f.In>0 && f.dIn_deta>0)) break;
    const double residual=std::log(norm*f.In/ne);
    const double step=residual*f.In/f.dIn_deta;
    if(std::abs(residual)<2e-13 && std::abs(step)<2e-12*(1+std::abs(eta))) {
      const double theta=f.dIn_deta/f.In;
      return {eta,theta,(f.d2In_detadlnb-f.d2In_deta2*f.dIn_dlnb/f.dIn_deta)/f.In,
        f.d2In_deta2/f.dIn_deta-theta};
    }
    eta-=std::clamp(step,-4*(1+std::abs(eta)),4*(1+std::abs(eta)));
  }
  throw std::runtime_error("PPChains: electron susceptibility inversion failed");
}
} // namespace

ThermonuclearRate pp_bare_rate(double T,PPReaction which,PPRates prescription) {
  if(!std::isfinite(T) || T<=0) throw std::domain_error("pp_bare_rate: invalid temperature");
  const auto r=reaction(which);
  if(prescription==PPRates::solar_fusion_ii && T>2e7)
    throw std::domain_error("PPChains: Solar Fusion II low-energy expansion limited to T<=2e7 K");
  if(T<1e5) return {};
  if(prescription==PPRates::legacy) {
    const auto v=which==PPReaction::pp?pp(T*1e-9):which==PPReaction::he3_he3?he3he3(T*1e-9):he3he4(T*1e-9);
    return {v.v,v.dlnv_dlnT};
  }
  // Bare SFII quadrature depends only on temperature and reaction, whereas
  // burning Newton iterations change abundances at fixed thermal states.
  // Exact, bounded, thread-local memoization leaves screening and every
  // composition response fully state-dependent and preserves result bits.
  using CachedRates=std::array<std::optional<ThermonuclearRate>,3>;
  thread_local std::unordered_map<double,CachedRates> cache;
  const auto index=static_cast<std::size_t>(which); // reaction() validated it above
  const auto found=cache.find(T);
  if(found!=cache.end() && found->second[index])return *found->second[index];
  static const Quadrature q;
  const double kt=kB*T, eg=gamow_energy(r), peak=std::cbrt(eg*kt*kt/4);
  const double mu=r.m1*r.m2/(r.m1+r.m2);
  double integral=0,moment=0;
  for(double mid:{-2.,2.}) for(int i=0;i<q.n;++i) {
    const double energy=peak*std::exp(mid+2*q.x[i]), E=energy/mev;
    const double S=(r.s0+E*(r.s1+.5*E*r.s2))*mev*1e-24;
    const double term=2*q.w[i]*energy*S*std::exp(-energy/kt-std::sqrt(eg/energy));
    integral+=term;moment+=term*energy/kt;
  }
  const double rate=NA*std::sqrt(8/(M_PI*mu))/std::pow(kt,1.5)*integral;
  const ThermonuclearRate result{rate,integral>0?-1.5+moment/integral:0};
  if(found==cache.end() && cache.size()>=8192)cache.clear();
  cache[T][index]=result;
  return result;
}

ScreeningState pp_screening(double T,double rho,const Composition& comp,PPReaction which,PPScreening model) {
  validate(T,rho,comp);
  const auto r=reaction(which);
  using D=detail::Differential<NSPEC+2>;
  using detail::exp;using detail::log;
  const auto temp=exp(D::variable(std::log(T),0)),density=exp(D::variable(std::log(rho),1));
  D ye,ions;
  for(std::size_t j=0;j<NSPEC;++j) {
    if(j>=3 && comp.metal_inventory==MetalInventory::gs98) {
      const D metal=D::variable(comp.X[j],j+2);
      ye=ye+metal*comp.metal_ion_moment(1);ions=ions+metal*comp.metal_ion_moment(2);
      continue;
    }
    const D y=D::variable(comp.X[j],j+2)/comp.abundance_weight(j);
    ye=ye+y*nuclides[j].Z;ions=ions+y*nuclides[j].Z*nuclides[j].Z;
  }
  ScreeningState out;
  D theta(1);
  if(model!=PPScreening::legacy_weak) {
    const auto f=electrons(T,rho*NA*ye.value);
    out.electron_eta=f.eta;theta.value=f.theta;
    theta.d[0]=f.dtheta_dlnT;theta.d[1]=f.dtheta_dlnne;
    for(std::size_t j=0;j<NSPEC;++j) theta.d[j+2]=f.dtheta_dlnne*ye.d[j+2]/ye.value;
  } else out.electron_eta=std::numeric_limits<double>::quiet_NaN();
  out.electron_susceptibility=theta.value;
  const auto ge=e2/(kB*temp)*exp(log(4*M_PI*NA*density*ye/3)/3);
  out.gamma_e=ge.value;
  const double g12=ge.value*2*r.z1*r.z2/(std::cbrt(r.z1)+std::cbrt(r.z2));
  out.zeta=g12/std::cbrt(gamow_energy(r)/(4*kB*T));
  // Explicit classical-ion thermonuclear domain, not a cap or an extrapolation.
  if(model!=PPScreening::legacy_weak && out.zeta>.2)
    throw std::domain_error("PPChains: classical-ion screening requires zeta<=0.2; quantum burning unavailable");
  const auto weak=r.z1*r.z2*e2/(kB*temp)*exp(.5*log(4*M_PI*e2*NA*density*(ions+theta*ye)/(kB*temp)));
  D exponent=weak;
  if(model==PPScreening::legacy_weak && weak.value>=2) exponent=D(2);
  if(model==PPScreening::salpeter_van_horn) {
    const auto strong=.9*ge*(std::pow(r.z1+r.z2,5./3)-std::pow(r.z1,5./3)-std::pow(r.z2,5./3));
    exponent=weak*strong/exp(.5*log(weak*weak+strong*strong));
  }
  out.log_factor=exponent.value;out.dlog_dlnT=exponent.d[0];out.dlog_dlnRho=exponent.d[1];
  for(std::size_t j=0;j<NSPEC;++j) out.dlog_dX[j]=exponent.d[j+2];
  return out;
}

NuclearResponse PPChains::composition_response(double T,double rho,const Composition& comp) const {
  validate(T,rho,comp);
  NuclearResponse result;
  if(T<1e5) return result; // retained explicit negligible-burning cutoff
  const std::array rates{pp_bare_rate(T,PPReaction::pp,rates_),pp_bare_rate(T,PPReaction::he3_he3,rates_),
    pp_bare_rate(T,PPReaction::he3_he4,rates_)};
  // Classical screening depends on charge, so both helium reactions share it.
  const auto f1=pp_screening(T,rho,comp,PPReaction::pp,screening_);
  const auto f2=pp_screening(T,rho,comp,PPReaction::he3_he3,screening_);
  const std::array screens{f1,f2,f2};
  const std::array w{comp.abundance_weight(0),comp.abundance_weight(1),comp.abundance_weight(2)};
  const std::array y{comp.X[0]/w[0],comp.X[1]/w[1],comp.X[2]/w[2]};
  const double m1=nuclides[0].A,m3=nuclides[1].A,m4=nuclides[2].A;
  const std::array total_q{(3*m1-m3)*c_light*c_light,(2*m3-2*m1-m4)*c_light*c_light,
    (m3+m1-m4)*c_light*c_light};
  // Retained reduced ppII neutrino approximation, qualified in docs/NUCLEAR.md.
  const std::array nu_q{.265*mev*NA,0.,.861*mev*NA};
  constexpr std::array<std::array<double,3>,3> stoich{{{-3,1,0},{2,-2,1},{-1,-1,1}}};
  constexpr std::array a{0,1,1},b{0,1,2};
  auto& s=result.state;
  double epsT=0,epsR=0;
  for(std::size_t k=0;k<3;++k) {
    const double coefficient=(k<2?.5:1.)*rho*rates[k].molar_rate*std::exp(screens[k].log_factor);
    // Moles of reactions / g / s, not individual reactions / g / s.
    const double rate=coefficient*y[a[k]]*y[b[k]], heat=total_q[k]-nu_q[k];
    s.eps+=rate*heat;s.eps_neutrino+=rate*nu_q[k];
    epsT+=rate*heat*(rates[k].dlnrate_dlnT+screens[k].dlog_dlnT);
    epsR+=rate*heat*(1+screens[k].dlog_dlnRho);
    for(std::size_t i=0;i<3;++i) s.dXdt[i]+=stoich[k][i]*rate*w[i];
    for(std::size_t j=0;j<NSPEC;++j) {
      const double derivative=rate*screens[k].dlog_dX[j]
        +(j==static_cast<std::size_t>(a[k])?coefficient*y[b[k]]/w[a[k]]:0)
        +(j==static_cast<std::size_t>(b[k])?coefficient*y[a[k]]/w[b[k]]:0);
      result.deps_dX[j]+=heat*derivative;
      for(std::size_t i=0;i<3;++i) result.d_dXdt_dX[i][j]+=stoich[k][i]*w[i]*derivative;
    }
  }
  if(s.eps>0) {s.dlneps_dlnT=epsT/s.eps;s.dlneps_dlnRho=epsR/s.eps;}
  return result;
}
NuclearState PPChains::eval(double T,double rho,const Composition& comp) const {
  return composition_response(T,rho,comp).state;
}
const char* PPChains::name() const {
  if(rates_==PPRates::legacy) {
    if(screening_==PPScreening::legacy_weak) return "legacy pp fits; capped classical weak screening";
    if(screening_==PPScreening::debye_fermi) return "legacy pp fits; finite-degeneracy Debye screening";
    return "legacy pp fits; finite-degeneracy Salpeter--Van Horn screening";
  }
  if(screening_==PPScreening::legacy_weak) return "Solar Fusion II quadrature; capped classical weak screening";
  if(screening_==PPScreening::debye_fermi) return "Solar Fusion II quadrature; finite-degeneracy Debye screening";
  return "Solar Fusion II quadrature; finite-degeneracy Salpeter--Van Horn screening";
}
} // namespace ember
