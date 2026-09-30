#include "ember/detail/ion_free_energy.hpp"
#include "ember/nuclear.hpp"
#include "ember/deuterium.hpp"
#include "ember/constants.hpp"
#include "ember/detail/differential.hpp"
#include "fermi.hpp"
#include <algorithm>
#include <atomic>
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
Reaction reaction(PPReaction r,PPRates prescription=PPRates::solar_fusion_ii) {
  // Kinematic nuclear masses exclude electrons (electronic binding neglected).
  const double m1=nuclides[0].A*amu-me, m3=nuclides[1].A*amu-2*me;
  const double m4=nuclides[2].A*amu-2*me;
  if(r==PPReaction::deuterium_p)
    return {1,1,m1,deuterium_atomic_mass*amu-me,0,0,0};
  if(prescription==PPRates::solar_fusion_iii) {
    // Acharya et al. (2025), equations 8--9 and section V.C.
    // The 34 reaction uses the full equation 14 below, not a Taylor fit.
    switch(r) {
    case PPReaction::pp: return {1,1,m1,m1,4.09e-25,4.09e-25*11.0,4.09e-25*242.};
    case PPReaction::he3_he3: return {2,2,m3,m3,5.21,-4.9,22.42};
    case PPReaction::he3_he4: return {2,2,m3,m4,.0005610,0,0};
    case PPReaction::deuterium_p: break; // handled above
    }
  }
  switch(r) {
  case PPReaction::pp: return {1,1,m1,m1,4.01e-25,4.01e-25*11.2,0};
  case PPReaction::he3_he3: return {2,2,m3,m3,5.21,-4.9,22.};
  case PPReaction::he3_he4: return {2,2,m3,m4,.00056,-.00036,.000151};
  case PPReaction::deuterium_p: break; // handled above
  }
  throw std::invalid_argument("PPReaction: unknown reaction");
}
double gamow_energy(const Reaction& r) {
  const double hbar=h/(2*M_PI), mu=r.m1*r.m2/(r.m1+r.m2);
  return 2*mu*std::pow(M_PI*r.z1*r.z2*e2/hbar,2);
}
Reaction cn_reaction(PPRates rates,CNReaction which=CNReaction::n14_p) {
  if(rates!=PPRates::solar_fusion_ii && rates!=PPRates::solar_fusion_iii)
    throw std::invalid_argument("CNCycle: Solar Fusion II or III required");
  const double mp=nuclides[0].A*amu-me;
  // SFIII Table I; SFII Table XII. s2 is the second derivative, not
  // the coefficient of E^2. N14 retains SFII's derivatives with SFIII S(0).
  switch(which) {
  case CNReaction::c12_p:
    return rates==PPRates::solar_fusion_iii
      ?Reaction{1,6,mp,nuclides[3].A*amu-6*me,1.44e-3,2.71e-3,3.74e-2}
      :Reaction{1,6,mp,nuclides[3].A*amu-6*me,1.34e-3,2.6e-3,8.3e-2};
  case CNReaction::c13_p:
    return rates==PPRates::solar_fusion_iii
      ?Reaction{1,6,mp,nuclides[4].A*amu-6*me,6.1e-3,1.04e-2,9.20e-2}
      :Reaction{1,6,mp,nuclides[4].A*amu-6*me,7.6e-3,-7.83e-3,7.29e-1};
  case CNReaction::n14_p:
    return {1,7,mp,nuclides[5].A*amu-7*me,
      rates==PPRates::solar_fusion_iii?1.68e-3:1.66e-3,-3.3e-3,4.4e-2};
  }
  throw std::invalid_argument("CNReaction: unknown reaction");
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
Susceptibility electrons_exact(double T,double ne) {
  // pp and CN captures at one state need the same electron inversion.
  // Reuse its result only for identical temperature and electron density;
  // reaction charges and composition chain-rule factors remain independent.
  struct Key {
    double temperature,density;
    bool operator==(const Key&) const = default;
  };
  struct Hash {
    std::size_t operator()(const Key& key) const {
      const auto a=std::hash<double>{}(key.temperature),b=std::hash<double>{}(key.density);
      return a^(b+0x9e3779b9+(a<<6)+(a>>2));
    }
  };
  thread_local std::unordered_map<Key,Susceptibility,Hash> cache;
  const Key key{T,ne};
  if(const auto found=cache.find(key);found!=cache.end())return found->second;
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
      const Susceptibility result{eta,theta,(f.d2In_detadlnb-f.d2In_deta2*f.dIn_dlnb/f.dIn_deta)/f.In,
        f.d2In_deta2/f.dIn_deta-theta};
      if(cache.size()>=8192)cache.clear();
      cache.emplace(key,result);return result;
    }
    eta-=std::clamp(step,-4*(1+std::abs(eta)),4*(1+std::abs(eta)));
  }
  throw std::runtime_error("PPChains: electron susceptibility inversion failed");
}
std::atomic<double> screening_reuse_spacing{0};
// Optional first-order reuse: exact at the nearest point of a grid of spacing
// h in (ln T, ln ne), extended with its exact derivatives. The anchor depends
// only on the state, so results do not depend on thread or call order. eta is
// the anchor's (diagnostic only).
Susceptibility electrons(double T,double ne) {
  const double h=screening_reuse_spacing.load(std::memory_order_relaxed);
  if(!(h>0))return electrons_exact(T,ne);
  const double lt=std::log(T),ln=std::log(ne);
  const double at=std::nearbyint(lt/h)*h,an=std::nearbyint(ln/h)*h;
  auto s=electrons_exact(std::exp(at),std::exp(an));
  s.theta+=s.dtheta_dlnT*(lt-at)+s.dtheta_dlnne*(ln-an);
  return s;
}
} // namespace

std::atomic<double> quantum_screening_zeta_max{0};
std::atomic<double> quantum_burning_fuel_limit{0};
void set_quantum_burning_fuel_limit(double max_fraction) {
  if(!std::isfinite(max_fraction) || max_fraction<0 || max_fraction>1e-6)
    throw std::invalid_argument("quantum burning fuel limit must lie in [0,1e-6]");
  quantum_burning_fuel_limit=max_fraction;
}
void set_quantum_screening(double zeta_max) {
  if(!std::isfinite(zeta_max) || zeta_max<0 || zeta_max>1.6)
    throw std::invalid_argument("quantum screening zeta_max must lie in [0,1.6] (mean-field WKB domain)");
  quantum_screening_zeta_max=zeta_max;
}

void set_screening_reuse(double h) {
  if(!std::isfinite(h) || h<0 || h>.01)throw std::invalid_argument("screening reuse spacing must lie in [0,0.01]");
  screening_reuse_spacing=h;
}

ThermonuclearRate pp_bare_rate(double T,PPReaction which,PPRates prescription) {
  if(which==PPReaction::deuterium_p) {
    if(prescription!=PPRates::solar_fusion_iii)
      throw std::invalid_argument("deuterium capture requires Solar Fusion III");
    return deuterium_bare_rate(T);
  }
  if(!std::isfinite(T) || T<=0) throw std::domain_error("pp_bare_rate: invalid temperature");
  if(prescription!=PPRates::legacy && prescription!=PPRates::solar_fusion_ii
      && prescription!=PPRates::solar_fusion_iii)
    throw std::invalid_argument("PPRates: unknown prescription");
  const auto r=reaction(which,prescription);
  if(prescription!=PPRates::legacy && T>2e7)
    throw std::domain_error("PPChains: Solar Fusion low-energy rates limited to T<=2e7 K");
  if(T<1e5) return {};
  if(prescription==PPRates::legacy) {
    const auto v=which==PPReaction::pp?pp(T*1e-9):which==PPReaction::he3_he3?he3he3(T*1e-9):he3he4(T*1e-9);
    return {v.v,v.dlnv_dlnT};
  }
  // Bare quadrature depends on temperature, reaction and rate prescription;
  // burning Newton iterations change abundances at fixed thermal states.
  // Exact, bounded, thread-local memoization leaves screening and every
  // composition response fully state-dependent and preserves result bits.
  using CachedRates=std::array<std::optional<ThermonuclearRate>,6>;
  thread_local std::unordered_map<double,CachedRates> cache;
  const auto index=static_cast<std::size_t>(which)+(prescription==PPRates::solar_fusion_iii?3:0);
  const auto found=cache.find(T);
  if(found!=cache.end() && found->second[index])return *found->second[index];
  static const Quadrature q;
  const double kt=kB*T, eg=gamow_energy(r), peak=std::cbrt(eg*kt*kt/4);
  const double mu=r.m1*r.m2/(r.m1+r.m2);
  double integral=0,moment=0;
  for(double mid:{-2.,2.}) for(int i=0;i<q.n;++i) {
    const double energy=peak*std::exp(mid+2*q.x[i]), E=energy/mev;
    double S;
    if(prescription==PPRates::solar_fusion_iii && which==PPReaction::he3_he4) {
      // Equation 14 is specified only through 1.6 MeV. At T<=20 MK the
      // omitted Maxwell tail starts beyond 928 kT; do not extrapolate it.
      S=E>1.6?0:r.s0*std::exp(-.5374*E)*(1+E*E*(-.4829+E*(.6310-.1527*E)))*mev*1e-24;
    } else S=(r.s0+E*(r.s1+.5*E*r.s2))*mev*1e-24;
    const double term=2*q.w[i]*energy*S*std::exp(-energy/kt-std::sqrt(eg/energy));
    integral+=term;moment+=term*energy/kt;
  }
  const double rate=NA*std::sqrt(8/(M_PI*mu))/std::pow(kt,1.5)*integral;
  const ThermonuclearRate result{rate,integral>0?-1.5+moment/integral:0};
  if(found==cache.end() && cache.size()>=8192)cache.clear();
  cache[T][index]=result;
  return result;
}

static ScreeningState screening_response(double T,double rho,const Composition& comp,const Reaction& r,PPScreening model,
                                         std::optional<Susceptibility>* shared_electrons) {
  validate(T,rho,comp);
  using D=detail::Differential<NSPEC+2>;
  using detail::exp;using detail::log;
  const auto temp=exp(D::variable(std::log(T),0)),density=exp(D::variable(std::log(rho),1));
  D ye,ions;
  for(std::size_t j=0;j<NSPEC;++j) {
    if(is_metal_species(j) && comp.metal_inventory==MetalInventory::gs98) {
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
    if(shared_electrons && !*shared_electrons)*shared_electrons=electrons(T,rho*NA*ye.value);
    const auto f=shared_electrons?**shared_electrons:electrons(T,rho*NA*ye.value);
    out.electron_eta=f.eta;theta.value=f.theta;
    theta.d[0]=f.dtheta_dlnT;theta.d[1]=f.dtheta_dlnne;
    for(std::size_t j=0;j<NSPEC;++j) theta.d[j+2]=f.dtheta_dlnne*ye.d[j+2]/ye.value;
  } else out.electron_eta=std::numeric_limits<double>::quiet_NaN();
  out.electron_susceptibility=theta.value;
  const auto ge=e2/(kB*temp)*exp(log(4*M_PI*NA*density*ye/3)/3);
  out.gamma_e=ge.value;
  const double g12=ge.value*2*r.z1*r.z2/(std::cbrt(r.z1)+std::cbrt(r.z2));
  out.zeta=g12/std::cbrt(gamow_energy(r)/(4*kB*T));
  const double zeta_max=quantum_screening_zeta_max.load(std::memory_order_relaxed);
  const bool quantum=zeta_max>0 && model==PPScreening::salpeter_van_horn;
  // Explicit classical-ion thermonuclear domain, not a cap or an extrapolation.
  if(model!=PPScreening::legacy_weak && !quantum && out.zeta>.2)
    throw std::domain_error("PPChains: classical-ion screening requires zeta<=0.2; quantum burning unavailable");
  if(quantum && (out.zeta>zeta_max || g12>200)) {
    const double limit=quantum_burning_fuel_limit.load(std::memory_order_relaxed);
    const double fuel=comp[Species::H1]+comp[Species::H2]+comp[Species::He3];
    if(limit>0 && fuel<=limit) {
      out.reaction_omitted=true;
      out.log_factor=-std::numeric_limits<double>::infinity();
      return out;
    }
    throw std::domain_error("PPChains: beyond the mean-field quantum screening domain (thermo-pycnonuclear); no extrapolation");
  }
  const auto weak=r.z1*r.z2*e2/(kB*temp)*exp(.5*log(4*M_PI*e2*NA*density*(ions+theta*ye)/(kB*temp)));
  D exponent=weak;
  if(model==PPScreening::legacy_weak && weak.value>=2) exponent=D(2);
  if(model==PPScreening::salpeter_van_horn) {
    const auto strong=.9*ge*(std::pow(r.z1+r.z2,5./3)-std::pow(r.z1,5./3)-std::pow(r.z2,5./3));
    exponent=weak*strong/exp(.5*log(weak*weak+strong*strong));
  }
  if(quantum) {
    // CD09 eqs. 23-25 with the A4 weak-coupling interpolation: only the finite-zeta
    // tunnelling change h(Gamma,zeta)-h(Gamma,0) is added, so the established
    // classical SVH exponent (with electron polarization) is retained as zeta->0.
    const auto g12d=ge*(2*r.z1*r.z2/(std::cbrt(r.z1)+std::cbrt(r.z2)));
    const auto zeta=g12d/detail::cbrt(D(gamow_energy(r)/(4*kB))/temp);
    const double y=4*r.z1*r.z2/((r.z1+r.z2)*(r.z1+r.z2));
    const auto t=detail::cbrt(1+.013*y*y*zeta+.406*std::pow(y,.14)*zeta*zeta
                              +(.062*std::pow(y,.19)+1.8/g12d)*zeta*zeta*zeta);
    const double a=std::pow(r.z1,5./3),b=std::pow(r.z2,5./3),c=std::pow(r.z1+r.z2,5./3);
    auto mixing=[&](const D& scale){return detail::classical_ocp_free_energy(ge*a/scale)+detail::classical_ocp_free_energy(ge*b/scale)-detail::classical_ocp_free_energy(ge*c/scale);};
    const auto moment=ions/ye;   // <Z^2>/<Z> of the ion mixture
    const auto cfac=3*r.z1*r.z2*detail::sqrt(moment)
                    /(std::pow(r.z1+r.z2,2.5)-std::pow(r.z1,2.5)-std::pow(r.z2,2.5));
    const auto weight=(cfac+g12d*g12d)/(1+g12d*g12d);
    exponent=exponent+weight*(mixing(t)-mixing(D(1)));
  }
  out.log_factor=exponent.value;out.dlog_dlnT=exponent.d[0];out.dlog_dlnRho=exponent.d[1];
  for(std::size_t j=0;j<NSPEC;++j) out.dlog_dX[j]=exponent.d[j+2];
  return out;
}

ScreeningState pp_screening(double T,double rho,const Composition& comp,PPReaction which,PPScreening model) {
  return screening_response(T,rho,comp,reaction(which),model,nullptr);
}

NuclearResponse PPChains::composition_response(double T,double rho,const Composition& comp) const {
  validate(T,rho,comp);
  NuclearResponse result;
  if(T<1e5) return result; // retained explicit negligible-burning cutoff
  const std::array rates{pp_bare_rate(T,PPReaction::pp,rates_),pp_bare_rate(T,PPReaction::he3_he3,rates_),
    pp_bare_rate(T,PPReaction::he3_he4,rates_)};
  // The classical helium factors agree; quantum tunnelling also depends on
  // reduced mass, so compute the two factors separately when selected.
  std::optional<Susceptibility> shared_electrons;
  const auto f1=screening_response(T,rho,comp,reaction(PPReaction::pp),screening_,&shared_electrons);
  const auto f2=screening_response(T,rho,comp,reaction(PPReaction::he3_he3),screening_,&shared_electrons);
  const auto f3=quantum_screening_zeta_max.load(std::memory_order_relaxed)>0
    && screening_==PPScreening::salpeter_van_horn
    ?screening_response(T,rho,comp,reaction(PPReaction::he3_he4),screening_,&shared_electrons):f2;
  const std::array screens{f1,f2,f3};
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
  if(rates_==PPRates::solar_fusion_iii) {
    if(screening_==PPScreening::legacy_weak) return "Solar Fusion III rates; capped classical weak screening";
    if(screening_==PPScreening::debye_fermi) return "Solar Fusion III rates; finite-degeneracy Debye screening";
    return "Solar Fusion III rates; finite-degeneracy Salpeter--Van Horn screening";
  }
  if(screening_==PPScreening::legacy_weak) return "Solar Fusion II quadrature; capped classical weak screening";
  if(screening_==PPScreening::debye_fermi) return "Solar Fusion II quadrature; finite-degeneracy Debye screening";
  return "Solar Fusion II quadrature; finite-degeneracy Salpeter--Van Horn screening";
}
ThermonuclearRate cn_bare_rate(double T,PPRates prescription) {
  return cn_bare_rate(T,CNReaction::n14_p,prescription);
}
ThermonuclearRate cn_bare_rate(double T,CNReaction which,PPRates prescription) {
  const auto r=cn_reaction(prescription,which);
  if(!std::isfinite(T) || T<=0 || T>2e7)
    throw std::domain_error("CNCycle: positive finite T<=2e7 K required for low-energy rates");
  if(T<1e5)return {};
  using Entry=std::array<std::optional<ThermonuclearRate>,6>;
  thread_local std::unordered_map<double,Entry> cache;
  const auto index=2*static_cast<std::size_t>(which)+(prescription==PPRates::solar_fusion_iii?1:0);
  const auto found=cache.find(T);
  if(found!=cache.end() && found->second[index])return *found->second[index];
  static const Quadrature q;
  const double kt=kB*T,eg=gamow_energy(r),peak=std::cbrt(eg*kt*kt/4);
  const double mu=r.m1*r.m2/(r.m1+r.m2);
  double integral=0,moment=0;
  for(double mid:{-2.,2.})for(int i=0;i<q.n;++i) {
    const double energy=peak*std::exp(mid+2*q.x[i]),E=energy/mev;
    const double S=(r.s0+E*(r.s1+.5*E*r.s2))*mev*1e-24;
    const double term=2*q.w[i]*energy*S*std::exp(-energy/kt-std::sqrt(eg/energy));
    integral+=term;moment+=term*energy/kt;
  }
  const ThermonuclearRate result{NA*std::sqrt(8/(M_PI*mu))/std::pow(kt,1.5)*integral,
    integral>0?-1.5+moment/integral:0};
  if(found==cache.end() && cache.size()>=8192)cache.clear();
  cache[T][index]=result;
  return result;
}

ScreeningState cn_screening(double T,double rho,const Composition& c,PPScreening model) {
  return cn_screening(T,rho,c,CNReaction::n14_p,model);
}
ScreeningState cn_screening(double T,double rho,const Composition& c,CNReaction which,PPScreening model) {
  return screening_response(T,rho,c,cn_reaction(PPRates::solar_fusion_iii,which),model,nullptr);
}

CNCycle::CNCycle(PPRates rates,PPScreening screening,double converted_carbon)
    :rates_(rates),screening_(screening),converted_carbon_(converted_carbon) {
  (void)cn_reaction(rates);
  if(screening!=PPScreening::debye_fermi && screening!=PPScreening::salpeter_van_horn)
    throw std::invalid_argument("CNCycle: finite-degeneracy screening required");
  if(!std::isfinite(converted_carbon) || converted_carbon<0 || converted_carbon>1)
    throw std::invalid_argument("CNCycle: converted carbon fraction must be in [0,1]");
}

NuclearResponse CNCycle::composition_response(double T,double rho,const Composition& c) const {
  validate(T,rho,c);
  if(c.metal_inventory!=MetalInventory::gs98)
    throw std::domain_error("CNCycle: fixed GS98 catalyst approximation requires GS98 inventory");
  const auto bare=cn_bare_rate(T,rates_);
  NuclearResponse out;
  if(T<1e5)return out;
  const auto scr=cn_screening(T,rho,c,screening_);
  const double catalyst_per_Z=(converted_carbon_*gs98_metals[0].fraction/12
    +gs98_metals[1].fraction/14)/(c.basis==AbundanceBasis::baryon_mass?1.:gs98_atomic_mass_scale());
  const double wH=c.abundance_weight(0),wHe=c.abundance_weight(2);
  const double hydrogen=c.X[0]/wH,catalyst=c.Z()*catalyst_per_Z;
  const double coefficient=rho*bare.molar_rate*std::exp(scr.log_factor);
  const double rate=coefficient*hydrogen*catalyst; // mol of cycles / g / s
  const double q=(4*nuclides[0].A-nuclides[2].A)*c_light*c_light;
  const double nu=(.706+.996)*mev*NA,heat=q-nu;
  auto& s=out.state;
  s.eps=rate*heat;s.eps_neutrino=rate*nu;
  s.dXdt[0]=-4*wH*rate;s.dXdt[2]=wHe*rate;
  if(s.eps>0) {
    s.dlneps_dlnT=bare.dlnrate_dlnT+scr.dlog_dlnT;
    s.dlneps_dlnRho=1+scr.dlog_dlnRho;
  }
  for(std::size_t j=0;j<NSPEC;++j) {
    const double dr=rate*scr.dlog_dX[j]+(j==0?coefficient*catalyst/wH:0)
      +(is_metal_species(j)?coefficient*hydrogen*catalyst_per_Z:0);
    out.deps_dX[j]=heat*dr;
    out.d_dXdt_dX[0][j]=-4*wH*dr;out.d_dXdt_dX[2][j]=wHe*dr;
  }
  return out;
}
NuclearState CNCycle::eval(double T,double rho,const Composition& c) const {
  return composition_response(T,rho,c).state;
}
NuclearResponse PPCNO::composition_response(double T,double rho,const Composition& c) const {
  auto out=pp_.composition_response(T,rho,c);
  const auto cn=cn_.composition_response(T,rho,c);
  auto& s=out.state;
  const double eps=s.eps+cn.state.eps;
  if(eps>0) {
    s.dlneps_dlnT=(s.eps*s.dlneps_dlnT+cn.state.eps*cn.state.dlneps_dlnT)/eps;
    s.dlneps_dlnRho=(s.eps*s.dlneps_dlnRho+cn.state.eps*cn.state.dlneps_dlnRho)/eps;
  }
  s.eps=eps;s.eps_neutrino+=cn.state.eps_neutrino;
  for(std::size_t i=0;i<NSPEC;++i) {
    s.dXdt[i]+=cn.state.dXdt[i];out.deps_dX[i]+=cn.deps_dX[i];
    for(std::size_t j=0;j<NSPEC;++j)out.d_dXdt_dX[i][j]+=cn.d_dXdt_dX[i][j];
  }
  return out;
}
NuclearState PPCNO::eval(double T,double rho,const Composition& c) const {
  return composition_response(T,rho,c).state;
}
} // namespace ember
