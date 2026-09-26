#include "ember/deuterium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ember {
namespace {
using namespace constants;
constexpr double mev = 1e6*eV;
// Acharya et al. 2025, RMP 97, 035002, Table IV. Energy in MeV;
// S in 10^-7 MeV barn. Do not interpret the table units as keV barn.
constexpr std::array energies{0.,.010,.020,.040,.080,.091,.100,.120};
constexpr std::array sfactors{2.028,2.644,3.276,4.579,7.31,8.11,8.77,10.24};
double s_factor(double E) {
  const auto hi=std::upper_bound(energies.begin(),energies.end(),E);
  const auto i=std::clamp<std::ptrdiff_t>(hi-energies.begin()-1,0,energies.size()-2);
  return std::lerp(sfactors[i],sfactors[i+1],
      (E-energies[i])/(energies[i+1]-energies[i]))*1e-7*mev*1e-24;
}
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
}

ThermonuclearRate deuterium_bare_rate(double T) {
  if(!std::isfinite(T) || T<=0 || T>2e7)
    throw std::domain_error("deuterium capture: positive temperature up to 20 MK required");
  if(T<1e4) return {}; // Explicit negligible thermal burning cutoff in cool PMS layers.
  static const Quadrature q;
  const double mp=nuclides[0].A*amu-me,md=deuterium_atomic_mass*amu-me;
  const double mu=mp*md/(mp+md),kt=kB*T,hbar=h/(2*M_PI);
  constexpr double e2=4.803204673e-10*4.803204673e-10;
  const double eg=2*mu*std::pow(M_PI*e2/hbar,2);
  const double peak=std::cbrt(eg*kt*kt/4);
  // Split at the peak and every S-table join. The low boundary is 2981 times
  // below the Gamow peak; its Coulomb suppression is negligible throughout
  // the declared thermal domain. No S-factor extrapolation above 120 keV.
  std::vector<double> cuts{std::log(peak)-8,std::log(peak)};
  for(std::size_t i=1;i<energies.size();++i) cuts.push_back(std::log(energies[i]*mev));
  std::sort(cuts.begin(),cuts.end());
  double integral=0,moment=0;
  for(std::size_t j=1;j<cuts.size();++j) {
    const double mid=(cuts[j]+cuts[j-1])/2,half=(cuts[j]-cuts[j-1])/2;
    for(int i=0;i<q.n;++i) {
      const double E=std::exp(mid+half*q.x[i]);
      const double term=half*q.w[i]*E*s_factor(E/mev)*std::exp(-E/kt-std::sqrt(eg/E));
      integral+=term;moment+=term*E/kt;
    }
  }
  return {NA*std::sqrt(8/(M_PI*mu))/std::pow(kt,1.5)*integral,
          integral>0?-1.5+moment/integral:0};
}

DeuteriumCapture deuterium_capture(double T,double rho,const Composition& other,
    double xd,PPScreening screening) {
  if(other.basis!=AbundanceBasis::baryon_mass || !std::isfinite(xd) || xd<0 ||
      other[Species::H2]!=0 ||
      !std::isfinite(rho) || rho<=0 || std::abs(other.sum()+xd-1)>1e-10)
    throw std::domain_error("deuterium capture: normalized baryon inventory required");
  for(const double x:other.X) if(!std::isfinite(x) || x<0)
    throw std::domain_error("deuterium capture: negative/nonfinite inventory");
  if(!std::isfinite(T) || T<=0)
    throw std::domain_error("deuterium capture: positive finite temperature required");
  DeuteriumCapture out;
  out.heat_per_mole=(nuclides[0].A+deuterium_atomic_mass-nuclides[1].A)*c_light*c_light;
  if(xd==0 || other.X[0]==0) return out;
  const auto bare=deuterium_bare_rate(T);
  // Do not evaluate a fully ionized screening model in cold molecular layers
  // where the explicitly omitted thermal capture rate is negligible.
  if(bare.molar_rate==0) return out;
  // Fully ionized D has exactly the same charge moments as half as much H1
  // baryon mass. This temporary mixture is ONLY for Coulomb screening, never
  // for material lookup or the conserved abundance array.
  auto screening_counts=other;screening_counts.X[0]+=.5*xd;
  const auto s=pp_screening(T,rho,screening_counts,PPReaction::deuterium_p,screening);
  const double r=rho*bare.molar_rate*std::exp(s.log_factor)*other.X[0]*xd/2;
  out.molar_reactions_per_gram_second=r;
  out.dX_deuterium_dt=-2*r;
  out.source.dXdt[0]=-r;out.source.dXdt[1]=3*r;
  out.source.eps=r*out.heat_per_mole;
  out.source.dlneps_dlnT=bare.dlnrate_dlnT+s.dlog_dlnT;
  out.source.dlneps_dlnRho=1+s.dlog_dlnRho;
  // Radiative capture produces no escaping nuclear neutrino.
  return out;
}

NuclearResponse PPDeuterium::composition_response(double T,double rho,const Composition& c) const {
  if(c.basis!=AbundanceBasis::baryon_mass)
    throw std::invalid_argument("PMS network: baryon abundances required");
  // This is only the light-isotope source. A separate CN source can use the
  // same composition; its catalysts are neither consumed nor reset here.
  auto out=pp_.composition_response(T,rho,c);
  const auto bare=deuterium_bare_rate(T);
  if(bare.molar_rate==0)return out;
  constexpr auto d=static_cast<std::size_t>(Species::H2);
  const auto screen=pp_screening(T,rho,c,PPReaction::deuterium_p,screening_);
  const double coefficient=rho*bare.molar_rate*std::exp(screen.log_factor);
  const double rate=coefficient*c.X[0]*c.X[d]/2;
  const double Q=(nuclides[0].A+deuterium_atomic_mass-nuclides[1].A)*c_light*c_light;
  const double heat=rate*Q,pp_heat=out.state.eps;
  out.state.eps+=heat;
  if(out.state.eps>0) {
    out.state.dlneps_dlnT=(pp_heat*out.state.dlneps_dlnT+
        heat*(bare.dlnrate_dlnT+screen.dlog_dlnT))/out.state.eps;
    out.state.dlneps_dlnRho=(pp_heat*out.state.dlneps_dlnRho+
        heat*(1+screen.dlog_dlnRho))/out.state.eps;
  }
  out.state.dXdt[0]-=rate;out.state.dXdt[1]+=3*rate;out.state.dXdt[d]-=2*rate;
  for(std::size_t j=0;j<NSPEC;++j) {
    const double dr=rate*screen.dlog_dX[j]+coefficient*
        ((j==0?c.X[d]/2:0)+(j==d?c.X[0]/2:0));
    out.deps_dX[j]+=dr*Q;
    out.d_dXdt_dX[0][j]-=dr;out.d_dXdt_dX[1][j]+=3*dr;out.d_dXdt_dX[d][j]-=2*dr;
  }
  return out;
}
NuclearState PPDeuterium::eval(double T,double rho,const Composition& c) const {
  return composition_response(T,rho,c).state;
}
} // namespace ember
