#include "ember/viscosity.hpp"
#include "ember/constants.hpp"
#include "ember/detail/differential.hpp"
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace ember {
namespace {
using D=detail::Differential<4>;
constexpr double pi=std::numbers::pi, electron_charge=4.803204712570263e-10;
constexpr double hbar=constants::h/(2*pi), e2=electron_charge*electron_charge;
constexpr double alpha=e2/(hbar*constants::c);
struct Quadrature {
  std::array<double,32> x{},w{};
  Quadrature() {
    for(int i=0;i<16;++i) {
      double z=std::cos(pi*(i+.75)/32.5),derivative=0;
      for(int iteration=0;iteration<20;++iteration) {
        double a=1,b=z;
        for(int j=2;j<=32;++j) {const double n=((2*j-1)*z*b-(j-1)*a)/j;a=b;b=n;}
        derivative=32*(z*b-a)/(z*z-1);
        const double correction=b/derivative;z-=correction;
        if(std::abs(correction)<2e-15)break;
      }
      x[i]=.5*(1-z);x[31-i]=.5*(1+z);
      w[i]=w[31-i]=1/((1-z*z)*derivative*derivative);
    }
  }
};
D expm1(const D& x) {
  D y(std::expm1(x.value));
  for(std::size_t i=0;i<4;++i)y.d[i]=std::exp(x.value)*x.d[i];
  return y;
}
// Equation 12 with the effective potential (20). Integrate in
// u=ln(1+t/s), t=(q/2kF)^2, to resolve the screened forward peak.
D coulomb_integral(const D& s,const D& w,const D& beta) {
  static const Quadrature rule;
  const auto upper=detail::log1p(1/s);D sum;
  for(std::size_t j=0;j<rule.x.size();++j) {
    const auto u=upper*rule.x[j], t=s*expm1(u);
    sum+=rule.w[j]*.5*(1-detail::exp(-u))*(1-t)*(1-beta*beta*t)*(-expm1(-w*t));
  }
  return upper*sum;
}
}

LiquidViscosityResponse ocp_liquid_viscosity(double temperature,double density,
                                          double charge,double mass_number) {
  const std::array<double,4> input{temperature,density,charge,mass_number};
  for(double v:input)if(!std::isfinite(v)||v<=0)
    throw std::invalid_argument("liquid viscosity: require positive finite T, rho, Z, A");
  if(charge>mass_number)
    throw std::invalid_argument("liquid viscosity: ion charge exceeds mass number");
  D T(temperature),rho(density),Z(charge),A(mass_number);
  T.d[0]=temperature;rho.d[1]=density;Z.d[2]=charge;A.d[3]=mass_number;
  const auto ne=rho*Z/(A*constants::amu),kf=detail::cbrt(3*pi*pi*ne);
  const auto pf=hbar*kf,x=pf/(constants::me*constants::c),rel=detail::sqrt(1+x*x);
  const auto mstar=constants::me*rel,v=pf/mstar,beta=v/constants::c;
  const auto TF=constants::me*constants::c*constants::c/constants::kB*x*x/(rel+1);
  const auto ai=detail::cbrt(3*Z/(4*pi*ne)),g=Z*Z*e2/(ai*constants::kB*T);
  const auto omega=detail::sqrt(4*pi*ne*Z*e2/(A*constants::amu));
  const auto theta=T*constants::kB/(hbar*omega);
  if((T/TF).value>.05 || x.value>.5 || g.value>175)
    throw std::domain_error("liquid viscosity: outside nonrelativistic degenerate liquid domain");
  const auto sd=3*g/(4*kf*kf*ai*ai),bz=pi*alpha*Z*beta;
  const auto s=(sd*(1+.06*g)*detail::exp(-detail::sqrt(g))+alpha/(pi*beta))*detail::exp(-bz);
  const auto w=13/sd*(1+bz/3);
  const auto G=(1+.122*bz*bz)/detail::sqrt(1+.0361*detail::pow(Z,-1./3)/(theta*theta));
  const auto debye=detail::exp(-.42*2.8*detail::sqrt(x/(A*Z))*detail::exp(-9.1*theta));
  const auto log=coulomb_integral(s,w,beta)*G*debye;
  // e^4 is required here; the printed e^2 in one version of Eq. 11 is a typo
  // (its equivalent alpha^2 expression and dimensions both require e^4).
  const auto ei=12*pi*Z*e2*e2*log*ne/(pf*pf*v);
  const auto ee=5*pi*pi*alpha*alpha*detail::pow(constants::kB*T,2)
      /(2*mstar*constants::c*constants::c*hbar)
      /(2*detail::sqrt(alpha/(pi*beta)))*(1+6/(5*x*x)+2/(5*detail::pow(x,4)));
  const auto electron=ne*pf*v/(5*(ei+ee));
  // Daligault et al., table 4, in units rho*a^2*omega_p.
  const auto fit=.794811/(detail::pow(g,2.5)*detail::log1p(.862151/detail::pow(g,1.5)))
      *(1+.0425698*g+.00205782*g*g+7.03658e-5*detail::pow(g,3))
      /(1+.0429942*g-.000270798*g*g+3.25441e-6*detail::pow(g,3)-1.15019e-8*detail::pow(g,4));
  const auto ion=rho*ai*ai*omega*fit;
  if(!(electron.value>0 && ion.value>0) || !std::isfinite(electron.value+ion.value))
    throw std::domain_error("liquid viscosity: invalid coefficient");
  LiquidViscosityResponse result;
  result.value={electron.value,ion.value,ei.value,ee.value,log.value,g.value,
                1/theta.value,(T/TF).value};
  result.partials={electron.d,ion.d};
  return result;
}
} // namespace ember
