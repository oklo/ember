#include "ember/fingering.hpp"
#include "differential.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
struct Mode {
  double lambda{}, derivative_numerator{};
};

Mode mode(double q, double pr, double tau, double inverse_R) {
  const double b2=1+pr+tau, b1=tau*pr+pr+tau;
  const double a2=q*b2, a1=q*q*b1+pr*(1-inverse_R);
  const double driving=pr*(inverse_R-tau);
  const double a0=q*(q*q*tau*pr-driving);
  if(a0>=0 || q==0)
    return {0, 3*q*q*tau*pr-driving};

  // On R0>1 this cubic is strictly increasing for lambda>=0 and has a
  // unique positive root. Each bound drops other positive polynomial terms.
  // A bracketed Newton solve avoids cancellation in closed-form cubics.
  double lo=0, hi=std::min({-a0/a1,std::sqrt(-a0/a2),std::cbrt(-a0)});
  double x=hi;
  constexpr double tolerance=32*std::numeric_limits<double>::epsilon();
  bool converged=false;
  for(int i=0;i<100;++i) {
    const double positive=((x+a2)*x+a1)*x;
    const double f=positive+a0;
    if(std::abs(f)<=tolerance*(positive-a0)) {converged=true;break;}
    if(f>0)hi=x;else lo=x;
    const double slope=(3*x+2*a2)*x+a1;
    const double candidate=x-f/slope;
    x=(candidate>lo && candidate<hi)?candidate:.5*(lo+hi);
  }
  if(!converged || !(x>0) || !std::isfinite(x))
    throw std::runtime_error("fingering: positive growth root did not converge");
  return {x, b2*x*x+2*q*b1*x+3*q*q*tau*pr-driving};
}
} // namespace

FingeringFlux brown_fingering_flux(double pr, double tau, double R) {
  if(!std::isfinite(pr) || !std::isfinite(tau) || !std::isfinite(R)
      || pr<=0 || tau<=0 || tau>=1 || R<=0)
    throw std::invalid_argument("fingering: require finite Pr>0, 0<tau<1, R0>0");
  FingeringFlux result;
  result.within_calibrated_ratio=tau<=pr;
  if(R<=1) {result.regime=FingeringRegime::overturning;return result;}
  const double inverse_R=1/R;
  if(inverse_R<=tau)return result;
  const double critical_q=std::sqrt((inverse_R-tau)/tau);
  if(!std::isfinite(critical_q) || critical_q<=0)
    throw std::domain_error("fingering: parameters exceed resolved numerical range");

  // At the fastest growing mode d(lambda)/dq=0, so the numerator of the
  // implicit derivative vanishes (equation 20). Bracket in q=l^2; do not
  // perturb physical inputs to make a Newton iteration converge.
  double lo=0, hi=std::min(1.,critical_q);
  while(mode(hi,pr,tau,inverse_R).derivative_numerator<0) {
    const double next=std::min(2*hi,critical_q);
    if(next<=hi)throw std::runtime_error("fingering: cannot bracket fastest mode");
    hi=next;
  }
  for(int i=0;i<80;++i) {
    const double q=.5*(lo+hi);
    if(mode(q,pr,tau,inverse_R).derivative_numerator<0)lo=q;else hi=q;
    if(hi-lo<=2e-12*(hi+lo))break;
  }
  const double q=.5*(lo+hi), lambda=mode(q,pr,tau,inverse_R).lambda;
  result.regime=FingeringRegime::fingering;
  result.growth_rate=lambda;result.wavenumber_squared=q;
  result.thermal_nusselt_excess=49*(lambda/q)*(lambda/(lambda+q));
  result.chemical_nusselt_excess=49*(lambda/q)*(lambda/(lambda+tau*q))/tau;
  if(!(lambda>0 && result.thermal_nusselt_excess>0
      && result.chemical_nusselt_excess>0)
      || !std::isfinite(result.thermal_nusselt_excess)
      || !std::isfinite(result.chemical_nusselt_excess))
    throw std::runtime_error("fingering: unresolved or nonfinite transport");
  return result;
}

FingeringResponse brown_fingering_response(double pr, double tau, double R) {
  FingeringResponse result;
  result.value=brown_fingering_flux(pr,tau,R);
  if(result.value.regime!=FingeringRegime::fingering)return result;
  using D=detail::Differential<3>;
  D p(pr),t(tau),inv(1/R);
  p.d[0]=pr;t.d[1]=tau;inv.d[2]=-1/R;
  const double lambda=result.value.growth_rate,q=result.value.wavenumber_squared;
  const D b2=1+p+t,b1=t*p+p+t,driving=p*(inv-t);
  // F(lambda,q)=0 and F_q(lambda,q)=0. At the stationary mode,
  // F_q vanishes, so the first equation gives lambda's parameter derivative
  // without differencing nearby eigenvalues or subtracting growing roots.
  const D F=lambda*lambda*lambda+q*b2*lambda*lambda
      +(q*q*b1+p*(1-inv))*lambda+q*(q*q*t*p-driving);
  const D Fq=b2*lambda*lambda+2*q*b1*lambda+3*q*q*t*p-driving;
  const double Flambda=3*lambda*lambda+2*q*b2.value*lambda+q*q*b1.value+pr*(1-1/R);
  const double Fqlambda=2*b2.value*lambda+2*q*b1.value;
  const double Fqq=2*b1.value*lambda+6*q*tau*pr;
  if(!(Flambda>0 && Fqq>0))throw std::runtime_error("fingering: singular response");
  D growth(lambda),wavenumber(q);
  for(std::size_t j=0;j<3;++j) {
    growth.d[j]=-F.d[j]/Flambda;
    wavenumber.d[j]=-(Fq.d[j]+Fqlambda*growth.d[j])/Fqq;
  }
  const D heat=49*(growth/wavenumber)*(growth/(growth+wavenumber));
  const D chemical=49*(growth/wavenumber)*(growth/(growth+t*wavenumber))/t;
  result.partials={growth.d,wavenumber.d,heat.d,chemical.d};
  for(const auto& row:result.partials)for(double x:row)
    if(!std::isfinite(x))throw std::runtime_error("fingering: nonfinite response");
  result.derivatives_defined=true;
  return result;
}
} // namespace ember
