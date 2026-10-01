#include "ember/fingering.hpp"
#include "ember/detail/differential.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>
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

TwoCompositionFingering two_composition_fingering(double pr,
    std::array<double,2> tau,std::array<double,2> driving) {
  if(!(pr>0 && tau[0]>0 && tau[1]>=tau[0] && tau[1]<1)
      || !std::isfinite(pr+tau[0]+tau[1]+driving[0]+driving[1]))
    throw std::invalid_argument("two-composition fingering: invalid diffusivity or driving");
  TwoCompositionFingering out;
  const double net=driving[0]+driving[1];
  if(net>=1){out.regime=FingeringRegime::overturning;return out;}
  if(driving[0]<=0 && driving[1]<=0)return out;
  if(tau[0]==tau[1]) {
    if(net<=0)return out;
    const auto f=brown_fingering_flux(pr,tau[0],1/net);
    out.regime=f.regime;out.growth_rate=f.growth_rate;
    out.wavenumber_squared=f.wavenumber_squared;out.thermal_nusselt_excess=f.thermal_nusselt_excess;
    out.mixing_over_thermal.fill(tau[0]*f.chemical_nusselt_excess);return out;
  }
  // The quartic follows by multiplying
  // lambda+Pr*q+Pr/(lambda+q)-Pr*sum(gamma_i/(lambda+tau_i*q))=0
  // by its three positive denominators. Net driving <1 makes its second
  // derivative positive for lambda>=0: there is at most one positive minimum.
  // Solve on the rising side of that minimum, without complex root formulas.
  const double sum=tau[0]+tau[1],product=tau[0]*tau[1];
  const double cross=driving[0]*tau[1]+driving[1]*tau[0];
  const double A=1+pr+sum,B=(1+pr)*sum+pr+product;
  const double C=(1+pr)*product+pr*sum,D=sum-net-cross,E=product-cross;
  const double upper=std::sqrt(pr*(std::max(0.,driving[0])+std::max(0.,driving[1])));
  auto at=[&](double q) -> Mode {
    const double a3=q*A,a2=q*q*B+pr*(1-net),a1=q*(q*q*C+pr*D);
    const double a0=pr*q*q*(product*q*q+E);
    auto polynomial=[&](double x){return (((x+a3)*x+a2)*x+a1)*x+a0;};
    auto slope=[&](double x){return ((4*x+3*a3)*x+2*a2)*x+a1;};
    double lo=0;
    if(a1<0) {
      double a=0,b=upper;
      for(int i=0;i<55;++i) {const double x=.5*(a+b);if(slope(x)<0)a=x;else b=x;}
      lo=.5*(a+b);
    }
    if(polynomial(lo)>=0)return {};
    double hi=upper,x=.5*(lo+hi);
    if(polynomial(hi)<0)throw std::runtime_error("two-composition fingering: unresolved root bound");
    for(int i=0;i<100;++i) {
      const double f=polynomial(x),df=slope(x);
      if(f>0)hi=x;else lo=x;
      if(hi-lo<=8*std::numeric_limits<double>::epsilon()*(hi+lo))break;
      const double next=x-f/df;
      x=next>lo && next<hi?next:.5*(lo+hi);
    }
    x=.5*(lo+hi);
    const double dq=A*x*x*x+2*q*B*x*x+(3*q*q*C+pr*D)*x
        +pr*(4*product*q*q*q+2*E*q);
    return {x,dq};
  };
  double best=0,bestq=0;
  if(driving[0]<0) {
    // Opposite buoyancy signs do not identify the fastest mode: stationary
    // and oscillatory branches can coexist. Routh-Hurwitz applied after a
    // real shift gives the spectral abscissa without complex quartic roots.
    const long double h2=A*B-C,h0=A*pr*(1-net)-pr*D;
    const long double aa=C*h2-A*A*pr*product;
    const long double bb=C*h0+pr*D*h2-A*A*pr*E,cc=pr*D*h0;
    if(D>=0 && E>=0 && h2>=0 && h0>=0 && aa>0
        && (bb>=0?cc>=0:4*aa*cc>=bb*bb))return out;
    if(!(h2>0 && aa>0))
      throw std::runtime_error("two-composition fingering: unresolved stability bound");
    // Above these four determinant roots every mode is damped. The last
    // determinant is quadratic in q^2; use a cancellation-safe upper root.
    long double end2=std::max<long double>({0.L,-pr*D/C,-E/product,-h0/h2});
    const long double discriminant=bb*bb-4*aa*cc;
    if(discriminant>=0) {
      const long double root=std::sqrt(discriminant);
      const long double largest=bb<0?(-bb+root)/(2*aa):
          (bb+root>0?-2*cc/(bb+root):0);
      end2=std::max(end2,largest);
    }
    if(!(end2>0))
      throw std::runtime_error("two-composition fingering: unresolved unstable range");
    auto growth=[&](double q) {
      const long double a3=q*A,a2=q*q*B+pr*(1-net),a1=q*(q*q*C+pr*D);
      const long double a0=pr*q*q*(product*q*q+E);
      auto stable_shift=[&](long double x) {
        const long double b3=a3+4*x,b2=a2+3*a3*x+6*x*x;
        const long double b1=a1+x*(2*a2+x*(3*a3+4*x));
        const long double b0=a0+x*(a1+x*(a2+x*(a3+x)));
        return b1>0 && b0>0 && b3*b2>b1 &&
            b1*(b3*b2-b1)>b3*b3*b0;
      };
      if(stable_shift(0))return 0.;
      const double stationary=at(q).lambda;
      if(stationary>0 && stable_shift(stationary*(1+1e-10)))return stationary;
      double lo=stationary,hi=upper;
      if(!stable_shift(hi))
        throw std::runtime_error("two-composition fingering: unresolved spectral bound");
      for(int j=0;j<70;++j) {
        const double x=.5*(lo+hi);
        if(stable_shift(x))hi=x;else lo=x;
        if(hi-lo<=2e-13*(hi+lo))break;
      }
      return .5*(lo+hi);
    };
    const double qmax=std::sqrt(double(end2)),qmin=std::min(1.,qmax)*1e-10;
    constexpr int samples=128;
    std::array<double,samples+1> logs{},rates{};
    for(int i=0;i<=samples;++i) {
      logs[i]=std::log(qmin)+double(i)/samples*std::log(qmax/qmin);
      rates[i]=growth(std::exp(logs[i]));
    }
    for(int i=1;i<samples;++i)if(rates[i]>0 && rates[i]>=rates[i-1] && rates[i]>=rates[i+1]) {
      double a=logs[i-1],b=logs[i+1];
      constexpr double ratio=.6180339887498948482;
      double x=b-ratio*(b-a),y=a+ratio*(b-a),gx=growth(std::exp(x)),gy=growth(std::exp(y));
      for(int j=0;j<70 && b-a>2e-10;++j) {
        if(gx<gy) {a=x;x=y;gx=gy;y=a+ratio*(b-a);gy=growth(std::exp(y));}
        else {b=y;y=x;gy=gx;x=b-ratio*(b-a);gx=growth(std::exp(x));}
      }
      const double q=std::exp(.5*(a+b)),value=growth(q);
      if(value>best){best=value;bestq=q;}
    }
    if(!(best>0))throw std::runtime_error("two-composition fingering: spectral maximum not resolved");
    const double stationary=at(bestq).lambda;
    if(stationary<best*(1-1e-7)) {
      std::ostringstream message;
      message<<std::setprecision(17)<<"two-composition fingering: fastest mode is oscillatory; closure unavailable"
          <<" Pr="<<pr<<" tau0="<<tau[0]<<" tau1="<<tau[1]
          <<" g0="<<driving[0]<<" g1="<<driving[1]<<" growth="<<best<<" q="<<bestq;
      throw std::domain_error(message.str());
    }
    // A value maximum alone leaves O(sqrt(epsilon)) wavelength noise. Refine
    // F_q=0 so composition finite differences see smooth transport coefficients.
    double a=bestq*.999,b=bestq*1.001;
    if(!(at(a).derivative_numerator<0 && at(b).derivative_numerator>0))
      throw std::runtime_error("two-composition fingering: stationary maximum not bracketed");
    for(int j=0;j<60;++j) {
      const double q=.5*(a+b);
      if(at(q).derivative_numerator<0)a=q;else b=q;
      if(b-a<=2e-12*(a+b))break;
    }
    bestq=.5*(a+b);best=at(bestq).lambda;
  } else {
    const double drive_weight=driving[0]/tau[0]+driving[1]/tau[1];
    if(drive_weight<=1)return out;
    const double qmax=std::sqrt(drive_weight-1);
    // Resolve every sampled local maximum, rather than assuming the two
    // diffusivities create a single extremum. Refine by the analytic F_q=0.
    constexpr int samples=96;
    const double qmin=std::min(1.,qmax)*1e-8;
    double previous_q=qmin;auto previous=at(previous_q);
    for(int i=1;i<=samples;++i) {
      const double q=std::exp(std::log(qmin)+(std::log(qmax)-std::log(qmin))*double(i)/samples);
      const auto current=at(q);
      if(previous.lambda>0 && previous.derivative_numerator<0
          && (current.lambda==0 || current.derivative_numerator>=0)) {
        double a=previous_q,b=q;
        for(int j=0;j<60;++j) {
          const double middle=.5*(a+b);const auto m=at(middle);
          if(m.lambda>0 && m.derivative_numerator<0)a=middle;else b=middle;
          if(b-a<=2e-12*(a+b))break;
        }
        const double peakq=.5*(a+b),growth=at(peakq).lambda;
        if(growth>best){best=growth;bestq=peakq;}
      }
      previous=current;previous_q=q;
    }
    if(!(best>0 && bestq>0))throw std::runtime_error("two-composition fingering: fastest mode not resolved");
  }
  out.regime=FingeringRegime::fingering;out.growth_rate=best;out.wavenumber_squared=bestq;
  const double velocity_squared=49*best*best/bestq;
  out.thermal_nusselt_excess=velocity_squared/(best+bestq);
  for(std::size_t j=0;j<2;++j)out.mixing_over_thermal[j]=velocity_squared/(best+tau[j]*bestq);
  return out;
}
} // namespace ember
