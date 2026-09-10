#include "ember/losses.hpp"
#include "differential.hpp"

namespace ember {
LossState PlasmaNeutrinoLosses::eval(double T,double rho,const Composition& comp) const {
  const double ye=comp.mu_elec_inv();
  if(!std::isfinite(T) || !std::isfinite(rho) || !std::isfinite(ye) || T<=0 || rho<=0 || ye<=0)
    throw std::domain_error("PlasmaNeutrinoLosses: positive finite state and electron count required");
  for(double x:comp.X) if(!std::isfinite(x) || x<0)
    throw std::domain_error("PlasmaNeutrinoLosses: invalid abundance");
  using D=detail::Differential<2>;
  using detail::exp;using detail::log;
  const auto lt=D::variable(std::log(T),0),lr=D::variable(std::log(rho),1);
  const auto lry=lr+std::log(ye);
  const auto lg=.5*(std::log(1.1095e11)+lry-2*lt
      -.5*log(1+exp((2./3)*(std::log(1.019e-6)+lry))));
  const auto gamma=exp(lg);
  const auto ft=2.4+.6*exp(.5*lg)+.51*gamma+1.25*exp(1.5*lg);
  const auto fl=(8.6*exp(2*lg)+1.35*exp(3.5*lg))/(225-17*gamma+exp(2*lg));
  const auto x=(17.5+(std::log(2.)+lry-3*lt)/std::log(10.))/6;
  const auto y=(-24.5+(std::log(2.)+lry+3*lt)/std::log(10.))/6;
  D fxy(1);
  if(std::abs(x.value)<=.7 && y.value>=0) {
    const auto angle=4.5*x;
    D sine(std::sin(angle.value));
    for(std::size_t i=0;i<2;++i)sine.d[i]=std::cos(angle.value)*angle.d[i];
    const auto shift=angle+.9;
    auto numerator=y-1.6+1.25*x;
    if(numerator.value>0)numerator=D(0);
    const auto argument=numerator/(.57-.25*x);
    fxy=1.05+(.39-1.25*x-.35*sine-.3*exp(-shift*shift))*exp(-argument*argument);
  }
  // Sum C_V^2 for three flavours with sin^2(theta_W)=.23 is .9248.
  // Evaluate the exponential in log space to retain very small cooling rates.
  const auto logarithm=std::log(.9248*3e21)+9*(std::log(1.686e-10)+lt)
      +6*lg-gamma+log(ft+fl)+log(fxy)-lr;
  const double epsilon=std::exp(logarithm.value);
  if(!std::isfinite(epsilon) || !std::isfinite(logarithm.d[0]) || !std::isfinite(logarithm.d[1]))
    throw std::domain_error("PlasmaNeutrinoLosses: nonfinite analytic fit");
  return {epsilon,epsilon>0?logarithm.d[0]:0.,epsilon>0?logarithm.d[1]:0.};
}
} // namespace ember
