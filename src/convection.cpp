#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {

CompositionBuoyancy composition_buoyancy(const Eos& eos,double T,double P,double delta,
    double contrast,const Composition& lo,const Composition& hi,double guess) {
  if(!std::isfinite(T+P+delta+contrast) || T<=0 || P<=0 || delta<=0 || lo.basis!=hi.basis
      || lo.metal_inventory!=hi.metal_inventory)
    throw std::domain_error("composition_buoyancy: invalid state or abundance basis");
  CompositionBuoyancy out{};
  if(lo.X==hi.X)return out;
  if(std::abs(contrast)<32*std::numeric_limits<double>::epsilon())
    throw std::domain_error("composition_buoyancy: unresolved pressure contrast across a composition gradient");
  const double rlo=eos.rho_from_PT(T,P,lo,guess),rhi=eos.rho_from_PT(T,P,hi,guess);
  const auto a=eos.eval(T,rlo,lo),b=eos.eval(T,rhi,hi);
  const double factor=1/(delta*contrast);
  out.B=std::log1p((rhi-rlo)/rlo)*factor;
  out.dB_dlnT=(a.delta-b.delta)*factor;
  out.dB_dlnP=(1/b.chiRho-1/a.chiRho)*factor;
  out.dB_ddelta=-out.B/delta;
  out.dB_dpressure_contrast=-out.B/contrast;
  return out;
}

ConvectionState ledoux_mixing_length_gradient(double rad,double ad,double B,double U) {
  if(!std::isfinite(B))throw std::domain_error("ledoux_mixing_length_gradient: invalid buoyancy");
  auto result=mixing_length_gradient(rad-B,ad,U);
  if(!result.unstable)result.grad=result.grad_element=rad;
  else {result.grad+=B;result.grad_element+=B;result.superadiabaticity+=B;}
  return result;
}

ConvectionState mixing_length_gradient(double grad_rad, double grad_ad, double U) {
  if (!std::isfinite(grad_rad) || !std::isfinite(grad_ad) || grad_ad < 0.0
      || !std::isfinite(U) || !(U > 0.0))
    throw std::domain_error("mixing_length_gradient: invalid gradient or U");

  ConvectionState s{};
  s.grad = s.grad_element = grad_rad;
  if (grad_rad <= grad_ad) return s;
  s.unstable = true;

  // Put q = sqrt(grad - grad_element), W = grad_rad - grad_ad. The
  // element's cooling and total flux conservation give
  //
  //   W = q^2 + 2*U*q + 9*q^3/(8*U).
  //
  // Solve for q directly: xi=sqrt(grad-grad_ad+U^2) loses xi-U at large U.
  // Scale q by the smallest positive root of each individual term = W.
  // Then A*t^2 + B*t + C*t^3 = 1, with A,B,C <= 1 and a root in [1/3,1].
  // Logarithms avoid overflow in U^2 or U*W at either efficiency extreme.
  const double W = grad_rad - grad_ad;
  const double logW = std::log(W), logU = std::log(U);
  const double logq2 = 0.5 * logW;
  const double logq1 = logW - std::log(2.0) - logU;
  const double logq3 = (std::log(8.0 / 9.0) + logU + logW) / 3.0;
  const double logq = std::min({logq2, logq1, logq3});
  // One coefficient is exactly one, so roundoff cannot move the root above
  // the upper bracket when the other two terms are vanishingly small.
  const double A = std::exp(2.0 * (logq - logq2));
  const double B = std::exp(logq - logq1);
  const double C = std::exp(3.0 * (logq - logq3));

  double lo = 0.0, hi = 1.0, t = 0.5;
  bool converged = false;
  for (int it = 0; it < 80; ++it) {
    const double f = ((C * t + A) * t + B) * t - 1.0;
    if (std::abs(f) <= 8.0 * std::numeric_limits<double>::epsilon()) {
      converged = true;
      break;
    }
    if (f > 0.0) hi = t; else lo = t;
    const double slope = (3.0 * C * t + 2.0 * A) * t + B;
    const double next = t - f / slope;
    t = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
  }
  if (!converged)
    throw std::runtime_error("mixing_length_gradient: cubic did not converge");

  // Fractions of W in the element contrast, element's cooling, and
  // convective flux. Normalize away the root solve's last rounding error.
  const double total = ((C * t + A) * t + B) * t;
  const double a = A * t * t / total;
  const double b = B * t / total;
  const double c = C * t * t * t / total;
  s.element_contrast = W * a;
  s.superadiabaticity = W * (a + b);
  s.convective_excess = W * c;
  s.grad = std::clamp(grad_ad + s.superadiabaticity, grad_ad, grad_rad);
  s.grad_element = std::clamp(grad_ad + W * b, grad_ad, s.grad);

  // Implicit differentiation of the same positive-term cubic. These forms
  // avoid differences of nearly equal quantities in both limiting regimes.
  const double denom = 2.0 * a + b + 3.0 * c;
  s.dgrad_dgrad_rad = (2.0 * a + b) / denom;
  s.dgrad_dgrad_ad = 3.0 * c / denom;
  s.dgrad_dlnU = W * c * (2.0 * a + 4.0 * b) / denom;
  return s;
}

double mixing_length_U(double T, double rho, double kappa, double gravity,
                       const EosState& eos, double alpha) {
  for (double x : {T, rho, kappa, gravity, eos.P, eos.cp, eos.delta, alpha})
    if (!std::isfinite(x) || !(x > 0.0))
      throw std::domain_error("mixing_length_U: inputs must be positive and finite");

  const double logHp = std::log(eos.P) - std::log(rho) - std::log(gravity);
  const double logU = std::log(3.0 * constants::a_rad * constants::c)
      + 3.0 * std::log(T) - std::log(eos.cp) - 2.0 * std::log(rho)
      - std::log(kappa) - 2.0 * std::log(alpha) - 1.5 * logHp
      + 0.5 * (std::log(8.0) - std::log(gravity) - std::log(eos.delta));
  const double U = std::exp(logU);
  if (!std::isfinite(U) || !(U > 0.0))
    throw std::domain_error("mixing_length_U: U outside representable range");
  return U;
}

} // namespace ember
