#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

namespace ember {

namespace {
// L=ln(rho_hi/rho_lo) at (T,P) and its exact first derivatives.
struct Contrast {double L{},dlnT{},dlnP{},rho_lo{},rho_hi{},delta_lo{},delta_hi{},chi_lo{},chi_hi{};};
Contrast exact_contrast(const Eos& eos,double T,double P,const Composition& lo,const Composition& hi,double guess) {
  const double rlo=eos.rho_from_PT(T,P,lo,guess),rhi=eos.rho_from_PT(T,P,hi,guess);
  const auto a=eos.eval(T,rlo,lo),b=eos.eval(T,rhi,hi);
  return {std::log1p((rhi-rlo)/rlo),a.delta-b.delta,1/b.chiRho-1/a.chiRho,rlo,rhi,a.delta,b.delta,a.chiRho,b.chiRho};
}
struct Key {
  std::array<std::uint64_t,2*NSPEC+12> bits{};
  bool operator==(const Key&) const=default;
};
struct KeyHash {
  std::size_t operator()(const Key& k) const {
    std::uint64_t h=1469598103934665603ull;
    for(auto b:k.bits){h^=b;h*=1099511628211ull;h^=h>>29;}
    return static_cast<std::size_t>(h);
  }
};
struct Shard {std::mutex mutex;std::unordered_map<Key,Contrast,KeyHash> anchors;};
constexpr std::size_t shard_count=64,shard_capacity=4096;
std::atomic<double> reuse_spacing{0};
std::atomic<bool> reuse_verify{false};
std::atomic<std::size_t> reuse_hits{0},reuse_exact{0},reuse_verified{0};
std::mutex worst_mutex;double worst_error=0;
Shard& shard(std::size_t i) {static Shard shards[shard_count];return shards[i%shard_count];}
std::uint64_t bits(double x) {return std::bit_cast<std::uint64_t>(x);}
Key make_key(const Eos& eos,const Composition& lo,const Composition& hi,std::int64_t iT,std::int64_t iP) {
  Key k;std::size_t n=0;
  for(double x:lo.X)k.bits[n++]=bits(x);
  for(double x:hi.X)k.bits[n++]=bits(x);
  std::uint64_t cn=0;
  for(const auto* c:{&lo,&hi})if(c->cn_molality)for(double x:*c->cn_molality)cn=(cn^bits(x))*1099511628211ull;
  k.bits[n++]=cn;
  k.bits[n++]=static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(&eos));
  k.bits[n++]=static_cast<std::uint64_t>(iT);k.bits[n++]=static_cast<std::uint64_t>(iP);
  k.bits[n++]=static_cast<std::uint64_t>(lo.basis)|static_cast<std::uint64_t>(lo.metal_inventory)<<8
      |static_cast<std::uint64_t>(lo.cn_mass_convention)<<16;
  k.bits[n++]=bits(reuse_spacing.load(std::memory_order_relaxed));
  k.bits[n++]=static_cast<std::uint64_t>(hi.basis)|static_cast<std::uint64_t>(hi.metal_inventory)<<8
      |static_cast<std::uint64_t>(hi.cn_mass_convention)<<16;
  return k;
}
} // namespace

void set_composition_buoyancy_reuse(double h,bool verify) {
  if(!std::isfinite(h) || h<0 || h>.01)throw std::invalid_argument("composition buoyancy reuse spacing must lie in [0,0.01]");
  // Configure before starting evaluations; also clears EOS pointer identities
  // from a previous calculation in this process.
  for(std::size_t i=0;i<shard_count;++i) {auto& s=shard(i);std::lock_guard lock(s.mutex);s.anchors.clear();}
  reuse_hits=0;reuse_exact=0;reuse_verified=0;
  {std::lock_guard lock(worst_mutex);worst_error=0;}
  reuse_spacing=h;reuse_verify=verify;
}
CompositionBuoyancyReuse composition_buoyancy_reuse_statistics() {
  std::lock_guard lock(worst_mutex);
  return {reuse_hits.load(),reuse_exact.load(),reuse_verified.load(),worst_error};
}

CompositionBuoyancy composition_buoyancy(const Eos& eos,double T,double P,double delta,
    double contrast,const Composition& lo,const Composition& hi,double guess) {
  if(!std::isfinite(T+P+delta+contrast) || T<=0 || P<=0 || delta<=0 || lo.basis!=hi.basis
      || lo.metal_inventory!=hi.metal_inventory)
    throw std::domain_error("composition_buoyancy: invalid state or abundance basis");
  CompositionBuoyancy out{};
  if(lo.X==hi.X)return out;
  if(std::abs(contrast)<32*std::numeric_limits<double>::epsilon())
    throw std::domain_error("composition_buoyancy: unresolved pressure contrast across a composition gradient");
  const double h=reuse_spacing.load(std::memory_order_relaxed);
  Contrast c;
  if(h>0) {
    const double lnT=std::log(T),lnP=std::log(P);
    const auto iT=static_cast<std::int64_t>(std::llround(lnT/h)),iP=static_cast<std::int64_t>(std::llround(lnP/h));
    const double lnTa=double(iT)*h,lnPa=double(iP)*h;
    const auto key=make_key(eos,lo,hi,iT,iP);
    auto& s=shard(KeyHash{}(key));
    Contrast anchor;bool found=false;
    {
      std::lock_guard lock(s.mutex);
      if(auto it=s.anchors.find(key);it!=s.anchors.end()){anchor=it->second;found=true;}
    }
    if(!found) {
      try {anchor=exact_contrast(eos,std::exp(lnTa),std::exp(lnPa),lo,hi,guess);}
      catch(const std::domain_error&) {
        const auto e=exact_contrast(eos,T,P,lo,hi,guess);
        const double f=1/(delta*contrast);
        return {e.L*f,e.dlnT*f,e.dlnP*f,-e.L*f/delta,-e.L*f/contrast};
      }
      std::lock_guard lock(s.mutex);
      if(s.anchors.size()>=shard_capacity)s.anchors.clear();
      s.anchors.emplace(key,anchor);
      reuse_exact.fetch_add(1,std::memory_order_relaxed);
    } else reuse_hits.fetch_add(1,std::memory_order_relaxed);
    // Preserve actual-temperature source bounds on reuse. Near a density
    // edge use the exact inversion, with a generous one-percent margin
    // relative to the <=1e-4 logarithmic response step.
    const auto interior=[&](const Composition& composition,double rho,double expansion,double chi) {
      const double predicted=std::log(rho)-expansion*(lnT-lnTa)+(lnP-lnPa)/chi;
      try {
        const auto range=eos.density_range_near(T,composition,std::exp(predicted));
        if(!range)return true;
        return predicted>std::log(range->min)+.01 && predicted<std::log(range->max)-.01;
      } catch(const std::domain_error&) {
        return false; // Let the exact inversion check the actual state.
      }
    };
    if(interior(lo,anchor.rho_lo,anchor.delta_lo,anchor.chi_lo)
        && interior(hi,anchor.rho_hi,anchor.delta_hi,anchor.chi_hi))
      c={anchor.L+anchor.dlnT*(lnT-lnTa)+anchor.dlnP*(lnP-lnPa),anchor.dlnT,anchor.dlnP};
    else c=exact_contrast(eos,T,P,lo,hi,guess);
    if(reuse_verify.load(std::memory_order_relaxed)) {
      const auto e=exact_contrast(eos,T,P,lo,hi,guess);
      // Absolute error in B itself; B enters as grad_ad+B, with grad_ad ~ 0.4.
      const double err=std::abs(c.L-e.L)/std::abs(delta*contrast);
      reuse_verified.fetch_add(1,std::memory_order_relaxed);
      std::lock_guard lock(worst_mutex);worst_error=std::max(worst_error,err);
    }
  } else c=exact_contrast(eos,T,P,lo,hi,guess);
  const double factor=1/(delta*contrast);
  out.B=c.L*factor;
  out.dB_dlnT=c.dlnT*factor;
  out.dB_dlnP=c.dlnP*factor;
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
