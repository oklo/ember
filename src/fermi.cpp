#include "fermi.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace ember::fermi {
namespace {

// Gauss-Legendre nodes, generated once by Newton on the Legendre polynomial so
// there is no table to mistype.
constexpr int NQ = 64;
struct GL {
  std::array<double, NQ> x{}, w{};
  GL() {
    for (int i = 0; i < NQ; ++i) {
      double z = std::cos(M_PI * (i + 0.75) / (NQ + 0.5));
      for (int it = 0; it < 100; ++it) {
        double p0 = 1.0, p1 = 0.0;
        for (int j = 0; j < NQ; ++j) { const double p2 = p1; p1 = p0;
          p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1); }
        const double dp = NQ * (z * p0 - p1) / (z * z - 1.0);
        const double dz = p0 / dp; z -= dz;
        if (std::abs(dz) < 1e-15) break;
      }
      double p0 = 1.0, p1 = 0.0;
      for (int j = 0; j < NQ; ++j) { const double p2 = p1; p1 = p0;
        p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1); }
      const double dp = NQ * (z * p0 - p1) / (z * z - 1.0);
      x[i] = z; w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
  }
};
const GL& gl() { static const GL g; return g; }

} // namespace

// One Gauss-Legendre panel in t, where the kinetic energy is g = t^2.  The
// square-root substitution removes the dx/dg singularity at the origin, and
// panelling at the Fermi surface is what makes strong degeneracy tractable:
// there the occupancy falls from one to zero across a shell of width ~beta in
// g, and at beta ~ 1e-4 a single panel spread over the whole range cannot see
// it at all - the derivative integrals come back as zero and the Newton solve
// divides by them.
static void accumulate(Integrals& R, double t0, double t1,
                       double eta, double beta) {
  if (!(t1 > t0)) return;
  const GL& q = gl();
  const double h = 0.5 * (t1 - t0), m = 0.5 * (t1 + t0);
  for (int i = 0; i < NQ; ++i) {
    const double t  = m + h * q.x[i];
    const double wt = q.w[i] * h;
    const double t2 = t * t;
    const double g  = t2;                       // kinetic energy / (m c^2)
    const double sr = std::sqrt(t2 + 2.0);      // x = t * sr
    const double e  = g / beta - eta;

    const double z = std::exp(-std::abs(e));
    const double f = e > 0.0 ? z / (1.0 + z) : 1.0 / (1.0 + z);
    const double fp = z / ((1.0 + z) * (1.0 + z));

    // Jacobians already folded in (dg = 2t dt, dx = (1+g)/x dg):
    const double kn = wt * 2.0 * t2 * sr * (1.0 + g);
    const double kp = wt * 2.0 * t2 * t2 * sr * (t2 + 2.0);
    const double ku = kn * g;

    R.In += kn * f;  R.Ip += kp * f;  R.Iu += ku * f;
    R.dIn_deta += kn * fp;  R.dIp_deta += kp * fp;  R.dIu_deta += ku * fp;
    const double gb = g / beta;
    const double s = e >= 0.0 ? (1.0-z)/(1.0+z) : (z-1.0)/(1.0+z);
    R.d2In_deta2 += kn * fp * s;
    R.d2In_detadlnb += kn * fp * s * gb;
    R.dIn_dlnb += kn * fp * gb;
    R.dIp_dlnb += kp * fp * gb;
    R.dIu_dlnb += ku * fp * gb;
  }
}

Integrals evaluate(double eta, double beta) {
  Integrals R{};
  // Fermi kinetic energy, and a window wide enough that the occupancy at its
  // edges is below 1e-13.
  const double gF   = (eta > 0.0) ? beta * eta : 0.0;
  const double wide = 30.0 * beta;
  const double glo  = (gF > wide) ? gF - wide : 0.0;
  const double ghi  = gF + wide;
  accumulate(R, 0.0, std::sqrt(glo), eta, beta);   // deep interior, f = 1
  // Split the Fermi shell at its center as well. A single 64-node shell was
  // adequate for bulk P, but its derivative error is amplified by the second
  // derivatives needed for convection's response to density.
  accumulate(R, std::sqrt(glo), std::sqrt(gF), eta, beta);
  accumulate(R, std::sqrt(gF), std::sqrt(ghi), eta, beta);
  return R;
}

double entropy(double eta, double beta, const Integrals& f) {
  if(eta<=100.) return (f.Iu+f.Ip/3.)/beta-eta*f.In;
  // The occupation entropy is even in e=(g-gF)/beta and decays
  // exponentially. Pair both sides of the Fermi surface in e, retaining
  // the exact relativistic density of states. eta>100 keeps gF-64 beta>0.
  const double gf=eta*beta;
  const auto& q=gl();
  double sum=0;
  for(int i=0;i<NQ;++i) {
    const double e=32*(1+q.x[i]),z=std::exp(-e);
    const double occupation_entropy=std::log1p(z)+e*z/(1+z);
    const double gp=gf+beta*e,gm=gf-beta*e;
    const double density=(1+gp)*std::sqrt(gp*(gp+2))
                        +(1+gm)*std::sqrt(gm*(gm+2));
    sum+=32*q.w[i]*beta*density*occupation_entropy;
  }
  return sum;
}

// In a cold Fermi gas, the small thermal response comes from paired points
// on either side of the Fermi surface. Pair in e=(g-gF)/beta, and evaluate
// kernel sums/differences algebraically. Subtracting independently rounded
// g coordinates or odd quadrature sums loses the mixed density derivative.
static DensityResponse degenerate_density_response(double eta,double beta,
    const Integrals& integrals,bool second) {
  const double gf=eta*beta,vf=1+gf;
  struct Node {double e,s,sum,difference,dp_sum,dp_difference;};
  std::array<Node,NQ> nodes;
  double weight=0,moment=0;
  const auto& q=gl();
  for(int i=0;i<NQ;++i) {
    const double e=32*(1+q.x[i]),du=beta*e;
    const double gp=gf+du,gm=gf-du;
    const double xp=std::sqrt(gp*(gp+2)),xm=std::sqrt(gm*(gm+2));
    const double z=std::exp(-e),fp=z/((1+z)*(1+z));
    const double w=32*q.w[i]*beta*fp;
    const double sum=w*((vf+du)*xp+(vf-du)*xm);
    const double difference=w*2*du*((xp+xm)/2+2*vf*vf/(xp+xm));
    const double product=vf*vf-du*du;
    const double dp_sum=-2*du*du/(vf*product);
    const double dp_difference=2*du*(1+1/product);
    nodes[i]={e,(1-z)/(1+z),sum,difference,dp_sum,dp_difference};
    weight+=sum;moment+=e*difference;
  }
  // eta+eta_T = -<e>, evaluated without subtracting two large eta terms.
  const double shift=-moment/weight,eta_rho=integrals.In/weight;
  double B=0,C=0,D=0;
  if(second) {
    for(const auto& n:nodes) {
      B+=n.s*((n.e*n.e+shift*shift)*n.difference+2*n.e*shift*n.sum);
      C+=n.s*(n.e*n.sum+shift*n.difference);
      D+=n.s*n.difference;
    }
    B/=weight;C/=weight;D/=weight;
  }
  DensityResponse out{};
  if(second)out.d2Ip_dlnRho2=gf*(gf+2)/(gf+1)*integrals.In;
  for(const auto& n:nodes) {
    const double u_sum=beta*n.e*n.difference,u_difference=beta*n.e*n.sum;
    const double p_sum=(n.sum*n.dp_sum+n.difference*n.dp_difference)/2;
    const double p_difference=(n.sum*n.dp_difference+n.difference*n.dp_sum)/2;
    out.dIu_dlnT+=n.e*u_difference+shift*u_sum;
    out.dIp_dlnT+=n.e*p_difference+shift*p_sum;
    if(!second)continue;
    const double odd=n.s*(n.e*n.e+shift*shift)-n.e;
    const double even=2*n.s*n.e*shift-shift-B;
    out.d2Iu_dlnT2+=odd*u_difference+even*u_sum;
    out.d2Ip_dlnT2+=odd*p_difference+even*p_sum;
    out.d2Iu_dlnTdlnRho+=eta_rho*((n.s*n.e-C)*u_sum+n.s*shift*u_difference);
    out.d2Ip_dlnTdlnRho+=eta_rho*((n.s*n.e-C)*p_sum+n.s*shift*p_difference);
    out.d2Ip_dlnRho2+=eta_rho*((1-eta_rho*D)*p_sum+eta_rho*n.s*p_difference);
  }
  return out;
}

DensityResponse density_response(double eta, double beta, const Integrals& integrals, bool second) {
  if(eta>100.)return degenerate_density_response(eta,beta,integrals,second);
  const double gf = beta * std::max(eta, 0.0);
  const double tlo = std::sqrt(std::max(gf - 30.0 * beta, 0.0));
  const double thi = std::sqrt(gf + 30.0 * beta);
  const double eta_T = -integrals.dIn_dlnb / integrals.dIn_deta;
  const double eta_rho = integrals.In / integrals.dIn_deta;

  // With h=g/beta+eta_T and s=1-2f, density conservation gives the
  // occupation responses fp*(s*h^2-h-B), fp*eta_rho*(s*h-C), and
  // fp*(eta_rho + eta_rho^2*(s-D)). B,C,D are fp-weighted number moments.
  auto each_node = [&](const auto& visit) {
    const auto& quadrature = gl();
    for (const auto interval : {std::array{0.0, tlo}, std::array{tlo, std::sqrt(gf)},
                               std::array{std::sqrt(gf), thi}}) {
      const double half = 0.5 * (interval[1] - interval[0]);
      if (!(half > 0.0)) continue;
      const double mid = 0.5 * (interval[1] + interval[0]);
      for (int i = 0; i < NQ; ++i) {
        const double t = mid + half * quadrature.x[i], g = t * t;
        const double e = g / beta - eta;
        const double z = std::exp(-std::abs(e));
        const double fp = z / ((1.0 + z) * (1.0 + z));
        const double s = e >= 0.0 ? (1.0 - z) / (1.0 + z) : (z - 1.0) / (1.0 + z);
        const double hT = eta > 0.0 ? (g - gf) / beta + (eta + eta_T) : g / beta + eta_T;
        const double weight = quadrature.w[i] * half * 2.0 * g * std::sqrt(g + 2.0) * (1.0 + g) * fp;
        visit(g, weight, s, hT);
      }
    }
  };
  double B = 0.0, C = 0.0, D = 0.0, weight_sum = 0.0;
  if(second) {
    each_node([&](double, double w, double s, double hT) {
      weight_sum += w;
      B += w * s * hT * hT; C += w * s * hT; D += w * s;
    });
    B /= weight_sum; C /= weight_sum; D /= weight_sum;
  }

  DensityResponse out{};
  // Subtract a constant multiple of the number kernel. Its T, TT and T-rho
  // derivatives vanish. In rho-rho its derivative is n, restored below.
  // The pressure-kernel difference has this cancellation-free expression.
  const double pressure_ref = gf * (gf + 2.0) / (gf + 1.0);
  out.d2Ip_dlnRho2 = pressure_ref * integrals.In;
  each_node([&](double g, double w, double s, double hT) {
    const double pressure = (g - gf) * (1.0 + 1.0 / ((1.0 + g) * (1.0 + gf)));
    const double energy = g - gf;
    // The constant number-kernel contribution vanishes at fixed density.
    // Integrate the thermal response directly instead of subtracting two
    // O(eta) terms to recover an O(1/eta) cold-electron heat capacity.
    out.dIp_dlnT += w * pressure * hT;
    out.dIu_dlnT += w * energy * hT;
    if(!second) return;
    const double tt = s * hT * hT - hT - B;
    const double tr = eta_rho * (s * hT - C);
    const double rr = eta_rho * (1.0 + eta_rho * (s - D));
    out.d2Ip_dlnT2 += w * pressure * tt;
    out.d2Ip_dlnTdlnRho += w * pressure * tr;
    out.d2Ip_dlnRho2 += w * pressure * rr;
    out.d2Iu_dlnT2 += w * energy * tt;
    out.d2Iu_dlnTdlnRho += w * energy * tr;
  });
  return out;
}


} // namespace ember::fermi
