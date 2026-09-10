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

DensityResponse density_response(double eta, double beta, const Integrals& integrals) {
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
  each_node([&](double, double w, double s, double hT) {
    weight_sum += w;
    B += w * s * hT * hT; C += w * s * hT; D += w * s;
  });
  B /= weight_sum; C /= weight_sum; D /= weight_sum;

  DensityResponse out{};
  // Subtract a constant multiple of the number kernel. Its T, TT and T-rho
  // derivatives vanish. In rho-rho its derivative is n, restored below.
  // The pressure-kernel difference has this cancellation-free expression.
  const double pressure_ref = gf * (gf + 2.0) / (gf + 1.0);
  out.d2Ip_dlnRho2 = pressure_ref * integrals.In;
  each_node([&](double g, double w, double s, double hT) {
    const double pressure = (g - gf) * (1.0 + 1.0 / ((1.0 + g) * (1.0 + gf)));
    const double energy = g - gf;
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
