#include "fermi.hpp"
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

    double f;
    if (e > 0.0) { const double z = std::exp(-e); f = z / (1.0 + z); }
    else         { const double z = std::exp(e);  f = 1.0 / (1.0 + z); }
    const double fp = f * (1.0 - f);

    // Jacobians already folded in (dg = 2t dt, dx = (1+g)/x dg):
    const double kn = wt * 2.0 * t2 * sr * (1.0 + g);
    const double kp = wt * 2.0 * t2 * t2 * sr * (t2 + 2.0);
    const double ku = kn * g;

    R.In += kn * f;  R.Ip += kp * f;  R.Iu += ku * f;
    R.dIn_deta += kn * fp;  R.dIp_deta += kp * fp;  R.dIu_deta += ku * fp;
    const double gb = g / beta;
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
  accumulate(R, std::sqrt(glo), std::sqrt(ghi), eta, beta); // the Fermi shell
  return R;
}


} // namespace ember::fermi
